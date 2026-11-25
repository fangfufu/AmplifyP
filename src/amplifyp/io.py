# Copyright (C) 2026 AmplifyP Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Primer I/O, validation, reaction input, and YAML state serialisation."""

from __future__ import annotations

import csv
import functools
import io
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import yaml

from .dimer import PrimerDimerGenerator
from .dna import (
    DNA,
    DNAType,
    Nucleotides,
    Primer,
    clean_sequence,
    format_sequence,
)
from .pcr import PCR
from .settings import (
    GLOBAL_PRIMER_DIMER_SETTINGS,
    GLOBAL_REPLICATION_SETTINGS,
    PrimerDimerSettings,
    ReplicationSettings,
    TMSettings,
    coerce_bool,
)

__all__ = [
    "ProjectState",
    "ReactionInput",
    "count_active_primer_names_and_sequences",
    "get_duplicate_primer_indices",
    "load_state_yaml",
    "parse_primer_clipboard_text",
    "parse_primers_from_text",
    "parse_state_yaml",
    "serialise_primers_to_tsv",
    "serialise_state",
    "validate_primer",
    "validate_primers",
]

_PRIMER_VALID_SET = frozenset(Nucleotides.PRIMER)


@functools.lru_cache(maxsize=4096)
def validate_primer(
    name: str, seq: str, show_empty_errors: bool = False
) -> tuple[str | None, str | None]:
    """Validate a single primer name and DNA sequence.

    Args:
        name: The primer name to validate.
        seq: The primer DNA sequence to validate.
        show_empty_errors: If True, reports empty name or sequence as errors.

    Returns:
        A tuple of ``(name_error, seq_error)`` where each element is an error
        message string or ``None`` if valid.
    """
    name_err: str | None = None
    seq_err: str | None = None

    if not name.strip() and show_empty_errors:
        name_err = "Name cannot be empty"

    cleaned_seq = clean_sequence(seq)
    if not cleaned_seq:
        if show_empty_errors:
            seq_err = "Sequence cannot be empty"
    else:
        try:
            Primer(sequence=cleaned_seq, name=name)
        except ValueError as ex:
            seq_err = str(ex)

    return name_err, seq_err


def count_active_primer_names_and_sequences(
    primers: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, int], dict[str, int]]:
    """Count active primer names (case-sensitive) and sequences (uppercased)."""
    active_names_count: dict[str, int] = {}
    active_seqs_count: dict[str, int] = {}

    for p in primers:
        if not p.get("active", True):
            continue
        name_key = str(p.get("name") or "").strip()
        seq_key = clean_sequence(str(p.get("seq") or "")).upper()

        if name_key:
            active_names_count[name_key] = (
                active_names_count.get(name_key, 0) + 1
            )
        if seq_key:
            active_seqs_count[seq_key] = active_seqs_count.get(seq_key, 0) + 1

    return active_names_count, active_seqs_count


def get_duplicate_primer_indices(
    primers: Sequence[Mapping[str, Any]],
) -> set[int]:
    """Find indices of duplicate active primers by name or sequence."""
    active_names_count, active_seqs_count = (
        count_active_primer_names_and_sequences(primers)
    )

    dup_indices: set[int] = set()
    for idx, p in enumerate(primers):
        if not p.get("active", True):
            continue
        name_key = str(p.get("name") or "").strip()
        seq_key = clean_sequence(str(p.get("seq") or "")).upper()

        is_name_dup = bool(name_key) and active_names_count.get(name_key, 0) > 1
        is_seq_dup = bool(seq_key) and active_seqs_count.get(seq_key, 0) > 1
        if is_name_dup or is_seq_dup:
            dup_indices.add(idx)
    return dup_indices


def validate_primers(
    primers: Sequence[Mapping[str, Any]],
    check_empty: bool = True,
    check_duplicates: bool = True,
) -> list[dict[str, str | None]]:
    """Validate a sequence of primer dicts for format and duplicate errors."""
    active_names_count, active_seqs_count = (
        count_active_primer_names_and_sequences(primers)
    )

    errors: list[dict[str, str | None]] = []
    for p in primers:
        name_val = str(p.get("name") or "").strip()
        seq_val = clean_sequence(str(p.get("seq") or ""))

        is_active = bool(p.get("active", True))
        show_empty_errors = check_empty and (
            is_active or bool(p.get("show_empty_errors", False))
        )
        name_err, seq_err = validate_primer(
            name_val, seq_val, show_empty_errors=show_empty_errors
        )

        seq_key = seq_val.upper()

        if check_duplicates and is_active and not seq_err and seq_key:
            if active_seqs_count.get(seq_key, 0) > 1:
                seq_err = "Duplicate primer sequence"

        if check_duplicates and is_active and not name_err and name_val:
            if active_names_count.get(name_val, 0) > 1:
                name_err = "Duplicate primer name"

        errors.append({"name": name_err, "seq": seq_err})
    return errors


def parse_primers_from_text(
    content: str,
    csv_reader: Callable[..., Iterable[list[str]]] = csv.reader,
) -> list[dict[str, Any]]:
    """Parse primers from CSV or TSV text.

    Expects rows ordered as ``sequence, name[, extra_description]``, matching
    :func:`serialise_primers_to_tsv`. Supports both comma-separated and
    tab-separated formats, skips blank lines and lines starting with ``#``,
    and appends an optional third column to the primer name.

    Args:
        content: Raw CSV/TSV text content.
        csv_reader: Optional CSV reader callable (defaults to ``csv.reader``).

    Returns:
        List of primer dicts with ``name``, ``seq``, and ``active`` keys.
    """
    parsed_primers: list[dict[str, Any]] = []
    for raw_line in content.splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        line = raw_line if "\t" in raw_line else stripped
        delimiter = "\t" if "\t" in line else ","
        reader = csv_reader([line], delimiter=delimiter)
        try:
            row = next(iter(reader))
        except (StopIteration, csv.Error):
            continue

        if not row or len(row) < 2:
            continue

        seq = clean_sequence(row[0])
        name = row[1].strip()

        if len(row) >= 3:
            extra = row[2].strip()
            if extra:
                name = f"{name} - {extra}"

        if not name and not seq:
            continue

        name_err, seq_err = validate_primer(name, seq, show_empty_errors=False)
        is_valid_primer = name_err is None and seq_err is None

        parsed_primers.append(
            {
                "name": name,
                "seq": seq,
                "active": is_valid_primer,
            }
        )
    return parsed_primers


def serialise_primers_to_tsv(primers: Sequence[Mapping[str, Any]]) -> str:
    r"""Serialise primer dicts to a tab-separated ``sequence\tname`` string."""
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    for p in primers:
        seq = clean_sequence(str(p.get("seq", "")))
        name = " ".join(str(p.get("name", "")).splitlines()).strip()
        writer.writerow([seq, name])
    tsv_content = output.getvalue()
    output.close()
    return tsv_content


def parse_primer_clipboard_text(text: str) -> list[dict[str, Any]]:
    r"""Parse pasted clipboard text containing primers.

    Supports ``name\tsequence`` (tab-separated), ``name,sequence``
    (comma-separated), single-sequence, or single-name lines. Compare with
    :func:`parse_primers_from_text`, which reads ``sequence, name`` file rows.

    Args:
        text: The raw clipboard or pasted text.

    Returns:
        A list of dicts with keys ``name``, ``seq``, and ``active``.
    """
    parsed: list[dict[str, Any]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            name = parts[0].strip()
            seq = clean_sequence(parts[1])
        elif len(parts) == 1:
            subparts = line.split(",")
            if len(subparts) >= 2:
                name = subparts[0].strip()
                seq = clean_sequence(subparts[1])
            else:
                val = line.strip()
                cleaned = clean_sequence(val)
                is_seq = bool(cleaned) and all(
                    c in _PRIMER_VALID_SET for c in cleaned.upper()
                )
                if is_seq:
                    name = ""
                    seq = cleaned
                else:
                    name = val
                    seq = ""
        else:
            continue
        parsed.append({"name": name, "seq": seq, "active": False})
    return parsed


class ReactionInput:
    """Encapsulates DNA template and primers list for PCR/dimer workflows."""

    def __init__(
        self,
        template: str = "",
        template_circular: bool = False,
        primers: list[dict[str, Any]] | None = None,
    ) -> None:
        """Initialise ReactionInput."""
        self.template: str = template
        self.template_circular: bool = template_circular
        self.primers: list[dict[str, Any]] = (
            primers if primers is not None else []
        )

    def get_active_primers(self) -> list[dict[str, Any]]:
        """Get the active primers with non-empty names and sequences."""
        return [
            p
            for p in self.primers
            if p.get("active", True)
            and str(p.get("name") or "").strip()
            and clean_sequence(str(p.get("seq") or "")).strip()
        ]

    def _parse_primers(
        self, primers: Sequence[Mapping[str, Any]]
    ) -> list[dict[str, Any]]:
        """Validate and normalise primer entries without mutating self."""
        primers_obj = cast("object", primers)
        if not isinstance(primers_obj, Sequence) or isinstance(
            primers_obj, (str, bytes)
        ):
            raise ValueError("Invalid primers format: expected a sequence.")
        parsed_primers: list[dict[str, Any]] = []
        for raw_p in cast("Sequence[object]", primers_obj):
            if not isinstance(raw_p, Mapping):
                raise ValueError("Invalid primer entry: expected a mapping.")
            p = cast("Mapping[str, Any]", raw_p)
            p_seq = p.get("seq", "")
            if isinstance(p_seq, list):
                seq_items = cast("Sequence[object]", p_seq)
                seq_str = "".join(str(x) for x in seq_items)
            else:
                seq_str = str(p_seq)
            parsed_primers.append(
                {
                    "name": str(p.get("name", "")),
                    "seq": clean_sequence(seq_str),
                    "active": bool(p.get("active", True)),
                }
            )
        return parsed_primers

    def to_dict(self) -> dict[str, Any]:
        """Convert input data to a dictionary for serialisation."""
        return {
            "template": format_sequence(self.template),
            "template_circular": self.template_circular,
            "primers": [
                {
                    "name": str(p.get("name", "")),
                    "seq": format_sequence(str(p.get("seq", ""))),
                    "active": bool(p.get("active", True)),
                }
                for p in self.primers
            ],
        }

    def from_dict(self, state_dict: Mapping[str, Any]) -> None:
        """Update input data from a dictionary."""
        state_obj = cast("object", state_dict)
        if not isinstance(state_obj, Mapping):
            raise ValueError("Invalid input state: expected a mapping.")
        state_map = cast("Mapping[str, Any]", state_obj)
        new_template = self.template
        new_circular = self.template_circular
        new_primers = self.primers
        if "template" in state_map:
            raw_tpl = state_map["template"]
            new_template = (
                "" if raw_tpl is None else clean_sequence(str(raw_tpl))
            )
        if "template_circular" in state_map:
            new_circular = coerce_bool(state_map["template_circular"], False)
        if "primers" in state_map:
            raw_primers = state_map["primers"]
            new_primers = self._parse_primers(
                [] if raw_primers is None else raw_primers
            )
        self.template = new_template
        self.template_circular = new_circular
        self.primers = new_primers

    def activate_only_primers(
        self, candidates: Sequence[tuple[str, str]]
    ) -> None:
        """Activate only the specified ``(name, seq)`` candidate primers.

        Deactivates all existing entries, then activates matching sequences
        (case-insensitively updating their name) or appends new entries for
        unmatched sequences.
        """
        for p in self.primers:
            p["active"] = False

        claimed_indices: set[int] = set()
        for name, seq in candidates:
            clean_seq = clean_sequence(seq)
            clean_seq_upper = clean_seq.upper()
            if not clean_seq:
                continue
            matched = False
            for idx, p in enumerate(self.primers):
                if idx in claimed_indices:
                    continue
                if (
                    clean_sequence(str(p.get("seq", ""))).upper()
                    == clean_seq_upper
                ):
                    p["active"] = True
                    p["name"] = name
                    claimed_indices.add(idx)
                    matched = True
                    break
            if not matched:
                self.primers.append(
                    {"name": name, "seq": clean_seq, "active": True}
                )
                claimed_indices.add(len(self.primers) - 1)

    def to_template_dna(self) -> DNA:
        """Construct a core ``DNA`` template from the stored sequence."""
        clean_tpl = clean_sequence(self.template)
        t_type = DNAType.CIRCULAR if self.template_circular else DNAType.LINEAR
        return DNA(clean_tpl, dna_type=t_type)

    def to_pcr(
        self, settings: ReplicationSettings = GLOBAL_REPLICATION_SETTINGS
    ) -> PCR:
        """Construct a ``PCR`` reaction populated with active primers."""
        pcr = PCR(self.to_template_dna(), settings=settings)
        for p in self.get_active_primers():
            name = str(p["name"])
            seq = clean_sequence(str(p["seq"]))
            pcr.add_primer(Primer(sequence=seq, name=name))
        return pcr

    def to_dimer_generator(
        self, settings: PrimerDimerSettings = GLOBAL_PRIMER_DIMER_SETTINGS
    ) -> PrimerDimerGenerator:
        """Construct a ``PrimerDimerGenerator`` from active primers."""
        generator = PrimerDimerGenerator(settings=settings)
        for p in self.get_active_primers():
            name = str(p["name"])
            seq = clean_sequence(str(p["seq"]))
            generator.add_primer(Primer(sequence=seq, name=name))
        return generator


def serialise_state(state: Mapping[str, object]) -> str:
    """Serialise state dict to YAML string, handling multiline strings.

    Uses a custom YAML dumper that represents multiline strings using
    the '|' block style for better readability.

    Args:
        state: The state dictionary to serialise.

    Returns:
        A YAML string representation of the state.
    """

    def multiline_presenter(
        dumper: yaml.SafeDumper, data: str
    ) -> yaml.ScalarNode:
        """Represent multiline strings using the '|' block style in YAML."""
        represent = cast(
            "Callable[..., yaml.ScalarNode]",
            dumper.represent_scalar,  # pyright: ignore[reportUnknownMemberType]
        )
        text = str(data)
        if "\n" in text:
            return represent("tag:yaml.org,2002:str", text, style="|")
        return represent("tag:yaml.org,2002:str", text)

    class _StateDumper(yaml.SafeDumper):
        """Custom SafeDumper using multiline representer for str/StrEnum."""

    _StateDumper.add_representer(str, multiline_presenter)
    _StateDumper.add_multi_representer(str, multiline_presenter)
    return yaml.dump(dict(state), Dumper=_StateDumper, sort_keys=False)


def _empty_settings_dict() -> dict[str, Any]:
    """Return an empty typed settings dictionary."""
    return {}


_TYPED_SETTING_PREFIXES = (
    "bp_score_",
    "rep_pd_score_",
    "match_weight_",
    "run_weight_",
    "pd_score_",
    "tm_enthalpy_",
    "tm_entropy_",
)


@dataclass(slots=True)
class ProjectState:
    """Represents a complete serialisable AmplifyP project state."""

    input: ReactionInput = field(default_factory=ReactionInput)
    replication_settings: ReplicationSettings = field(
        default_factory=ReplicationSettings
    )
    primer_dimer_settings: PrimerDimerSettings = field(
        default_factory=PrimerDimerSettings
    )
    tm_settings: TMSettings = field(default_factory=TMSettings)
    raw_settings: dict[str, Any] = field(default_factory=_empty_settings_dict)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ProjectState:
        """Construct a ProjectState from a parsed YAML/JSON mapping."""
        reaction_input = ReactionInput()
        if "input" in data:
            raw_input = data["input"]
            if isinstance(raw_input, Mapping):
                reaction_input.from_dict(cast("Mapping[str, Any]", raw_input))
            elif raw_input is not None:
                raise ValueError(
                    "Invalid state format: 'input' must be a mapping."
                )
        else:
            reaction_input.from_dict(data)

        raw_settings_obj = data.get("settings")
        if raw_settings_obj is not None and not isinstance(
            raw_settings_obj, Mapping
        ):
            raise ValueError(
                "Invalid state format: 'settings' must be a mapping."
            )
        raw_settings: dict[str, Any] = (
            dict(cast("Mapping[str, Any]", raw_settings_obj))
            if isinstance(raw_settings_obj, Mapping)
            else {}
        )
        return cls(
            input=reaction_input,
            replication_settings=ReplicationSettings.from_dict(raw_settings),
            primer_dimer_settings=PrimerDimerSettings.from_dict(raw_settings),
            tm_settings=TMSettings.from_dict(raw_settings),
            raw_settings=raw_settings,
        )

    def to_dict(self, include_settings: bool = True) -> dict[str, Any]:
        """Serialise the project state to a dictionary."""
        result: dict[str, Any] = {"input": self.input.to_dict()}
        if include_settings:
            combined_settings = {
                k: v
                for k, v in self.raw_settings.items()
                if not (
                    isinstance(cast("object", k), str)
                    and k.startswith(_TYPED_SETTING_PREFIXES)
                )
            }
            combined_settings.update(self.replication_settings.to_dict())
            combined_settings.update(self.primer_dimer_settings.to_dict())
            combined_settings.update(self.tm_settings.to_dict())
            if "amp4_compat" in combined_settings:
                combined_settings["amp4_compat"] = (
                    self.replication_settings.amplify4_compatibility_mode
                )
            result["settings"] = combined_settings
        return result

    def to_yaml(self, include_settings: bool = True) -> str:
        """Serialise the project state to a YAML string."""
        return serialise_state(self.to_dict(include_settings=include_settings))


def parse_state_yaml(content: str) -> ProjectState:
    """Parse a YAML state string into a ``ProjectState``.

    Args:
        content: Raw YAML state string.

    Returns:
        A parsed ``ProjectState`` instance.

    Raises:
        yaml.YAMLError: If the YAML syntax is invalid.
        ValueError: If the root value is empty/comment-only, not a mapping,
            or contains malformed ``input`` / ``settings`` fields.
    """
    parsed = yaml.safe_load(content)
    if not isinstance(parsed, dict):
        raise ValueError("Invalid state YAML format: expected a mapping.")
    return ProjectState.from_dict(cast("dict[str, Any]", parsed))


def load_state_yaml(path: str | Path) -> ProjectState:
    """Load a YAML state file from disk into a ``ProjectState``.

    Args:
        path: Path to the YAML state file.

    Returns:
        A parsed ``ProjectState`` instance.

    Raises:
        OSError: If the file cannot be read.
        yaml.YAMLError: If the YAML syntax is invalid.
        ValueError: If the parsed content is not a valid state mapping.
    """
    with open(path, encoding="utf-8") as f:
        return parse_state_yaml(f.read())
