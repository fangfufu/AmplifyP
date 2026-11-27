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

"""Replication configuration-related classes for AmplifyP."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import cached_property, lru_cache

from .dir_idx import DirIdx as DirIdx
from .dir_idx import DirIdxDb as DirIdxDb
from .dna import (
    DNA,
    GLOBAL_COMPLEMENT_TABLE,
    DNADirection,
    DNAType,
    Nucleotides,
    Primer,
)
from .origin import ReplicationOrigin
from .settings import (
    GLOBAL_REPLICATION_SETTINGS,
    ReplicationSettings,
)

_GAP_BYTE = Nucleotides.GAP.encode("ascii")


@dataclass(slots=True, frozen=True)
class OriginContextAlignment:
    """Formatted alignment of a primer binding site with template context.

    Attributes:
        coordinates_line: 1-based start and end genomic coordinate header.
        pipe_line: Vertical markers pointing to binding boundaries.
        arrow_line: Arrow markers pointing to binding boundaries.
        primer_line: Formatted primer label and sequence line.
        bonds_line: Base-pairing bond symbol line.
        complement_line: Optional 3'-5' template complement line.
        template_line: Formatted 5'-3' template context sequence line.
    """

    coordinates_line: str
    pipe_line: str
    arrow_line: str
    primer_line: str
    bonds_line: str
    complement_line: str | None
    template_line: str

    @property
    def top_block(self) -> str:
        """Return the combined coordinate, pipe, and arrow header lines."""
        return f"{self.coordinates_line}\n{self.pipe_line}\n{self.arrow_line}"

    @property
    def bottom_block(self) -> str:
        """Return the bonds, optional complement, and template lines."""
        if self.complement_line is not None:
            return (
                f"{self.bonds_line}\n"
                f"{self.complement_line}\n"
                f"{self.template_line}"
            )
        return f"{self.bonds_line}\n{self.template_line}"

    def to_text(self) -> str:
        """Render the full origin context alignment as a multi-line string."""
        return f"{self.top_block}\n{self.primer_line}\n{self.bottom_block}"


class _SettingsWrapper:
    """Wrapper to make settings hashable by content for lru_cache."""

    __slots__ = ("_key", "settings")

    def __init__(self, settings: ReplicationSettings) -> None:
        """Initialise the wrapper with an immutable snapshot of scoring tables.

        Args:
            settings (ReplicationSettings): The settings object to wrap.
        """
        self.settings = settings
        self._key = (
            settings.base_pair_scores.cache_key(),
            settings.match_weight.cache_key(),
            settings.run_weights.cache_key(),
            settings.primability_cutoff,
        )

    def __hash__(self) -> int:
        """Return the hash of the settings key.

        Returns:
            int: Hash of the settings tuple.
        """
        return hash(self._key)

    def __eq__(self, other: object) -> bool:
        """Check equality based on the immutable settings snapshot.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the wrapped settings have identical scoring values.
        """
        if isinstance(other, _SettingsWrapper):
            return self._key == other._key
        return NotImplemented


@lru_cache(maxsize=128)
def _get_scoring_tables(
    primer_rev_seq: str, settings_wrapper: _SettingsWrapper
) -> tuple[
    float,
    float,
    tuple[tuple[tuple[tuple[float, float], ...], float], ...],
    tuple[float, ...],
    float,
    bytes | None,
    int,
]:
    """Calculate scoring tables, denominators, and pruning metadata.

    This function caches the expensive table construction for primer scoring.

    Args:
        primer_rev_seq (str): The reversed primer sequence (uppercase).
        settings_wrapper (_SettingsWrapper): Wrapped replication settings.

    Returns:
        tuple: (prim_denom, stab_denom, score_pairs, run_weights, prim_req,
            k0_mask, good_acgt_count).
    """
    settings = settings_wrapper.settings
    m = settings.match_weight
    S = settings.base_pair_scores
    r = settings.run_weights

    prim_denom = 0.0
    stab_denom_base = 0.0

    lookup_list: list[tuple[tuple[float, float], ...]] = []
    step_max_list: list[float] = []

    # Prepare set of characters for lookup keys
    lookup_keys = set(Nucleotides.TEMPLATE)
    lookup_keys.update(c.lower() for c in tuple(lookup_keys))

    for k, base_p in enumerate(primer_rev_seq):
        row_max = S.row_max(base_p)
        prim_denom += m[k] * row_max
        stab_denom_base += row_max

        k_lookup = [(0.0, 0.0)] * 256
        for base_t in lookup_keys:
            score = S[base_p, base_t]
            k_lookup[ord(base_t)] = (m[k] * score, score)

        k_tuple = tuple(k_lookup)
        lookup_list.append(k_tuple)
        step_max_list.append(max(pv for pv, _ in k_tuple))

    length = len(primer_rev_seq)
    rem_max_list = [0.0] * length
    acc = 0.0
    for k in range(length - 1, -1, -1):
        rem_max_list[k] = acc
        acc += step_max_list[k]

    score_pairs = tuple(zip(lookup_list, rem_max_list, strict=True))
    run_weights = tuple(r[k] for k in range(length))
    stab_denom = stab_denom_base * (
        run_weights[length - 1] if length > 0 else r[0]
    )

    prim_cutoff = settings.primability_cutoff
    if length > 0 and prim_denom > 0:
        raw_req = prim_cutoff * prim_denom
        prim_req = raw_req - 1e-12 * abs(raw_req) - 1e-12
        rem0 = rem_max_list[0]
        k0_table = lookup_list[0]
        mask_arr = bytearray(256)
        for b in range(256):
            if k0_table[b][0] + rem0 >= prim_req:
                mask_arr[b] = 1
        k0_mask: bytes | None = bytes(mask_arr)
        good_acgt_count = sum(mask_arr[c] for c in (65, 67, 71, 84))
    else:
        prim_req = float("-inf")
        k0_mask = None
        good_acgt_count = 4

    return (
        prim_denom,
        stab_denom,
        score_pairs,
        run_weights,
        prim_req,
        k0_mask,
        good_acgt_count,
    )


class Repliconf:
    """A class representing a Replication Configuration.

    A 'Repliconf' encapsulates the setup required to search for replication
    origins given a specific Primer and Template DNA. It prepares the template
    sequences (including padding and complementary strands) and manages the
    search process.

    Attributes:
        padding_len (int): The length of padding added to the template,
            typically equal to the primer length.
        primer (Primer): The primer sequence being analysed.
        template (DNA): The original template DNA sequence.
        template_seq (dict[DNADirection, str]): A dictionary mapping direction
            (FWD/REV) to the processed (padded/complemented) template strings.
        settings (ReplicationSettings): The configuration settings (weights,
            cutoffs) used for scoring.
        origin_db (DirIdxDb): The results database storing found origin indices.
    """

    def __init__(
        self,
        template: DNA,
        primer: Primer,
        settings: ReplicationSettings = GLOBAL_REPLICATION_SETTINGS,
    ) -> None:
        """Initialise a Repliconf object.

        Args:
            template (DNA): The template DNA sequence to be searched.
            primer (Primer): The primer sequence to search for.
            settings (Settings): The settings used for evaluating origins.

        Raises:
            ValueError: If uppercase normalisation changes the primer length.
        """
        if len(primer.seq_upper) != len(primer):
            msg = (
                "Primer sequence length must be preserved under uppercase "
                "normalisation."
            )
            raise ValueError(msg)

        self.padding_len = len(primer)

        self.primer = primer
        self.template = template

        self.settings = settings
        self.origin_db = DirIdxDb([], [], False)

        # Pre-calculate reversed primer sequence
        self._rev_primer_seq = self.primer.seq[::-1]

    @cached_property
    def template_seq(self) -> dict[DNADirection, str]:
        """Return padded template strings indexed by direction."""
        seqs: dict[DNADirection, str] = {}
        template = self.template
        # Add padding to the 5' end of the template
        if template.type in (DNAType.LINEAR, DNAType.PRIMER):
            gap = Nucleotides.GAP * self.padding_len
            seqs[DNADirection.FWD] = gap + template.seq
            seqs[DNADirection.REV] = template.complement_seq + gap
        elif template.type == DNAType.CIRCULAR:
            # Matches existing behaviour including 0 case where padding is empty
            padding = (
                template.seq[-self.padding_len :]
                if self.padding_len > 0
                else ""
            )
            comp = template.complement_seq
            seqs[DNADirection.FWD] = padding + template.seq
            seqs[DNADirection.REV] = comp + comp[: self.padding_len]
        return seqs

    @cached_property
    def rev_comp_primer(self) -> DNA:
        """Return the reverse complement primer for amplicon generation."""
        return self.primer.reverse_complement()

    def range(self) -> range:
        """Return the range of valid starting indices for search.

        Returns:
            range: A range object covering all valid start positions in the
                padded template.
        """
        if "template_seq" in self.__dict__:
            padded_len = len(self.template_seq[DNADirection.FWD])
        else:
            tpl = self.template
            if tpl.type in (DNAType.LINEAR, DNAType.PRIMER):
                padded_len = self.padding_len + len(tpl)
            elif tpl.type == DNAType.CIRCULAR:
                padded_len = min(self.padding_len, len(tpl)) + len(tpl)
            else:
                padded_len = len(self.template_seq[DNADirection.FWD])
        return range(padded_len - len(self.primer) + 1)

    def slice(self, i: int) -> slice:
        """Create a slice object for extracting a target segment.

        Args:
            i (int): The start index.

        Returns:
            slice: A slice of length equal to the primer length.
        """
        return slice(i, i + len(self.primer))

    def origin(self, var: DirIdx) -> ReplicationOrigin:
        """Construct a ReplicationOrigin object for a specific location.

        This method extracts the target sequence from the pre-processed template
        at the specified location and direction, and creates a ReplicationOrigin
        object to evaluate it.

        Args:
            var (DirIdx): A `DirIdx` object containing both direction and index.

        Returns:
            ReplicationOrigin: An evaluated replication origin object.
        """
        direction = var.direction
        i = var.index

        if direction == DNADirection.FWD:
            # Optimised reverse slicing: template_seq[end-1:start-1:-1]
            # This creates only 1 string object instead of 2 (slice then
            # reverse)
            end = i + len(self.primer)
            target = (
                self.template_seq[direction][end - 1 : i - 1 : -1]
                if i > 0
                else self.template_seq[direction][end - 1 :: -1]
            )
        else:
            target = self.template_seq[direction][self.slice(i)]

        return ReplicationOrigin(
            target,
            self._rev_primer_seq,
            self.settings,
            _primability=var.primability,
            _stability=var.stability,
        )

    def origin_from_db(
        self, direction: DNADirection, i: int
    ) -> ReplicationOrigin:
        """Retrieve a ReplicationOrigin based on database index.

        Args:
            direction (DNADirection): The direction to look up in the database.
            i (int): The index *within the list of found origins* (not the
                genomic index).

        Returns:
            ReplicationOrigin: The replication origin object.
        """
        return self.origin(self.origin_db[direction, i])

    def _get_template_bytes(self, direction: DNADirection) -> bytes:
        """Return uppercase ASCII template bytes for the given direction."""
        if "template_seq" in self.__dict__:
            return self.template_seq[direction].upper().encode("ascii")
        tpl = self.template
        pad_len = self.padding_len
        if tpl.type in (DNAType.LINEAR, DNAType.PRIMER):
            gap_bytes = _GAP_BYTE * pad_len
            if direction == DNADirection.FWD:
                return gap_bytes + tpl.seq_upper_bytes
            return tpl.comp_upper_bytes + gap_bytes
        if tpl.type == DNAType.CIRCULAR:
            if direction == DNADirection.FWD:
                tpl_bytes = tpl.seq_upper_bytes
                padding = tpl_bytes[-pad_len:] if pad_len > 0 else b""
                return padding + tpl_bytes
            comp_bytes = tpl.comp_upper_bytes
            return comp_bytes + comp_bytes[:pad_len]
        return self.template_seq[direction].upper().encode("ascii")

    def search(self) -> None:
        """Perform a search for valid replication origins.

        Scans the entire template in both forward and reverse directions.
        Locations where the primability and stability scores exceed the
        configured cutoffs are stored in `origin_db`.

        Raises:
            ValueError: If uppercase normalisation changes the primer length.
        """
        primer_len = len(self.primer)
        primer_rev = self.primer.seq_upper[::-1]
        if len(primer_rev) != primer_len:
            msg = (
                "Primer sequence length must be preserved under uppercase "
                "normalisation."
            )
            raise ValueError(msg)

        self.origin_db.clear()

        (
            prim_denom,
            stab_denom,
            score_pairs,
            run_weights,
            prim_req,
            k0_mask,
            good_acgt_count,
        ) = _get_scoring_tables(primer_rev, _SettingsWrapper(self.settings))

        prim_cutoff = self.settings.primability_cutoff
        stab_cutoff = self.settings.stability_cutoff

        search_range = self.range()
        n_positions = len(search_range)
        if n_positions <= 0:
            self.origin_db.searched = True
            return

        has_any_good_k0 = k0_mask is not None and b"\x01" in k0_mask
        use_mask = has_any_good_k0 and good_acgt_count <= 2

        for direction in (DNADirection.FWD, DNADirection.REV):
            # FWD direction:
            # - Searches `template_seq[FWD]`, which is the Sense strand (5'-3').
            # - Matches the primer sequence against the Sense strand (reversed).
            # - Effectively finds binding sites for Forward Primers (which are
            #   identical to the Sense strand) by checking for identity.
            #
            # REV direction:
            # - Searches `template_seq[REV]`, which is the Antisense strand
            #   (complement of Sense, read 3'-5' in physical space, but stored
            #   as a string).
            # - Matches the primer sequence against the Antisense strand.
            # - Effectively finds binding sites for Reverse Primers (which are
            #   complementary to the Sense strand) by checking for identity
            #   against the Antisense strand.
            if k0_mask is not None and not has_any_good_k0:
                continue

            seq_bytes = self._get_template_bytes(direction)
            origin_list = (
                self.origin_db.fwd
                if direction == DNADirection.FWD
                else self.origin_db.rev
            )

            if direction == DNADirection.FWD:
                start_offset = primer_len - 1
                step = -1
                k0_offset = primer_len - 1
            else:
                start_offset = 0
                step = 1
                k0_offset = 0

            candidate_indices: Iterable[int]
            if use_mask and k0_mask is not None:
                anchored = seq_bytes[
                    k0_offset : k0_offset + n_positions
                ].translate(k0_mask)
                hits: list[int] = []
                append_hit = hits.append
                pos = anchored.find(b"\x01")
                while pos != -1:
                    append_hit(pos)
                    pos = anchored.find(b"\x01", pos + 1)
                candidate_indices = hits
            else:
                candidate_indices = search_range

            for i in candidate_indices:
                prim_num = 0.0
                stab_num = 0.0
                run_len = 0
                run_score = 0.0
                seq_idx = i + start_offset

                for k_lookup, r_max in score_pairs:
                    prim_val, val = k_lookup[seq_bytes[seq_idx]]
                    prim_num += prim_val
                    if prim_num + r_max < prim_req:
                        break
                    if val > 0:
                        run_len += 1
                        run_score += val
                    elif run_len > 0:
                        stab_num += run_weights[run_len - 1] * run_score
                        run_len = 0
                        run_score = 0.0
                    seq_idx += step
                else:
                    if run_len > 0:
                        stab_num += run_weights[run_len - 1] * run_score

                    primability = (
                        prim_num / prim_denom if prim_denom != 0 else 0.0
                    )
                    stability = (
                        stab_num / stab_denom if stab_denom != 0 else 0.0
                    )

                    if primability > prim_cutoff and stability > stab_cutoff:
                        origin_list.append(
                            DirIdx(direction, i, primability, stability)
                        )

        self.origin_db.searched = True

    @property
    def searched(self) -> bool:
        """Check if a search has been performed.

        Returns:
            bool: True if `search()` has been called, False otherwise.
        """
        return self.origin_db.searched

    def __eq__(self, other: object) -> bool:
        """Check equality with another Repliconf.

        Equality is based on having the same primer and template.

        Args:
            other (object): The object to compare.

        Returns:
            bool: True if equal, False otherwise.
        """
        if not isinstance(other, Repliconf):
            return NotImplemented
        return (
            self.primer == other.primer
            and self.template == other.template
            and self.settings == other.settings
        )

    def __hash__(self) -> int:
        """Compute the hash of the Repliconf.

        Returns:
            int: The hash based on primer and template.
        """
        return hash((self.primer, self.template))

    def __str__(self) -> str:
        """Return a string representation.

        Returns:
            str: Description of the configuration.
        """
        return (
            f"ReplicationConfig: Primer: {self.primer}, Target: {self.template}"
        )

    def format_context_alignment(
        self,
        var: DirIdx,
        primer_name: str | None = None,
        origin: ReplicationOrigin | None = None,
        improved_visualisation: bool = False,
        separator: str = " ",
        context_bases: int = 20,
    ) -> OriginContextAlignment:
        """Format a multi-line alignment showing a binding site in context.

        Args:
            var: Directed index specifying the binding site and strand.
            primer_name: Optional override for the displayed primer name.
            origin: Optional precomputed `ReplicationOrigin` for bond symbols.
            improved_visualisation: If True, includes the complementary strand
                for forward primer alignments.
            separator: Separator between 5'/3' markers and sequences
                (``" "`` or ``"-"``).
            context_bases: Number of flanking template bases on each side.

        Returns:
            OriginContextAlignment containing structured alignment lines.

        Raises:
            ValueError: If ``separator`` is invalid, ``context_bases < 0``,
                or the template/primer sequence is empty.
        """
        if separator not in (" ", "-"):
            raise ValueError("separator must be ' ' or '-'")
        if context_bases < 0:
            raise ValueError("context_bases must be non-negative.")
        sep = separator
        resolved_name = self.primer.name if primer_name is None else primer_name
        resolved_origin = (
            Repliconf.origin(self, var) if origin is None else origin
        )
        primer_len = len(self.primer)
        template_len = len(self.template)
        if template_len == 0 or primer_len == 0:
            raise ValueError(
                "Cannot format context alignment for empty template or primer."
            )

        is_circular = self.template.type == DNAType.CIRCULAR
        if var.direction == DNADirection.FWD:
            raw_start = var.index - primer_len
            start_genomic = (
                raw_start % template_len if is_circular else raw_start
            )
            primer_display_seq = self.primer.seq
            primer_label = resolved_name
            primer_ends = ("5'", "3'")
            strength_display = resolved_origin.binding_strength_str[::-1]
        else:
            start_genomic = (
                var.index % template_len if is_circular else var.index
            )
            primer_display_seq = self.primer.seq[::-1]
            primer_label = f"{resolved_name} (Reverse)"
            primer_ends = ("3'", "5'")
            strength_display = resolved_origin.binding_strength_str

        end_genomic = start_genomic + primer_len
        if is_circular:
            start_num_str = str((start_genomic % template_len) + 1)
            end_num_str = str(((end_genomic - 1) % template_len) + 1)
        else:
            start_num_str = str(
                start_genomic + 1 if start_genomic >= 0 else start_genomic
            )
            end_pos = end_genomic - 1
            end_num_str = str(end_pos + 1 if end_pos >= 0 else end_pos)

        base_label_width = 9 + context_bases
        label_width = max(base_label_width, len(primer_label))
        extra_spaces = label_width - base_label_width
        bind_offset = 12 + context_bases + extra_spaces
        prefix_pad = " " * bind_offset

        primer_line = (
            f"{primer_label:<{label_width}}"
            f"{primer_ends[0]}{sep}{primer_display_seq}{sep}{primer_ends[1]}"
        )
        bonds_line = f"{prefix_pad}{strength_display}"

        if primer_len == 1:
            pipe_line = f"{prefix_pad}|"
            arrow_line = f"{prefix_pad}V"
        else:
            gap_spaces = " " * (primer_len - 2)
            pipe_line = f"{prefix_pad}|{gap_spaces}|"
            arrow_line = f"{prefix_pad}V{gap_spaces}V"

        col1 = max(0, bind_offset - len(start_num_str) // 2)
        col2 = max(0, (bind_offset + primer_len - 1) - len(end_num_str) // 2)
        if primer_len > 1 and col2 < col1 + len(start_num_str) + 1:
            col2 = col1 + len(start_num_str) + 1
        line_len = max(
            bind_offset + primer_len + context_bases,
            col1 + len(start_num_str),
            col2 + len(end_num_str) if primer_len > 1 else 0,
        )
        num_line_chars = [" "] * line_len
        for idx, char in enumerate(start_num_str):
            num_line_chars[col1 + idx] = char
        if primer_len > 1:
            for idx, char in enumerate(end_num_str):
                num_line_chars[col2 + idx] = char
        coordinates_line = "".join(num_line_chars).rstrip()

        upstream_seq = self.template.padded_slice(
            start_genomic - context_bases, context_bases
        )
        binding_seq = self.template.padded_slice(start_genomic, primer_len)
        downstream_seq = self.template.padded_slice(
            start_genomic + primer_len, context_bases
        )

        context_line_prefix = "Context  "
        complement_line: str | None = None
        if improved_visualisation and var.direction == DNADirection.FWD:
            comp_upstream = upstream_seq.translate(GLOBAL_COMPLEMENT_TABLE)
            comp_binding = binding_seq.translate(GLOBAL_COMPLEMENT_TABLE)
            comp_downstream = downstream_seq.translate(GLOBAL_COMPLEMENT_TABLE)
            complement_line = (
                f"{' ' * (9 + extra_spaces)}3'{sep}"
                f"{comp_upstream}{comp_binding}{comp_downstream}{sep}5'"
            )

        template_line = (
            f"{context_line_prefix}{' ' * extra_spaces}"
            f"5'{sep}{upstream_seq}{binding_seq}{downstream_seq}{sep}3'"
        )

        return OriginContextAlignment(
            coordinates_line=coordinates_line,
            pipe_line=pipe_line,
            arrow_line=arrow_line,
            primer_line=primer_line,
            bonds_line=bonds_line,
            complement_line=complement_line,
            template_line=template_line,
        )
