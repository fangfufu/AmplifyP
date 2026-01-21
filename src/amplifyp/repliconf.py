"""Replication configuration-related classes for AmplifyP."""

from __future__ import annotations

import logging
from dataclasses import dataclass

from .dna import DNA, DNADirection, DNAType, Nucleotides, Primer
from .origin import ReplicationOrigin
from .settings import Settings

COMPLEMENT_TABLE = str.maketrans("ACGTMKRYBDHVacgtmkrybdhv", "TGCAKMYRVHDBtgcakmyrvhdb")


@dataclass(slots=True, frozen=True)
class DirIdx:
    """A class for storing the index of a specific direction."""

    direction: DNADirection
    index: int

    def __int__(self) -> int:
        """Return the index as an integer."""
        return self.index

    def __index__(self) -> int:
        """Return the index as an index."""
        return self.index

    def __add__(self, other: object) -> DirIdx:
        """Add an integer to the index."""
        if isinstance(other, int):
            return DirIdx(self.direction, self.index + other)
        if not isinstance(other, DirIdx):
            return NotImplemented
        return DirIdx(self.direction, self.index + other.index)

    def __sub__(self, other: object) -> DirIdx:
        """Subtract an integer from the index."""
        if isinstance(other, int):
            return DirIdx(self.direction, self.index - other)
        if not isinstance(other, DirIdx):
            return NotImplemented
        return DirIdx(self.direction, self.index - other.index)

    def __eq__(self, other: object) -> bool:
        """Check if two DirectionIndex objects are equal."""
        if isinstance(other, int):
            return self.index == other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.direction == other.direction and self.index == other.index

    def __lt__(self, other: object) -> bool:
        """Check if one DirectionIndex is less than another."""
        if isinstance(other, int):
            return self.index < other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index < other.index

    def __gt__(self, other: object) -> bool:
        """Check if one DirectionIndex is greater than another."""
        if isinstance(other, int):
            return self.index > other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index > other.index

    def __le__(self, other: object) -> bool:
        """Check if one DirectionIndex is less than or equal to another."""
        if isinstance(other, int):
            return self.index <= other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index <= other.index

    def __ge__(self, other: object) -> bool:
        """Check if one DirectionIndex is greater than or equal to another."""
        if isinstance(other, int):
            return self.index >= other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index >= other.index

    def __str__(self) -> str:
        """Return the index as a string."""
        return f"{self.index}"


@dataclass(slots=True)
class DirIdxDb:
    """A class for storing the locations of replication origins.

    This class holds the indices of forward and reverse replication origins
    found within a DNA template. It also tracks whether a search has been
    performed.

    Attributes:
        fwd (List[DirectionIndex]): A list of DirectionIndex for forward
        replication origins.
        rev (List[DirectionIndex]): A list of DirectionIndex for reverse
        replication origins.
        searched (bool): A flag indicating whether a search for origins has
                         been completed.
    """

    fwd: list[DirIdx]
    rev: list[DirIdx]
    searched: bool = False

    def clear(self) -> None:
        """Clear the origin index, resetting all lists and flags."""
        self.fwd.clear()
        self.rev.clear()
        self.searched = False

    def __getitem__(self, key: tuple[DNADirection, int]) -> DirIdx:
        """Get the DirectionIndex at the given index."""
        direction, i = key
        if direction == DNADirection.FWD:
            return self.fwd[i]
        else:
            return self.rev[i]


class Repliconf:
    """A class representing a replication configuration.

    A "replication configuration" is defined as a combination of a single
    primer and a padded target DNA sequence. This class is instantiated when
    a new primer is used, as the padding of the target sequence depends on the
    primer length. It handles both the forward and reverse complement strands
    of the target DNA.

    Attributes:
        padding_len (int): The length of the padding, equal to the primer length.
        primer (Primer): The primer used in this configuration.
        template (DNA): The template DNA sequence.
        template_seq (Dict[DNADirection, str]): A dictionary containing the
                                                padded forward and reverse
                                                template sequences.
        settings (Settings): The settings for replication analysis.
        origin_db (DirIdxDb): For storing the indices of replication origins.
    """

    def __init__(
        self,
        template: DNA,
        primer: Primer,
        settings: Settings,
    ) -> None:
        """Construct the Repliconf object.

        Args:
            template (DNA): The template DNA sequence.
            primer (Primer): The primer for this configuration.
            settings (Settings): The replication analysis settings.
        """
        self.padding_len = len(primer)

        self.primer = primer
        self.template = template

        self.template_seq: dict[DNADirection, str] = {}
        # Add padding the 5' end of the template
        if template.type == DNAType.LINEAR:
            self.template_seq[DNADirection.FWD] = (
                Nucleotides.GAP * self.padding_len + template.seq
            )
            self.template_seq[DNADirection.REV] = (
                template.seq + Nucleotides.GAP * self.padding_len
            ).translate(COMPLEMENT_TABLE)
        elif template.type == DNAType.CIRCULAR:
            # Matches existing behavior including 0 case where padding is empty
            padding = template.seq[-self.padding_len :] if self.padding_len > 0 else ""
            self.template_seq[DNADirection.FWD] = padding + template.seq
            self.template_seq[DNADirection.REV] = (
                template.seq + template.seq[: self.padding_len]
            ).translate(COMPLEMENT_TABLE)
        else:
            raise TypeError("Invalid DNA type for padding operation.")

        logging.debug(
            f"Repliconf.__init__(): FWD: {self.template_seq[DNADirection.FWD]}"
        )
        logging.debug(
            f"Repliconf.__init__(): REV: {self.template_seq[DNADirection.REV]}"
        )

        self.settings = settings
        self.origin_db = DirIdxDb([], [], False)

        # Pre-calculate reversed primer sequence
        self._rev_primer_seq = self.primer.seq[::-1]

    def range(self) -> range:
        """Return the valid search range for replication origins.

        The search range is determined by the length of the padded template
        sequence and the length of the primer.

        Returns:
            range: The valid search range.
        """
        return range(len(self.template_seq[DNADirection.FWD]) - len(self.primer) + 1)

    def slice(self, i: int) -> slice:
        """Return a slice object for constructing a ReplicationOrigin.

        This slice is used to extract the target sequence for a potential
        replication origin at a given index.

        Args:
            i (int): The starting index of the slice.

        Returns:
            slice: A slice object.
        """
        return slice(i, i + len(self.primer))

    def origin(self, var: DirIdx | DNADirection, *idx: int) -> ReplicationOrigin:
        """Return the ReplicationOrigin at the given index.

        Args:
            var (DirectionIndex | DNADirection): The direction of the DNA strand.
            *idx (int): The index of the replication origin.

        Returns:
            ReplicationOrigin: The replication origin object.
        """
        if isinstance(var, DirIdx):
            direction = var.direction
            i = var.index
        elif isinstance(var, DNADirection):
            direction = var
            i = idx[0]
        else:
            raise TypeError("var must be DirectionIndex or DNADirection")

        if direction:
            # Optimized reverse slicing: template_seq[end-1:start-1:-1]
            # This creates only 1 string object instead of 2 (slice then reverse)
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
        )

    def origin_from_db(self, direction: DNADirection, i: int) -> ReplicationOrigin:
        """Return the ReplicationOrigin at a given index from the origin database.

        Args:
            direction (DNADirection): The direction of the DNA strand.
            i (int): The index for the replication origin in the origin database.

        Returns:
            ReplicationOrigin: The replication origin object.
        """
        return self.origin(self.origin_db[direction, i])

    def search(self) -> None:
        """Search for valid replication origins in both directions."""
        self.origin_db.clear()

        # Optimization: Pre-calculate constants and lookup tables
        m = self.settings.match_weight
        S = self.settings.base_pair_scores
        r = self.settings.run_weights

        primer_rev = self.primer.seq[::-1]
        L = len(primer_rev)

        prim_denom = 0.0
        stab_denom_base = 0.0

        prim_score_lookup: list[dict[str, float]] = []
        stab_score_lookup: list[dict[str, float]] = []

        # Prepare set of characters for lookup keys
        lookup_keys = set(Nucleotides.LINEAR)
        lookup_keys.update(c.lower() for c in list(lookup_keys))

        for k, base_p in enumerate(primer_rev):
            row_max = S.row_max(base_p)
            prim_denom += m[k] * row_max
            stab_denom_base += row_max

            p_dict = {}
            s_dict = {}

            for base_t in lookup_keys:
                score = S[base_p, base_t]
                p_dict[base_t] = m[k] * score
                s_dict[base_t] = score

            prim_score_lookup.append(p_dict)
            stab_score_lookup.append(s_dict)

        stab_denom = stab_denom_base * r[int(max(0, L - 1))]

        prim_cutoff = self.settings.primability_cutoff
        stab_cutoff = self.settings.stability_cutoff

        for direction in [DNADirection.FWD, DNADirection.REV]:
            logging.debug(f"Repliconf.search(): {direction}")

            seq = self.template_seq[direction]
            search_range = self.range()

            if direction == DNADirection.FWD:
                # FWD: target is slice reversed.
                # target[k] = seq[i + L - 1 - k]
                for i in search_range:
                    prim_num = 0.0
                    stab_num = 0.0
                    run_len = 0
                    run_score = 0.0

                    for k in range(L):
                        base_t = seq[i + L - 1 - k]

                        # Primability
                        prim_num += prim_score_lookup[k][base_t]

                        # Stability
                        val = stab_score_lookup[k][base_t]
                        if val > 0:
                            run_len += 1
                            run_score += val
                        else:
                            if run_len > 0:
                                idx = run_len - 1
                                stab_num += r[idx] * run_score
                                run_len = 0
                                run_score = 0.0

                    # Finish stability run
                    if run_len > 0:
                        idx = run_len - 1
                        stab_num += r[idx] * run_score

                    primability = prim_num / prim_denom if prim_denom != 0 else 0.0
                    stability = stab_num / stab_denom if stab_denom != 0 else 0.0

                    if primability > prim_cutoff and stability > stab_cutoff:
                        origin = self.origin(direction, i)
                        logging.debug(
                            f"Repliconf.search(): adding [{direction}, {i}]: {origin}"
                        )
                        self.origin_db.fwd.append(DirIdx(direction, i))

            else:
                # REV: target is slice (not reversed).
                # target[k] = seq[i + k]
                for i in search_range:
                    prim_num = 0.0
                    stab_num = 0.0
                    run_len = 0
                    run_score = 0.0

                    for k in range(L):
                        base_t = seq[i + k]

                        prim_num += prim_score_lookup[k][base_t]

                        val = stab_score_lookup[k][base_t]
                        if val > 0:
                            run_len += 1
                            run_score += val
                        else:
                            if run_len > 0:
                                idx = run_len - 1
                                stab_num += r[idx] * run_score
                                run_len = 0
                                run_score = 0.0

                    if run_len > 0:
                        idx = run_len - 1
                        stab_num += r[idx] * run_score

                    primability = prim_num / prim_denom if prim_denom != 0 else 0.0
                    stability = stab_num / stab_denom if stab_denom != 0 else 0.0

                    if primability > prim_cutoff and stability > stab_cutoff:
                        origin = self.origin(direction, i)
                        logging.debug(
                            f"Repliconf.search(): adding [{direction}, {i}]: {origin}"
                        )
                        self.origin_db.rev.append(DirIdx(direction, i))

        self.origin_db.searched = True

    @property
    def searched(self) -> bool:
        """Return whether we have searched for the replication origins."""
        return self.origin_db.searched

    def __eq__(self, other: object) -> bool:
        """Check if two Repliconf objects are equal.

        Two Repliconf objects are considered equal if they have the same primer
        and template DNA.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if not isinstance(other, Repliconf):
            return NotImplemented
        return self.primer == other.primer and self.template == other.template

    def __hash__(self) -> int:
        """Return a hash value for the Repliconf object.

        The hash is based on the primer and template DNA.

        Returns:
            int: The hash value.
        """
        return hash((self.primer, self.template))

    def __str__(self) -> str:
        """Return the string representation of the Repliconf object.

        Returns:
            str: A string representation of the object.
        """
        return f"ReplicationConfig: Primer: {self.primer}, Target: {self.template}"
