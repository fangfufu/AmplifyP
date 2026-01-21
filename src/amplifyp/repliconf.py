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
    """A class representing a directed index.

    This class encapsulates a specific position (index) on a DNA strand,
    along with the direction (forward or reverse) of that strand.

    Attributes:
        direction (DNADirection): The direction of the DNA strand.
        index (int): The integer index position.
    """

    direction: DNADirection
    index: int

    def __int__(self) -> int:
        """Return the index value as an integer.

        Returns:
            int: The index.
        """
        return self.index

    def __index__(self) -> int:
        """Return the index value for slicing/indexing usage.

        Returns:
            int: The index.
        """
        return self.index

    def __add__(self, other: object) -> DirIdx:
        """Add an integer or another DirIdx to this index.

        Args:
            other (int | DirIdx): The value to add.

        Returns:
            DirIdx: A new DirIdx with the updated index value. The direction remains
                unchanged.
        """
        if isinstance(other, int):
            return DirIdx(self.direction, self.index + other)
        if not isinstance(other, DirIdx):
            return NotImplemented
        return DirIdx(self.direction, self.index + other.index)

    def __sub__(self, other: object) -> DirIdx:
        """Subtract an integer or another DirIdx from this index.

        Args:
            other (int | DirIdx): The value to subtract.

        Returns:
            DirIdx: A new DirIdx with the updated index value. The direction remains
                unchanged.
        """
        if isinstance(other, int):
            return DirIdx(self.direction, self.index - other)
        if not isinstance(other, DirIdx):
            return NotImplemented
        return DirIdx(self.direction, self.index - other.index)

    def __eq__(self, other: object) -> bool:
        """Check equality with an integer or another DirIdx.

        Args:
            other (object): The object to compare.

        Returns:
            bool: True if equal, False otherwise. If comparing with int, checks
                index only.
        """
        if isinstance(other, int):
            return self.index == other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.direction == other.direction and self.index == other.index

    def __lt__(self, other: object) -> bool:
        """Check if index is less than another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if less than.
        """
        if isinstance(other, int):
            return self.index < other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index < other.index

    def __gt__(self, other: object) -> bool:
        """Check if index is greater than another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if greater than.
        """
        if isinstance(other, int):
            return self.index > other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index > other.index

    def __le__(self, other: object) -> bool:
        """Check if index is less than or equal to another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if less than or equal.
        """
        if isinstance(other, int):
            return self.index <= other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index <= other.index

    def __ge__(self, other: object) -> bool:
        """Check if index is greater than or equal to another.

        Args:
            other (object): Value to compare.

        Returns:
            bool: True if greater than or equal.
        """
        if isinstance(other, int):
            return self.index >= other
        if not isinstance(other, DirIdx):
            return NotImplemented
        return self.index >= other.index

    def __str__(self) -> str:
        """Return the string representation of the index.

        Returns:
            str: The index value as a string.
        """
        return f"{self.index}"


@dataclass(slots=True)
class DirIdxDb:
    """A database for storing valid replication origin locations.

    This container holds lists of indices where valid replication origins were found,
    segregated by strand direction (forward and reverse). It also tracks whether
    a search operation has been performed.

    Attributes:
        fwd (list[DirIdx]): A list of locations (DirIdx) for origins on the forward
            strand.
        rev (list[DirIdx]): A list of locations (DirIdx) for origins on the reverse
            strand.
        searched (bool): Indicates if the search has been executed. Defaults to False.
    """

    fwd: list[DirIdx]
    rev: list[DirIdx]
    searched: bool = False

    def clear(self) -> None:
        """Clear the database, removing stored indices and resetting search flag."""
        self.fwd.clear()
        self.rev.clear()
        self.searched = False

    def __getitem__(self, key: tuple[DNADirection, int]) -> DirIdx:
        """Retrieve a stored DirIdx by direction and list index.

        Args:
            key (tuple[DNADirection, int]): A tuple specifying the direction (FWD/REV)
                and the integer index within that list.

        Returns:
            DirIdx: The requested directed index.
        """
        direction, i = key
        if direction == DNADirection.FWD:
            return self.fwd[i]
        else:
            return self.rev[i]


class Repliconf:
    """A class representing a Replication Configuration.

    A 'Repliconf' encapsulates the setup required to search for replication origins
    given a specific Primer and Template DNA. It prepares the template sequences
    (including padding and complementary strands) and manages the search process.

    Attributes:
        padding_len (int): The length of padding added to the template, typically
            equal to the primer length.
        primer (Primer): The primer sequence being analyzed.
        template (DNA): The original template DNA sequence.
        template_seq (dict[DNADirection, str]): A dictionary mapping direction
            (FWD/REV) to the processed (padded/complemented) template strings.
        settings (Settings): The configuration settings (weights, cutoffs) used
            for scoring.
        origin_db (DirIdxDb): The results database storing found origin indices.
    """

    def __init__(
        self,
        template: DNA,
        primer: Primer,
        settings: Settings,
    ) -> None:
        """Initialize a Repliconf object.

        Args:
            template (DNA): The template DNA sequence to be searched.
            primer (Primer): The primer sequence to search for.
            settings (Settings): The settings used for evaluating origins.
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
        """Return the range of valid starting indices for search.

        Returns:
            range: A range object covering all valid start positions in the padded
                template.
        """
        return range(len(self.template_seq[DNADirection.FWD]) - len(self.primer) + 1)

    def slice(self, i: int) -> slice:
        """Create a slice object for extracting a target segment.

        Args:
            i (int): The start index.

        Returns:
            slice: A slice of length equal to the primer length.
        """
        return slice(i, i + len(self.primer))

    def origin(self, var: DirIdx | DNADirection, *idx: int) -> ReplicationOrigin:
        """Construct a ReplicationOrigin object for a specific location.

        This method extracts the target sequence from the pre-processed template
        at the specified location and direction, and creates a ReplicationOrigin
        object to evaluate it.

        Args:
            var (DirIdx | DNADirection): Either a `DirIdx` object containing both
                direction and index, or a `DNADirection` enum.
            *idx (int): If `var` is a `DNADirection`, this argument must be provided
                as the integer index.

        Returns:
            ReplicationOrigin: An evaluated replication origin object.

        Raises:
            TypeError: If arguments are not of the expected types.
        """
        if isinstance(var, DirIdx):
            direction = var.direction
            i = var.index
        elif isinstance(var, DNADirection):
            direction = var
            i = idx[0]
        else:
            raise TypeError("var must be DirIdx or DNADirection")

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
        """Retrieve a ReplicationOrigin based on an index stored in the database.

        Args:
            direction (DNADirection): The direction to look up in the database.
            i (int): The index *within the list of found origins* (not the genomic
                index).

        Returns:
            ReplicationOrigin: The replication origin object.
        """
        return self.origin(self.origin_db[direction, i])

    def search(self) -> None:
        """Perform a search for valid replication origins.

        Scans the entire template in both forward and reverse directions.
        Locations where the primability and stability scores exceed the configured
        cutoffs are stored in `origin_db`.
        """
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
        return self.primer == other.primer and self.template == other.template

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
        return f"ReplicationConfig: Primer: {self.primer}, Target: {self.template}"
