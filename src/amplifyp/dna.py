"""DNA-related classes and enumerations for AmplifyP."""

from enum import Flag, IntEnum, StrEnum


class Nucleotides(StrEnum):
    """Enumeration of valid nucleotide characters for use in AmplifyP.

    This class defines groups of nucleotide characters that are valid in
    different contexts within the AmplifyP application. It includes single
    nucleotides, ambiguous (double and triple) nucleotides, wildcards, and gaps.
    It also defines composite groups for circular, linear, and primer DNA types.
    """

    SINGLE = "GATC"
    DOUBLE = "MRWSYK"
    TRIPLE = "VHDB"
    WILDCARD = "N"
    GAP = "-"

    CIRCULAR = SINGLE + WILDCARD
    LINEAR = CIRCULAR + GAP
    PRIMER = SINGLE + DOUBLE + TRIPLE + WILDCARD


class DNAType(IntEnum):
    """An enumeration representing the type of DNA.

    This enumeration is used to specify whether a DNA sequence is linear,
    circular, or a primer. This distinction is important for certain
    operations, such as padding and rotation.
    """

    LINEAR = 1
    CIRCULAR = 2
    PRIMER = 3


class DNADirection(Flag):
    """An enumeration representing the direction of a DNA strand.

    This enumeration defines the two possible directions for a DNA strand:
    forward (5' to 3') and reverse (3' to 5'). It is used to track the
    orientation of DNA sequences.
    """

    FWD = True
    REV = False


class DNA:
    """A class representing a DNA sequence.

    This class encapsulates a DNA sequence and its properties, such as its type
    (linear, circular, or primer), name, and direction. It provides methods for
    manipulating and analyzing the DNA sequence, such as complementing,
    reversing, and padding.

    Attributes:
        seq (str): The DNA sequence.
        type (DNAType): The type of the DNA sequence.
        name (str): The name of the DNA sequence.
        direction (DNADirection): The direction of the DNA sequence.
    """

    def __init__(
        self,
        seq: str,
        dna_type: DNAType = DNAType.LINEAR,
        name: str | None = None,
        direction: bool | DNADirection = DNADirection.FWD,
    ) -> None:
        """Initializes a DNA object.

        Args:
            seq (str): The DNA sequence.
            dna_type (DNAType, optional): The type of DNA. Defaults to
              DNAType.LINEAR.
            name (str, optional): The name of the DNA sequence. Defaults to None.
            direction (bool | DNADirection, optional): The direction of the DNA
              sequence. Defaults to DNADirection.FWD.

        Raises:
            TypeError: If the DNA type is invalid.
            ValueError: If the DNA sequence contains invalid characters.
        """
        self.__seq: str = "".join(seq.split())
        self.__type: DNAType = dna_type

        if name is None:
            self.__name = seq
        else:
            self.__name = name.strip()
        self.__direction: bool | DNADirection = direction

        if dna_type == DNAType.LINEAR:
            check_str = Nucleotides.LINEAR
        elif dna_type == DNAType.CIRCULAR:
            check_str = Nucleotides.CIRCULAR
        elif dna_type == DNAType.PRIMER:
            check_str = Nucleotides.PRIMER
        else:
            raise TypeError("Invalid DNA type.")

        if not set(self.__seq.upper()) <= set(check_str):
            raise ValueError("The DNA sequence contains invalid characters.")

    @property
    def seq(self) -> str:
        """The DNA sequence as a string."""
        return self.__seq

    @property
    def type(self) -> DNAType:
        """The DNA type."""
        return self.__type

    @property
    def name(self) -> str:
        """The name of the DNA sequence."""
        return self.__name

    @name.setter
    def name(self, value: str) -> None:
        """Set the name of the DNA sequence."""
        self.__name = value.strip()

    @property
    def direction(self) -> bool | DNADirection:
        """The direction of the DNA sequence."""
        return self.__direction

    def lower(self) -> "DNA":
        """Return a new DNA object with the sequence in lower case.

        Returns:
            DNA: A new DNA object with the sequence in lower case.
        """
        return DNA(self.seq.lower(), self.type, self.name, self.direction)

    def upper(self) -> "DNA":
        """Return a new DNA object with the sequence in upper case.

        Returns:
            DNA: A new DNA object with the sequence in upper case.
        """
        return DNA(self.seq.upper(), self.type, self.name, self.direction)

    def complement(self) -> "DNA":
        """Return the complement of the DNA sequence.

        The complement is created by swapping A with T, C with G, and vice versa.
        The direction of the new DNA object is also inverted.

        Returns:
            DNA: A new DNA object representing the complement of the original
                 sequence.
        """
        return DNA(
            self.seq.translate(
                str.maketrans("ACGTMKRYBDHVacgtmkrybdhv", "TGCAKMYRVHDBtgcakmyrvhdb")
            ),
            self.type,
            self.name,
            not self.direction,
        )

    def reverse(self) -> "DNA":
        """Return the reverse of the DNA sequence.

        The direction of the new DNA object is also inverted.

        Returns:
            DNA: A new DNA object representing the reverse of the original
                 sequence.
        """
        return DNA(self.seq[::-1], self.type, self.name, not self.direction)

    def reverse_complement(self) -> "DNA":
        """Return the reverse complement of the DNA sequence.

        The direction of the new DNA object is also inverted.

        Returns:
            DNA: A new DNA object representing the reverse complement of the
                 original sequence.
        """
        return self.reverse().complement()

    def __eq__(self, other: object) -> bool:
        """Check if two DNA objects are identical.

        DNA objects are considered identical if they have the same sequence,
        direction, and type. The comparison is case-insensitive.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if the DNA sequences are identical, False otherwise.
        """
        if not isinstance(other, DNA):
            return NotImplemented
        return (
            self.seq.upper() == other.seq.upper()
            and self.direction == other.direction
            and self.type == other.type
        )

    def __hash__(self) -> int:
        """Return a hash value for the DNA object.

        The hash value is calculated based on the DNA sequence (in upper case),
        direction, and type.

        Returns:
            int: The hash value for the DNA object.
        """
        return hash((self.seq.upper(), self.direction, self.type))

    def is_complement_of(self, other: "DNA") -> bool:
        """Check if another DNA object is the complement of this one.

        This method checks if the sequence of the other DNA object is the
        complement of this object's sequence and if their directions are
        opposite.

        Args:
            other (DNA): The other DNA object to compare with.

        Returns:
            bool: True if the other DNA object is the complement of this one,
                  False otherwise.
        """
        return (
            self.seq.upper() == other.complement().seq.upper()
            and self.direction != other.direction
        )

    def __len__(self) -> int:
        """Return the length of the DNA sequence.

        Returns:
            int: The length of the DNA sequence.
        """
        return len(self.seq)

    def pad(self, i: int) -> "DNA":
        """Add padding to the beginning of the DNA sequence string.

        For linear DNA, the padding consists of gap characters ('-'). For
        circular DNA, the padding is taken from the end of the sequence,
        simulating a circular molecule.

        Args:
            i (int): The amount of padding to add.

        Returns:
            DNA: A new DNA object with the padded sequence.

        Raises:
            TypeError: If padding is attempted on a DNA type that does not
                       support it.
        """
        if self.type == DNAType.LINEAR:
            padding_str = Nucleotides.GAP * i
        elif self.type == DNAType.CIRCULAR:
            padding_str = self.seq[-i::]
        else:
            raise TypeError("Invalid DNA type for padding operation.")

        new_str = padding_str + self.seq

        return DNA(
            new_str,
            self.type,
            self.name,
            self.direction,
        )

    def rot(self, i: int) -> "DNA":
        """Rotate the DNA sequence by a specified number of bases.

        This operation is only supported for circular DNA.

        Args:
            i (int): The number of bases to rotate the sequence by.

        Returns:
            DNA: A new DNA object with the rotated sequence.

        Raises:
            TypeError: If rotation is attempted on non-circular DNA.
        """
        if self.type == DNAType.LINEAR:
            raise TypeError("Rotation is unsupported for linear DNA.")
        padded_dna = self.pad(i)

        return padded_dna[0 : len(self.seq)]

    def __getitem__(self, key: slice) -> "DNA":
        """Return a slice of the DNA sequence.

        This method allows for slicing of the DNA object, returning a new DNA
        object representing the specified slice.

        Args:
            key (slice): The slice to be taken from the DNA sequence.

        Returns:
            DNA: A new DNA object representing the slice.
        """
        return DNA(self.seq[key], self.type, self.name, self.direction)

    def __str__(self) -> str:
        """Return a string representation of the DNA object.

        The string includes the DNA's name, type, and direction.

        Returns:
            str: A string representation of the DNA object.
        """
        return (
            f"DNA: {self.name}, {self.type.name}, {DNADirection(self.direction).name}"
        )

    def __add__(self, other: "DNA") -> "DNA":
        """Add two DNA sequences.

        Args:
            other (DNA): The other DNA sequence to add.

        Returns:
            DNA: A new DNA object representing the sum of the two sequences.
        """
        return DNA(self.seq + other.seq, DNAType.LINEAR, self.name, self.direction)


class Primer(DNA):
    """A class representing a primer sequence.

    A primer is a short single-stranded DNA sequence used as a starting point
    for DNA synthesis. This class is a specialized subclass of `DNA` with the
    type set to `DNAType.PRIMER` and the direction set to `DNADirection.FWD`.

    Attributes:
        seq (str): The primer sequence.
        name (str): The name of the primer.
    """

    def __init__(
        self,
        sequence: str,
        name: str | None = None,
    ) -> None:
        """Initializes a Primer object.

        Args:
            sequence (str): The primer sequence.
            name (str, optional): The name of the primer. Defaults to None.
        """
        super().__init__(sequence, DNAType.PRIMER, name, DNADirection.FWD)
