# pylint: disable=W0613

"""Cythonized oligonucleotide functions."""

def reverse(sequence: str) -> str:
    """Reverse a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        (str) Reverse a string.

    Examples:
        .. code-block:: python

            reverse("ATATAT") == "TATATA"
            reverse("AATATA") == "ATATAA"

    """

def complement(sequence: str, dna: bool = ...) -> str:
    """Complement a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (str) Complement of a nucleotide sequence string.

    Examples:
        .. code-block:: python

            complement("ATGC", True) == "TACG"
            complement("ATGC", False) == "UACG"

    """

def reverse_complement(sequence: str, dna: bool = ...) -> str:
    """Reverse complement a nucleotide sequence.

    Args:
        sequence (str): Nucelotide sequence string.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (str) Reverse complement of sequence string.

    Examples:
        .. code-block:: python

            reverse_complement("ATGC", True) == "GCAT"
            reverse_complement("ATGC", False) == "GCAU"

    """

def palindrome(sequence: str, dna: bool = ...) -> str:
    """Find the longest palindromic substring within a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (str): longest palindromic subsequence within sequence.

    Examples:
        .. code-block:: python

            palindrome("ATAT") == "ATAT"
            palindrome("GATATG") == "ATAT"

    Notes:
        * Uses a modified center expansion method (Manacher's algorithm) to identify the
        longest substring that is palindromic.
        * If a sequence contains two or more palindromic substrings of equal size, the
        first leftmost palindrome is prioritized.

    """

def stretch(sequence: str) -> int:
    """Return the maximum length of a single letter (nucleotide) repeat in a string.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        (int): Length of maximum run of a single letter.

    Examples:
        .. code-block:: python

            stretch("ATATAT") == 0  # True
            stretch("AATATA") == 1  # True

    """

def nrepeats(sequence: str, n: int) -> int:
    """Calculate the maximum observed repeats of composite pattern size n characters.

    Args:
        sequence (str): Nucleotide sequence string.
        n (int): stretch of k-mers to observe.

    Returns:
        (int) The longest tandem run of nucleotides comprised of a composite pattern
        of length n characters.

    Raises:
        ValueError: if value of n is less than 1.

    Examples:
        .. code-block:: python

            nrepeats("AAAA", 1) == 3  #  True
            nrepeats("AAAA", 2) == 1  #  True
            nrepeats("ACAACAACA", 3) == 2  #  True

    """
