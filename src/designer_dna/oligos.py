# BSD 3-Clause License
#
# Copyright (c) 2025, Spill-Tea
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Common utility functions to work with and analyze oligonucleotide sequences."""

from ._oligonucleotides import manacher
from ._oligos import (
    complement,
    nrepeats,
    palindrome,
    reverse,
    reverse_complement,
    stretch,
)


__all__ = [
    "complement",
    "complement_py",
    "manacher",
    "nrepeats",
    "nrepeats_py",
    "palindrome",
    "palindrome_py",
    "reverse",
    "reverse_complement",
    "reverse_complement_py",
    "reverse_py",
    "stretch",
    "stretch_py",
]

BASEPAIRS_DNA: dict[str, str] = dict(
    zip(
        "AGTCNURYSWKMBVDH-.",
        "TCAGNAYRSWMKVBHD-.",
        strict=True,
    )
)
BASEPAIRS_DNA.update({k.lower(): v for k, v in BASEPAIRS_DNA.items()})
BASEPAIRS_RNA: dict[str, str] = BASEPAIRS_DNA.copy()
BASEPAIRS_RNA.update({"A": "U", "a": "U"})
DEFAULT_ENCODING: str = "utf_8"


def _make_translation(
    mapping: dict[str, str],
    encoding: str = DEFAULT_ENCODING,
) -> bytes:
    """Construct a string translation table from a dictionary mapping.

    Args:
        mapping (dict): a dictionary of complements.
        encoding (str): Valid and supported encoding schema.

    Returns:
        (bytes): translation table of the provided mapping, prepared as specified by the
        encoding.

    """
    keys: bytes = "".join(mapping.keys()).encode(encoding)
    values: bytes = "".join(mapping.values()).encode(encoding)

    return bytes.maketrans(keys, values)


DNA: bytes = _make_translation(BASEPAIRS_DNA, DEFAULT_ENCODING)
RNA: bytes = _make_translation(BASEPAIRS_RNA, DEFAULT_ENCODING)


def reverse_py(sequence: str) -> str:
    """Reverse a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        (str) Reverse a string.

    Examples:
        .. code-block:: python

            reverse_py("ATATAT") == "TATATA"
            reverse_py("AATATA") == "ATATAA"

    """
    return sequence[::-1]


def complement_py(sequence: str, dna: bool = True) -> str:
    """Return the complement of a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): If true, treat sequence as DNA, otherwise treat as RNA

    Returns:
        (str): Complement of input sequence.

    Examples:
        .. code-block:: python

            complement_py("ATGC", True) == "TACG"
            complement_py("ATGC", False) == "UACG"

    """
    return (
        sequence.encode(DEFAULT_ENCODING)
        .translate(DNA if dna else RNA)
        .decode(DEFAULT_ENCODING)
    )


def reverse_complement_py(sequence: str, dna: bool = True) -> str:
    """Reverse complement a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): sequence is dna, else rna.

    Returns:
        (str) Reverse complement of sequence string.

    Examples:
        .. code-block:: python

            reverse_complement_py("ATGC", True) == "GCAT"
            reverse_complement_py("ATGC", False) == "GCAU"

    """
    return complement_py(reverse_py(sequence), dna)


def _center_expansion(
    s: str,
    c: str,
    left: int,
    right: int,
    length: int,
) -> tuple[int, int]:
    while left > -1 and right < length and s[left] == c[right] and s[right] == c[left]:
        left -= 1
        right += 1

    return left + 1, right


def palindrome_py(sequence: str, dna: bool = True) -> str:
    """Find the longest substring palindrome within a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): If true, treat sequence as DNA, otherwise treat as RNA

    Returns:
        (str): longest palindromic subsequence within sequence.

    Examples:
        .. code-block:: python

            palindrome_py("ATAT") == "ATAT"
            palindrome_py("GATATG") == "ATAT"

    Notes:
        * Algorithmic time complexity is O(N).
        * If a sequence contains two or more palindromic substrings of equal size, the
          first leftmost palindrome is prioritized.

    """
    seq_length: int = len(sequence)
    comp: str = complement_py(sequence, dna)
    best_left: int = 0
    best_right: int = 0
    left: int = 0
    right: int = 0
    length: int = 0

    for i in range(seq_length - 1):
        # If we only consider ATGC based sequences, then Palindromic nucleotides are
        # only even length, reducing search space in half
        left, right = _center_expansion(sequence, comp, i, i + 1, seq_length)
        current: int = right - left
        if current > length:
            length = current
            best_left = left
            best_right = right

        # Handle Degenerate bases which equals it's complement (e.g. N, S, W)
        # We can have odd length palindromes if a degenerate base is at the center
        if sequence[i] != comp[i]:
            continue

        left, right = _center_expansion(sequence, comp, i - 1, i + 1, seq_length)
        current = right - left
        if current > length:
            length = current
            best_left = left
            best_right = right

    return sequence[best_left:best_right]


def stretch_py(sequence: str) -> int:
    """Calculate the maximum stretch of a single character in a string.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        (int): maximum length observed within sequence of a repeated character.

    Examples:
        .. code-block:: python

            stretch_py("AAAA") == 3
            stretch_py("AATT") == 1

    """
    if not sequence:
        return 0

    longest: int = 0
    current: int = 0
    last: str = sequence[0]
    char: str

    for char in sequence[1:]:
        if char == last:
            current += 1
            if current > longest:
                longest = current
        else:  # reset
            current = 0
            last = char

    return longest


def nrepeats_py(sequence: str, n: int) -> int:
    """Calculate the longest substring of n repeating characters.

    Args:
        sequence (str): Nucleotide string or Series of string
        n (int): stretch of k-mer to observe

    Returns:
        (int) The longest run of repeating n-length characters.

    Raises:
        ValueError: when n < 1

    Examples:
        .. code-block:: python

            nrepeats_py("AAAA", 1) == 3  #  True
            nrepeats_py("AAAA", 2) == 1  #  True
            nrepeats_py("ACAACAACA", 3) == 2  #  True

    """
    max_val: int = 0
    length: int = len(sequence)

    for k in range(n):
        previous: str = sequence[k : n + k]
        current: int = 0
        for j in range(n, length, n):
            phase: str = sequence[j + k : j + k + n]
            if phase == previous:
                current += 1
                if current > max_val:
                    max_val = current
            else:
                current = 0
                previous = phase

    return max_val
