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

from array import array

def m_reverse(sequence: array[int]) -> None:
    """Reverse a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.

    Returns:
        (void) Reverse a sequence in place.

    """

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

def m_complement(sequence: array[int], dna: bool = ...) -> None:
    """Complement a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (void) Complement nucleotide sequence in place.

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

def m_reverse_complement(sequence: array[int], dna: bool = ...) -> None:
    """Reverse complement a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (void) Reverse complement nucleotide sequence in place.

    """

def reverse_complement(sequence: str, dna: bool = ...) -> str:
    """Reverse complement a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (str) Reverse complement of sequence string.

    Examples:
        .. code-block:: python

            reverse_complement("ATGC", True) == "GCAT"
            reverse_complement("ATGC", False) == "GCAU"

    """

def m_palindrome(sequence: array[int], dna: bool = ...) -> tuple[int, int]:
    """Find the longest palindromic substring within a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (int, int) start and end indices denoting the longest found
        palindromic subsequence within sequence.

    """

def palindrome(sequence: str, dna: bool = ...) -> str:
    """Find the longest palindromic substring within a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (str) longest palindromic subsequence within sequence.

    Examples:
        .. code-block:: python

            palindrome("ATAT") == "ATAT"
            palindrome("GATATG") == "ATAT"
            palindrome("ANT") == "ANT" # Handles degenerate bases

    Notes:
        * If a sequence contains two or more palindromic substrings of equal size, the
          first leftmost palindrome is prioritized.

    """

def m_stretch(sequence: array[int]) -> int:
    """Return the maximum length of a single letter (nucleotide) repeat in a string.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.

    Returns:
        (int) Length of maximum run of a single letter.

    """

def stretch(sequence: str) -> int:
    """Return the maximum length of a single letter (nucleotide) repeat in a string.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        (int) Length of maximum run of a single letter.

    Examples:
        .. code-block:: python

            stretch("ATATAT") == 0  # True
            stretch("AATATA") == 1  # True

    """

def m_nrepeats(sequence: array[int], n: int) -> int:
    """Calculate the maximum observed repeats of composite pattern size n characters.

    Args:
        sequence (uchar[]): Nucleotide sequence string.
        n (int): Size of k-mers (composite pattern) to observe.

    Returns:
        (int) The longest tandem run of nucleotides comprised of a composite pattern
        of length n characters.

    Raises:
        ZeroDivisionError: if value of n is 0.

    """

def nrepeats(sequence: str, n: int) -> int:
    """Calculate the maximum observed repeats of composite pattern size n characters.

    Args:
        sequence (str): Nucleotide sequence string.
        n (int): Size of k-mers (composite pattern) to observe.

    Returns:
        (int) The longest tandem run of nucleotides comprised of a composite pattern
        of length n characters.

    Raises:
        ZeroDivisionError: if value of n is 0.

    Examples:
        .. code-block:: python

            nrepeats("AAAA", 1) == 3  #  True
            nrepeats("AAAA", 2) == 1  #  True
            nrepeats("ACAACAACA", 3) == 2  #  True

    """
