# BSD 3-Clause License

# Copyright (c) 2025, Spill-Tea

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

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

# cython: boundscheck=False, wraparound=False, nonecheck=False
"""Cythonized oligonucleotide functions."""

from libc.stdlib cimport free

cdef extern from "Python.h":
    Py_ssize_t PyUnicode_GET_LENGTH(object)
    bytes PyUnicode_AsUTF8String(object)
    Py_ssize_t PyBytes_GET_SIZE(object)

cimport common

cdef extern from "oligos.h":
    const unsigned char DNA[0x100]
    const unsigned char RNA[0x100]


cdef void c_reverse(char* seq, Py_ssize_t length) noexcept:
    """Reverse a C string in place.

    Args:
        seq (char*): buffer sequence.
        length (Py_ssize_t): length of seq.

    """
    cdef Py_ssize_t start, end, x = length // 2

    for start in range(x):
        end = length - start - 1
        seq[start], seq[end] = seq[end], seq[start]


cpdef str reverse(str sequence):
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
    return sequence[::-1]


cdef void c_complement(char* sequence, Py_ssize_t length, unsigned char[] table):
    """Complement sequence C string in place.

    Args:
        seq (char*): buffer sequence.
        length (Py_ssize_t): length of seq.
        table (char[]): translation table.

    """
    cdef:
        Py_ssize_t j, end, idx = length // 2

    for j in range(idx):
        end = (length - 1) - j
        sequence[j] = table[<unsigned char> sequence[j]]
        sequence[end] = table[<unsigned char> sequence[end]]

    if length % 2:
        sequence[idx] = table[<unsigned char> sequence[idx]]


cpdef str complement(str sequence, bint dna = True):
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
    cdef common.StringView view = common.str_to_view(sequence)

    if dna:
        c_complement(view.ptr, view.size, DNA)
    else:
        c_complement(view.ptr, view.size, RNA)

    return common.to_str(view)


cdef void c_reverse_complement(
    char* sequence,
    Py_ssize_t length,
    unsigned char[] table
):
    """Reverse complement sequence C string in place.

    Args:
        sequence (char*): buffer pointer to nucleotide char sequence.
        length (Py_ssize_t): length of seq.
        table (char[]): translation table.

    """
    cdef:
        char* end_ptr = sequence + (length - 1)

    while end_ptr > sequence:
        sequence[0], end_ptr[0] = (
            table[<unsigned char> end_ptr[0]],
            table[<unsigned char> sequence[0]]
        )
        end_ptr -= 1
        sequence += 1

    if length % 2:
        sequence[0] = table[<unsigned char> sequence[0]]


cpdef str reverse_complement(str sequence, bint dna = True):
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
    cdef common.StringView view = common.str_to_view(sequence)

    if dna:
        c_reverse_complement(view.ptr, view.size, DNA)
    else:
        c_reverse_complement(view.ptr, view.size, RNA)

    return common.to_str(view)


cdef bytes _expand_from_center(
    bytes seq,
    bytes comp,
    Py_ssize_t left,
    Py_ssize_t right,
    Py_ssize_t length,
):
    while (
        left > -1
        and right < length
        and seq[left] == comp[right]
        and seq[right] == comp[left]  # required to detect dna to rna based complements
    ):
        left -= 1
        right += 1

    return seq[left + 1 : right]


cpdef str palindrome(str sequence, bint dna = True):
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
    cdef:
        bytes temp = PyUnicode_AsUTF8String(sequence)
        bytes comp = PyUnicode_AsUTF8String(complement(sequence, dna))
        bytes even, pal = b""
        Py_ssize_t seq_length = PyUnicode_GET_LENGTH(sequence)
        Py_ssize_t i, current, length = 0

    if seq_length < 2:  # noqa: PLR2004
        return ""

    for i in range(seq_length - 1):
        # NOTE: Palindromic nucleotides are only even length, halving search space
        even = _expand_from_center(temp, comp, i, i + 1, seq_length)
        current = PyBytes_GET_SIZE(even)
        if current > length:
            pal = even
            length = current

    return pal.decode("utf8")


cpdef int stretch(str sequence):
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
    cdef:
        common.StringView view = common.str_to_view(sequence)
        Py_ssize_t j
        int longest = 0, current = 0
        char prev = view.ptr[0]

    for j in range(1, view.size):
        if view.ptr[j] == prev:
            current += 1
            if current > longest:
                longest = current
        else:
            current = 0
            prev = view.ptr[j]
    free(view.ptr)

    return longest


cpdef int nrepeats(str sequence, int n):
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
    cdef:
        common.StringView view = common.str_to_view(sequence)
        Py_ssize_t t = <Py_ssize_t> n
        Py_ssize_t v = view.size // t
        Py_ssize_t i, j, k
        int max_val = 0
        list[char] previous = [view.ptr[i : t + i] for i in range(t)]
        list[int] current = [0 for i in range(t)]
        bytes phase

    for j in range(1, v):
        for k in range(t):
            phase = view.ptr[j * t + k : j * t + k + t]
            if phase == previous[k]:
                current[k] += 1
                if current[k] > max_val:
                    max_val = current[k]
            else:
                current[k] = 0
                previous[k] = phase
    free(view.ptr)

    return max_val
