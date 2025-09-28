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

# cython: boundscheck=False, wraparound=False, nonecheck=False
"""Cythonized oligonucleotide functions."""

from libc.string cimport memcpy
from libc.stdlib cimport free, malloc

from common cimport (
    StringView,
    str_to_view,
    to_str
)

cdef extern from "oligos.h":
    const unsigned char DNA[0x100]
    const unsigned char RNA[0x100]


cdef inline void c_reverse(unsigned char* sequence, Py_ssize_t length) noexcept:
    """Reverse a C string in place.

    Args:
        sequence (uchar*): Buffer sequence.
        length (Py_ssize_t): Length of sequence.

    Returns:
        (void) Reverse sequence in place.

    """
    cdef Py_ssize_t start, end, x = length // 2

    for start in range(x):
        end = length - start - 1
        sequence[start], sequence[end] = sequence[end], sequence[start]


cdef inline void v_reverse(StringView* view) noexcept:
    """Handle reverse in place on StringView directly.

    Args:
        view (StringView*): Nucleotide sequence view.

    Returns:
        (void) Reverse char in place.

    """
    c_reverse(view[0].ptr, view[0].size)


cpdef void m_reverse(unsigned char[:] sequence):
    """Reverse a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.

    Returns:
        (void) Reverse a sequence in place.

    """
    cdef Py_ssize_t length = sequence.shape[0]

    c_reverse(&sequence[0], length)


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


cdef inline void _c_complement(
    unsigned char* sequence,
    Py_ssize_t length,
    unsigned char* table
) noexcept:
    """Complement sequence C string in place.

    Args:
        sequence (uchar*):  Nucleotide sequence pointer.
        length (Py_ssize_t): Length of seq.
        table (uchar*): translation table.

    Returns:
        (void) Complement seq in place.

    """
    cdef:
        Py_ssize_t j, end, idx = length // 2

    for j in range(idx):
        end = (length - 1) - j
        sequence[j] = table[<ssize_t> sequence[j]]
        sequence[end] = table[<ssize_t> sequence[end]]

    if length % 2:
        sequence[idx] = table[<ssize_t> sequence[idx]]


cdef void c_complement(
    unsigned char* sequence,
    Py_ssize_t length,
    bint dna
) noexcept:
    """Complement sequence C string in place.

    Args:
        sequence (uchar*):  Nucleotide sequence pointer.
        length (Py_ssize_t): Length of seq.
        dna (bint): DNA else RNA.

    Returns:
        (void) Complement sequence in place.

    """
    if dna:
        _c_complement(sequence, length, &DNA[0])
    else:
        _c_complement(sequence, length, &RNA[0])


cdef inline void v_complement(StringView* view, bint dna) noexcept:
    """Handle complement in place on StringView directly.

    Args:
        view (StringView*): Nucleotide sequence view.

    Returns:
        (void) Complement char in place.

    """
    c_complement(view[0].ptr, view[0].size, dna)


cpdef void m_complement(unsigned char[:] sequence, bint dna = True):
    """Complement a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (void) Complement nucleotide sequence in place.

    """
    cdef:
        Py_ssize_t length = sequence.shape[0]

    c_complement(&sequence[0], length, dna)


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
    cdef StringView view = str_to_view(sequence)
    v_complement(&view, dna)

    return to_str(view)


cdef inline void _c_reverse_complement(
    unsigned char* sequence,
    Py_ssize_t length,
    unsigned char* table
) noexcept:
    """Reverse complement sequence C string in place.

    Args:
        sequence (uchar*): Nucleotide sequence pointer.
        length (Py_ssize_t): Length of sequence.
        table (uchar*): Translation table.

    Returns:
        (void) Reverse complement sequence in place.

    """
    cdef:
        unsigned char* end_ptr = sequence + (length - 1)

    while end_ptr > sequence:
        sequence[0], end_ptr[0] = (
            table[<ssize_t> end_ptr[0]],
            table[<ssize_t> sequence[0]]
        )
        end_ptr -= 1
        sequence += 1

    if length % 2:
        sequence[0] = table[<ssize_t> sequence[0]]


cdef void c_reverse_complement(
    unsigned char* sequence,
    Py_ssize_t length,
    bint dna
) noexcept:
    """Reverse complement sequence C string in place.

    Args:
        sequence (uchar*): Nucleotide sequence pointer.
        length (Py_ssize_t): Length of sequence.
        dna (bint): Sequence is DNA, else RNA.

    Returns:
        (void) Complement sequence in place.

    """
    if dna:
        _c_reverse_complement(sequence, length, &DNA[0])
    else:
        _c_reverse_complement(sequence, length, &RNA[0])


cdef inline void v_reverse_complement(StringView* view, bint dna) noexcept:
    """Handle reverse complement in place on StringView directly.

    Args:
        view (StringView*): Nucleotide sequence view.

    Returns:
        (void) Reverse complement char in place.

    """
    c_reverse_complement(view[0].ptr, view[0].size, dna)


cpdef void m_reverse_complement(unsigned char[:] sequence, bint dna = True):
    """Reverse complement a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (void) Reverse complement nucleotide sequence in place.

    """
    cdef:
        Py_ssize_t length = sequence.shape[0]

    c_reverse_complement(&sequence[0], length, dna)


cpdef str reverse_complement(str sequence, bint dna = True):
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
    cdef StringView view = str_to_view(sequence)
    v_reverse_complement(&view, dna)

    return to_str(view)


cdef inline void _center(
    unsigned char* seq,
    unsigned char* comp,
    Py_ssize_t* left,
    Py_ssize_t* right,
    Py_ssize_t length,
) noexcept:
    while (left[0] > -1 and right[0] < length):
        if seq[left[0]] != comp[right[0]] or seq[right[0]] != comp[left[0]]:
            break
        left[0] -= 1
        right[0] += 1
    left[0] += 1


cdef inline void _update_bounds(
    Py_ssize_t left,
    Py_ssize_t right,
    Py_ssize_t* current,
    Py_ssize_t* length,
    Py_ssize_t* start,
    Py_ssize_t* end
) noexcept:
    current[0] = right - left
    if current[0] > length[0]:
        length[0] = current[0]
        start[0] = left
        end[0] = right


cdef inline (Py_ssize_t, Py_ssize_t) c_palindrome(
    unsigned char* seq,
    Py_ssize_t size,
    bint dna
) noexcept:
    """Find the longest palindromic substring within a nucleotide sequence.

    Args:
        seq (uchar*): Nucleotide sequence pointer.
        size (Py_ssize_t): Length of seq.
        dna (bint): Sequence is DNA, else RNA.

    Returns:
        (tuple[Py_ssize_t, Py_ssize_t]) longest palindromic subsequence within sequence,
        described by start (inclusive) and end (exclusive) index.

    """
    cdef:
        unsigned char* com = <unsigned char*> malloc((size + 1) * sizeof(unsigned char))
        Py_ssize_t i, left, right, current, length = 0, start = 0, end = 0

    memcpy(com, seq, size)
    com[size] = "\0"
    c_complement(com, size, dna)

    for i in range(size - 1):
        # Check even length palindromes first (more common for ATGC based sequences)
        left = i
        right = i + 1
        _center(seq, com, &left, &right, size)
        _update_bounds(left, right, &current, &length, &start, &end)

        # Only check odd length palindromes in case of (center) degenerate bases
        if seq[i] != com[i]:
            continue

        left = i - 1
        right = i + 1
        _center(seq, com, &left, &right, size)
        _update_bounds(left, right, &current, &length, &start, &end)

    free(com)

    return start, end


cpdef (int, int) m_palindrome(
    unsigned char[:] sequence,
    bint dna = True
):
    """Find the longest palindromic substring within a nucleotide sequence.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (int, int) start and end indices denoting the longest found
        palindromic subsequence within sequence.

    """
    cdef:
        int start, end
        Py_ssize_t length = sequence.shape[0]

    start, end = c_palindrome(&sequence[0], length, dna)

    return start, end


cpdef str palindrome(str sequence, bint dna = True):
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
    cdef:
        StringView seq = str_to_view(sequence)
        Py_ssize_t start, end

    start, end = c_palindrome(seq.ptr, seq.size, dna)
    free(seq.ptr)

    return sequence[start: end]


cdef inline int c_stretch(
    unsigned char* sequence,
    Py_ssize_t length,
) noexcept:
    """Compute nucleotide stretch of a single character repeat.

    Args:
        sequence (uchar*): buffer pointer to nucleotide char sequence.
        length (Py_ssize_t): length of sequence.

    Returns:
        (int) length of longest single character repeat subsequence.

    """
    cdef:
        unsigned char prev = sequence[0]
        Py_ssize_t j
        int current = 0, longest = 0

    for j in range(1, length):
        if sequence[j] == prev:
            current += 1
            if current > longest:
                longest = current
        else:
            prev = sequence[j]
            current = 0

    return longest


cpdef int m_stretch(unsigned char[:] sequence):
    """Return the maximum length of a single letter (nucleotide) repeat in a string.

    Args:
        sequence (uchar[]): Nucleotide sequence writeable memory view.

    Returns:
        (int) Length of maximum run of a single letter.

    """
    cdef:
        Py_ssize_t length = sequence.shape[0]

    return c_stretch(&sequence[0], length)


cpdef int stretch(str sequence):
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
    cdef:
        int longest
        StringView view = str_to_view(sequence)

    longest = c_stretch(view.ptr, view.size)
    free(view.ptr)

    return longest


cdef inline bint _compare(
    unsigned char* p,
    unsigned char* q,
    Py_ssize_t start,
    Py_ssize_t end
) noexcept:
    """Awkward slice comparison between two different size chars."""
    cdef:
        Py_ssize_t j, count = 0

    for j in range(start, end):
        if p[j] != q[count]:
            return False
        count += 1

    return True


cdef inline void _assign(
    unsigned char* src,
    unsigned char* dest,
    Py_ssize_t start,
    Py_ssize_t end
):
    """Overcome slice assignment between two char variables of different sizes."""
    cdef:
        Py_ssize_t j, count = 0

    for j in range(start, end):
        dest[count] = src[j]
        count += 1


cdef inline int c_nrepeats(
    unsigned char* sequence,
    int length,
    int n
) noexcept:
    """Calculate the maximum observed repeats of composite pattern size n characters.

    Args:
        sequence (str): Nucleotide sequence string.
        n (int): Size of k-mers (composite pattern) to observe.

    Returns:
        (int) The longest tandem run of nucleotides comprised of a composite pattern
        of length n characters.

    """
    cdef:
        Py_ssize_t t = <Py_ssize_t> n
        Py_ssize_t v = length // t
        Py_ssize_t j, k
        int current, max_val = 0
        unsigned char* previous = <unsigned char *> malloc(
            (t + 1) * sizeof(unsigned char)
        )

    for k in range(t):
        _assign(sequence, previous, k, t + k)
        current = 0
        for j in range(1, v):
            if _compare(sequence, previous, j * t + k, j * t + k + t):
                current += 1
                if current > max_val:
                    max_val = current
            else:
                current = 0
                _assign(sequence, previous, j * t + k, j * t + k + t)

    free(previous)

    return max_val


cpdef int m_nrepeats(unsigned char[:] sequence, int n):
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
    return c_nrepeats(&sequence[0], sequence.shape[0], n)


cpdef int nrepeats(str sequence, int n):
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
    cdef:
        StringView view = str_to_view(sequence)
        int max_val

    max_val = c_nrepeats(view.ptr, view.size, n)
    free(view.ptr)

    return max_val
