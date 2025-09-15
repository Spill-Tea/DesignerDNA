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

# distutils: language = c++

"""Oligonucleotide functions with the help of C++."""

from libc.stdlib cimport free

cdef extern from "Python.h":
    str PyUnicode_Join(str, str)

from common cimport StringView, str_to_view
from narray cimport NumericArray

from designer_dna._oligos cimport v_complement


cdef inline void _compute(
    unsigned char* s,
    unsigned char* c,
    NumericArray[int]* arr,
    ssize_t n,
):
    """Primary computation behind manacher's algorithm.

    Args:
        s (uchar*): nucleotide sequence
        c (uchar*): complement of nucleotide sequence
        arr (NumericArray[int]*): an array of integers
        n (ssize_t): length of input sequence, s.

    Returns:
        (void) relevant data saved in place to NumericArray

    """
    cdef:
        ssize_t mirror, a, b, i, stemp, center = 0, radius = 0
        int temp, zero = 0

    for i in range(1, n - 1):
        # parity: 0 1 0 1 0 1 0 1 0
        # string: ^ # A # N # T # $
        #  index: 0 1 2 3 4 5 6 7 8

        # skip odd length palindromes (when index is centered on an even position)
        # if character from ref seq and complement do not match
        if i % 2 == 0 and s[i] != c[i]:
            arr[0][i] = zero
            continue

        # Look ahead at mirror position
        mirror = 2 * center - i
        if i < radius:
            temp = <int> (radius - i)
            arr[0][i] = min(temp, arr[0][mirror])
        else:
            arr[0][i] = zero

        # Center expansion method
        stemp = <ssize_t> arr[0][i]
        a = i + 1 + stemp
        b = i - 1 - stemp
        while s[a] == c[b] and s[b] == c[a]:
            arr[0][i] += 1
            a += 1
            b -= 1

        stemp = <ssize_t> arr[0][i]
        if i + stemp > radius:
            center = i
            radius = i + stemp


cpdef str manacher(str sequence, bint dna = True):
    """Find the longest palindromic substring within a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.
        dna (bool): Sequence is DNA, else RNA.

    Returns:
        (str) Longest palindromic substring within a sequence.

    Notes:
        * This is a cython/c++ implementation of the O(n) Manacher's algorithm.

    """
    cdef:
        str k = PyUnicode_Join("#", f"^{sequence}$")
        StringView ref = str_to_view(k)
        StringView com = str_to_view(k)
        NumericArray[int]* arr
        ssize_t i, center = 0
        int radius = 0

    arr = new NumericArray[int](ref.size)
    arr.fill(radius)
    v_complement(&com, dna)

    _compute(ref.ptr, com.ptr, arr, <ssize_t> ref.size)
    free(ref.ptr)
    free(com.ptr)

    # Enumerate, capturing index (center) at value of max (radius)
    for i in range(1, ref.size - 1):
        if arr[0][i] > radius:
            radius = arr[0][i]
            center = i
    del arr

    # By nature, a palindrome is symmetrical around center (+/- radius)
    return sequence[(center - radius + 1) // 2 - 1: (center + radius) // 2]
