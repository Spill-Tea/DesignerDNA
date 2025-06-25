# cython: boundscheck=False, wraparound=False, nonecheck=False
"""Cythonized oligonucleotide functions."""

from libc.string cimport memcpy
from libc.stdlib cimport malloc, free

cdef extern from "Python.h":
    str PyUnicode_FromStringAndSize(char*, Py_ssize_t)
    Py_ssize_t PyUnicode_GET_LENGTH(object)
    # bint PyBytes_Check(object)
    # char* PyBytes_AS_STRING(object)
    # Py_ssize_t PyBytes_GET_SIZE(object)

cdef extern from "oligos.h":
    const unsigned char DNA[0x100]
    const unsigned char RNA[0x100]

# ctypedef fused StrT:
#     str
#     bytes


cdef struct StringView:
    char* ptr
    Py_ssize_t size


cdef inline StringView to_view(str sequence):
    """Construct StringView, using Cpython C-API to construct a c char string."""
    cdef:
        Py_ssize_t length = PyUnicode_GET_LENGTH(sequence)
        bytes temp = sequence.encode("utf8")
        char* buffer = temp
        StringView view

    view.ptr = <char *> malloc((length + 1) * sizeof(char))
    memcpy(view.ptr, buffer, length + 1)
    view.ptr[length] = "\0"  # c string terminator
    view.size = length

    return view


cdef inline str to_str(StringView view):
    """Convert StringView back into a python string object, safely releasing memory."""
    cdef str obj = PyUnicode_FromStringAndSize(view.ptr, view.size)
    free(view.ptr)

    return obj


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
    cdef StringView view = to_view(sequence)

    if dna:
        c_complement(view.ptr, view.size, DNA)
    else:
        c_complement(view.ptr, view.size, RNA)

    return to_str(view)


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
    cdef StringView view = to_view(sequence)

    if dna:
        c_reverse_complement(view.ptr, view.size, DNA)
    else:
        c_reverse_complement(view.ptr, view.size, RNA)

    return to_str(view)


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
        bytes temp = sequence.encode("utf8")
        bytes comp = complement(sequence, dna).encode("utf8")
        bytes even, pal = b""
        Py_ssize_t i, current, seq_length = len(sequence), length = 0

    if seq_length < 2:  # noqa: PLR2004
        return ""

    for i in range(seq_length - 1):
        # NOTE: Palindromic nucleotides are only even length, halving search space
        even = _expand_from_center(temp, comp, i, i + 1, seq_length)
        current = len(even)
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
        bytes temp = sequence.encode("utf8")
        char* buffer = temp
        int longest = 0, current = 0
        char c, prev = buffer[0]

    for c in buffer[1:]:
        if c == prev:
            current += 1
            if current > longest:
                longest = current
        else:
            current = 0
            prev = c

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
    if n < 1:
        raise ValueError("n must be greater than 0.")
    if n == 1:
        return stretch(sequence)

    cdef:
        bytes phase, temp = sequence.encode("utf8")
        char* buffer = temp
        int i, j, k, max_val = 0
        list[char] previous = [buffer[i : n + i] for i in range(n)]
        list[int] current = [0 for _ in range(n)]

    for j in range(n, len(sequence), n):
        for k in range(n):
            phase = buffer[j + k : j + k + n]
            if phase == previous[k]:
                current[k] += 1
                if current[k] > max_val:
                    max_val = current[k]
            else:
                current[k] = 0
                previous[k] = phase

    return max_val
