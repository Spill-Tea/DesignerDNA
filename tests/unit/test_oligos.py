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

"""Unit test oligos module."""

from functools import wraps
from typing import Any, Callable

import pytest

from designer_dna import _oligos, oligos


def wrapper(f: Callable[[str], str]) -> Callable[[str], str]:
    """Wrap a function which acts inplace on string."""

    @wraps(f)
    def inner(s: str, *args) -> str:
        b: bytes = s.encode("utf8")
        ba: bytearray = bytearray(b)
        m: memoryview = memoryview(ba)

        f(m, *args)

        return ba.decode("utf8")

    return inner


def wrap_out(f: Callable[[str], Any]) -> Callable[[str], Any]:
    """Wrap a function which uses a string, but provides a different output."""

    @wraps(f)
    def inner(s: str, *args) -> Any:
        b: bytes = s.encode("utf8")
        ba: bytearray = bytearray(b)
        m: memoryview = memoryview(ba)

        return f(m, *args)

    return inner


@pytest.mark.parametrize(
    "function",
    [
        oligos.reverse,
        wrapper(_oligos.m_reverse),
        oligos.reverse_py,
    ],
)
@pytest.mark.parametrize(
    ["seq", "expected"],
    [
        ("", ""),
        ("A", "A"),
        ("AT", "TA"),
        ("TAT", "TAT"),
        ("GGC", "CGG"),
        ("GATC", "CTAG"),
    ],
)
def test_reverse(seq: str, expected: str, function: Callable[[str], str]) -> None:
    """Test reversing a nucleotide sequence."""
    result: str = function(seq)
    assert result == expected, "Unexpected reverse seq."


# 0, and 1 length complements --> important for cython woes
complements = [
    ("", True, ""),
    ("", False, ""),
]

for _k, _v in oligos.BASEPAIRS_DNA.items():
    for _j in (True, False):
        complements.append((_k, _j, _v if _j else oligos.BASEPAIRS_RNA[_k]))


@pytest.mark.parametrize(
    "function",
    [
        oligos.complement,
        wrapper(_oligos.m_complement),
        oligos.complement_py,
    ],
)
@pytest.mark.parametrize(
    ["seq", "dna", "expected"],
    [
        *complements,
        ("AA", True, "TT"),
        ("AA", False, "UU"),
        ("AT", True, "TA"),
        ("TAT", True, "ATA"),
        ("GGC", True, "CCG"),
        ("GATC", True, "CTAG"),
        ("GATC", False, "CUAG"),
        ("AGTCNURYSWKMBVDH-.", True, "TCAGNAYRSWMKVBHD-."),
        ("AGTCNURYSWKMBVDH-.", False, "UCAGNAYRSWMKVBHD-."),
        ("agtcnuryswkmbvdh-.", True, "TCAGNAYRSWMKVBHD-."),
        ("agtcnuryswkmbvdh-.", False, "UCAGNAYRSWMKVBHD-."),
        (
            "".join(oligos.BASEPAIRS_DNA.keys()),
            True,
            "".join(oligos.BASEPAIRS_DNA.values()),
        ),
        (
            "".join(oligos.BASEPAIRS_RNA.keys()),
            False,
            "".join(oligos.BASEPAIRS_RNA.values()),
        ),
    ],
)
def test_complement(
    seq: str,
    dna: bool,
    expected: str,
    function: Callable[[str, bool], str],
) -> None:
    """Test complement of a nucleotide sequence."""
    result: str = function(seq, dna)
    assert result == expected, f"Unexpected complement seq: {result}"


@pytest.mark.parametrize(
    "function",
    [
        oligos.reverse_complement,
        wrapper(_oligos.m_reverse_complement),
        oligos.reverse_complement_py,
    ],
)
@pytest.mark.parametrize(
    ["seq", "dna", "expected"],
    [
        *complements,
        ("AT", True, "AT"),
        ("AT", False, "AU"),
        ("TAT", True, "ATA"),
        ("GGC", True, "GCC"),
        ("GATC", True, "GATC"),
        ("GATC", False, "GAUC"),
        ("AGTCNURYSWKMBVDH-.", True, ".-DHBVKMWSRYANGACT"),
        ("AGTCNURYSWKMBVDH-.", False, ".-DHBVKMWSRYANGACU"),
        ("agtcnuryswkmbvdh-.", True, ".-DHBVKMWSRYANGACT"),
        ("agtcnuryswkmbvdh-.", False, ".-DHBVKMWSRYANGACU"),
    ],
)
def test_reverse_complement(
    seq: str,
    dna: bool,
    expected: str,
    function: Callable[[str, bool], str],
) -> None:
    """Test reverse complement of a nucleotide sequence."""
    result: str = function(seq, dna)
    assert result == expected, "Unexpected reverse complement."


@pytest.mark.parametrize(
    "function",
    [
        oligos.stretch,
        wrap_out(_oligos.m_stretch),
        oligos.stretch_py,
    ],
)
@pytest.mark.parametrize(
    ["seq", "expected"],
    [
        ("", 0),
        ("A", 0),
        ("AA", 1),
        ("TAA", 1),
        ("AAT", 1),
        ("ATGC", 0),
        ("AAAAACCCCCCGGGGGGG", 6),
    ],
)
def test_stretch(seq, expected: int, function: Callable[[str], int]) -> None:
    """Test calculation of longest observed run of a single nucleotide within a seq."""
    result: int = function(seq)
    assert result == expected, f"Unexpected stretch calculation: {result}"


@pytest.mark.parametrize("function", [oligos.nrepeats, oligos.nrepeats_py])
@pytest.mark.parametrize(
    ["seq", "n", "expected"],
    [
        ("", 1, 0),
        ("A", 1, 0),
        ("ATGC", 1, 0),
        ("AAAAACCCCCCGGGGGGG", 1, 6),
        ("ACACAC", 2, 2),
        ("ATC" * 4, 3, 3),
        ("A" + "ATC" * 4, 3, 3),
        ("AG" + "ATC" * 4, 3, 3),
        ("AG" + "ATC" * 4 + "C", 3, 3),
        ("AG" + "ATC" * 4 + "CG", 3, 3),
        ("G" + "ATC" * 4 + "C", 3, 3),
        ("G" + "ATC" * 4 + "CG", 3, 3),
    ],
)
def test_nrepeats(
    seq,
    n: int,
    expected: int,
    function: Callable[[str, int], int],
) -> None:
    """Test calculation of longest observed run of n characters within a seq."""
    result: int = function(seq, n)
    assert result == expected, f"Unexpected stretch calculation: {result}"


@pytest.mark.parametrize(
    "function",
    [
        oligos.palindrome,
        oligos.palindrome_py,
        oligos.manacher,
    ],
)
@pytest.mark.parametrize(
    ["seq", "dna", "expected"],
    [
        ("", True, ""),
        ("", False, ""),
        ("A", True, ""),
        ("A", False, ""),
        ("AAAA", True, ""),
        ("AAAA", False, ""),
        ("ATATATATATAT", True, "ATATATATATAT"),
        ("ATATATATATAT", False, ""),
        ("TGGATCCA", True, "TGGATCCA"),
        ("ATGGATCCA", True, "TGGATCCA"),
        ("AATGGATCCA", True, "TGGATCCA"),
        ("TGGATCCAT", True, "TGGATCCA"),
        ("TGGATCCATT", True, "TGGATCCA"),
        ("GAATTC", True, "GAATTC"),
        ("ATGAATTC", True, "GAATTC"),
        ("CTTAAG", True, "CTTAAG"),
        ("ANT", True, "ANT"),
        ("AANT", True, "ANT"),
        ("AWSNSWT", True, "AWSNSWT"),
        ("AWSSWT", True, "AWSSWT"),
    ],
)
def test_palindromes(
    seq: str,
    dna: bool,
    expected: str,
    function: Callable[[str, int], int],
) -> None:
    """Test detection of longest palindrome within a sequence."""
    result = function(seq, dna)
    assert result == expected, f"Unexpected palindrome: {result}"
    if result:
        assert result == oligos.reverse_complement(result)
