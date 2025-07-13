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

"""Customized python lexer."""

from collections import deque
from collections.abc import Iterator
from typing import ClassVar

from pygments.lexer import bygroups, include
from pygments.lexers.python import PythonLexer
from pygments.token import (
    Comment,
    Keyword,
    Name,
    Number,
    Punctuation,
    String,
    Text,
    Whitespace,
    _TokenType,
)
from utils import get_bracket_level


def _find(it, obj, key=lambda a, b: a == b) -> int:
    for n, j in enumerate(it):
        if key(j, obj):
            return n
    raise IndexError("Unable to find object.")


def _get_index(n: int):
    def inner(a, b) -> bool:
        return a[n] == b

    return inner


root: list = [
    (r"\n", Whitespace),
    (  # single line docstrings (edge case)
        r'^(\s*)([rRuUbB]{,2})("""(?:.)*?""")',
        bygroups(Whitespace, String.Affix, String.Doc),
    ),
    (  # Modfied triple double quote docstrings to highlight docstring titles
        r'^(\s*)([rRuUbB]{,2})(""")',
        bygroups(Whitespace, String.Affix, String.Doc),
        "docstring-double",
    ),
    (  # Intentionally treat text encapsulated within single triple quotes as String
        r"^(\s*)([rRuUbB]{,2})('''(?:.|\n)*?''')",
        bygroups(Whitespace, String.Affix, String),
    ),
    (r"\A#!.+$", Comment.Hashbang),
    (
        # Format Special Common Keyword Comments
        # NOTE: Must come before Comment.Single token in order to be matched.
        r"(#\s*)(TODO|FIXME|NOTE|BUG|HACK|XXX)(:?)(.*$)",
        bygroups(Comment.Single, Comment.Special, Comment.Special, Comment.Single),
    ),
    (r"#.*$", Comment.Single),
    (r"\\\n", Text),
    (r"\\", Text),
    include("keywords"),
    include("soft-keywords"),
    (
        r"(def)((?:\s|\\\s)+)",
        bygroups(Keyword.Declare, Whitespace),
        "funcname",
    ),
    (
        r"(class)((?:\s|\\\s)+)",
        bygroups(Keyword.Declare, Whitespace),
        "classname",
    ),
    (
        r"(from)((?:\s|\\\s)+)",
        bygroups(Keyword.Namespace, Whitespace),
        "fromimport",
    ),
    (
        r"(import)((?:\s|\\\s)+)",
        bygroups(Keyword.Namespace, Whitespace),
        "import",
    ),
    include("expr"),
]


python_tokens: dict[str, list] = PythonLexer.tokens.copy()
python_tokens["root"] = root
python_tokens["docstring-double"] = [
    (
        r"(?<=\n)(\s*)(Args|Attributes|Returns|Raises|"
        r"Examples|Yields|References|Notes|Equations)(:)(\s*)",
        bygroups(Whitespace, String.Doc.Title, String.Doc, Whitespace),
    ),
    (r'^\s*(?:""")', String.Doc, "#pop"),
    (r".+[\r\n]*", String.Doc),
]

# Tokenize function names when used (i.e. function calls)
# NOTE: Must be inserted before general `Name` token but after `Name.Builtins` token
# NOTE: Implementation limitations -> we cannot distinguish between class and function
#       calls using regex based parsing alone (i.e without semantic analysis).
python_tokens["name"].insert(
    _find(python_tokens["name"], Name, _get_index(1)),
    (r"\b([a-zA-Z_]\w*)(?=\s*\()", Name.Function),
)

python_tokens["numbers"] = [
    (
        r"(\d(?:_?\d)*\.(?:\d(?:_?\d)*)?|(?:\d(?:_?\d)*)?\.\d(?:_?\d)*)"
        r"([eE][+-]?\d(?:_?\d)*)?([jJ]?)",
        bygroups(Number.Float, Number.Float, Number.Other),
    ),
    (r"(\d(?:_?\d)*[eE][+-]?\d(?:_?\d)*)([jJ]?)", bygroups(Number.Float, Number.Other)),
    (r"(0[oO])((?:_?[0-7])+)", bygroups(Number.Other, Number.Oct)),
    (r"(0[bB])((?:_?[01])+)", bygroups(Number.Other, Number.Bin)),
    (r"(0[xX])((?:_?[a-fA-F0-9])+)", bygroups(Number.Other, Number.Hex)),
    (r"(\d(?:_?\d)*)([jJ]?)", bygroups(Number.Integer, Number.Other)),
]


class CustomPythonLexer(PythonLexer):
    """Enhanced regex-based python Lexer.

    Notes:
        1. Implemented a simple stack based rainbow bracket colorizer.
            * limitation: Only detects errors that close more brackets than opens.
        2. Highlight Docstring titles (assumes google docstring format)
        3. Improved highlighting function calls (with limitations)
        4. Modify display of number components which indicate a different base number.

    """

    n_brackets: int
    _stack: deque[int]
    tokens: ClassVar[dict[str, list]] = python_tokens

    def __init__(self, **options) -> None:
        super().__init__(**options)
        self._stack = deque[int]()
        self.n_brackets = int(options.get("n_brackets", 4))

    def _enter(self) -> _TokenType:
        """Retrieve next token in cycle."""
        idx = len(self._stack) % self.n_brackets
        self._stack.append(idx)

        return get_bracket_level(idx)

    def _exit(self) -> _TokenType:
        """Remove element from stack and return token."""
        try:
            idx: int = self._stack.pop()
            return get_bracket_level(idx)

        except IndexError:
            return Punctuation.Error

    def get_tokens_unprocessed(
        self,
        text,
        stack=("root",),
    ) -> Iterator[tuple[int, _TokenType, str]]:
        _token: _TokenType
        for idx, token, value in super().get_tokens_unprocessed(text, stack):
            _token = token
            if token is Name and value.isupper():
                _token = Name.Constant

            elif token is Punctuation:
                match value:
                    case "(" | "[" | "{":
                        _token = self._enter()
                    case "}" | "]" | ")":
                        _token = self._exit()
                    case _:
                        ...

            yield idx, _token, value
