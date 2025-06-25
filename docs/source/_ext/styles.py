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

"""Custom Pygment styles."""

from typing import ClassVar

from pygments.style import Style
from pygments.token import (
    Comment,
    Error,
    Keyword,
    Name,
    Number,
    Operator,
    Punctuation,
    String,
    _TokenType,
)


class VSCodeDarkPlus(Style):
    """VSCode Dark+ Style."""

    background_color: str = "#1E1E1E"

    styles: ClassVar[dict[_TokenType, str]] = {
        Number: "#B6CEA9",
        Operator: "#D4D4D4",
        Operator.Word: "#C586C0",
        Comment: "#6D9957",
        Comment.Preproc: "#639BD4",
        Keyword.Namespace: "#C287A0",
        # Keyword.Reserved: "#C287A0",
        Keyword.Reserved: "#639BD4",
        Keyword.Type: "#61C8B0",
        Keyword.Constant: "#4FC1FF",
        # Keyword: "#639BD4",
        Keyword: "#C586C0",
        Name: "#7FD0FD",
        Name.Class: "#61C8B0",
        Name.Namespace: "#61C8B0",
        Name.Function: "#DCDCAA",
        # Name.Builtin: "#DCDCAA",
        Name.Builtin: "#4EC9B0",
        Name.Type: "#4EC9B0",
        Name.Builtin.Pseudo: "#9CDCFE",
        Name.Variable: "#9CDCFE",
        Name.Variable.Class: "#61C8B0",
        Name.Variable.Magic: "#DCDCAA",
        Name.Exception: "#61C8B0",
        Error: "#61C8B0",
        String: "#C9937A",
        Punctuation: "#F9C922",
    }
