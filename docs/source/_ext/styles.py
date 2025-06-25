"""pygment styles."""

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
