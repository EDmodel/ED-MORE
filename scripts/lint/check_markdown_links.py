# SPDX-License-Identifier: Apache-2.0
"""Every relative link and image in a Markdown file points at a file or folder the repository
tracks.

On GitHub a link to a file the repository does not hold is broken. The usual cause is a figure that
.gitignore kept out without anyone noticing, or a page that moved. Web links and anchors within a
page are not checked.
"""
import posixpath
import re
from urllib.parse import unquote

from lintlib import Problem

INLINE = re.compile(r"!?\[[^\]]*\]\(\s*(<[^>]*>|[^)\s]+)")
REFERENCE = re.compile(r"^\s{0,3}\[(?!\^)[^\]]+\]:\s*(\S+)")   # a footnote, [^1]:, is not a link
HTML = re.compile(r"\b(?:src|href)\s*=\s*([\"'])(.+?)\1")
CODE_SPAN = re.compile(r"(`+).+?\1")
FENCE = re.compile(r"^\s*(`{3,}|~{3,})(.*)$")
EXTERNAL = re.compile(r"^(?:[a-zA-Z][a-zA-Z0-9+.-]*:|//|#)")


def link_lines(text: str):
    """Each line that can hold a link, as (number, line), leaving out fenced code and HTML comments.
    A fence may be indented, as under a list item, and closes on a bare fence at least as long."""
    fence = ""
    in_comment = False
    for no, line in enumerate(text.split("\n"), 1):
        marker = FENCE.match(line)
        if fence:
            if marker and marker.group(1)[0] == fence[0] and len(marker.group(1)) >= len(fence) \
                    and not marker.group(2).strip():
                fence = ""
            continue
        if marker and not in_comment:
            fence = marker.group(1)
            continue
        kept = ""
        while line:
            if in_comment:
                end = line.find("-->")
                if end < 0:
                    line = ""
                else:
                    line, in_comment = line[end + 3:], False
            else:
                start = line.find("<!--")
                if start < 0:
                    kept, line = kept + line, ""
                else:
                    kept, line, in_comment = kept + line[:start], line[start + 4:], True
        yield no, kept


def targets(line: str):
    """The link targets on a line. A reference definition must start the line, so it is matched
    before code spans are removed: `x` [unit]: text is prose, not a definition."""
    for match in REFERENCE.finditer(line):
        yield match.group(1).strip("<>")
    prose = CODE_SPAN.sub("", line)
    for match in INLINE.finditer(prose):
        yield match.group(1).strip("<>")
    for match in HTML.finditer(prose):
        yield match.group(2)


def check(repo):
    for path in repo.text_files((".md",)):
        folder = posixpath.dirname(path)
        for no, line in link_lines(repo.text(path)):
            for target in targets(line):
                if not target or EXTERNAL.match(target):
                    continue
                rel = unquote(target.split("#", 1)[0].split("?", 1)[0])
                if not rel:
                    continue
                rel = posixpath.normpath(rel.lstrip("/") if rel.startswith("/") else posixpath.join(folder, rel))
                if rel.startswith(".."):
                    yield Problem(path, no, f"link '{target}' points outside the repository")
                elif not repo.exists(rel):
                    yield Problem(path, no, f"link '{target}' points at '{rel}', which the repository does not hold")
