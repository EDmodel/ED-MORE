# SPDX-License-Identifier: Apache-2.0
"""Every source file states its licence, `SPDX-License-Identifier: Apache-2.0`, in a comment on its
first line, or on the second line when the first is a `#!` line.

The identifier lets anyone who copies a single file, and the tools that check licences, see the
terms without the rest of the repository.
"""
from lintlib import Problem

TAG = "SPDX-License-Identifier: Apache-2.0"
SOURCE_SUFFIXES = (".f90", ".F90", ".fpp", ".py", ".sh", ".cmake", ".c", ".h")
SOURCE_NAMES = ("CMakeLists.txt",)


def check(repo):
    for path in repo.text_files(SOURCE_SUFFIXES, SOURCE_NAMES):
        head = repo.text(path).split("\n", 2)[:2]
        line = head[1] if head[0].startswith("#!") and len(head) > 1 else head[0]
        if TAG not in line:
            yield Problem(path, 0, f"no '{TAG}' comment at the top")
