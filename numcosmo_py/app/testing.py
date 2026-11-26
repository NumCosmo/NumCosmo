#
# testing.py
#
# Thu Nov 26 16:23:29 2026
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Helpers for tests of the numcosmo command-line application."""

import re


def plain_output(output: str) -> str:
    """Return the text of a command's output for substring matching.

    The ANSI codes (CI forces colors) are removed, the Rich box-drawing characters
    become spaces, and every run of whitespace, line breaks included, becomes one
    space, so a message matches whatever the terminal width.

    :param output: the output of the command
    :return: the output as one line of text
    """
    output = re.sub(r"\x1b\[[0-9;]*m", "", output)
    output = re.sub(r"[─-╿]", " ", output)
    return " ".join(output.split())
