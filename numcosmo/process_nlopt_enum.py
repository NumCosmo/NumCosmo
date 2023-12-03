#!/usr/bin/env python3
#
# process_nlopt_enum.py
#
# Wed Nov 22 21:51:44 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# process_nlopt_enum.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Process generated nlopt enum header file."""

import re
import sys


def process_file(input_file, output_file):
    """Process the input removing the comments and replacing the enum names."""

    with open(input_file, "r", encoding="utf-8") as file:
        lines = file.readlines()

    lines = [line for line in lines if not line.startswith("#")]
    lines = [re.sub("nlopt_algorithm", "NcmFitNloptAlgorithm", line) for line in lines]
    lines = [re.sub("nlopt_result", "NcmFitNloptResult", line) for line in lines]

    with open(output_file, "w", encoding="utf-8") as file:
        file.writelines(lines)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    process_file(input_file_path, output_file_path)
