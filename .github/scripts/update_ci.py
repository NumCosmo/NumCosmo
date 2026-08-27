"""Update the CI environment file."""

import re
import sys

python_re = re.compile(
    r"^(?P<prefix>\s*-\s*python)(?P<rest>(\s*[=<>!].*)$)", re.IGNORECASE
)
python_version = sys.argv[1]

with open("environment.yml", encoding="utf-8") as fh:
    lines = fh.readlines()
found = False
for i, line in enumerate(lines):
    m = python_re.match(line)
    if not m:
        continue
    prefix = m.group("prefix")
    rest = m.group("rest")
    lines[i] = f"{prefix} ={python_version}\n"
    found = True
    break

if found:
    with open("environment-ci.yml", "w", encoding="utf-8") as f:
        f.writelines(lines)
    sys.exit(0)

sys.exit(1)
