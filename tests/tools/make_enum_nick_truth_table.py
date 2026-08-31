#!/usr/bin/env python
"""Record the introspected nick of every enum and flags member.

A member's nick is what introspection turns into its Python name, what
ncm_cfg_get_enum_by_id_name_nick() looks up, and what value_nick prints. When an
enum leaves the stripped prefix implicit, glib-mkenums infers it from whatever
the members happen to share -- so adding one member can silently rename every
other member of that enum. Every enum here pins it with /*< prefix=... >*/ on
its typedef; this table is what proves they all still match.

Regenerate it only when a nick is meant to change:

    python tests/tools/make_enum_nick_truth_table.py
"""

import importlib.util
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "data", "truth_tables", "enums", "enum_nicks.json")

# The collection lives with the test that enforces it, so the two can never
# disagree about what a nick is.
_TEST = os.path.join(HERE, "..", "python", "ncm", "core", "test_enum_nicks.py")
_spec = importlib.util.spec_from_file_location("test_enum_nicks", _TEST)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

collect_enum_nicks = _mod.collect_enum_nicks


def main():
    """Write the truth table."""
    nicks = collect_enum_nicks()

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(nicks, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"{len(nicks)} enum/flags members recorded in {OUT}")


if __name__ == "__main__":
    main()
