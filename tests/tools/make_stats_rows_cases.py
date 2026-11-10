#!/usr/bin/env python3
#
# make_stats_rows_cases.py
#
# Fri Oct 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_stats_rows_cases.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Print the reference table of tests/c/ncm/stats/test_ncm_stats_rows.c.

The sample is an AR(1) series per column, 1500 rows and 3 columns, drawn from NcmRNG
seed 12345 exactly as the C test draws it: coefficients 0.9, 0.5 and 0, started at 30, 15
and 6, shifted by 1, 2 and 3. The table holds the results of the six row analyses on it,
through the NcmStatsVec saved-row API, at full double precision.

The committed table was printed at commit 701a5314, the last one with the original
NcmStatsVec implementation of the analyses; the test requires the moved functions to
reproduce it to its tolerance. Printing it from a later version only shows what the current code
gives. The output replaces the block between the BEGIN and END markers of the test:

    python3 tests/tools/make_stats_rows_cases.py
"""

from numcosmo_py import Ncm

NROWS = 1500
NCOLS = 3
SEED = 12345
RHO = (0.9, 0.5, 0.0)
START = (30.0, 15.0, 6.0)
VISUAL_FIRST_INDEX = 7


def sample() -> Ncm.StatsVec:
    """Draw the sample into a NcmStatsVec with saved rows."""
    rng = Ncm.RNG.seeded_new(None, SEED)
    sv = Ncm.StatsVec.new(NCOLS, Ncm.StatsVecType.COV, True)
    x = list(START)
    for _ in range(NROWS):
        v = Ncm.Vector.new(NCOLS)
        for p in range(NCOLS):
            x[p] = RHO[p] * x[p] + rng.gaussian_gen(0.0, 1.0)
            v.set(p, x[p] + (p + 1.0))
        sv.append(v, True)
    return sv


def c_array(values: list[float]) -> str:
    """Format a C initializer list."""
    return "{" + ", ".join(repr(v) for v in values) + "}"


def main() -> None:
    """Print the table."""
    Ncm.cfg_init()
    sv = sample()

    ess, ess_bindex, ess_wp, ess_wp_order, _ = sv.max_ess_time(0)
    pvals, hd_bindex, hd_wp, hd_wp_order, _ = sv.heidel_diag(0, 0.0)
    cs, mean, var = sv.visual_heidel_diag(0, VISUAL_FIRST_INDEX)
    cs = cs.dup_array()
    const_break = [sv.estimate_const_break(p) for p in range(NCOLS)]
    cov_diag = sv.compute_cov_robust_diag()
    cov_ogk = sv.compute_cov_robust_ogk()

    print(f"static const gdouble ref_ess[TEST_NCOLS] = {c_array(ess.dup_array())};")
    print(f"static const gint ref_ess_bindex = {ess_bindex};")
    print(f"static const guint ref_ess_wp = {ess_wp};")
    print(f"static const guint ref_ess_wp_order = {ess_wp_order};")
    print(f"static const gdouble ref_pvals[TEST_NCOLS] = {c_array(pvals.dup_array())};")
    print(f"static const gint ref_hd_bindex = {hd_bindex};")
    print(f"static const guint ref_hd_wp = {hd_wp};")
    print(f"static const guint ref_hd_wp_order = {hd_wp_order};")
    print(f"static const gdouble ref_cs_first[3] = {c_array(cs[:3])};")
    print(f"static const gdouble ref_cs_last[3] = {c_array(cs[-3:])};")
    print(f"static const gdouble ref_visual_mean = {mean!r};")
    print(f"static const gdouble ref_visual_var = {var!r};")
    print(f"static const gdouble ref_const_break[TEST_NCOLS] = {c_array(const_break)};")
    print(
        "static const gdouble ref_cov_diag[TEST_NCOLS] = "
        f"{c_array([cov_diag.get(i, i) for i in range(NCOLS)])};"
    )
    rows = ",\n  ".join(
        c_array([cov_ogk.get(i, j) for j in range(NCOLS)]) for i in range(NCOLS)
    )
    print(
        "static const gdouble ref_cov_ogk[TEST_NCOLS][TEST_NCOLS] = {\n  "
        + rows
        + ",\n};"
    )


if __name__ == "__main__":
    main()
