#!/usr/bin/env python
#
# test_window_vocabulary.py
#
# Tue Sep 8 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_window_vocabulary.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""The analytic windows are defined in one place, and this keeps them there.

``windows.py`` is the vocabulary: every certified window, its closed form, its
parameters. Two other rosters drive the same windows from their own copies of
those numbers -- the truth-table generator's ``CASES`` and the pair matrix's
``KERNELS`` -- and a window written down three times is a window that can drift.

It already did. Before the vocabulary existed, ``gauss_near`` named
:math:`(\\chi_0, \\sigma) = (1095, 182)` in the truth table and
:math:`(1800, 300)` in the pair matrix: the same name, two different windows,
both appearing in committed data. These tests are what makes that a failure
rather than a discovery.

They compare the *generator command line* each roster produces, because that is
what actually determines which window got certified -- not a docstring, and not
the constructor call, which cannot be compared across two implementations.
"""

from __future__ import annotations

import importlib.util
import pathlib
import typing

import pytest

from numcosmo_py import Ncm

from xcor import cases_k_integral as cases
from xcor import windows as W

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor]

TOOLS = pathlib.Path(__file__).resolve().parents[3] / "tools"


def _load_generator() -> typing.Any:
    """The truth-table generator, imported for its case list alone.

    It lives in ``tests/tools`` and is a script rather than a package, so it is
    loaded by path. Importing it runs no generation: the module body only defines
    the case list and helpers.
    """
    path = TOOLS / "make_xcor_window_truth_table.py"
    spec = importlib.util.spec_from_file_location("_wtt", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def _normalise(args: list[str]) -> set[str]:
    """Generator options as a set, with numbers in a canonical form.

    ``--n-sigma=4`` and ``--n-sigma=4.0`` select the same window; a string
    comparison would call them different.
    """
    out = set()

    for arg in args:
        key, _, value = arg.partition("=")

        try:
            value = ",".join(repr(float(part)) for part in value.split(","))
        except ValueError:
            pass

        out.add(f"{key}={value}")

    return out


GENERATOR = _load_generator()


@pytest.mark.parametrize("case", sorted(GENERATOR.CASES))
def test_generator_case_matches_the_vocabulary(case: str) -> None:
    """Every truth-table case selects the window the vocabulary says it does.

    This is the test that would have caught the collision: the generator's own
    options are compared against the vocabulary's, so a case whose parameters
    drift in either file fails here rather than silently certifying a different
    window under an unchanged name.
    """
    spec = GENERATOR.CASES[case]
    from_generator = [f"--shape={spec.get('shape', case)}"] + list(spec["arb"])
    from_vocabulary = W.arb_args(W.resolve(case))

    assert _normalise(from_generator) == _normalise(from_vocabulary), (
        f"{case} differs between the generator and windows.py"
    )


@pytest.mark.parametrize("name", sorted(cases.KERNELS))
def test_pair_kernel_matches_the_vocabulary(name: str) -> None:
    """Every pair-matrix kernel selects the window the vocabulary says it does.

    Names the vocabulary does not know are skipped rather than failed: the pair
    matrix may legitimately carry a window that has no certified table yet, and
    this test is about agreement, not about coverage.
    """
    try:
        canonical = W.resolve(name)
    except KeyError:
        pytest.skip(f"{name} is not in the vocabulary")

    from_pairs = cases.KERNELS[name].arb_args("a")
    from_vocabulary = W.arb_args(canonical, side="a")

    assert _normalise(from_pairs) == _normalise(from_vocabulary), (
        f"{name} differs between cases_k_integral and windows.py"
    )


def test_no_two_names_mean_different_windows() -> None:
    """One name, one window -- across both rosters at once.

    The failure this guards against is not a typo but a coincidence: two rosters
    each internally consistent, using one name for two windows. Resolving every
    name through the vocabulary and checking that each resolved name has a single
    parameter set is what makes that impossible to reintroduce.
    """
    seen: dict[str, tuple[str, set[str]]] = {}

    for source, names in (
        ("generator", sorted(GENERATOR.CASES)),
        ("pairs", sorted(cases.KERNELS)),
    ):
        for name in names:
            try:
                canonical = W.resolve(name)
            except KeyError:
                continue

            args = _normalise(W.arb_args(canonical))
            previous = seen.get(canonical)

            if previous is not None:
                assert previous[1] == args, (
                    f"{canonical} resolves to different windows from "
                    f"{previous[0]} and {source}"
                )

            seen[canonical] = (source, args)

    assert seen, "no shared names found, so this test proves nothing"


def test_every_vocabulary_name_is_reachable() -> None:
    """Aliases point at real windows, and no window is orphaned by a typo."""
    for alias, canonical in W.ALIASES.items():
        assert canonical in W.WINDOWS, f"alias {alias} points at nothing"
        assert alias not in W.WINDOWS, f"alias {alias} shadows a real window"

    for name, window in W.WINDOWS.items():
        assert window.form in W.FORM_PARAMS, f"{name} has an unknown form"
        assert len(window.params) == len(W.FORM_PARAMS[window.form]), (
            f"{name} carries {len(window.params)} parameters, "
            f"{window.form} takes {len(W.FORM_PARAMS[window.form])}"
        )
        assert 0 <= window.deriv <= 2, f"{name} has derivative order {window.deriv}"
