#!/usr/bin/env python
#
# windows.py
#
# Mon Sep 8 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# windows.py
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

"""The analytic windows, named once and used by every consumer.

Two rosters used to carry these parameters independently -- the truth-table
generator in ``tests/tools`` and the pair matrix in ``cases_k_integral`` -- which
meant the same window was typed twice under two names, and one name,
``gauss_near``, meant *different windows* in the two files. This module is the
single source: parameters, the closed form each case uses, and what regime it is
there to probe.

**The name is three axes, not one.** A case is a closed *form* at a *place* --
a centre and width -- optionally carrying a Bessel-derivative *weight*:

    <form>_<place>[_d1|_d2|_kdep]

Forms are the certified closed forms (``gauss``, ``tophat``, ``tophat_smooth``,
``student_t``, ``power_exp``, ``lensing``, ``multi``). Places are a small fixed
vocabulary, so a reader can tell from the name where a window sits:

    near        chi ~ 1.1 Gpc, an LSST-SRD Y10 lens bin at the low end
    mid         chi ~ 1.5 Gpc, the historical baseline of this suite
    over        chi ~ 1.8 Gpc, overlapping `mid` rather than displaced from it
    low, high   chi ~ 0.6 and 3.5 Gpc, chosen to be far apart
    far         chi ~ 3.8 Gpc, an SRD Y10 lens bin at the high end
    broad       chi ~ 4.5 Gpc and wide, the scale of an SRD Y1 source bin
    thin        sigma = 50 Mpc, the thin-bin regime
    shell_*     a hard shell, narrow or wide
    skew_far    skewed and high, the shape a source bin takes
    isw         broad and low-weighted, the ISW range
    gal, cmb    lensing sources: a galaxy bin, or the last-scattering surface

Weights: ``_d1`` and ``_d2`` integrate against $j_\\ell'$ and $j_\\ell''$, the
latter being what a redshift-space distortion term carries; ``_kdep`` carries a
scale-dependent growth so the window is not separable in chi and k.

Locations and widths come from the LSST-SRD bins, but only as *regimes*: the
survey spans centres from 1.1 to 5.4 Gpc with sigma/chi from 0.04 to 0.36, and
these cases sit across that range. They are not fits to its dn/dz, which would
track a data release and could not be certified in closed form anyway.
"""

from __future__ import annotations

import typing

# form -> the parameters that form takes, in constructor order. Documented here
# so a case's tuple can be read without opening the generator.
FORM_PARAMS: typing.Final[dict[str, tuple[str, ...]]] = {
    "gauss": ("chi_mean", "chi_sigma", "n_sigma"),
    "tophat": ("chi_lower", "chi_upper"),
    "tophat_smooth": ("chi_lower", "chi_upper", "chi_sigma", "n_sigma"),
    "student_t": ("chi_mean", "chi_scale", "nu", "n_scale"),
    "power_exp": ("chi_scale", "alpha", "beta", "chi_lower", "chi_upper"),
    "lensing": ("chi_lower", "chi_source_lower", "chi_source_upper"),
    "multi": ("chi_mean", "chi_sigma", "weight", "n_sigma"),
}


# The generator's option names, where they differ from the constructor's
# parameter names. Only `multi` does: it takes the bump arrays as --mu/--sigma,
# while the library property is chi-mean/chi-sigma. Keeping the two vocabularies
# explicit is what makes arb_args() reproduce the committed tables exactly.
FORM_ARB_FLAGS: typing.Final[dict[str, tuple[str, ...]]] = {
    "multi": ("mu", "sigma", "weight", "n-sigma"),
}


class Window(typing.NamedTuple):
    """One analytic window: a closed form, its parameters, and its regime.

    :ivar form: the closed form, a key of :data:`FORM_PARAMS`.
    :ivar params: parameters in the order that form takes them.
    :ivar probes: what this case exists to exercise.
    :ivar deriv: Bessel-derivative order of the weight, 0, 1 or 2.
    :ivar kdep: the scale-dependent growth's parameters -- amplitude, transition
        wavenumber and reference distance -- or empty for none. They belong here
        rather than in a flag because they select the window as much as its width
        does, and the generator has to be handed them.
    """

    form: str
    params: tuple
    probes: str
    deriv: int = 0
    kdep: tuple = ()


WINDOWS: typing.Final[dict[str, Window]] = {
    # --- Gaussians, across the SRD's range of centres and widths -------------
    "gauss_mid": Window("gauss", (1500.0, 300.0, 4.0), "baseline"),
    "gauss_near": Window("gauss", (1095.0, 182.0, 4.0), "SRD Y10 lens, low end"),
    "gauss_far": Window("gauss", (3800.0, 165.0, 4.0), "SRD Y10 lens, high end; sigma/chi = 0.043"),
    "gauss_broad": Window("gauss", (4520.0, 715.0, 4.0), "source-bin width, distant"),
    "gauss_over": Window("gauss", (1800.0, 300.0, 4.0), "overlaps gauss_mid"),
    "gauss_low": Window("gauss", (600.0, 100.0, 4.0), "near bin, for far-separated pairs"),
    "gauss_high": Window("gauss", (3500.0, 150.0, 4.0), "far bin, for far-separated pairs"),
    "gauss_thin": Window("gauss", (1500.0, 50.0, 4.0), "thin bin"),
    "gauss_thin_shift": Window("gauss", (1650.0, 50.0, 4.0), "thin bin, displaced"),
    # --- Top-hats, paired with the Gaussians at equal variance ---------------
    "tophat_mid": Window("tophat", (500.0, 2500.0), "hard edges, wide"),
    "tophat_near": Window("tophat", (780.0, 1410.0), "pairs gauss_near at equal variance"),
    "tophat_far": Window("tophat", (3514.0, 4086.0), "pairs gauss_far at equal variance"),
    "tophat_broad": Window("tophat", (3282.0, 5758.0), "pairs gauss_broad at equal variance"),
    "shell_narrow": Window("tophat", (1000.0, 1056.0), "narrow hard shell"),
    "shell_wide": Window("tophat", (1100.0, 1500.0), "wide hard shell"),
    # --- The remaining forms -------------------------------------------------
    "tophat_smooth_mid": Window(
        "tophat_smooth", (1000.0, 2000.0, 150.0, 6.0), "smoothed edge"
    ),
    "student_t_mid": Window("student_t", (1500.0, 200.0, 2.0, 6.0), "power-law tail"),
    "power_exp_mid": Window(
        "power_exp", (1200.0, 2.0, 1.5, 50.0, 4000.0), "skewed and broad"
    ),
    "power_exp_skew_far": Window(
        "power_exp", (2400.0, 6.0, 2.0, 3500.0, 7200.0), "skewed at high chi, source-bin shape"
    ),
    "power_exp_isw": Window(
        "power_exp", (2500.0, 1.0, 1.5, 50.0, 8000.0), "broad and low-weighted, ISW range"
    ),
    "lensing_gal": Window("lensing", (50.0, 2000.0, 3000.0), "galaxy sources"),
    "lensing_cmb": Window(
        "lensing", (50.0, 13900.0, 14146.0), "sources at last scattering; support to 14 Gpc"
    ),
    "multi_mid": Window(
        "multi", ((1000.0, 1600.0), (300.0, 300.0), (1.0, 0.6), 4.0), "overlapping bumps"
    ),
    "multi_disjoint": Window(
        "multi", ((600.0, 2600.0), (100.0, 150.0), (1.0, 1.0), 4.0), "disconnected support"
    ),
}

# Derivative and scale-dependent variants, generated rather than retyped: the
# weight is orthogonal to the window, and writing them out by hand is how the two
# rosters drifted apart in the first place.
for _base, _deriv in (
    ("gauss_mid", 1),
    ("gauss_mid", 2),
    ("gauss_far", 2),
    ("tophat_mid", 2),
    ("tophat_smooth_mid", 2),
    ("student_t_mid", 2),
    ("power_exp_mid", 2),
    ("lensing_gal", 2),
    ("multi_disjoint", 2),
    ("shell_narrow", 2),
    ("shell_wide", 2),
):
    _w = WINDOWS[_base]
    WINDOWS[f"{_base}_d{_deriv}"] = _w._replace(
        probes=f"{_w.probes}; weighted by j_ell{chr(39) * _deriv}", deriv=_deriv
    )

WINDOWS["gauss_mid_kdep"] = WINDOWS["gauss_mid"]._replace(
    probes="non-separable W: sqrt(P) no longer factors out",
    kdep=(0.3, 0.05, 3000.0),
)

KDEP_FLAGS: typing.Final[tuple[str, ...]] = (
    "kdep-amplitude",
    "kdep-k-transition",
    "kdep-chi-ref",
)


# Names used before the vocabulary was unified, mapped to it. Consumers resolve
# through this so a committed table written under the old names still reads, and
# the data migration can happen whenever the generators are idle rather than
# being forced by a code change.
ALIASES: typing.Final[dict[str, str]] = {
    "gauss": "gauss_mid",
    "tophat": "tophat_mid",
    "tophat_smooth": "tophat_smooth_mid",
    "student_t": "student_t_mid",
    "power_exp": "power_exp_mid",
    "power_exp_srd4": "power_exp_skew_far",
    "power_exp_broad": "power_exp_isw",
    "lensing": "lensing_gal",
    "multi": "multi_mid",
    "srd_lens0": "gauss_near",
    "srd_lens9": "gauss_far",
    "srd_source4": "power_exp_skew_far",
    "srd_lens9_rsd": "gauss_far_d2",
    "gauss_rsd": "gauss_mid_d2",
    "gauss_deriv1": "gauss_mid_d1",
    "tophat_rsd": "tophat_mid_d2",
    "tophat_smooth_rsd": "tophat_smooth_mid_d2",
    "student_t_rsd": "student_t_mid_d2",
    "power_exp_rsd": "power_exp_mid_d2",
    "lensing_rsd": "lensing_gal_d2",
    "multi_disjoint_rsd": "multi_disjoint_d2",
    "shell_narrow_rsd": "shell_narrow_d2",
    "shell_wide_rsd": "shell_wide_d2",
    "kdep": "gauss_mid_kdep",
}

# The one genuine collision the unification had to resolve: `gauss_near` named
# (1800, 300) in the pair roster and (1095, 182) in the truth table. The pair
# roster's window is now `gauss_over`, which is what it always was.
COLLISIONS: typing.Final[dict[str, str]] = {"gauss_near@pairs": "gauss_over"}


def resolve(name: str) -> str:
    """The current name of a case, given a current or historical name.

    :param name: a case name, possibly one retired by the unification.
    :return: the current name.
    :raises KeyError: if the name is neither current nor a known alias.
    """
    if name in WINDOWS:
        return name

    if name in ALIASES:
        return ALIASES[name]

    raise KeyError(f"unknown window {name!r}")


def get(name: str) -> Window:
    """The window a name refers to, resolving historical names.

    :param name: a case name.
    :return: its :class:`Window`.
    """
    return WINDOWS[resolve(name)]


def arb_args(name: str, side: str | None = None) -> list[str]:
    """The generator arguments for a case.

    :param name: a case name.
    :param side: ``"a"`` or ``"b"`` for the two-window generator, or ``None`` for
        the single-window one, which takes bare options.
    :return: the command line options selecting this window.
    """
    window = get(name)
    prefix = "--" if side is None else f"--{side}:"
    out = [f"{prefix}shape={window.form}"]

    flags = FORM_ARB_FLAGS.get(
        window.form, tuple(k.replace("_", "-") for k in FORM_PARAMS[window.form])
    )

    for flag, value in zip(flags, window.params):

        if isinstance(value, tuple):
            value = ",".join(repr(float(v)) for v in value)

        out.append(f"{prefix}{flag}={value}")

    if window.deriv:
        out.append(f"{prefix}bessel-deriv={window.deriv}")

    for flag, value in zip(KDEP_FLAGS, window.kdep):
        out.append(f"{prefix}{flag}={value}")

    return out
