#!/usr/bin/env python
#
# make_mset_presets.py
#
# Thu Nov 20 12:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_mset_presets.py
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

"""Write the NcmMSetPreset files of data/mset_presets from published chains.

Each preset sets, for the parameters NumCosmo samples, the best fit as value, the
marginalized standard deviation as scale and the fit flag; and, for the parameters
the publication fixes, the value with the fit flag off. No preset sets bounds: the
model bounds give the range where the code computes, and each analysis narrows them.

Inputs, downloaded once:

* Planck 2018: the directory ``base`` extracted from the PLA archive
  ``COM_CosmoParams_base-plikHM_R3.01.zip``, with the ``.minimum`` (best fit) and
  ``dist/*.margestats`` (marginalized statistics) of each data combination.
* DESI DR2: the directory with the four cobaya chains ``cobaya_chain.[1-4].txt`` and
  the iminuit best fit ``iminuit_bestfit.minimum`` of the LCDM DESI+CMB run, from
  https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params.

Usage::

    python tools/make_mset_presets.py PLANCK_BASE DESI_DIR data/mset_presets
"""

import argparse
from dataclasses import dataclass, field
from pathlib import Path

PLANCK_REF = (
    "Planck 2018 results VI, A&A 641, A6 (2020), arXiv:1807.06209; "
    "PLA chains COM_CosmoParams_base-plikHM_R3.01"
)
DESI_REF = (
    "DESI DR2 results II, Phys. Rev. D 112, 083515 (2025), arXiv:2503.14738; "
    "chains data.desi.lbl.gov/public/papers/y3/bao-cosmo-params"
)
CMB_REPARAM = {"NcHICosmo": "NcHICosmoDEReparamCMB"}

# Mass contribution to N_eff of one massive neutrino at NumCosmo's default temperature
# T_nu = 0.71611 T_gamma (NC_HICOSMO_DE_DEFAULT_NU_T, as CLASS): (T_nu / (4/11)^(1/3))^4.
MASSIVE_NU_NEFF = 1.0132

# NumCosmo full name and the source name of each parameter NumCosmo samples. H0 and
# z_re are derived parameters in the sources (which sample theta_MC and tau).
COSMO_PARAMS = [
    ("NcHICosmo:H0", "H0"),
    ("NcHICosmo:omegac", "omegach2"),
    ("NcHICosmo:omegab", "omegabh2"),
    ("NcHIPrim:ln10e10ASA", "logA"),
    ("NcHIPrim:n_SA", "ns"),
    ("NcHIReion:z_re", "zrei"),
]
DESI_NAMES = {"omegach2": "omch2", "omegabh2": "ombh2"}

PLANCK_TT_NUISANCE = [
    ("NcPlanckFI:A_planck", "calPlanck"),
    ("NcPlanckFI:A_cib_217", "acib217"),
    ("NcPlanckFI:xi_sz_cib", "xi"),
    ("NcPlanckFI:A_sz", "asz143"),
    ("NcPlanckFI:ps_A_100_100", "aps100"),
    ("NcPlanckFI:ps_A_143_143", "aps143"),
    ("NcPlanckFI:ps_A_143_217", "aps143217"),
    ("NcPlanckFI:ps_A_217_217", "aps217"),
    ("NcPlanckFI:ksz_norm", "aksz"),
    ("NcPlanckFI:gal545_A_100", "kgal100"),
    ("NcPlanckFI:gal545_A_143", "kgal143"),
    ("NcPlanckFI:gal545_A_143_217", "kgal143217"),
    ("NcPlanckFI:gal545_A_217", "kgal217"),
    ("NcPlanckFI:calib_100T", "cal0"),
    ("NcPlanckFI:calib_217T", "cal2"),
]
PLANCK_TE_NUISANCE = [
    ("NcPlanckFI:galf_TE_A_100", "galfTE100"),
    ("NcPlanckFI:galf_TE_A_100_143", "galfTE100143"),
    ("NcPlanckFI:galf_TE_A_100_217", "galfTE100217"),
    ("NcPlanckFI:galf_TE_A_143", "galfTE143"),
    ("NcPlanckFI:galf_TE_A_143_217", "galfTE143217"),
    ("NcPlanckFI:galf_TE_A_217", "galfTE217"),
]


@dataclass
class Entry:
    """One parameter of a preset."""

    fullname: str
    value: float
    fit: bool
    scale: float | None = None
    comment: str = ""


@dataclass
class Preset:
    """One preset file."""

    name: str
    desc: str
    data_combination: str
    header: list[str]
    reparam: dict[str, str] = field(default_factory=dict)
    entries: list[Entry] = field(default_factory=list)


def fixed_lcdm(n_eff: float, source: str, t_cmb_source: str) -> list[Entry]:
    """Return the parameters a flat LCDM analysis with one 0.06 eV neutrino fixes.

    @source names the file stating Omega_k, w, the neutrino mass and N_eff, and
    @t_cmb_source the one stating T_CMB.
    """
    n_ur = round(n_eff - MASSIVE_NU_NEFF, 4)
    return [
        Entry("NcHICosmo:Omegak", 0.0, False, comment=f"omegak, {source}"),
        Entry("NcHICosmo:w", -1.0, False, comment=f"w, {source}"),
        Entry("NcHICosmo:Tgamma0", 2.7255, False, comment=f"T_CMB [K], {t_cmb_source}"),
        Entry("NcHICosmo:massnu_0", 0.06, False, comment=f"mnu [eV], {source}"),
        Entry(
            "NcHICosmo:ENnu",
            n_ur,
            False,
            comment=f"nnu = {n_eff} minus {MASSIVE_NU_NEFF} of the massive nu, {source}",
        ),
    ]


def read_minimum(path: Path) -> dict[str, float]:
    """Read a GetDist .minimum file: source name to best-fit value."""
    out = {}
    for line in path.read_text().splitlines():
        tok = line.split()
        if len(tok) >= 3 and tok[0].isdigit():
            out[tok[2].rstrip("*")] = float(tok[1])
    return out


def read_margestats(path: Path) -> dict[str, float]:
    """Read a GetDist .margestats file: source name to standard deviation."""
    out = {}
    for line in path.read_text().splitlines()[3:]:
        tok = line.split()
        if len(tok) > 5:
            out[tok[0].rstrip("*")] = float(tok[2])
    return out


def sampled(
    params: list[tuple[str, str]],
    best: dict[str, float],
    sd: dict[str, float],
    names: dict[str, str] | None = None,
) -> list[Entry]:
    """Return the entries of sampled parameters: best fit, standard deviation, free."""
    names = names or {}
    out = []
    for fullname, src in params:
        key = names.get(src, src)
        out.append(Entry(fullname, best[key], True, sd[key], comment=key))
    return out


def planck_presets(base: Path) -> list[Preset]:
    """Return the cosmology and nuisance presets of the three Planck combinations."""
    combos = [
        ("tt_lowl_lowe", "plikHM_TT_lowl_lowE", "TT+lowE", "2.1", False),
        ("ttteee_lowl_lowe", "plikHM_TTTEEE_lowl_lowE", "TT,TE,EE+lowE", "2.5", True),
        (
            "ttteee_lowl_lowe_lensing",
            "plikHM_TTTEEE_lowl_lowE_lensing",
            "TT,TE,EE+lowE+lensing",
            "2.17",
            True,
        ),
    ]
    presets = []
    for tag, chain, label, section, te in combos:
        best = read_minimum(base / chain / f"base_{chain}.minimum")
        sd = read_margestats(base / chain / "dist" / f"base_{chain}.margestats")
        header = [
            f"Planck 2018 {label}, plik high-ell likelihood, base LCDM.",
            f"value: best fit of base_{chain}.minimum; scale: sddev of its .margestats.",
            f"Checked against the PLA 68% tables (2018, v2), section {section}, and",
            "Planck 2018 VI Tables 1 and 2.",
        ]
        cosmo = Preset(
            f"planck2018_{tag}",
            f"{PLANCK_REF}, base_{chain}",
            f"Planck 2018 {label}",
            header,
            dict(CMB_REPARAM),
        )
        cosmo.entries = sampled(COSMO_PARAMS, best, sd) + fixed_lcdm(
            3.046, f"base_{chain}.ranges", "Planck 2018 VI (Fixsen 2009)"
        )
        nuisance_params = PLANCK_TT_NUISANCE + (PLANCK_TE_NUISANCE if te else [])
        nuisance = Preset(
            f"planck2018_{tag}_nuisance",
            f"{PLANCK_REF}, base_{chain}",
            f"Planck 2018 {label}",
            header + ["Foreground and calibration parameters of the plik likelihood."],
        )
        nuisance.entries = sampled(nuisance_params, best, sd)
        presets += [cosmo, nuisance]
    return presets


def desi_preset(desi_dir: Path) -> Preset:
    """Return the LCDM DESI DR2 BAO+CMB cosmology preset."""
    from getdist import loadMCSamples

    best = read_minimum(desi_dir / "iminuit_bestfit.minimum")
    samples = loadMCSamples(
        str(desi_dir / "cobaya_chain"), settings={"ignore_rows": 0.3}
    )
    stats = samples.getMargeStats()
    sd = {}
    for _, src in COSMO_PARAMS:
        key = DESI_NAMES.get(src, src)
        sd[key] = float(stats.parWithName(key).err)

    header = [
        "DESI DR2 BAO + CMB (Planck PR3 low-ell TT and EE, Planck PR4 CamSpec TTTEEE,",
        "Planck PR4 + ACT DR6 lensing), base LCDM. Cosmological parameters only: the",
        "CamSpec nuisance parameters have no NumCosmo likelihood.",
        "value: best fit (posterior maximum) of iminuit/.../bestfit.minimum; scale:",
        "standard deviation of the cobaya chains (GetDist, ignore_rows 0.3).",
        "Checked against the paper's Table tab:cosmo_constraints and chain.margestats.",
        "N_eff = 3.044 (DESI), not the Planck 2018 3.046; T_CMB is the CAMB default.",
    ]
    preset = Preset(
        "desi_dr2_bao_cmb",
        f"{DESI_REF}, base/desi-bao-all_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik"
        "_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing",
        "DESI DR2 BAO + CMB",
        header,
        dict(CMB_REPARAM),
    )
    preset.entries = sampled(COSMO_PARAMS, best, sd, DESI_NAMES) + fixed_lcdm(
        3.044, "cobaya chain.updated.yaml", "CAMB default"
    )
    return preset


def yaml_float(x: float) -> str:
    """Return the shortest round-trip form of @x, always with a decimal point."""
    s = repr(float(x))
    if "." not in s and "inf" not in s and "nan" not in s:
        mant, _, exp = s.partition("e")
        s = f"{mant}.0e{exp}" if exp else f"{mant}.0"
    return s


def render(preset: Preset) -> str:
    """Return the YAML text of @preset, in the form of ncm_serialize_to_yaml()."""
    lines = [f"# {h}" for h in preset.header]
    lines += [
        "NcmMSetPreset:",
        f"  desc: '{preset.desc}'",
        f"  data-combination: '{preset.data_combination}'",
        "  statistic: 'best fit'",
    ]
    if preset.reparam:
        lines.append("  reparam:")
        lines += [f"    {k}: '{v}'" for k, v in preset.reparam.items()]
    lines.append("  value:")
    lines += [
        f"    {e.fullname}: {yaml_float(e.value)}  # {e.comment}"
        for e in preset.entries
    ]
    lines.append("  scale:")
    lines += [
        f"    {e.fullname}: {yaml_float(e.scale)}"
        for e in preset.entries
        if e.scale is not None
    ]
    lines.append("  fit:")
    lines += [
        f"    {e.fullname}: {'true' if e.fit else 'false'}" for e in preset.entries
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    """Write every preset file."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("planck_base", type=Path)
    parser.add_argument("desi_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    presets = planck_presets(args.planck_base) + [desi_preset(args.desi_dir)]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for preset in presets:
        path = args.output_dir / f"{preset.name}.yaml"
        path.write_text(render(preset))
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
