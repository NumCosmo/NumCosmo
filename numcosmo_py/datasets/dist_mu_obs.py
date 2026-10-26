#
# dist_mu_obs.py
#
# Sun Oct 18 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# dist_mu_obs.py
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

"""Catalogs of distance moduli with sky positions (NcDistMuObs).

Builds :class:`Nc.DistMuObs` catalogs from external tables and names the catalogs
shipped with NumCosmo.
"""

from enum import StrEnum, auto
from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy.table import Table

from numcosmo_py import GEnum, Nc, Ncm

SPEED_OF_LIGHT_KM_S = 299792.458

EDD_CF4_SNIA_URL = "https://edd.ifa.hawaii.edu/dfirst.php"
EDD_CF4_SNIA_TABLE = "All CF4 SNIa Samples (kcf4snia): Cosmicflows-4 Table 2, column 6"
EDD_CF4_SNIA_REFERENCE = "Tully et al. 2023, ApJ 944, 94"
EDD_CF4_SNIA_RENAMED = ("RA", "DE", "DM_zp", "eDM")


class DistMuObsCatalogID(GEnum):
    """Catalogs of the NumCosmo data release."""

    # pylint: disable=no-member
    CF4_SNIA = Nc.DistMuObsCatalogId.CF4_SNIA


def dist_mu_obs_from_table(
    table: Table, meta: dict[str, str] | None = None
) -> Nc.DistMuObs:
    """Build a :class:`Nc.DistMuObs` from an astropy table.

    The table must have the columns ``ra``, ``dec``, ``z``, ``mu`` and ``sigma_mu``;
    every other numeric column is kept, with its integer or boolean type, and string
    columns are dropped. Empty cells of a float column become NaN.

    :param table: the table.
    :param meta: string metadata stored in the catalog, see
        :meth:`Ncm.Catalog.set_meta`.
    :return: the catalog.
    """
    columns: list[str] = []
    col_types: list[Ncm.CatalogColType] = []
    arrays: list[np.ndarray] = []

    for col in table.colnames:
        column = table[col]
        if np.issubdtype(column.dtype, np.bool_):
            col_types.append(Ncm.CatalogColType.BOOL)
        elif np.issubdtype(column.dtype, np.integer):
            col_types.append(Ncm.CatalogColType.INT)
        elif np.issubdtype(column.dtype, np.floating):
            col_types.append(Ncm.CatalogColType.DOUBLE)
        else:
            continue
        columns.append(col)
        values = np.asarray(column, dtype=np.float64)
        if hasattr(column, "mask"):
            values = np.where(np.asarray(column.mask), np.nan, values)
        arrays.append(values)

    nrows = len(table)
    obs = Nc.DistMuObs.new_full(nrows, columns, col_types)
    data = np.column_stack(arrays)
    obs.peek_data().set_from_array(data.reshape(-1).tolist())

    if meta is not None:
        var_dict = Ncm.VarDict.new()
        for key, value in meta.items():
            var_dict.set_string(key, value)
        obs.set_meta(var_dict)

    return obs


def edd_cf4_snia_to_dist_mu_obs(path: Path | str, downloaded: str) -> Nc.DistMuObs:
    """Build the CF4 SNIa catalog from an EDD text export of All CF4 SNIa Samples.

    The export is the ``|``-separated text of the EDD query with every column, whose
    first five lines are the table title, the column names, the C formats, the units
    and the column descriptions. The catalog columns are ``ra`` and ``dec`` from
    ``RA`` and ``DE``, ``z`` as ``Vcmb`` divided by the speed of light, ``mu`` as
    ``DM_zp`` and ``sigma_mu`` as ``eDM``; every other numeric column of the export is
    kept under its EDD name, so ``DM_av`` and the per-sample moduli are available as
    alternatives. ``eDM`` is the smallest ``eDM_X`` of the object, and the ``eDM_X``
    columns differ in whether they contain a peculiar velocity term.

    :param path: the export file.
    :param downloaded: the download date, stored in the metadata.
    :return: the catalog.
    """
    table = ascii.read(
        str(path),
        format="basic",
        delimiter="|",
        header_start=1,
        data_start=5,
        fill_values=[("", "nan")],
    )
    out = Table()
    out["ra"] = np.asarray(table["RA"], dtype=np.float64)
    out["dec"] = np.asarray(table["DE"], dtype=np.float64)
    out["z"] = np.asarray(table["Vcmb"], dtype=np.float64) / SPEED_OF_LIGHT_KM_S
    out["mu"] = np.asarray(table["DM_zp"], dtype=np.float64)
    out["sigma_mu"] = np.asarray(table["eDM"], dtype=np.float64)
    for col in table.colnames:
        if col not in EDD_CF4_SNIA_RENAMED and col not in out.colnames:
            out[col] = table[col]

    meta = {
        "source": EDD_CF4_SNIA_TABLE,
        "url": EDD_CF4_SNIA_URL,
        "reference": EDD_CF4_SNIA_REFERENCE,
        "downloaded": downloaded,
        "ra": "RA, right ascension of the host, J2000, degrees",
        "dec": "DE, declination of the host, J2000, degrees",
        "z": "Vcmb / 299792.458, Vcmb the radial velocity (cz) in the CMB frame",
        "mu": "DM_zp, weighted average distance modulus on the absolute scale set by the calibrators",
        "sigma_mu": "eDM, minimum error of all measured moduli",
    }
    return dist_mu_obs_from_table(out, meta)


class SampleRule(StrEnum):
    """Rules assigning a sample index to each object of a catalog."""

    NONE = auto()
    SMALLEST_ERROR = auto()
    HEMISPHERE = auto()


def survey_names(obs: Nc.DistMuObs) -> list[str]:
    """Surveys of a catalog: the ``X`` of every ``DM_X``/``eDM_X`` column pair.

    :param obs: the catalog.
    :return: the survey names, sorted.
    """
    columns = list(obs.peek_columns())
    return sorted(
        c[3:] for c in columns if c.startswith("DM_") and f"eDM_{c[3:]}" in columns
    )


def samples_by_rule(
    obs: Nc.DistMuObs, rule: SampleRule, selected: np.ndarray | None = None
) -> tuple[list[float], list[str]]:
    """Sample index of each object by a rule.

    ``SMALLEST_ERROR`` assigns each object to the survey with the smallest error among
    those that measured it, indexed into :func:`survey_names`, -1 when none did. The
    group of a survey holds only the objects where its error is the smallest, so it is
    not the sample of that survey, and offsets fitted on these groups are differences
    between subsets of the catalog, not zero points of the surveys;
    ``HEMISPHERE`` assigns 0 to dec >= 0 and 1 to dec < 0; ``NONE`` assigns -1.

    :param obs: the catalog.
    :param rule: the rule.
    :param selected: row indices to label, all rows when ``None``.
    :return: the sample index of each selected row and the sample names.
    """
    rows = np.arange(obs.len()) if selected is None else np.asarray(selected)
    columns = list(obs.peek_columns())
    data = obs.peek_data().to_numpy()
    if rule == SampleRule.NONE:
        return [-1.0] * len(rows), []
    if rule == SampleRule.HEMISPHERE:
        dec = data[rows, columns.index("dec")]
        return [0.0 if d >= 0.0 else 1.0 for d in dec], ["north", "south"]
    names = survey_names(obs)
    errors = np.column_stack([data[rows, columns.index(f"eDM_{n}")] for n in names])
    errors = np.where(np.isfinite(errors), errors, np.inf)
    best = np.argmin(errors, axis=1)
    none = ~np.isfinite(errors).any(axis=1)
    return [-1.0 if n else float(b) for b, n in zip(best, none)], names
