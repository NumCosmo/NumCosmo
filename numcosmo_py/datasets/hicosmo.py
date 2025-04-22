#
# hicosmo.py
#
# Sun Mar 5 11:35:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# hicosmo.py
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


"""Factory functions to create background cosmology likelihoods.

This module contains factory functions to create likelihoods for cosmology observables
that do not involve perturbations.
"""

from enum import StrEnum, auto

from numcosmo_py import Nc
from numcosmo_py import Ncm
from numcosmo_py import GEnum


class SNIaID(GEnum):
    """Possible SNIa data sets ids."""

    # pylint: disable=no-member
    COV_PANTHEON_PLUS_SH0ES_SYS_STAT = Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT
    SIMPLE_UNION2_1 = Nc.DataSNIAId.SIMPLE_UNION2_1
    COV_DES_Y5_STAT_SYS = Nc.DataSNIAId.COV_DES_Y5_STAT_SYS
    COV_DES_Y5_STATONLY = Nc.DataSNIAId.COV_DES_Y5_STATONLY


class BAOID(StrEnum):
    """Possible BAO data sets ids."""

    SDSS_ALL_COMBINED = auto()
    ALL_COMBINED_JAN_2023 = SDSS_ALL_COMBINED


class HID(StrEnum):
    """Possible Hubble data sets ids."""

    ALL_COMBINED_JAN_2023 = auto()
    ALL_COMBINED_APR_2025 = auto()


def add_snia_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    dist: Nc.Distance,
    snia_id: SNIaID = SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
) -> None:
    """Generate a likelihood for SNIa data."""
    cosmo = modelset.peek(Nc.HICosmo.id())
    assert cosmo is not None

    snia_data = Nc.DataSNIACov.new_from_cat_id(snia_id.genum, False)
    if snia_id == SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT:
        snia0_data = snia_data.apply_filter_sh0es_z(0.01, True)
        snia_data = snia0_data

    snia_model = Nc.SNIADistCov.new_by_id(dist, snia_id.genum)
    if modelset.peek(snia_model.id()) is None:
        modelset.set(snia_model)

    dataset.append_data(snia_data)


def add_bao_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    dist: Nc.Distance,
    bao_id: BAOID = BAOID.ALL_COMBINED_JAN_2023,
) -> None:
    """Generate a likelihood for BAO data."""
    assert modelset.peek(Nc.HICosmo.id()) is not None

    match bao_id:
        case BAOID.ALL_COMBINED_JAN_2023:
            bao_enums = [
                Nc.DataBaoId.RDV_BEUTLER2011,
                Nc.DataBaoId.EMPIRICAL_FIT_ROSS2015,
                Nc.DataBaoId.DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE,
                Nc.DataBaoId.DTR_DHR_SDSS_DR16_LRG_2021,
                Nc.DataBaoId.DTR_DHR_SDSS_DR16_QSO_2021,
                Nc.DataBaoId.EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021,
                Nc.DataBaoId.EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021,
                Nc.DataBaoId.EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021,
            ]
            for bao_enum in bao_enums:
                bao_likelihood = Nc.data_bao_create(dist, bao_enum)
                dataset.append_data(bao_likelihood)
        case _:
            raise ValueError(f"Unknown BAO data set id: {bao_id}")


def add_h_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    h_id: HID = HID.ALL_COMBINED_JAN_2023,
) -> None:
    """Generate a likelihood for Hubble data."""
    assert modelset.peek(Nc.HICosmo.id()) is not None

    match h_id:
        case HID.ALL_COMBINED_JAN_2023:
            h_enums = [Nc.DataHubbleId.GOMEZ_VALENT_COMP2018]
            for h_enum in h_enums:
                h_likelihood = Nc.DataHubble.new_from_id(h_enum)
                dataset.append_data(h_likelihood)
        case HID.ALL_COMBINED_APR_2025:
            h_enums = [
                Nc.DataHubbleId.GOMEZ_VALENT_COMP2018,
                Nc.DataHubbleId.TOMASETTI2023,
                Nc.DataHubbleId.JIMENEZ2023,
            ]
        case _:
            raise ValueError(f"Unknown Hubble data set id: {h_id}")
