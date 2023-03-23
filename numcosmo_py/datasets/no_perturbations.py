#
# no_perturbations.py
#
# Sun Mar 5 11:35:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# no_perturbations.py
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


"""Factory functions to generate cosmology likelihoods that do
not depend on perturbations.
"""

from enum import Enum

from numcosmo_py import Nc
from numcosmo_py import Ncm
from numcosmo_py import GEnum


class SNIaID(GEnum):
    """Possible SNIa data sets ids."""

    # pylint: disable=no-member
    COV_PANTHEON_PLUS_SH0ES_SYS_STAT = Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT


class BAOID(str, Enum):
    """Possible BAO data sets ids."""

    # pylint: disable=no-member
    ALL_COMBINED_JAN_2023 = "ALL_COMBINED_JAN_2023"


class HID(str, Enum):
    """Possible Hubble data sets ids."""

    # pylint: disable=no-member
    ALL_COMBINED_JAN_2023 = "ALL_COMBINED_JAN_2023"


def gen_snia_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    dist: Nc.Distance,
    snia_id: SNIaID = SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
) -> None:
    """Generate a likelihood for SNIa data."""

    assert modelset.peek(Nc.HICosmo.id()) is not None

    snia_data = Nc.DataSNIACov.new_from_cat_id(snia_id.genum, False)
    if snia_id == SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT:
        snia0_data = snia_data.apply_filter_sh0es_z(0.01, True)
        snia_data = snia0_data

    snia_model = Nc.SNIADistCov.new_by_id(dist, snia_id.genum)
    if modelset.peek(snia_model.id()) is None:
        modelset.set(snia_model)

    dataset.append_data(snia_data)


def gen_bao_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    dist: Nc.Distance,
    bao_id: BAOID = BAOID.ALL_COMBINED_JAN_2023,
) -> None:
    """Generate a likelihood for BAO data."""

    assert modelset.peek(Nc.HICosmo.id()) is not None

    if bao_id == BAOID.ALL_COMBINED_JAN_2023:
        bao_ids = [
            Nc.DataBaoId.RDV_BEUTLER2011,
            Nc.DataBaoId.EMPIRICAL_FIT_ROSS2015,
            Nc.DataBaoId.DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE,
            Nc.DataBaoId.DTR_DHR_SDSS_DR16_LRG_2021,
            Nc.DataBaoId.DTR_DHR_SDSS_DR16_QSO_2021,
            Nc.DataBaoId.EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021,
            Nc.DataBaoId.EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021,
            Nc.DataBaoId.EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021,
        ]
        for bao_id in bao_ids:
            bao_likelihood = Nc.data_bao_create(dist, bao_id)
            dataset.append_data(bao_likelihood)
    else:
        raise ValueError(f"Unknown BAO data set id: {bao_id}")


def gen_h_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    h_id: HID = HID.ALL_COMBINED_JAN_2023,
) -> None:
    """Generate a likelihood for Hubble data."""

    assert modelset.peek(Nc.HICosmo.id()) is not None

    if h_id == HID.ALL_COMBINED_JAN_2023:
        h_ids = [Nc.DataHubbleId.GOMEZ_VALENT_COMP2018]
        for h_id in h_ids:
            h_likelihood = Nc.DataHubble.new_from_id(h_id)
            dataset.append_data(h_likelihood)
    else:
        raise ValueError(f"Unknown Hubble data set id: {h_id}")
