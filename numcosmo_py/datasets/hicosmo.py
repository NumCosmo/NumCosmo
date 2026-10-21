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

from numcosmo_py import GEnum, Nc, Ncm


class SNIaID(GEnum):
    """Possible SNIa data sets ids."""

    COV_PANTHEON_PLUS_SH0ES_SYS_STAT = Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT
    COV_PANTHEON_PLUS_SH0ES_HOYT26_SYS_STAT = (
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_HOYT26_SYS_STAT
    )
    SIMPLE_UNION2_1 = Nc.DataSNIAId.SIMPLE_UNION2_1
    COV_DES_Y5_STAT_SYS = Nc.DataSNIAId.COV_DES_Y5_STAT_SYS
    COV_DES_Y5_STATONLY = Nc.DataSNIAId.COV_DES_Y5_STATONLY
    COV_DES_DOVEKIE_STAT_SYS = Nc.DataSNIAId.COV_DES_DOVEKIE_STAT_SYS
    COV_DES_DOVEKIE_STATONLY = Nc.DataSNIAId.COV_DES_DOVEKIE_STATONLY


class OmegaBH2PriorID(StrEnum):
    """Gaussian priors on Omega_b h^2."""

    BBN_DESI_DR2_2025 = auto()


# Mean and standard deviation of each Omega_b h^2 prior.
OMEGA_B_H2_PRIORS: dict[OmegaBH2PriorID, tuple[float, float]] = {
    # DESI DR2 (arXiv:2503.14738, eq. 14): BBN with PRyMordial, marginalized over the
    # nuclear reaction rates.
    OmegaBH2PriorID.BBN_DESI_DR2_2025: (0.02218, 0.00055),
}


class BAOID(StrEnum):
    """Possible BAO data sets ids."""

    SDSS_ALL_COMBINED = auto()
    ALL_COMBINED_JAN_2023 = SDSS_ALL_COMBINED
    ALL_COMBINED_JUN_2025 = auto()
    DESI_DR2_2025 = auto()


class HID(StrEnum):
    """Possible Hubble data sets ids."""

    ALL_COMBINED_JAN_2023 = auto()
    ALL_COMBINED_APR_2025 = auto()
    ALL_COMBINED_OCT_2026 = auto()
    CC_MORESCO_2020 = auto()


def add_snia_likelihood(
    dataset: Ncm.Dataset,
    modelset: Ncm.MSet,
    dist: Nc.Distance,
    snia_id: SNIaID = SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
    sh0es_calibrators: bool = True,
) -> None:
    """Generate a likelihood for SNIa data.

    The Pantheon+SH0ES STAT+SYS catalogs keep z_HD > 0.01 and, when
    ``sh0es_calibrators`` is set, the SH0ES Cepheid calibrators; without them the
    selection is the one of the Pantheon+ SN-only likelihood.
    """
    cosmo = modelset.peek(Nc.HICosmo.id())
    assert cosmo is not None

    snia_data = Nc.DataSNIACov.new_from_cat_id(snia_id.genum, False)
    if snia_id in (
        SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
        SNIaID.COV_PANTHEON_PLUS_SH0ES_HOYT26_SYS_STAT,
    ):
        snia0_data = snia_data.apply_filter_sh0es_z(0.01, sh0es_calibrators)
        snia_data = snia0_data

    snia_model = Nc.SNIADistCov.new_by_id(dist, snia_id.genum)
    if modelset.peek(snia_model.id()) is None:
        modelset.set(snia_model)

    dataset.append_data(snia_data)


def add_omega_b_h2_prior(
    likelihood: Ncm.Likelihood,
    modelset: Ncm.MSet,
    prior_id: OmegaBH2PriorID = OmegaBH2PriorID.BBN_DESI_DR2_2025,
) -> None:
    """Free Omega_b and add a Gaussian prior on Omega_b h^2.

    The prior constrains Omega_b h^2, so Omega_b follows H0: its bounds are widened to
    [0.01, 0.15], which do not bind for 0.5 < h < 1 (the model default, [0.03, 0.05], does at
    h < 0.67). Call it before ``modelset.prepare_fparam_map()``, since it frees ``Omegab``.
    """
    cosmo = modelset.peek(Nc.HICosmo.id())
    assert cosmo is not None
    cosmo.param_set_desc(
        "Omegab", {"fit": True, "lower-bound": 0.01, "upper-bound": 0.15}
    )

    mean, sigma = OMEGA_B_H2_PRIORS[prior_id]
    func = Ncm.MSetFuncList.new("NcHICosmo:Omega_b0h2", None)
    likelihood.priors_add(Ncm.PriorGaussFunc.new(func, mean, sigma, 0.0))


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
                Nc.DataBaoId.RDV_KAZIN2014,
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
        case BAOID.ALL_COMBINED_JUN_2025:
            bao_enums = [
                Nc.DataBaoId.RDV_BEUTLER2011,
                Nc.DataBaoId.RDV_KAZIN2014,
                Nc.DataBaoId.EMPIRICAL_FIT_ROSS2015,
                Nc.DataBaoId.RDV_DESI_DR2_BGS_2025,
                Nc.DataBaoId.DVR_DTDH_DESI_DR2_2025,
            ]
            for bao_enum in bao_enums:
                bao_likelihood = Nc.data_bao_create(dist, bao_enum)
                dataset.append_data(bao_likelihood)
        case BAOID.DESI_DR2_2025:
            bao_enums = [
                Nc.DataBaoId.RDV_DESI_DR2_BGS_2025,
                Nc.DataBaoId.DVR_DTDH_DESI_DR2_2025,
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

    h_cov_enums: list[Nc.DataHubbleCovId] = []
    match h_id:
        case HID.ALL_COMBINED_JAN_2023:
            h_enums = [Nc.DataHubbleId.GOMEZ_VALENT_COMP2018]
        case HID.ALL_COMBINED_APR_2025:
            h_enums = [
                Nc.DataHubbleId.GOMEZ_VALENT_COMP2018,
                Nc.DataHubbleId.TOMASETTI2023,
                Nc.DataHubbleId.JIMENEZ2023,
            ]
        case HID.ALL_COMBINED_OCT_2026:
            # APR_2025 plus Loubser 2025 (DESI DR1), whose three points share a
            # systematic error and carry a full covariance.
            h_enums = [
                Nc.DataHubbleId.GOMEZ_VALENT_COMP2018,
                Nc.DataHubbleId.TOMASETTI2023,
                Nc.DataHubbleId.JIMENEZ2023,
            ]
            h_cov_enums = [Nc.DataHubbleCovId.LOUBSER2025]
        case HID.CC_MORESCO_2020:
            # The D4000n sample with the Moresco et al. (2020) covariance, as published.
            h_enums = []
            h_cov_enums = [Nc.DataHubbleCovId.MORESCO2020]
        case _:
            raise ValueError(f"Unknown Hubble data set id: {h_id}")

    for h_enum in h_enums:
        h_likelihood = Nc.DataHubble.new_from_id(h_enum)
        dataset.append_data(h_likelihood)

    for h_cov_enum in h_cov_enums:
        dataset.append_data(Nc.DataHubbleCov.new_from_id(h_cov_enum))
