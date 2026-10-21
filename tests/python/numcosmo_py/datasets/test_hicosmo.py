"""Tests for numcosmo_py.datasets.hicosmo."""

import math

import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.datasets.hicosmo import (
    BAOID,
    HID,
    OMEGA_B_H2_PRIORS,
    OmegaBH2PriorID,
    SNIaID,
    add_bao_likelihood,
    add_h_likelihood,
    add_omega_b_h2_prior,
    add_snia_likelihood,
)

Ncm.cfg_init()


@pytest.mark.parametrize(
    "h_id, n_diag, n_cov",
    [
        (HID.ALL_COMBINED_JAN_2023, 1, 0),
        (HID.ALL_COMBINED_APR_2025, 3, 0),
        (HID.ALL_COMBINED_OCT_2026, 3, 1),
        (HID.CC_MORESCO_2020, 0, 1),
    ],
)
def test_add_h_likelihood(h_id, n_diag, n_cov):
    """Each H(z) set appends its diagonal and full-covariance data objects."""
    cosmo = Nc.HICosmoDEXcdm.new()
    mset = Ncm.MSet.new_array([cosmo])
    dset = Ncm.Dataset.new()

    add_h_likelihood(dset, mset, h_id)

    data = [dset.get_data(i) for i in range(dset.get_length())]
    assert sum(isinstance(d, Nc.DataHubble) for d in data) == n_diag
    assert sum(isinstance(d, Nc.DataHubbleCov) for d in data) == n_cov
    assert len(data) == n_diag + n_cov

    lh = Ncm.Likelihood.new(dset)
    assert math.isfinite(lh.m2lnL_val(mset))


def test_add_bao_likelihood_desi_dr2():
    """The DESI DR2 set holds the BGS point and the table IV block only."""
    cosmo = Nc.HICosmoDEXcdm.new()
    mset = Ncm.MSet.new_array([cosmo])
    dset = Ncm.Dataset.new()

    add_bao_likelihood(dset, mset, Nc.Distance.new(3.0), BAOID.DESI_DR2_2025)

    data = [dset.get_data(i) for i in range(dset.get_length())]
    assert [type(d) for d in data] == [Nc.DataBaoRDV, Nc.DataBaoDvrDtDh]
    assert sum(d.get_length() for d in data) == 13


@pytest.mark.parametrize("sh0es_calibrators, n_snia", [(True, 1657), (False, 1590)])
def test_add_snia_likelihood_sh0es_calibrators(sh0es_calibrators, n_snia):
    """Pantheon+ keeps z_HD > 0.01, plus the SH0ES calibrators only on request."""
    cosmo = Nc.HICosmoDEXcdm.new()
    mset = Ncm.MSet.new_array([cosmo])
    dset = Ncm.Dataset.new()

    add_snia_likelihood(
        dset,
        mset,
        Nc.Distance.new(3.0),
        SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
        sh0es_calibrators,
    )

    assert dset.get_data(0).get_length() == n_snia


def test_add_omega_b_h2_prior():
    """The BBN prior frees Omega_b and adds -2 ln of a Gaussian in Omega_b h^2."""
    cosmo = Nc.HICosmoDEXcdm.new()
    mset = Ncm.MSet.new_array([cosmo])
    dset = Ncm.Dataset.new()
    add_bao_likelihood(dset, mset, Nc.Distance.new(3.0), BAOID.DESI_DR2_2025)
    likelihood = Ncm.Likelihood.new(dset)

    m2lnL_data = likelihood.m2lnL_val(mset)
    add_omega_b_h2_prior(likelihood, mset, OmegaBH2PriorID.BBN_DESI_DR2_2025)
    mset.prepare_fparam_map()

    mean, sigma = OMEGA_B_H2_PRIORS[OmegaBH2PriorID.BBN_DESI_DR2_2025]
    assert mset.fparam_get_pi_by_name("NcHICosmo:Omegab") is not None
    i = cosmo.param_index_from_name("Omegab")[1]
    assert (cosmo.param_get_lower_bound(i), cosmo.param_get_upper_bound(i)) == (
        0.01,
        0.15,
    )
    assert likelihood.m2lnL_val(mset) == pytest.approx(
        m2lnL_data + ((cosmo.Omega_b0h2() - mean) / sigma) ** 2, rel=1.0e-12
    )
