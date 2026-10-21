"""Tests for numcosmo_py.datasets.hicosmo."""

import math

import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.datasets.hicosmo import (
    BAOID,
    HID,
    SNIaID,
    add_bao_likelihood,
    add_h_likelihood,
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
