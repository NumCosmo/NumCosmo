"""Tests for numcosmo_py.datasets.hicosmo."""

import math

import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.datasets.hicosmo import HID, add_h_likelihood

Ncm.cfg_init()


@pytest.mark.parametrize(
    "h_id, n_diag, n_cov",
    [
        (HID.ALL_COMBINED_JAN_2023, 1, 0),
        (HID.ALL_COMBINED_APR_2025, 3, 0),
        (HID.ALL_COMBINED_OCT_2026, 3, 1),
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
