#!/usr/bin/env python
#
# test_py_prior.py
#
# Fri Jan 19 12:34:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_prior.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for NumCosmo Priors."""

from typing import List, Dict, Any
import re
import numpy as np
from numpy.testing import assert_allclose
import pytest

from numcosmo_py import Ncm, GLib
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


def test_py_prior_gauss() -> None:
    """Test NumCosmo Gaussian Priors on parameters."""
    pparams: List[Dict[str, Any]] = [
        {"model_ns": "NcmModelMVND", "parameter_name": "mu_0", "mu": 0.0, "sigma": 1.0},
        {"model_ns": "NcmModelMVND", "parameter_name": "mu_1", "mu": 2.0, "sigma": 3.0},
        {"model_ns": "NcmModelMVND", "parameter_name": "mu_2", "mu": 4.0, "sigma": 5.0},
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "nu",
            "mu": 0.0,
            "sigma": 1.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_0",
            "mu": 2.0,
            "sigma": 3.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_1",
            "mu": 4.0,
            "sigma": 5.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_2",
            "mu": 6.0,
            "sigma": 7.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_3",
            "mu": 8.0,
            "sigma": 9.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_4",
            "mu": 10.0,
            "sigma": 11.0,
        },
    ]
    priors: List[Ncm.PriorGaussParam] = []

    for pparam in pparams:
        priors.append(
            Ncm.PriorGaussParam(
                model_ns=pparam["model_ns"],
                parameter_name=pparam["parameter_name"],
                mu=pparam["mu"],
                sigma=pparam["sigma"],
            )
        )

    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(3), Ncm.ModelFunnel.new(5)])

    for prior, pparam in zip(priors, pparams):
        assert isinstance(prior, Ncm.Prior)
        assert isinstance(prior, Ncm.PriorGauss)
        assert isinstance(prior, Ncm.PriorGaussParam)
        assert prior.peek_model_ns() == pparam["model_ns"]
        assert prior.peek_param_name() == pparam["parameter_name"]
        assert prior.get_stack_pos() == 0
        assert prior.props.mu == pparam["mu"]
        assert prior.props.sigma == pparam["sigma"]
        assert prior.is_scalar()

        pi = mset.param_get_by_full_name(
            f"{pparam['model_ns']}:{pparam['parameter_name']}"
        )
        assert pi is not None
        mu = mset.param_get(pi.mid, pi.pid)
        prior_f = (mu - prior.props.mu) / prior.props.sigma

        assert_allclose(prior.eval0(mset), prior_f)


def test_py_prior_gauss_param_new_name() -> None:
    """Test NumCosmo Gaussian Priors on parameters."""
    prior = Ncm.PriorGaussParam.new_name(name="NcmModelMVND:mu_0", mu=0.0, sigma=1.0)
    assert isinstance(prior, Ncm.Prior)
    assert isinstance(prior, Ncm.PriorGauss)
    assert isinstance(prior, Ncm.PriorGaussParam)
    assert prior.peek_model_ns() == "NcmModelMVND"
    assert prior.peek_param_name() == "mu_0"
    assert prior.get_stack_pos() == 0
    assert prior.props.mu == 0.0
    assert prior.props.sigma == 1.0
    assert prior.is_scalar()


def test_py_prior_gauss_param_new_name_stackpos() -> None:
    """Test NumCosmo Gaussian Priors on parameters."""
    prior = Ncm.PriorGaussParam.new_name(name="NcmModelMVND:23:mu_0", mu=0.0, sigma=1.0)
    assert isinstance(prior, Ncm.Prior)
    assert isinstance(prior, Ncm.PriorGauss)
    assert isinstance(prior, Ncm.PriorGaussParam)
    assert prior.peek_model_ns() == "NcmModelMVND"
    assert prior.peek_param_name() == "mu_0"
    assert prior.get_stack_pos() == 23
    assert prior.props.mu == 0.0
    assert prior.props.sigma == 1.0
    assert prior.is_scalar()


def test_py_prior_gauss_param_new_name_invalid_stackpos() -> None:
    """Test NumCosmo Gaussian Priors on parameters."""
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_param_split_full_name: invalid stackpos "
            rf"number \(2301 \>\= 1000\). "
            rf"\({int(Ncm.MSetError.FULLNAME_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        _ = Ncm.PriorGaussParam.new_name(
            name="NcmModelMVND:2301:mu_0", mu=0.0, sigma=1.0
        )


def test_py_prior_gauss_param_new_name_invalid() -> None:
    """Test NumCosmo Gaussian Priors on parameters."""
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_prior_gauss_param_new_name: invalid parameter "
            rf"name `ncmModelMVND:23:mu_0'. "
            rf"\({int(Ncm.MSetError.FULLNAME_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        _ = Ncm.PriorGaussParam.new_name(name="ncmModelMVND:23:mu_0", mu=0.0, sigma=1.0)


def test_py_prior_flat() -> None:
    """Test NumCosmo Flat Priors on parameters."""
    pparams: List[Dict[str, Any]] = [
        {
            "model_ns": "NcmModelMVND",
            "parameter_name": "mu_0",
            "x_low": -1.0,
            "x_upp": 1.0,
        },
        {
            "model_ns": "NcmModelMVND",
            "parameter_name": "mu_1",
            "x_low": 0.0,
            "x_upp": 3.0,
        },
        {
            "model_ns": "NcmModelMVND",
            "parameter_name": "mu_2",
            "x_low": -0.0001,
            "x_upp": 5.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "nu",
            "x_low": -1.0,
            "x_upp": 1.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_0",
            "x_low": 0.0000001,
            "x_upp": 3.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_1",
            "x_low": 4.0,
            "x_upp": 5.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_2",
            "x_low": 6.0,
            "x_upp": 7.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_3",
            "x_low": -8.0,
            "x_upp": 9.0,
        },
        {
            "model_ns": "NcmModelFunnel",
            "parameter_name": "x_4",
            "x_low": 0.0,
            "x_upp": 11.0,
        },
    ]
    priors: List[Ncm.PriorFlatParam] = []

    for pparam in pparams:
        priors.append(
            Ncm.PriorFlatParam(
                model_ns=pparam["model_ns"],
                parameter_name=pparam["parameter_name"],
                x_low=pparam["x_low"],
                x_upp=pparam["x_upp"],
            )
        )

    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(3), Ncm.ModelFunnel.new(5)])

    for prior, pparam in zip(priors, pparams):
        assert isinstance(prior, Ncm.Prior)
        assert isinstance(prior, Ncm.PriorFlat)
        assert isinstance(prior, Ncm.PriorFlatParam)
        assert prior.peek_model_ns() == pparam["model_ns"]
        assert prior.peek_param_name() == pparam["parameter_name"]
        assert prior.get_stack_pos() == 0
        assert prior.props.x_low == pparam["x_low"]
        assert prior.props.x_upp == pparam["x_upp"]
        assert prior.is_scalar()

        pi = mset.param_get_by_full_name(
            f"{pparam['model_ns']}:{pparam['parameter_name']}"
        )
        assert pi is not None
        mu = mset.param_get(pi.mid, pi.pid)
        prior_f = 0.0 if (mu <= prior.props.x_low or mu >= prior.props.x_upp) else 1.0

        assert_allclose(np.exp(-0.5 * prior.eval0(mset) ** 2), prior_f)


def test_py_prior_flat_param_new_name() -> None:
    """Test NumCosmo Flat Priors on parameters."""
    prior = Ncm.PriorFlatParam.new_name(
        name="NcmModelMVND:mu_0", x_low=-1.0, x_upp=1.0, scale=1.0
    )
    assert isinstance(prior, Ncm.Prior)
    assert isinstance(prior, Ncm.PriorFlat)
    assert isinstance(prior, Ncm.PriorFlatParam)
    assert prior.peek_model_ns() == "NcmModelMVND"
    assert prior.peek_param_name() == "mu_0"
    assert prior.get_stack_pos() == 0
    assert prior.props.x_low == -1.0
    assert prior.props.x_upp == 1.0
    assert prior.is_scalar()


def test_py_prior_flat_param_new_name_stackpos() -> None:
    """Test NumCosmo Flat Priors on parameters."""
    prior = Ncm.PriorFlatParam.new_name(
        name="NcmModelMVND:23:mu_0", x_low=-1.0, x_upp=1.0, scale=1.0
    )
    assert isinstance(prior, Ncm.Prior)
    assert isinstance(prior, Ncm.PriorFlat)
    assert isinstance(prior, Ncm.PriorFlatParam)
    assert prior.peek_model_ns() == "NcmModelMVND"
    assert prior.peek_param_name() == "mu_0"
    assert prior.get_stack_pos() == 23
    assert prior.props.x_low == -1.0
    assert prior.props.x_upp == 1.0
    assert prior.is_scalar()


def test_py_prior_flat_param_new_name_invalid_stackpos() -> None:
    """Test NumCosmo Flat Priors on parameters."""
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_param_split_full_name: invalid stackpos "
            rf"number \(2300 \>\= 1000\). "
            rf"\({int(Ncm.MSetError.FULLNAME_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        _ = Ncm.PriorFlatParam.new_name(
            name="NcmModelMVND:2300:mu_0", x_low=-1.0, x_upp=1.0, scale=1.0
        )


def test_py_prior_flat_param_new_name_invalid() -> None:
    """Test NumCosmo Flat Priors on parameters."""
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_prior_flat_param_new_name: invalid parameter "
            rf"name `ncmModelMVND:23:mu_0'. "
            rf"\({int(Ncm.MSetError.FULLNAME_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        _ = Ncm.PriorFlatParam.new_name(
            name="ncmModelMVND:23:mu_0", x_low=-1.0, x_upp=1.0, scale=1.0
        )


@pytest.mark.parametrize("function_name", ["NcHICosmo:sigma8", "NcHICosmo:S8"])
def test_py_prior_mset_flist_psf(function_name):
    """Test NumCosmo priors on MSetFuncList."""
    cosmology = Cosmology.default()
    func = Ncm.MSetFuncList.new(function_name, cosmology.psf)
    assert isinstance(func, Ncm.MSetFuncList)
    assert isinstance(func, Ncm.MSetFunc)

    assert np.isfinite(func.eval0(cosmology.mset))

    mean = 0.8
    sd = 1.0
    var = 0.0

    prior = Ncm.PriorGaussFunc.new(func, mean, sd, var)
    assert isinstance(prior, Ncm.PriorGaussFunc)
    assert isinstance(prior, Ncm.PriorGauss)
    assert isinstance(prior, Ncm.Prior)

    prior_val = (func.eval0(cosmology.mset) - mean) / sd
    assert prior.eval0(cosmology.mset) == pytest.approx(prior_val, abs=1.0e-6)
