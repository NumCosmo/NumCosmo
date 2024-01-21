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
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm

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
