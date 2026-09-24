#!/usr/bin/env python
#
# test_sampling_helpers.py
#
# Sun Sep 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_sampling_helpers.py
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

"""Unit tests for the sampler builders in numcosmo_py.sampling.

These are the entry points a script or a notebook uses instead of the command
line: create_esmcmc() assembles an Ncm.FitESMCMC from a likelihood, and APES wraps
the same sampler behind an emcee-shaped interface. Both forward a list of proposal
settings to the walker, which is what is checked here -- a setting that stops being
forwarded is invisible until a run comes out different.
"""

import sys

import numpy as np
import pytest

from numcosmo_py import Ncm
from numcosmo_py.app.loading import register_firecrown
from numcosmo_py.experiments.gauss_constraint import (
    create_data_object,
    create_mset as create_gauss_constraint_mset,
)
from numcosmo_py.interpolation.stats_dist import (
    CrossValidationMethod,
    InterpolationKernel,
    InterpolationMethod,
)
from numcosmo_py.sampling.apes import APES
from numcosmo_py.sampling.esmcmc import WalkerTypes, create_esmcmc

Ncm.cfg_init()


@pytest.fixture(name="mvnd_likelihood")
def fixture_mvnd_likelihood():
    """A two-dimensional Gaussian likelihood and its model set."""
    rng = Ncm.RNG.seeded_new(None, 4321)
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    data = Ncm.DataGaussCovMVND.new_full(2, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data]))

    return likelihood, mset


def test_create_esmcmc_apes_forwards_settings(tmp_path, mvnd_likelihood):
    """Every APES option create_esmcmc() takes reaches the walker."""
    likelihood, mset = mvnd_likelihood
    esmcmc = create_esmcmc(
        likelihood,
        mset,
        (tmp_path / "helpers").as_posix(),
        verbose=False,
        sampler=WalkerTypes.APES,
        nwalkers=50,
        interpolation_method=InterpolationMethod.VKDE,
        interpolation_kernel=InterpolationKernel.GAUSS,
        cv_method=CrossValidationMethod.SPLIT_M2LNP,
        use_apes_center_shrink=True,
        apes_defensive_frac=0.05,
        apes_defensive_scale=6.0,
        apes_defensive_nu=5.0,
        apes_vkde_points_per_dim=12.0,
        apes_uniform_weights=True,
        auto_kernel=True,
        split_fraction=0.4,
        local_fraction=0.5,
        use_threads=False,
        use_apes_threads=False,
    )

    walker = esmcmc.peek_walker()
    assert isinstance(walker, Ncm.FitESMCMCWalkerAPES)
    assert walker.get_center_shrink()
    assert walker.get_defensive_frac() == 0.05
    assert walker.get_defensive_scale() == 6.0
    assert walker.get_defensive_nu() == 5.0
    assert walker.get_vkde_points_per_dim() == 12.0
    assert walker.get_uniform_weights()
    assert walker.get_auto_kernel()
    assert walker.get_split_frac() == 0.4
    assert walker.get_cv_type() == CrossValidationMethod.SPLIT_M2LNP.genum

    # The description a catalog carries has to name the same proposal.
    desc = walker.desc()
    assert "Shrink-" in desc
    assert "unif" in desc
    assert "ppd=12" in desc
    assert "cv=split-m2lnp" in desc


def test_create_esmcmc_unknown_sampler(tmp_path, mvnd_likelihood):
    """A sampler the builder does not know is refused rather than guessed."""
    likelihood, mset = mvnd_likelihood

    with pytest.raises(ValueError, match="Unknown sampler"):
        create_esmcmc(
            likelihood,
            mset,
            (tmp_path / "unknown").as_posix(),
            verbose=False,
            sampler="not-a-sampler",  # type: ignore[arg-type]
            nwalkers=50,
        )


def test_apes_wrapper_runs_from_an_initial_sample(tmp_path, monkeypatch):
    """The emcee-shaped wrapper seeds its catalog from the sample it is given."""
    monkeypatch.chdir(tmp_path)

    nwalkers = 50
    ndim = 2

    def log_prob(theta, _args):
        return -0.5 * float(np.dot(theta, theta))

    sampler = APES(
        nwalkers=nwalkers,
        ndim=ndim,
        model=None,
        log_prob=log_prob,
        interpolation_method=InterpolationMethod.VKDE,
        interpolation_kernel=InterpolationKernel.GAUSS,
        over_smooth=1.0,
        local_fraction=0.5,
        center_shrink=True,
        cv_method=CrossValidationMethod.SPLIT_M2LNP,
        split_fraction=0.4,
        auto_kernel=False,
    )

    rng = np.random.default_rng(1234)
    initial_sample = rng.normal(size=(nwalkers, ndim))

    with pytest.raises(ValueError, match="wrong shape"):
        sampler.run_mcmc(initial_sample[:-1], 2)

    sampler.run_mcmc(initial_sample, 3)

    catalog = sampler.get_catalog()
    assert catalog.get_mean().shape == (ndim,)
    assert catalog.get_covar().shape == (ndim, ndim)


def test_gauss_constraint_data_object_verbose(capsys):
    """The constrained-Gaussian data object reports its normalization when asked."""
    dim = 3
    mset, _ = create_gauss_constraint_mset(dim)

    # The covariance is drawn from the generator it is handed, so the two calls only
    # describe the same distribution when each starts from the same seed.
    quiet = create_data_object(mset, dim, Ncm.RNG.seeded_new(None, 7), verbose=False)
    assert capsys.readouterr().out == ""

    loud = create_data_object(mset, dim, Ncm.RNG.seeded_new(None, 7), verbose=True)
    assert "Constant normalization" in capsys.readouterr().out

    assert loud.get_log_norma(mset) == quiet.get_log_norma(mset)


def test_register_firecrown_without_firecrown(monkeypatch):
    """A missing connector is not an error, and the streams are left as they were."""
    import builtins  # pylint: disable=import-outside-toplevel

    real_import = builtins.__import__

    def refuse_firecrown(name, *args, **kwargs):
        if name.startswith("firecrown"):
            raise ImportError("no firecrown here")
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "firecrown", raising=False)
    monkeypatch.setattr(builtins, "__import__", refuse_firecrown)

    stdout, stderr = sys.stdout, sys.stderr
    register_firecrown()

    assert sys.stdout is stdout
    assert sys.stderr is stderr
