#
# bianchi.py
#
# Sun Oct 18 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# bianchi.py
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

"""Bianchi I experiments from catalogs of distance moduli with sky positions."""

import math

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.datasets.dist_mu_obs import (
    SPEED_OF_LIGHT_KM_S,
    SampleRule,
    samples_by_rule,
)


def build_bianchi_experiment(
    obs: Nc.DistMuObs,
    *,
    mu_col: str | None = None,
    sigma_col: str | None = None,
    z_min: float = 0.0,
    z_max: float = math.inf,
    h_min: float = 50.0,
    h_max: float = 100.0,
    omegam_mean: float | None = 0.3,
    omegam_sigma: float = 0.04,
    ordering: bool = True,
    wall_scale: float = 0.1,
    zmax_factor: float = 1.05,
    sigma_v: float = 0.0,
    fit_sigma_v: bool = False,
    sample_rule: SampleRule = SampleRule.NONE,
    sample_col: str | None = None,
    fit_offsets: bool = False,
    reference_sample: str | None = None,
    fundamental_domain: bool = True,
    fit_velocity: bool = False,
    desc: str = "",
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Build a Bianchi I experiment from a catalog.

    The model is :class:`Nc.HCosmoBI` with every parameter free, H1, H2 and H3 bounded
    by ``h_min`` and ``h_max``. The data is :class:`Nc.DataHCosmoBIDistMu` filled from
    ``obs`` with :meth:`Nc.DataHCosmoBIDistMu.set_obs`, and its
    :class:`Nc.HCosmoBIDist` reaches ``zmax_factor`` times the largest selected
    redshift. The priors are a Gaussian on Omegam when ``omegam_mean`` is given and,
    when ``ordering`` is set, flat priors with walls of width ``wall_scale`` keeping
    H2 - H1 and H3 - H2 in [0, h_max - h_min].

    :param obs: the catalog.
    :param mu_col: column of the distance moduli, ``None`` for ``mu``.
    :param sigma_col: column of their uncertainties, ``None`` for ``sigma_mu``.
    :param z_min: lower redshift limit of the selection.
    :param z_max: upper redshift limit of the selection.
    :param h_min: lower bound of H1, H2 and H3 [km/s/Mpc].
    :param h_max: upper bound of H1, H2 and H3 [km/s/Mpc].
    :param omegam_mean: mean of the Gaussian prior on Omegam, ``None`` for no prior.
    :param omegam_sigma: standard deviation of the Gaussian prior on Omegam.
    :param ordering: whether to add the H1 < H2 < H3 priors.
    :param wall_scale: width of the walls of the ordering priors [km/s/Mpc].
    :param zmax_factor: the distance calculator covers this factor times the largest
        selected redshift.
    :param sigma_v: peculiar velocity dispersion [km/s] of the :class:`Nc.DistMuNuisance`
        added to the model set when it is positive, ``fit_sigma_v`` or ``fit_offsets``.
    :param fit_sigma_v: whether ``sigma_v`` is a free parameter.
    :param sample_rule: rule assigning the sample of each object, see
        :func:`numcosmo_py.datasets.dist_mu_obs.samples_by_rule`.
    :param sample_col: integer column of the catalog with the sample indices; it
        overrides ``sample_rule``.
    :param fit_offsets: whether one free offset per sample is added. A common shift of
        every offset is degenerate with the scale of H1, H2 and H3, so one offset must
        stay fixed: that of ``reference_sample``.
    :param reference_sample: name of the sample whose offset stays fixed at zero, or the
        sample index as a string; required with ``fit_offsets``.
    :param fundamental_domain: whether the orientation is bounded to one fundamental domain
        of its symmetry. The moduli see each principal axis up to sign, so eight points of
        the (u1, u2, u3) box give the same likelihood: the four translations
        (u1, u2 + k/4, u3 - k/4) and their images under u1 -> 1 - u1. The domain
        u1 <= 1/2, u2 <= 1/4 holds one copy of each.
    :param fit_velocity: whether a :class:`Nc.ObserverVelocity` with its three components
        free is added to the model set, the velocity of the observer relative to the rest
        frame of the model. The rest-frame redshift of every object must stay positive, so
        ``z_min`` must exceed the largest velocity allowed over the speed of light; the
        distance calculator reaches that much beyond the largest redshift.
    :param desc: description stored in the data object.
    :return: the experiment (``likelihood`` and ``model-set``) and the derived
        functions H0, Omega_s0, Omega_L0, H2_minus_H1 and H3_minus_H2.
    """
    cosmo = Nc.HCosmoBI.new()
    for name in ("H1", "H2", "H3"):
        pid = cosmo.param_index_from_name(name)[1]
        cosmo.param_set_lower_bound(pid, h_min)
        cosmo.param_set_upper_bound(pid, h_max)
    h_mid = 0.5 * (h_min + h_max)
    cosmo["H1"] = 0.97 * h_mid
    cosmo["H2"] = h_mid
    cosmo["H3"] = 1.03 * h_mid
    for name in ("u1", "u2", "u3"):
        cosmo[name] = 0.5
    if fundamental_domain:
        for name, upper in (("u1", 0.5), ("u2", 0.25)):
            pid = cosmo.param_index_from_name(name)[1]
            cosmo.param_set_upper_bound(pid, upper)
            cosmo[name] = 0.5 * upper
    cosmo.params_set_default_ftype()

    data = Nc.DataHCosmoBIDistMu.new(None)
    data.set_obs(obs, mu_col, sigma_col, sample_col, z_min, z_max)
    z_sel = data.peek_z()
    if z_sel.len() == 0:
        raise ValueError("No object of the catalog lies in the redshift selection.")

    sample_names: list[str] = []
    if sample_col is None and sample_rule != SampleRule.NONE:
        # The selection of set_obs, redone here to label the selected rows.
        columns = list(obs.peek_columns())
        table = obs.peek_data().to_numpy()
        z_all = table[:, columns.index("z")]
        mu_all = table[:, columns.index(mu_col or "mu")]
        sigma_all = table[:, columns.index(sigma_col or "sigma_mu")]
        selected = np.nonzero(
            (z_all >= z_min)
            & (z_all <= z_max)
            & np.isfinite(mu_all)
            & np.isfinite(sigma_all)
        )[0]
        assert len(selected) == z_sel.len()
        samples, names = samples_by_rule(obs, sample_rule, selected)
        # Only the samples with objects get an offset, in the order of their names.
        used = sorted({int(x) for x in samples if x >= 0})
        remap = {old: new for new, old in enumerate(used)}
        samples = [float(remap[int(x)]) if x >= 0 else -1.0 for x in samples]
        sample_names = [names[i] for i in used]
        data.set_samples(Ncm.Vector.new_array(samples))
    nsamples = int(max(data.peek_samples().get(i) for i in range(z_sel.len()))) + 1
    z_top = max(z_sel.get(i) for i in range(z_sel.len()))
    beta_max = 0.0
    if fit_velocity:
        beta_max = Nc.OBSERVER_VELOCITY_MAX_V / SPEED_OF_LIGHT_KM_S
        z_low = min(z_sel.get(i) for i in range(z_sel.len()))
        if z_low <= beta_max * (1.0 + z_low) / (1.0 - beta_max):
            raise ValueError(
                f"The smallest selected redshift {z_low:.4g} allows a rest-frame redshift "
                f"below zero for velocities up to {Nc.OBSERVER_VELOCITY_MAX_V:g} km/s; "
                "raise z_min."
            )
    dist = Nc.HCosmoBIDist.new(zmax_factor * z_top + beta_max * (1.0 + z_top))
    data.set_dist(dist)
    data.set_desc(desc + (f", samples {sample_names}" if sample_names else ""))

    models: list[Ncm.Model] = [cosmo]
    if sigma_v > 0.0 or fit_sigma_v or fit_offsets:
        nuisance = Nc.DistMuNuisance.new(nsamples if fit_offsets else 0)
        nuisance["sigmav"] = sigma_v
        nuisance.param_set_ftype(
            nuisance.param_index_from_name("sigmav")[1],
            Ncm.ParamType.FREE if fit_sigma_v else Ncm.ParamType.FIXED,
        )
        if fit_offsets:
            if reference_sample is None:
                raise ValueError("fit_offsets needs a reference_sample.")
            if reference_sample in sample_names:
                reference = sample_names.index(reference_sample)
            elif reference_sample.isdigit() and int(reference_sample) < nsamples:
                reference = int(reference_sample)
            else:
                raise ValueError(
                    f"Reference sample {reference_sample} not among the samples "
                    f"{sample_names or list(range(nsamples))}."
                )
            for s in range(nsamples):
                nuisance.param_set_ftype(
                    nuisance.param_index_from_name(f"dmu_{s}")[1],
                    Ncm.ParamType.FIXED if s == reference else Ncm.ParamType.FREE,
                )
        models.append(nuisance)
    if fit_velocity:
        velocity = Nc.ObserverVelocity.new()
        for name in ("vx", "vy", "vz"):
            velocity.param_set_ftype(
                velocity.param_index_from_name(name)[1], Ncm.ParamType.FREE
            )
        models.append(velocity)
    mset = Ncm.MSet.new_array(models)
    mset.prepare_fparam_map()

    dset = Ncm.Dataset.new_array([data])
    likelihood = Ncm.Likelihood.new(dset)
    if omegam_mean is not None:
        likelihood.priors_add(
            Ncm.PriorGaussParam.new_name("NcHCosmoBI:Omegam", omegam_mean, omegam_sigma)
        )
    if ordering:
        for fname in ("H2_minus_H1", "H3_minus_H2"):
            likelihood.priors_add(
                Ncm.PriorFlatFunc.new(
                    Ncm.MSetFuncList.new(f"NcHCosmoBI:{fname}", None),
                    0.0,
                    h_max - h_min,
                    wall_scale,
                    0.0,
                )
            )

    experiment = Ncm.ObjDictStr.new()
    experiment.add("likelihood", likelihood)
    experiment.add("model-set", mset)

    functions = Ncm.ObjArray.new()
    for fname in ("H0", "Omega_s0", "Omega_L0", "H2_minus_H1", "H3_minus_H2"):
        functions.add(Ncm.MSetFuncList.new(f"NcHCosmoBI:{fname}", None))

    return experiment, functions
