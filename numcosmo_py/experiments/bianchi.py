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

from numcosmo_py import Nc, Ncm


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
    cosmo.params_set_default_ftype()

    mset = Ncm.MSet.new_array([cosmo])
    mset.prepare_fparam_map()

    data = Nc.DataHCosmoBIDistMu.new(None)
    data.set_obs(obs, mu_col, sigma_col, z_min, z_max)
    z_sel = data.peek_z()
    if z_sel.len() == 0:
        raise ValueError("No object of the catalog lies in the redshift selection.")
    dist = Nc.HCosmoBIDist.new(
        zmax_factor * max(z_sel.get(i) for i in range(z_sel.len()))
    )
    data.set_dist(dist)
    data.set_desc(desc)

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
