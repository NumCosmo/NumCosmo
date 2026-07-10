#
# cluster_richness_count.py
#
# Fri Jul 10 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cluster_richness_count.py
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

"""Factory function for a mock cluster mass-richness-count experiment.

This module builds a mock #NcDataClusterMassRichCount experiment: true masses
and redshifts are drawn uniformly, and the observed richness (galaxy count)
of each cluster is drawn from the Poisson-Lognormal distribution implied by a
#NcClusterMassAscaso mass-richness model, applying the model's selection cut.
"""

import numpy as np
from rich.console import Console
from rich.table import Table

from numcosmo_py import Ncm, Nc


def generate_cluster_richness_count(
    *,
    n_clusters: int,
    lnM_min: float,
    lnM_max: float,
    z_min: float,
    z_max: float,
    mup0: float,
    mup1: float,
    mup2: float,
    sigmap0: float,
    sigmap1: float,
    sigmap2: float,
    cut: float,
    seed: None | int,
    summary: bool,
) -> Ncm.ObjDictStr:
    """Generate a mock cluster mass-richness-count experiment.

    True log-masses and redshifts are drawn uniformly in
    [lnM_min, lnM_max] x [z_min, z_max]; the observed richness (galaxy count)
    of each cluster is a Poisson-Lognormal realization implied by the
    NcClusterMassAscaso model, keeping only clusters above the selection cut.
    """
    ascaso = Nc.ClusterMassAscaso(lnRichness_min=0.0, lnRichness_max=20.0)
    ascaso["mup0"] = mup0
    ascaso["mup1"] = mup1
    ascaso["mup2"] = mup2
    ascaso["sigmap0"] = sigmap0
    ascaso["sigmap1"] = sigmap1
    ascaso["sigmap2"] = sigmap2
    ascaso["cut"] = cut

    mset = Ncm.MSet.new_array([ascaso])

    if seed is not None:
        rng = Ncm.RNG.seeded_new("mt19937", seed)
    else:
        rng = Ncm.RNG.new("mt19937")
        rng.set_random_seed(True)

    lnM = np.array([rng.uniform_gen(lnM_min, lnM_max) for _ in range(n_clusters)])
    z = np.array([rng.uniform_gen(z_min, z_max) for _ in range(n_clusters)])

    lnM_v = Ncm.Vector.new_array(lnM)
    z_v = Ncm.Vector.new_array(z)
    N_v = Ncm.Vector.new(n_clusters)
    N_v.set_zero()

    dmrc = Nc.DataClusterMassRichCount.new()
    dmrc.set_data(lnM_v, z_v, N_v)
    dmrc.resample(mset, rng)

    dset = Ncm.Dataset.new_array([dmrc])
    likelihood = Ncm.Likelihood.new(dset)

    experiment = Ncm.ObjDictStr.new()
    experiment.add("likelihood", likelihood)
    experiment.add("model-set", mset)

    if summary:
        N = np.array(dmrc.peek_N().dup_array())

        console = Console()
        table = Table(title="Cluster Mass-Richness Count Mock")
        table.add_column("Quantity")
        table.add_column("Value")
        table.add_row("Clusters generated", f"{n_clusters}")
        table.add_row("Clusters kept (N >= cut)", f"{len(N)}")
        table.add_row("Mass range", f"[{np.exp(lnM_min):.2e}, {np.exp(lnM_max):.2e}]")
        table.add_row("Redshift range", f"[{z_min:.3f}, {z_max:.3f}]")
        if len(N) > 0:
            table.add_row("Richness range", f"[{N.min():.0f}, {N.max():.0f}]")
            table.add_row("Mean richness", f"{N.mean():.2f}")
        console.print(table)

    return experiment
