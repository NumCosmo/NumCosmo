/***************************************************************************
 *            nc_xcor_kernel_cluster.c
 *
 *  Mon April 21 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcXcorKernelCluster:
 *
 * Base class for cluster kernels in cross-correlations.
 *
 * This is a placeholder base class for cluster kernels. The full implementation of a
 * cluster kernel would potentially involve integration over halo mass function:
 *
 * \begin{equation}
 *   K_\alpha(z,k) = \frac{dV}{dz}(z) \left[ \int dM \, b(M,z) \, \frac{dn}{dM}(M,z) \, S_\alpha(M,z) \, p(M,z) \right] \sqrt{P(k,z)}
 * \end{equation}
 *
 * where:
 * - $\frac{dV}{dz}$ is the comoving volume element per steradian,
 * - $\frac{dn}{dM}(M,z)$ is the halo mass function,
 * - $b(M,z)$ is the halo bias,
 * - $S_\alpha(M,z)$ is the selection function of bin $\alpha$,
 * - $p(M,z)$ is the catalog purity,
 * - $P(k,z)$ is the matter power spectrum.
 *
 * Different implementations will provide different approximations or methods for
 * computing this kernel.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "nc/xcor/nc_xcor_kernel_cluster.h"

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcXcorKernelCluster, nc_xcor_kernel_cluster, NC_TYPE_XCOR_KERNEL)

static void
nc_xcor_kernel_cluster_init (NcXcorKernelCluster *xclkc)
{
}

static void
nc_xcor_kernel_cluster_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_XCOR_KERNEL_CLUSTER (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_cluster_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_XCOR_KERNEL_CLUSTER (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_cluster_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cluster_parent_class)->dispose (object);
}

static void
nc_xcor_kernel_cluster_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cluster_parent_class)->finalize (object);
}

static void
nc_xcor_kernel_cluster_class_init (NcXcorKernelClusterClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &nc_xcor_kernel_cluster_set_property;
  object_class->get_property = &nc_xcor_kernel_cluster_get_property;
  object_class->dispose      = &nc_xcor_kernel_cluster_dispose;
  object_class->finalize     = &nc_xcor_kernel_cluster_finalize;
}

