/***************************************************************************
 *            nc_xcor_kernel_cluster_tophat.c
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
 * NcXcorKernelClusterTophat:
 *
 * Implementation of cluster kernel using a top-hat redshift distribution.
 *
 * This implements a simplified cluster kernel using the thin-z approximation,
 * where all mass-dependent factors are assumed constant within the redshift bin
 * and are absorbed into an overall normalization.
 *
 * The kernel is given by:
 * \begin{equation}
 *   K(z,k) = \Theta(z; z_{\rm lower}, z_{\rm upper}) \sqrt{P(k,z)}
 * \end{equation}
 *
 * where $\Theta(z; z_{\rm lower}, z_{\rm upper})$ is a normalized top-hat function:
 * \begin{equation}
 *   \Theta(z; z_{\rm lower}, z_{\rm upper}) = \begin{cases}
 *     \frac{1}{z_{\rm upper} - z_{\rm lower}} & \text{if } z_{\rm lower} \le z \le z_{\rm upper} \\
 *     0 & \text{otherwise}
 *   \end{cases}
 * \end{equation}
 *
 * and $P(k,z)$ is the matter power spectrum.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor_kernel_cluster_tophat.h"
#include "xcor/nc_xcor_kernel_component.h"
#include "xcor/nc_xcor.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelClusterTophat
{
  /*< private >*/
  NcXcorKernelCluster parent_instance;

  NcDistance *dist;
  gdouble z_lower;
  gdouble z_upper;
  gdouble V_upper;
  gdouble V_lower;

  NcXcorKernelComponent *clustering_comp;
};

enum
{
  PROP_0,
  PROP_Z_LOWER,
  PROP_Z_UPPER,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelClusterTophat, nc_xcor_kernel_cluster_tophat, NC_TYPE_XCOR_KERNEL_CLUSTER)

/*
 * Clustering Component Definition
 * Implements: TopHat(z) * sqrt(P(k,z))
 */

typedef struct _ClusteringComponentData
{
  NcXcorKernelClusterTophat *xclkc;
  NcDistance *dist;
  NcmPowspec *ps;
} ClusteringComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_CLUSTER_TOPHAT_GET_DATA(comp) \
        ((ClusteringComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _clustering_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);
static gdouble _clustering_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _clustering_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);
static void _clustering_component_data_clear (ClusteringComponentData *data);
static NcXcorKernelComponent *_nc_xcor_kernel_component_cluster_tophat_new (NcXcorKernelClusterTophat *xclkc, NcDistance *dist, NcmPowspec *ps);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_CLUSTER_TOPHAT,
                                      NcXcorKernelComponentClusterTophat,
                                      nc_xcor_kernel_component_cluster_tophat,
                                      _clustering_component_eval_kernel,
                                      _clustering_component_eval_prefactor,
                                      _clustering_component_get_limits,
                                      ClusteringComponentData,
                                      _clustering_component_data_clear)

static void
nc_xcor_kernel_cluster_tophat_init (NcXcorKernelClusterTophat *xclkc)
{
  xclkc->dist            = NULL;
  xclkc->z_lower         = 0.0;
  xclkc->z_upper         = 1.0;
  xclkc->clustering_comp = NULL;
}

static void
nc_xcor_kernel_cluster_tophat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_CLUSTER_TOPHAT (object));

  switch (prop_id)
  {
    case PROP_Z_LOWER:
      xclkc->z_lower = g_value_get_double (value);
      break;
    case PROP_Z_UPPER:
      xclkc->z_upper = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_cluster_tophat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_CLUSTER_TOPHAT (object));

  switch (prop_id)
  {
    case PROP_Z_LOWER:
      g_value_set_double (value, xclkc->z_lower);
      break;
    case PROP_Z_UPPER:
      g_value_set_double (value, xclkc->z_upper);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_cluster_tophat_dispose (GObject *object)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (object);

  nc_xcor_kernel_component_clear (&xclkc->clustering_comp);
  nc_distance_clear (&xclkc->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cluster_tophat_parent_class)->dispose (object);
}

static void
nc_xcor_kernel_cluster_tophat_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cluster_tophat_parent_class)->finalize (object);
}

static gdouble _nc_xcor_kernel_cluster_tophat_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_cluster_tophat_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_cluster_tophat_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_cluster_tophat_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_cluster_tophat_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_cluster_tophat_obs_params_len (NcXcorKernel *xclk);
static void _nc_xcor_kernel_cluster_tophat_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static GPtrArray *_nc_xcor_kernel_cluster_tophat_get_component_list (NcXcorKernel *xclk);

static gdouble
_nc_xcor_kernel_cluster_tophat_dndz (NcXcorKernelClusterTophat *xclkc, gdouble z)
{
  /* The boundary should always be respected */
  return 1.0 / (xclkc->z_upper - xclkc->z_lower);
}

static void
nc_xcor_kernel_cluster_tophat_constructed (GObject *object)
{
  G_OBJECT_CLASS (nc_xcor_kernel_cluster_tophat_parent_class)->constructed (object);
  {
    NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (object);
    NcXcorKernel *xclk               = NC_XCOR_KERNEL (xclkc);
    NcDistance *dist                 = nc_xcor_kernel_peek_dist (xclk);
    NcmPowspec *ps                   = nc_xcor_kernel_peek_powspec (xclk);

    g_assert (xclkc->z_upper > xclkc->z_lower);
    xclkc->dist = nc_distance_ref (dist);

    /* Create clustering component */
    g_assert_null (xclkc->clustering_comp);
    xclkc->clustering_comp = _nc_xcor_kernel_component_cluster_tophat_new (xclkc, dist, ps);
  }
}

static void
nc_xcor_kernel_cluster_tophat_class_init (NcXcorKernelClusterTophatClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  object_class->constructed = &nc_xcor_kernel_cluster_tophat_constructed;
  object_class->dispose     = &nc_xcor_kernel_cluster_tophat_dispose;
  object_class->finalize    = &nc_xcor_kernel_cluster_tophat_finalize;
  model_class->set_property = &nc_xcor_kernel_cluster_tophat_set_property;
  model_class->get_property = &nc_xcor_kernel_cluster_tophat_get_property;

  ncm_model_class_set_name_nick (model_class, "Cluster kernel with top-hat redshift distribution", "ClusterTophat");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelClusterTophat:z-lower:
   *
   * Lower redshift bound of the top-hat distribution.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LOWER,
                                   g_param_spec_double ("z-lower",
                                                        NULL,
                                                        "Lower redshift bound",
                                                        0.0, 10.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelClusterTophat:z-upper:
   *
   * Upper redshift bound of the top-hat distribution.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_UPPER,
                                   g_param_spec_double ("z-upper",
                                                        NULL,
                                                        "Upper redshift bound",
                                                        0.0, 10.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval_limber_z           = &_nc_xcor_kernel_cluster_tophat_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_cluster_tophat_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_cluster_tophat_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_cluster_tophat_add_noise;
  parent_class->obs_len                 = &_nc_xcor_kernel_cluster_tophat_obs_len;
  parent_class->obs_params_len          = &_nc_xcor_kernel_cluster_tophat_obs_params_len;
  parent_class->get_z_range             = &_nc_xcor_kernel_cluster_tophat_get_z_range;
  parent_class->get_component_list      = &_nc_xcor_kernel_cluster_tophat_get_component_list;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_kernel_cluster_tophat_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @z_lower: lower redshift bound
 * @z_upper: upper redshift bound
 *
 * Creates a new #NcXcorKernelClusterTophat with a top-hat redshift distribution
 * between @z_lower and @z_upper.
 *
 * Returns: (transfer full): a new #NcXcorKernelClusterTophat
 */
NcXcorKernelClusterTophat *
nc_xcor_kernel_cluster_tophat_new (NcDistance *dist, NcmPowspec *ps, gdouble z_lower, gdouble z_upper)
{
  NcXcorKernelClusterTophat *xclkc = g_object_new (NC_TYPE_XCOR_KERNEL_CLUSTER_TOPHAT,
                                                   "distance", dist,
                                                   "powspec", ps,
                                                   "z-lower", z_lower,
                                                   "z-upper", z_upper,
                                                   NULL);

  return xclkc;
}

/*
 * Old Limber interface implementation
 */

static void
_nc_xcor_kernel_cluster_tophat_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (xclk);

  *zmin = xclkc->z_lower;
  *zmax = xclkc->z_upper;
  *zmid = 0.5 * (xclkc->z_lower + xclkc->z_upper);
}

static gdouble
_nc_xcor_kernel_cluster_tophat_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (xclk);
  const gdouble dndz               = _nc_xcor_kernel_cluster_tophat_dndz (xclkc, z);
  const gdouble xi_t               = nc_distance_transverse (xclkc->dist, cosmo, z);

  return xi_t * xi_t * dndz;
}

static gdouble
_nc_xcor_kernel_cluster_tophat_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (xclk);

  return 1.0 / (xclkc->V_upper - xclkc->V_lower);
}

static void
_nc_xcor_kernel_cluster_tophat_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (xclk);
  NcDistance *dist                 = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                   = nc_xcor_kernel_peek_powspec (xclk);

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  xclkc->V_lower = nc_distance_comoving_volume_radial_integral (dist, cosmo, xclkc->z_lower);
  xclkc->V_upper = nc_distance_comoving_volume_radial_integral (dist, cosmo, xclkc->z_upper);

  g_assert_nonnull (xclkc->clustering_comp);
  nc_xcor_kernel_component_prepare (xclkc->clustering_comp, cosmo);
}

static void
_nc_xcor_kernel_cluster_tophat_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  /* No noise for now */
}

static guint
_nc_xcor_kernel_cluster_tophat_obs_len (NcXcorKernel *xclk)
{
  return 1;
}

static guint
_nc_xcor_kernel_cluster_tophat_obs_params_len (NcXcorKernel *xclk)
{
  return 0;
}

static GPtrArray *
_nc_xcor_kernel_cluster_tophat_get_component_list (NcXcorKernel *xclk)
{
  NcXcorKernelClusterTophat *xclkc = NC_XCOR_KERNEL_CLUSTER_TOPHAT (xclk);
  GPtrArray *comp_list             = g_ptr_array_new_with_free_func (g_object_unref);

  g_ptr_array_add (comp_list, xclkc->clustering_comp);

  return comp_list;
}

/*
 * Component implementation
 */

static void
_clustering_component_data_clear (ClusteringComponentData *data)
{
  /* No need to clear, these are weak references from parent kernel */
}

static gdouble
_clustering_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k)
{
  ClusteringComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CLUSTER_TOPHAT_GET_DATA (comp);
  const gdouble z               = nc_distance_inv_comoving (data->dist, cosmo, xi);
  const gdouble powspec         = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble tophat          = _nc_xcor_kernel_cluster_tophat_dndz (data->xclkc, z);

  return xi * xi * tophat * sqrt (powspec);
}

static gdouble
_clustering_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  return 1.0;
}

static void
_clustering_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max)
{
  ClusteringComponentData *data    = _NC_XCOR_KERNEL_COMPONENT_CLUSTER_TOPHAT_GET_DATA (comp);
  NcDistance *dist                 = data->dist;
  NcmPowspec *ps                   = data->ps;
  NcXcorKernelClusterTophat *xclkc = data->xclkc;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  *xi_min = nc_distance_comoving (dist, cosmo, xclkc->z_lower);
  *xi_max = nc_distance_comoving (dist, cosmo, xclkc->z_upper);
  *k_min  = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
  *k_max  = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_cluster_tophat_new (NcXcorKernelClusterTophat *xclkc, NcDistance *dist, NcmPowspec *ps)
{
  NcXcorKernelComponent *comp   = g_object_new (nc_xcor_kernel_component_cluster_tophat_get_type (), NULL);
  ClusteringComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CLUSTER_TOPHAT_GET_DATA (comp);

  data->xclkc = xclkc;
  data->dist  = dist;
  data->ps    = ps;

  return comp;
}

