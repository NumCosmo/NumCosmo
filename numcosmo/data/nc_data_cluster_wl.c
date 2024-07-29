/***************************************************************************
 *            nc_data_cluster_wl.c
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 *  Tue Jun 15 16:00:13 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.c
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_data_cluster_wl
 * @title: NcDataClusterWL
 * @short_description: Cluster weak lensing likelihood.
 * @stability: Unstable
 *
 * This class implements the weak lensing likelihood for galaxy clusters.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_wl.h"
#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_stats_dist.h"
#include "math/ncm_stats_dist_kde.h"
#include "math/ncm_stats_dist_vkde.h"
#include "math/ncm_stats_dist_kernel_gauss.h"
#include "math/ncm_stats_dist_kernel_st.h"
#include "math/ncm_integral_nd.h"
#include "numcosmo/nc_enum_types.h"
#include <math.h>
#include <gsl/gsl_math.h>

#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"

struct _NcDataClusterWLPrivate
{
  NcGalaxyWLObs *obs;
  NcGalaxySDShape *s_dist;
  NcGalaxySDZProxy *zp_dist;
  NcGalaxySDPosition *rz_dist;
  NcmMatrix *data;
  NcmStatsDist *kde;
  gdouble cut_fraction;
  gdouble zp_min;
  gdouble zp_max;
  gdouble theta_min;
  gdouble theta_max;
  gint ndata;
  gdouble prec;
  gboolean constructed;
  guint len;
  gdouble z_cluster;
  gdouble ra_cluster;
  gdouble dec_cluster;
  gboolean use_kde;
};

enum
{
  PROP_0,
  PROP_OBS,
  PROP_S_DIST,
  PROP_ZP_DIST,
  PROP_RZ_DIST,
  PROP_THETA_MIN,
  PROP_THETA_MAX,
  PROP_ZP_MIN,
  PROP_ZP_MAX,
  PROP_NDATA,
  PROP_PREC,
  PROP_Z_CLUSTER,
  PROP_RA_CLUSTER,
  PROP_DEC_CLUSTER,
  PROP_USE_KDE,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterWL, nc_data_cluster_wl, NCM_TYPE_DATA);

static void
nc_data_cluster_wl_init (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = dcwl->priv = nc_data_cluster_wl_get_instance_private (dcwl);
  NcmStatsDistKernelGauss *kernel     = ncm_stats_dist_kernel_gauss_new (4);

  self->obs         = NULL;
  self->s_dist      = NULL;
  self->zp_dist     = NULL;
  self->rz_dist     = NULL;
  self->data        = NULL;
  self->kde         = NCM_STATS_DIST (ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (kernel), NCM_STATS_DIST_CV_NONE));
  self->len         = 0;
  self->theta_max   = 0.0;
  self->theta_min   = 0.0;
  self->ndata       = 0.0;
  self->prec        = 1.0e-11;
  self->constructed = FALSE;
  self->z_cluster   = 0.0;
  self->ra_cluster  = 0.0;
  self->dec_cluster = 0.0;
  self->use_kde     = FALSE;
}

static void
nc_data_cluster_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_wl_obs_clear (&self->obs);
      nc_data_cluster_wl_set_obs (dcwl, g_value_get_object (value));
      break;
    case PROP_S_DIST:
      self->s_dist = g_value_dup_object (value);
      break;
    case PROP_ZP_DIST:
      self->zp_dist = g_value_dup_object (value);
      break;
    case PROP_RZ_DIST:
      self->rz_dist = g_value_dup_object (value);
      break;
    case PROP_THETA_MIN:
      self->theta_min = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->theta_min, <, self->theta_max);

      break;
    case PROP_THETA_MAX:
      self->theta_max = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->theta_min, <, self->theta_max);

      break;
    case PROP_NDATA:
      nc_data_cluster_wl_set_ndata (dcwl, g_value_get_int (value));
      break;
    case PROP_PREC:
      nc_data_cluster_wl_set_prec (dcwl, g_value_get_double (value));
      break;
    case PROP_Z_CLUSTER:
      self->z_cluster = g_value_get_double (value);
      break;
    case PROP_RA_CLUSTER:
      self->ra_cluster = g_value_get_double (value);
      break;
    case PROP_DEC_CLUSTER:
      self->dec_cluster = g_value_get_double (value);
      break;
    case PROP_USE_KDE:
      nc_data_cluster_wl_set_use_kde (dcwl, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_cluster_wl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_data_cluster_wl_peek_obs (dcwl));
      break;
    case PROP_S_DIST:
      g_value_set_object (value, self->s_dist);
      break;
    case PROP_ZP_DIST:
      g_value_set_object (value, self->zp_dist);
      break;
    case PROP_RZ_DIST:
      g_value_set_object (value, self->rz_dist);
      break;
    case PROP_THETA_MIN:
      g_value_set_double (value, self->theta_min);
      break;
    case PROP_THETA_MAX:
      g_value_set_double (value, self->theta_max);
      break;
    case PROP_NDATA:
      g_value_set_int (value, self->ndata);
      break;
    case PROP_PREC:
      g_value_set_double (value, self->prec);
      break;
    case PROP_Z_CLUSTER:
      g_value_set_double (value, self->z_cluster);
      break;
    case PROP_RA_CLUSTER:
      g_value_set_double (value, self->ra_cluster);
      break;
    case PROP_DEC_CLUSTER:
      g_value_set_double (value, self->dec_cluster);
      break;
    case PROP_USE_KDE:
      g_value_set_boolean (value, self->use_kde);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_cluster_wl_dispose (GObject *object)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  nc_galaxy_wl_obs_clear (&self->obs);
  ncm_matrix_clear (&self->data);
  nc_galaxy_sd_shape_clear (&self->s_dist);
  nc_galaxy_sd_z_proxy_clear (&self->zp_dist);
  nc_galaxy_sd_position_clear (&self->rz_dist);

  ncm_stats_dist_clear (&self->kde);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->dispose (object);
}

static void
nc_data_cluster_wl_finalize (GObject *object)
{
  /*NcDataClusterWL *dcwl = NC_DATA_CLUSTER_WL (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->finalize (object);
}

static void
_nc_data_cluster_wl_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->constructed (object);
  {
    NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
    NcDataClusterWLPrivate * const self = dcwl->priv;

    g_assert_cmpfloat (self->theta_min, <, self->theta_max);

    self->constructed = TRUE;
  }
}

static void _nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_cluster_wl_get_len (NcmData *data);
static void _nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset);

static void
nc_data_cluster_wl_class_init (NcDataClusterWLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_wl_set_property;
  object_class->get_property = nc_data_cluster_wl_get_property;
  object_class->dispose      = nc_data_cluster_wl_dispose;
  object_class->finalize     = nc_data_cluster_wl_finalize;
  object_class->constructed  = _nc_data_cluster_wl_constructed;

  /**
   * NcGalaxyWLObs:obs:
   *
   * A #NcGalaxyWLObs object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Weak lensing observation data",
                                                        NC_TYPE_GALAXY_WL_OBS,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:s-dist:
   *
   * A #NcGalaxySDShape object.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_S_DIST,
                                   g_param_spec_object ("s-dist",
                                                        NULL,
                                                        "Galaxy sample shape distribution",
                                                        NC_TYPE_GALAXY_SD_SHAPE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:zp-dist:
   *
   * A #NcGalaxySDZProxy object.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_ZP_DIST,
                                   g_param_spec_object ("zp-dist",
                                                        NULL,
                                                        "Galaxy sample proxy redshift distribution",
                                                        NC_TYPE_GALAXY_SD_Z_PROXY,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:rz-dist:
   *
   * A #NcGalaxySDPosition object.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_RZ_DIST,
                                   g_param_spec_object ("rz-dist",
                                                        NULL,
                                                        "Galaxy sample position distribution",
                                                        NC_TYPE_GALAXY_SD_POSITION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:theta-min:
   *
   * Minimum radius of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_THETA_MIN,
                                   g_param_spec_double ("theta-min",
                                                        NULL,
                                                        "Minimum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:theta-max:
   *
   * Maximum radius of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_THETA_MAX,
                                   g_param_spec_double ("theta-max",
                                                        NULL,
                                                        "Maximum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:ndata:
   *
   * Number of data points to sample for KDE.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_NDATA,
                                   g_param_spec_int ("ndata",
                                                     NULL,
                                                     "Number of data points to sample for KDE",
                                                     0, G_MAXINT, 10000,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:prec:
   *
   * Precision for integral.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision for integral",
                                                        0.0, G_MAXDOUBLE, 1.e-11,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:z-cluster:
   *
   * Cluster (halo) redshift.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_double ("z-cluster",
                                                        NULL,
                                                        "Cluster (halo) redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:ra-cluster:
   *
   * Cluster (halo) right ascension in radians.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_RA_CLUSTER,
                                   g_param_spec_double ("ra-cluster",
                                                        NULL,
                                                        "Cluster (halo) right ascension in radians.",
                                                        -G_PI, G_PI, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:dec-cluster:
   *
   * Cluster (halo) declination in radians.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_DEC_CLUSTER,
                                   g_param_spec_double ("dec-cluster",
                                                        NULL,
                                                        "Cluster (halo) declination in radians",
                                                        -G_PI_2, G_PI_2, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:use-kde:
   *
   * Whether to use KDE method.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_USE_KDE,
                                   g_param_spec_boolean ("use-kde",
                                                         NULL,
                                                         "Whether to use KDE method",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  data_class->m2lnL_val  = &_nc_data_cluster_wl_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_wl_get_len;
  data_class->prepare    = &_nc_data_cluster_wl_prepare;
}

struct _NcDataClusterWLIntArg
{
  NcGalaxySDPosition *pos_dist;
  NcGalaxySDZProxy *zp_dist;
  NcGalaxySDShape *s_dist;
  NcHICosmo *cosmo;
  NcHaloDensityProfile *dp;
  NcWLSurfaceMassDensity *smd;
  gdouble z_cluster;
  gdouble zp;
  gdouble et;
  gdouble ex;
  gdouble theta;
};

static void nc_data_cluster_wl_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_data_cluster_wl_int_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, DATA_CLUSTER_WL_INT, NcDataClusterWLInt, nc_data_cluster_wl_integ, nc_data_cluster_wl_int_dim, nc_data_cluster_wl_integ, struct _NcDataClusterWLIntArg);

static void
nc_data_cluster_wl_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcDataClusterWLInt *lh_int = NC_DATA_CLUSTER_WL_INT (intnd);
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z = ncm_vector_get (x, i);
    gdouble res     = nc_galaxy_sd_position_integ (lh_int->data.pos_dist, lh_int->data.theta, z) *
                      nc_galaxy_sd_z_proxy_integ (lh_int->data.zp_dist, z, lh_int->data.zp) *
                      nc_galaxy_sd_shape_integ_optzs (lh_int->data.s_dist,
                                                      lh_int->data.cosmo,
                                                      lh_int->data.dp,
                                                      lh_int->data.smd,
                                                      lh_int->data.z_cluster,
                                                      z,
                                                      lh_int->data.et,
                                                      lh_int->data.ex);

    ncm_vector_set (fval, i, res);
  }
}

static void
nc_data_cluster_wl_int_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 1;
  *fdim = 1;
}

void
nc_data_cluster_wl_prepare_kde (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcmVector *sample                   = ncm_vector_new (4);
  NcmRNG *rng                         = ncm_rng_new (NULL);
  gdouble in_cut                      = 0;
  gdouble out_cut                     = 0;
  gint i                              = 0;

  ncm_stats_dist_reset (self->kde);

  while (i < self->ndata)
  {
    const gdouble theta = nc_galaxy_sd_position_gen_theta (self->rz_dist, rng);
    gdouble z           = 0.0;
    gdouble zp          = 0.0;

    while (TRUE)
    {
      z = nc_galaxy_sd_position_gen_z (self->rz_dist, rng);

      if (nc_galaxy_sd_z_proxy_gen (self->zp_dist, rng, z, &zp))
        break;
    }

    if ((self->theta_min <= theta) && (theta <= self->theta_max))
      in_cut++;
    else
      out_cut++;

    i++;

    ncm_vector_set (sample, 0, theta);
    ncm_vector_set (sample, 1, zp);
    nc_galaxy_sd_shape_gen (self->s_dist, cosmo, dp, smd, self->z_cluster, rng, theta, z, ncm_vector_ptr (sample, 2), ncm_vector_ptr (sample, 3));

    ncm_stats_dist_add_obs (self->kde, sample);
  }

  self->cut_fraction = (gdouble) in_cut / (gdouble) (in_cut + out_cut);

  ncm_stats_dist_prepare (self->kde);

  ncm_rng_clear (&rng);
  ncm_vector_clear (&sample);
}

/**
 * nc_data_cluster_wl_eval_m2lnP:
 * @dcwl: a #NcDataClusterWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @m2lnP_gal: (out) (optional): a #NcmVector
 *
 * Computes the observables probability given the theoretical modeling using
 * integration method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_data_cluster_wl_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *m2lnP_gal)
{
  NcDataClusterWLPrivate * const self     = dcwl->priv;
  NcDataClusterWLInt *likelihood_integral = g_object_new (nc_data_cluster_wl_integ_get_type (), NULL);
  NcmIntegralND *lh_int                   = NCM_INTEGRAL_ND (likelihood_integral);
  gdouble verr, vzpi, vzpf, vres, vtmp;
  NcmVector *err = ncm_vector_new_data_static (&verr, 1, 1);
  NcmVector *zpi = ncm_vector_new_data_static (&vzpi, 1, 1);
  NcmVector *zpf = ncm_vector_new_data_static (&vzpf, 1, 1);
  NcmVector *res = ncm_vector_new_data_static (&vres, 1, 1);
  NcmVector *tmp = ncm_vector_new_data_static (&vtmp, 1, 1);
  gdouble zp_i   = 1.0e-11;
  gdouble zp_f   = 100.0;
  gdouble result = 0;
  guint gal_i;

  ncm_vector_set (zpi, 0, zp_i);
  ncm_vector_set (zpf, 0, zp_f);

  ncm_integral_nd_set_reltol (lh_int, self->prec);
  ncm_integral_nd_set_abstol (lh_int, 0.0);
  ncm_integral_nd_set_method (lh_int, NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V);

  if ((self->obs != NULL) && (self->data != NULL))
    g_assert_cmpuint (ncm_matrix_ncols (self->data), ==, 4);
  else
    g_error ("Weak lensing observables matrix is not set.");

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    gdouble m2lnP_gal_i;
    gdouble z_min, z_max;

    likelihood_integral->data.pos_dist  = NC_GALAXY_SD_POSITION (self->rz_dist);
    likelihood_integral->data.zp_dist   = NC_GALAXY_SD_Z_PROXY (self->zp_dist);
    likelihood_integral->data.s_dist    = NC_GALAXY_SD_SHAPE (self->s_dist);
    likelihood_integral->data.cosmo     = cosmo;
    likelihood_integral->data.dp        = dp;
    likelihood_integral->data.smd       = smd;
    likelihood_integral->data.z_cluster = self->z_cluster;
    likelihood_integral->data.theta     = ncm_matrix_get (self->data, gal_i, 0);
    likelihood_integral->data.zp        = ncm_matrix_get (self->data, gal_i, 1);
    likelihood_integral->data.et        = ncm_matrix_get (self->data, gal_i, 2);
    likelihood_integral->data.ex        = ncm_matrix_get (self->data, gal_i, 3);

    nc_galaxy_sd_shape_integ_optzs_prep (self->s_dist, cosmo, dp, smd, self->z_cluster, likelihood_integral->data.theta);

    nc_galaxy_sd_z_proxy_get_true_z_lim (self->zp_dist, likelihood_integral->data.zp, &z_min, &z_max);

    if (z_min < self->z_cluster)
    {
      vzpi = z_min;
      vzpf = self->z_cluster;
      ncm_integral_nd_eval (lh_int, zpi, zpf, tmp, err);

      vzpi = self->z_cluster;
      vzpf = z_max;
      ncm_integral_nd_eval (lh_int, zpi, zpf, res, err);
      vres += vtmp;
    }
    else
    {
      vzpi = z_min;
      vzpf = z_max;
      ncm_integral_nd_eval (lh_int, zpi, zpf, res, err);
    }

    m2lnP_gal_i = -2.0 * log (vres);

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  ncm_vector_free (err);
  ncm_vector_free (zpi);
  ncm_vector_free (zpf);
  ncm_vector_free (res);
  ncm_vector_free (tmp);
  ncm_integral_nd_clear (&lh_int);

  return result;
}

/**
 * nc_data_cluster_wl_kde_eval_m2lnP:
 * @dcwl: a #NcDataClusterWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @m2lnP_gal: (out) (optional): a #NcmVector
 *
 * Computes the observables probability given the theoretical modeling using
 * kernel density estimation method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_data_cluster_wl_kde_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *m2lnP_gal)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcmVector *data_vec                 = ncm_vector_new (4);
  gdouble res                         = 0.0;
  glong in_cut                        = 0;
  guint gal_i;

  if ((self->obs != NULL) && (self->data != NULL))
    g_assert_cmpuint (ncm_matrix_ncols (self->data), ==, 4);
  else
    g_error ("Weak lensing observables matrix is not set.");

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  nc_data_cluster_wl_prepare_kde (dcwl, cosmo, dp, smd);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    const gdouble theta_i = ncm_matrix_get (self->data, gal_i, 0);
    const gdouble z_i     = ncm_matrix_get (self->data, gal_i, 1);
    const gdouble et_i    = ncm_matrix_get (self->data, gal_i, 2);
    const gdouble ex_i    = ncm_matrix_get (self->data, gal_i, 3);

    ncm_vector_set (data_vec, 0, theta_i);
    ncm_vector_set (data_vec, 1, z_i);
    ncm_vector_set (data_vec, 2, et_i);
    ncm_vector_set (data_vec, 3, ex_i);

    if ((self->theta_min <= theta_i) && (theta_i <= self->theta_max))
    {
      const gdouble m2lnP_gal_i = ncm_stats_dist_eval_m2lnp (self->kde, data_vec);

      if (m2lnP_gal != NULL)
        ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

      in_cut++;
      res += m2lnP_gal_i;
    }
  }

  res += 2.0 * in_cut * log (self->cut_fraction);

  ncm_vector_clear (&data_vec);

  return res;
}

static void
_nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcHICosmo *cosmo                    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd         = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  guint i;

  m2lnL[0] = 0.0;

  if (self->use_kde)
    m2lnL[0] += nc_data_cluster_wl_kde_eval_m2lnP (dcwl, cosmo, dp, smd, NULL);
  else
    m2lnL[0] += nc_data_cluster_wl_eval_m2lnP (dcwl, cosmo, dp, smd, NULL);

  return;
}

static guint
_nc_data_cluster_wl_get_len (NcmData *data)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  if (self->data == NULL)
    return 0;
  else
    return self->len;
}

static void
_nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcHICosmo *cosmo                    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd         = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));

  g_assert ((cosmo != NULL) && (smd != NULL) && (dp != NULL));

  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);
}

/**
 * nc_data_cluster_wl_new:
 * @s_dist: a #NcGalaxySDShape
 * @zp_dist: a #NcGalaxySDZProxy
 * @pos_dist: a #NcGalaxySDPosition
 * @z_cluster: cluster (halo) redshift
 * @ra_cluster: cluster (halo) right ascension in radians
 * @dec_cluster: cluster (halo) declination in radians
 * @coord: a #NcDataClusterWLCoord
 *
 * Creates a new galaxy weak lensing object.
 * Requires an instance of #NcGalaxySDShape, #NcGalaxySDZProxy, and #NcGalaxySDPosition.
 *
 * Returns: (transfer full): a new NcDataClusterWL.
 */
NcDataClusterWL *
nc_data_cluster_wl_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *pos_dist, gdouble z_cluster, gdouble ra_cluster, gdouble dec_cluster)
{
  NcDataClusterWL *dcwl = g_object_new (NC_TYPE_DATA_CLUSTER_WL,
                                        "s-dist", s_dist,
                                        "zp-dist", zp_dist,
                                        "rz-dist", pos_dist,
                                        "z-cluster", z_cluster,
                                        "ra-cluster", ra_cluster,
                                        "dec-cluster", dec_cluster,
                                        NULL);

  return dcwl;
}

/**
 * nc_data_cluster_wl_ref:
 * @dcwl: a #NcDataClusterWL
 *
 * Increases the reference count of @dcwl by one.
 *
 * Returns: (transfer full): @dcwl
 */
NcDataClusterWL *
nc_data_cluster_wl_ref (NcDataClusterWL *dcwl)
{
  return g_object_ref (dcwl);
}

/**
 * nc_data_cluster_wl_free:
 * @dcwl: a #NcDataClusterWL
 *
 * Atomically decrements the reference count of @dcwl by one. If the reference count drops to 0,
 * all memory allocated by @dcwl is released.
 *
 */
void
nc_data_cluster_wl_free (NcDataClusterWL *dcwl)
{
  g_object_unref (dcwl);
}

/**
 * nc_data_cluster_wl_clear:
 * @dcwl: a #NcDataClusterWL
 *
 * The reference count of @dcwl is decreased and the pointer is set to NULL.
 *
 */
void
nc_data_cluster_wl_clear (NcDataClusterWL **dcwl)
{
  g_clear_object (dcwl);
}

/**
 * nc_data_cluster_wl_set_use_kde:
 * @dcwl: a #NcDataClusterWL
 * @kde: whether to use KDE method
 *
 * The reference count of @dcwl is decreased and the pointer is set to NULL.
 *
 */
void
nc_data_cluster_wl_set_use_kde (NcDataClusterWL *data, gboolean use_kde)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  self->use_kde = use_kde;
}

/**
 * nc_data_cluster_wl_set_prec:
 * @dcwl: a #NcDataClusterWL
 * @prec: precision for integral
 *
 * Sets the number of samples ndata.
 *
 */
void
nc_data_cluster_wl_set_prec (NcDataClusterWL *dcwl, gdouble prec)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;

  self->prec = prec;
}

/**
 * nc_data_cluster_wl_set_ndata:
 * @dcwl: a #NcDataClusterWL
 * @ndata: number of samples to take for KDE
 *
 * Sets the number of samples ndata.
 *
 */
void
nc_data_cluster_wl_set_ndata (NcDataClusterWL *dcwl, gdouble ndata)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;

  self->ndata = ndata;
}

/**
 * nc_data_cluster_wl_set_obs:
 * @dcwl: a #NcDataClusterWL
 * @obs: a #NcGalaxyWLObs
 *
 * Sets the observables matrix @obs and observables matrix in polar coordinates @obs_polar. The observables matrix must contain the following columns:
 * - $ra$ (right ascension) in radians
 * - $dec$ (declination) in radians
 * - $z$ (redshift)
 * - $e_1$ (ellipticity component 1)
 * - $e_2$ (ellipticity component 2)
 *
 */
void
nc_data_cluster_wl_set_obs (NcDataClusterWL *dcwl, NcGalaxyWLObs *obs)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcDistance *dist                    = nc_distance_new (1100);
  GArray *arr                         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmMatrix *data;
  guint gal_i;

  g_assert_cmpuint (nc_galaxy_wl_obs_len (obs), >, 0);

  nc_galaxy_wl_obs_clear (&self->obs);

  for (gal_i = 0; gal_i < nc_galaxy_wl_obs_len (obs); gal_i++)
  {
    gdouble ra    = nc_galaxy_wl_obs_get (obs, gal_i, 0);
    gdouble dec   = nc_galaxy_wl_obs_get (obs, gal_i, 1);
    gdouble z     = nc_galaxy_wl_obs_get (obs, gal_i, 2);
    gdouble e1    = nc_galaxy_wl_obs_get (obs, gal_i, 3);
    gdouble e2    = nc_galaxy_wl_obs_get (obs, gal_i, 4);
    gdouble theta = acos (sin (dec) * sin (self->dec_cluster) + cos (dec) * cos (self->dec_cluster) * cos (ra - self->ra_cluster));
    gdouble phi;
    gdouble et;
    gdouble ex;

    if ((theta > self->theta_min) && (theta <= self->theta_max))
    {
      switch (nc_galaxy_wl_obs_get_coord (obs))
      {
        case NC_GALAXY_WL_OBS_COORD_EUCLIDEAN:
          phi2 = 2 * atan2 (dec - self->dec_cluster, -(ra - self->ra_cluster) * cos (self->dec_cluster));
          break;
        case NC_GALAXY_WL_OBS_COORD_CELESTIAL:
          phi2 = 2 * atan2 (dec - self->dec_cluster, (ra - self->ra_cluster) * cos (self->dec_cluster));
          break;
        default:
          g_error ("Undefined coordinate system.");
          break;
      }

      et = -(e1 * cos (2.0 * phi) + e2 * sin (2.0 * phi));

      ex = e1 * sin (2.0 * phi) - e2 * cos (2.0 * phi);

      g_array_append_val (arr, theta);
      g_array_append_val (arr, z);
      g_array_append_val (arr, et);
      g_array_append_val (arr, ex);
    }
  }

  data = ncm_matrix_new_array (arr, 4);

  self->len  = ncm_matrix_nrows (data);
  self->obs  = obs;
  self->data = data;
}

/**
 * nc_data_cluster_wl_set_cut:
 * @dcwl: a #NcDataClusterWL
 * @theta_min: minimum projected radius $\theta_\mathrm{min}$
 * @theta_max: maximum projected radius $\theta_\mathrm{max}$
 *
 * Sets the cut in the observables.
 *
 */
void
nc_data_cluster_wl_set_cut (NcDataClusterWL *dcwl, const gdouble theta_min, const gdouble theta_max)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;

  g_assert_cmpfloat (theta_min, <, theta_max);

  self->theta_min = theta_min;
  self->theta_max = theta_max;
}

/**
 * nc_data_cluster_wl_gen_data:
 * @dcwl: a #NcDataClusterWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @nobs: number of observables to generate
 * @rng: a #NcmRNG
 *
 * Generates @nobs observables.
 */
void
nc_data_cluster_wl_gen_data (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, guint nobs, NcmRNG *rng)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcmMatrix *obs                      = ncm_matrix_new (nobs, 4);
  NcmVector *sample                   = ncm_vector_new (4);
  guint i;

  for (i = 0; i < nobs; i++)
  {
    gdouble theta = 0.0;
    gdouble z     = 0.0;
    gdouble zp    = 0.0;

    do {
      theta = nc_galaxy_sd_position_gen_theta (self->rz_dist, rng);
    } while ((self->theta_min > theta) || (theta > self->theta_max));

    while (TRUE)
    {
      z = nc_galaxy_sd_position_gen_z (self->rz_dist, rng);

      if (nc_galaxy_sd_z_proxy_gen (self->zp_dist, rng, z, &zp))
        break;
    }

    ncm_vector_set (sample, 0, theta);
    ncm_vector_set (sample, 1, zp);
    nc_galaxy_sd_shape_gen (self->s_dist, cosmo, dp, smd, self->z_cluster, rng, theta, z, ncm_vector_ptr (sample, 2), ncm_vector_ptr (sample, 3));

    ncm_matrix_set_row (obs, i, sample);
  }

  self->data = obs;

  ncm_vector_free (sample);
}

/**
 * nc_data_cluster_wl_peek_obs:
 * @dcwl: a #NcDataClusterWL
 *
 * Gets the weak lensing observables object.
 *
 * Returns: (transfer none): the observables object.
 */
NcGalaxyWLObs *
nc_data_cluster_wl_peek_obs (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;

  return self->obs;
}

/**
 * nc_data_cluster_wl_peek_data:
 * @dcwl: a #NcDataClusterWL
 *
 * Gets the weak lensing observables object.
 *
 * Returns: (transfer none): the observables object.
 */
NcmMatrix *
nc_data_cluster_wl_peek_data (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;

  return self->data;
}

/**
 * nc_data_cluster_wl_peek_kde:
 * @dcwl: a #NcDataClusterWL
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmStatsDist *
nc_data_cluster_wl_peek_kde (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;

  return self->kde;
}

