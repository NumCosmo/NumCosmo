/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_likelihood.c
 *
 *  Mon May 08 16:12:03 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_likelihood.c
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
 * SECTION: nc_galaxy_wl_likelihood
 * @title: NcGalaxyWLLikelihood
 * @short_description: Class describing galaxy weak lensing distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy weak lensing distribution.
 * It is composed by three distributions: a shape distribution $P(s)$, a proxy redshift distribution $P(z_p)$, and a position distribution $P(z)P(r)$.
 * The shape distribution is defined by the abstract class #NcGalaxySDShape.
 * The proxy redshift distribution is defined by the abstract class #NcGalaxySDZProxy.
 * The position distribution is defined by the abstract class #NcGalaxySDPosition.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl_likelihood.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_stats_dist_kde.h"
#include "math/ncm_stats_dist_kernel_gauss.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLLikelihoodPrivate
{
  NcmMatrix *obs;
  NcGalaxySDShape *s_dist;
  NcGalaxySDZProxy *zp_dist;
  NcGalaxySDPosition *rz_dist;
  NcmStatsDistKDE *kde;
  gdouble cut_fraction;
  gdouble zp_min;
  gdouble zp_max;
  gdouble r_min;
  gdouble r_max;
  gint ndata;
  gboolean constructed;
  guint len;
};

enum
{
  PROP_0,
  PROP_OBS,
  PROP_S_DIST,
  PROP_ZP_DIST,
  PROP_RZ_DIST,
  PROP_R_MIN,
  PROP_R_MAX,
  PROP_ZP_MIN,
  PROP_ZP_MAX,
  PROP_NDATA,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLLikelihood, nc_galaxy_wl_likelihood, G_TYPE_OBJECT);

static void
nc_galaxy_wl_likelihood_init (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv = nc_galaxy_wl_likelihood_get_instance_private (gwl);

  NcmStatsDistKernelGauss *kernel = ncm_stats_dist_kernel_gauss_new (3);

  self->obs         = NULL;
  self->s_dist      = NULL;
  self->zp_dist     = NULL;
  self->rz_dist     = NULL;
  self->kde         = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (kernel), NCM_STATS_DIST_CV_NONE);
  self->len         = 0;
  self->r_max       = 0.0;
  self->r_min       = 0.0;
  self->zp_max      = 0.0;
  self->zp_min      = 0.0;
  self->ndata       = 100;
  self->constructed = FALSE;
}

static void
_nc_galaxy_wl_likelihood_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_wl_likelihood_set_obs (gwl, g_value_get_object (value));
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
    case PROP_R_MIN:
      self->r_min = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->r_min, <, self->r_max);

      break;
    case PROP_R_MAX:
      self->r_max = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->r_min, <, self->r_max);

      break;
    case PROP_ZP_MIN:
      self->zp_min = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->zp_min, <, self->zp_max);

      break;
    case PROP_ZP_MAX:
      self->zp_max = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->zp_min, <, self->zp_max);

      break;
    case PROP_NDATA:
      nc_galaxy_wl_likelihood_set_ndata (gwl, g_value_get_int (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_likelihood_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_likelihood_peek_obs (gwl));
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
    case PROP_R_MIN:
      g_value_set_double (value, self->r_min);
      break;
    case PROP_R_MAX:
      g_value_set_double (value, self->r_max);
      break;
    case PROP_ZP_MIN:
      g_value_set_double (value, self->zp_min);
      break;
    case PROP_ZP_MAX:
      g_value_set_double (value, self->zp_max);
      break;
    case PROP_NDATA:
      g_value_set_int (value, self->ndata);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_likelihood_dispose (GObject *object)
{
  NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  ncm_matrix_clear (&self->obs);
  nc_galaxy_sd_shape_clear (&self->s_dist);
  nc_galaxy_sd_z_proxy_clear (&self->zp_dist);
  nc_galaxy_sd_position_clear (&self->rz_dist);

  ncm_stats_dist_kde_clear (&self->kde);

  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_likelihood_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->finalize (object);
}

static void
_nc_galaxy_wl_likelihood_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->constructed (object);
  {
    NcGalaxyWLLikelihood *gwl                = NC_GALAXY_WL_LIKELIHOOD (object);
    NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

    g_assert_cmpfloat (self->r_min, <, self->r_max);
    g_assert_cmpfloat (self->zp_min, <, self->zp_max);

    self->constructed = TRUE;
  }
}

static void
nc_galaxy_wl_likelihood_class_init (NcGalaxyWLLikelihoodClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_likelihood_set_property;
  object_class->get_property = &_nc_galaxy_wl_likelihood_get_property;
  object_class->dispose      = &_nc_galaxy_wl_likelihood_dispose;
  object_class->finalize     = &_nc_galaxy_wl_likelihood_finalize;
  object_class->constructed  = &_nc_galaxy_wl_likelihood_constructed;

  /**
   * NcGalaxyWLLikelihood:obs:
   *
   * Galaxy weak lensing observables matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Galaxy weak lensing observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:s-dist:
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
   * NcGalaxyWLLikelihood:zp-dist:
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
   * NcGalaxyWLLikelihood:rz-dist:
   *
   * A #NcGalaxySDZPosition object.
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
   * NcGalaxyWLLikelihood:r-min:
   *
   * Minimum radius of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_R_MIN,
                                   g_param_spec_double ("r-min",
                                                        NULL,
                                                        "Minimum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:r-max:
   *
   * Maximum radius of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_R_MAX,
                                   g_param_spec_double ("r-max",
                                                        NULL,
                                                        "Maximum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:zp-min:
   *
   * Minimum redshift of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_ZP_MIN,
                                   g_param_spec_double ("zp-min",
                                                        NULL,
                                                        "Minimum redshift of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:zp-max:
   *
   * Maximum redshift of the weak lensing observables.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_ZP_MAX,
                                   g_param_spec_double ("zp-max",
                                                        NULL,
                                                        "Maximum redshift of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLLikelihood:ndata:
   *
   * Number of data points to sample for KDE.
   *
   */

  g_object_class_install_property (object_class,
                                   PROP_NDATA,
                                   g_param_spec_int ("ndata",
                                                     NULL,
                                                     "Number of data points to sample for KDE",
                                                     0.0, G_MAXINT, 100000.0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_galaxy_wl_likelihood_new:
 * @s_dist: a #NcGalaxySDShape
 * @zp_dist: a #NcGalaxySDZProxy
 * @rz_dist: a #NcGalaxySDPosition
 *
 * Creates a new galaxy weak lensing object.
 * Requires an instance of #NcGalaxySDShape, #NcGalaxySDZProxy, and #NcGalaxySDPosition.
 *
 * Returns: (transfer full): a new NcGalaxyWLLikelihood.
 */
NcGalaxyWLLikelihood *
nc_galaxy_wl_likelihood_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *rz_dist)
{
  NcGalaxyWLLikelihood *gwl = g_object_new (NC_TYPE_GALAXY_WL_LIKELIHOOD,
                                            "s-dist", s_dist,
                                            "zp-dist", zp_dist,
                                            "rz-dist", rz_dist,
                                            NULL);

  return gwl;
}

/**
 * nc_galaxy_wl_likelihood_ref:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Increase the reference of @gwl by one.
 *
 * Returns: (transfer full): @gwl.
 */
NcGalaxyWLLikelihood *
nc_galaxy_wl_likelihood_ref (NcGalaxyWLLikelihood *gwl)
{
  return g_object_ref (gwl);
}

/**
 * nc_galaxy_wl_likelihood_free:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Decrease the reference count of @gwl by one.
 *
 */
void
nc_galaxy_wl_likelihood_free (NcGalaxyWLLikelihood *gwl)
{
  g_object_unref (gwl);
}

/**
 * nc_galaxy_wl_likelihood_clear:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Decrease the reference count of @gwl by one, and sets the pointer *@gwl to
 * NULL.
 *
 */
void
nc_galaxy_wl_likelihood_clear (NcGalaxyWLLikelihood **gwl)
{
  g_clear_object (gwl);
}

/**
 * nc_galaxy_wl_likelihood_set_obs:
 * @gwl: a #NcGalaxyWLLikelihood
 * @obs: a #NcmMatrix
 *
 * Sets the observables matrix @obs.
 */
void
nc_galaxy_wl_likelihood_set_obs (NcGalaxyWLLikelihood *gwl, NcmMatrix *obs)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);

  ncm_matrix_clear (&self->obs);

  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_likelihood_peek_obs:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_likelihood_peek_obs (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  return self->obs;
}

/**
 * nc_galaxy_wl_likelihood_peek_kde:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmStatsDistKDE *
nc_galaxy_wl_likelihood_peek_kde (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  return self->kde;
}

void
nc_galaxy_wl_likelihood_prepare (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;
  NcmVector *sample                        = ncm_vector_new (3);
  NcmVector *pos                           = ncm_vector_new (2);
  NcmRNG *rng                              = ncm_rng_new (NULL);
  glong in_cut                             = 0;
  glong out_cut                            = 0;

  gint i;

  ncm_stats_dist_reset (NCM_STATS_DIST (self->kde));

  for (i = 0; i < self->ndata; i++)
  {
    while (TRUE)
    {
      gdouble zp = 0.0;
      gdouble z, r;

      nc_galaxy_sd_position_gen (self->rz_dist, pos, rng);

      z = ncm_vector_get (pos, 0);
      r = ncm_vector_get (pos, 1);

      if (nc_galaxy_sd_z_proxy_gen (self->zp_dist, rng, z, &zp))
      {
        const gdouble s = nc_galaxy_sd_shape_gen (self->s_dist, cosmo, dp, smd, z_cluster, rng, pos);

        ncm_vector_set (sample, 0, r);
        ncm_vector_set (sample, 1, zp);
        ncm_vector_set (sample, 2, s);

        if ((self->r_min <= r) && (r <= self->r_max) && (self->zp_min <= zp) && (zp <= self->zp_max))
          in_cut++;
        else
          out_cut++;

        ncm_stats_dist_add_obs (NCM_STATS_DIST (self->kde), sample);
        break;
      }
    }
  }

  self->cut_fraction = (gdouble) in_cut / (gdouble) (in_cut + out_cut);

  ncm_stats_dist_prepare (NCM_STATS_DIST (self->kde));

  ncm_rng_clear (&rng);
  ncm_vector_clear (&pos);
  ncm_vector_clear (&sample);
}

/**
 * nc_galaxy_wl_likelihood_eval_m2lnP:
 * @gwl: a #NcGalaxyWLLikelihood
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 *
 * Computes the observables probability given the theoretical modeling using
 * integration method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_galaxy_wl_likelihood_eval_m2lnP (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  return 0.0;
}

/**
 * nc_galaxy_wl_likelihood_kde_eval_m2lnP:
 * @gwl: a #NcGalaxyWLLikelihood
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 *
 * Computes the observables probability given the theoretical modeling using
 * kernel density estimation method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_galaxy_wl_likelihood_kde_eval_m2lnP (NcGalaxyWLLikelihood *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;
  NcmVector *data_vec                      = ncm_vector_new (3);
  gdouble res                              = 0.0;
  glong in_cut                             = 0;
  gint gal_i;

  nc_galaxy_wl_likelihood_prepare (gwl, cosmo, dp, smd, z_cluster);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    const gdouble r_i = ncm_matrix_get (self->obs, gal_i, 0);
    const gdouble z_i = ncm_matrix_get (self->obs, gal_i, 1);
    const gdouble s_i = ncm_matrix_get (self->obs, gal_i, 2);

    ncm_vector_set (data_vec, 0, r_i);
    ncm_vector_set (data_vec, 1, z_i);
    ncm_vector_set (data_vec, 2, s_i);

    if ((self->r_min <= r_i) && (r_i <= self->r_max) && (self->zp_min <= z_i) && (z_i <= self->zp_max))
    {
      in_cut++;
      res += ncm_stats_dist_eval_m2lnp (NCM_STATS_DIST (self->kde), data_vec);
    }
  }

  res += 2.0 * in_cut * log (self->cut_fraction);

  ncm_vector_clear (&data_vec);

  return res;
}

/**
 * nc_galaxy_wl_likelihood_len:
 * @gwll: a #NcGalaxyWL
 *
 * Returns: the number of galaxies in @gwl.
 */
guint
nc_galaxy_wl_likelihood_len (NcGalaxyWLLikelihood *gwll)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwll->priv;

  return self->len;
}

/**
 * nc_galaxy_wl_likelihood_set_cut:
 * @gwl: a #NcGalaxyWL
 * @zp_min: minimum redshift proxy $z_\mathrm{p,min}$
 * @zp_max: maximum redshift proxy $z_\mathrm{p,max}$
 * @r_min: minimum projected radius $r_\mathrm{min}$
 * @r_max: maximum projected radius $r_\mathrm{max}$
 *
 * Sets the cut in the observables.
 *
 */
void
nc_galaxy_wl_likelihood_set_cut (NcGalaxyWLLikelihood *gwl, const gdouble zp_min, const gdouble zp_max, const gdouble r_min, const gdouble r_max)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_assert_cmpfloat (zp_min, <, zp_max);
  g_assert_cmpfloat (r_min, <, r_max);

  self->zp_min = zp_min;
  self->zp_max = zp_max;
  self->r_min  = r_min;
  self->r_max  = r_max;
}

/**
 * nc_galaxy_wl_likelihood_set_ndata:
 * @gwl: a #NcGalaxyWL
 * @ndata: number of samples to take for KDE
 *
 * Sets the number of samples ndata.
 *
 */
void
nc_galaxy_wl_likelihood_set_ndata (NcGalaxyWLLikelihood *gwl, gdouble ndata)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  self->ndata = ndata;
}