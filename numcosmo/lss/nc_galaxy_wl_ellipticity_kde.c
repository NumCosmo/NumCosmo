/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_kde.c
 *
 *  Mon July 27 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_kde.c
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_galaxy_wl_ellipticity_kde
 * @title: NcGalaxyWLEllipticityKDE
 * @short_description: Class describing galaxy weak lensing ellipticity with Kernel Density Estimation
 * @stability: Unstable
 *
 *
 * Class defining a galaxy weak lensing ellipticity with Kernel Density Estimation.
 * probability distribution $P_\mathrm{wl}(g)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_ellipticity_kde.h"
#include "nc_enum_types.h"
#include "math/ncm_stats_dist1d_epdf.h"
#include "lss/nc_galaxy_redshift_spec.h"

#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLEllipticityKDEPrivate
{
  NcmMatrix *obs;
  NcmStatsDist1dEPDF *kde;
  NcmVector *e_vec;
  gdouble r;
  gdouble twolnN;
  guint len;
};

enum
{
  PROP_0,
  PROP_OBS,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLEllipticityKDE, nc_galaxy_wl_ellipticity_kde, NC_TYPE_GALAXY_WL_DIST);

static void
nc_galaxy_wl_ellipticity_kde_init (NcGalaxyWLEllipticityKDE *gekde)
{
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv = nc_galaxy_wl_ellipticity_kde_get_instance_private (gekde);

  self->obs    = NULL;
  self->kde    = ncm_stats_dist1d_epdf_new_full (20000, NCM_STATS_DIST1D_EPDF_BW_RoT, 0.1, 1.0e-2);
  self->e_vec  = NULL;
  self->r      = 0.0;
  self->twolnN = 0.0;
  self->len    = 0;

  ncm_stats_dist1d_set_compute_cdf (NCM_STATS_DIST1D (self->kde), FALSE);
}

static void
_nc_galaxy_wl_ellipticity_kde_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityKDE *gekde = NC_GALAXY_WL_ELLIPTICITY_KDE (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_KDE (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_wl_ellipticity_kde_set_obs (gekde, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_kde_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityKDE *gekde = NC_GALAXY_WL_ELLIPTICITY_KDE (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_KDE (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_ellipticity_kde_peek_obs (gekde));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_kde_dispose (GObject *object)
{
  NcGalaxyWLEllipticityKDE *gekde              = NC_GALAXY_WL_ELLIPTICITY_KDE (object);
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;

  ncm_matrix_clear (&self->obs);
  ncm_stats_dist1d_epdf_clear (&self->kde);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_kde_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_ellipticity_kde_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_kde_parent_class)->finalize (object);
}

static void _nc_galaxy_wl_ellipticity_kde_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
static gdouble _nc_galaxy_wl_ellipticity_kde_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
static gdouble _nc_galaxy_wl_ellipticity_kde_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
static guint _nc_galaxy_wl_ellipticity_kde_len (NcGalaxyWLDist *gwld);

static void
nc_galaxy_wl_ellipticity_kde_class_init (NcGalaxyWLEllipticityKDEClass *klass)
{
  NcGalaxyWLDistClass *wl_dist_class = NC_GALAXY_WL_DIST_CLASS (klass);
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_ellipticity_kde_set_property;
  object_class->get_property = &_nc_galaxy_wl_ellipticity_kde_get_property;
  object_class->dispose      = &_nc_galaxy_wl_ellipticity_kde_dispose;
  object_class->finalize     = &_nc_galaxy_wl_ellipticity_kde_finalize;

  /**
   * NcGalaxyWLEllipticityKDE:obs:
   *
   * Galaxy ellipticity observable matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Galaxy observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  wl_dist_class->m2lnP_initial_prep = &_nc_galaxy_wl_ellipticity_kde_m2lnP_initial_prep;
  wl_dist_class->m2lnP = &_nc_galaxy_wl_ellipticity_kde_m2lnP;
  wl_dist_class->gen   = &_nc_galaxy_wl_ellipticity_kde_gen;
  wl_dist_class->len   = &_nc_galaxy_wl_ellipticity_kde_len;
}

static void
_nc_galaxy_wl_ellipticity_kde_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLEllipticityKDE *gekde              = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;
  NcmVector *g_vec                             = ncm_vector_new (self->len);
  NcmVector *z_vec                             = nc_galaxy_redshift_spec_peek_z (NC_GALAXY_REDSHIFT_SPEC (gz));
  gdouble min_g_i                              = GSL_POSINF;
  gdouble max_g_i                              = GSL_NEGINF;
  int j                                        = 0;
  gint gal_i;

  ncm_stats_dist1d_epdf_reset (self->kde);
  ncm_stats_dist1d_epdf_set_bw_type (self->kde, NCM_STATS_DIST1D_EPDF_BW_RoT);
  ncm_vector_free (self->e_vec);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    const gdouble z_i = ncm_vector_get (z_vec, gal_i);
    const gdouble r_i = ncm_matrix_get (self->obs, gal_i, 0);
    const gdouble g_i = ncm_matrix_get (self->obs, gal_i, 1);
    const gdouble s_i = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r_i, z_i, z_cluster, z_cluster);

    max_g_i = MAX (max_g_i, g_i);
    min_g_i = MIN (min_g_i, g_i);

    ncm_vector_set (g_vec, j, g_i);
    ncm_stats_dist1d_epdf_add_obs (self->kde, s_i);
    j++;
  }

  ncm_stats_dist1d_epdf_set_max (self->kde, max_g_i * 1.01);

  if (min_g_i > 0)
    ncm_stats_dist1d_epdf_set_min (self->kde, min_g_i * 0.99);
  else
    ncm_stats_dist1d_epdf_set_min (self->kde, min_g_i * 1.01);

  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (self->kde));

  const gdouble h  = ncm_stats_dist1d_get_current_h (NCM_STATS_DIST1D (self->kde));
  const gdouble hp = gsl_hypot (h, ncm_matrix_get (self->obs, 0, 2));

  ncm_stats_dist1d_epdf_set_bw_type (self->kde, NCM_STATS_DIST1D_EPDF_BW_FIXED);
  self->kde->h_fixed = hp;
  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (self->kde));

  self->e_vec = g_vec;
}

static gdouble
_nc_galaxy_wl_ellipticity_kde_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  NcGalaxyWLEllipticityKDE *gekde              = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;
  const gdouble e_i                            = ncm_vector_get (self->e_vec, gal_i);
  const gdouble p                              = ncm_stats_dist1d_eval_p (NCM_STATS_DIST1D (self->kde), e_i);
  gdouble res                                  = 0.0;

  res += -2 * log (p);

  return res;
}

static gdouble
_nc_galaxy_wl_ellipticity_kde_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLEllipticityKDE *gekde              = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;

  return ncm_stats_dist1d_gen (NCM_STATS_DIST1D (self->kde), rng);
}

static guint
_nc_galaxy_wl_ellipticity_kde_len (NcGalaxyWLDist *gwld)
{
  NcGalaxyWLEllipticityKDE *gekde              = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;

  return self->len;
}

/**
 * nc_galaxy_wl_ellipticity_kde_new:
 *
 * Creates a new #NcGalaxyWLEllipticityKDE
 *
 * Returns: (transfer full): a new NcGalaxyWLEllipticityKDE.
 */
NcGalaxyWLEllipticityKDE *
nc_galaxy_wl_ellipticity_kde_new ()
{
  NcGalaxyWLEllipticityKDE *gekde = g_object_new (NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE,
                                                  NULL);

  return gekde;
}

/**
 * nc_galaxy_wl_ellipticity_kde_ref:
 * @gekde: a #NcGalaxyWLEllipticityKDE
 *
 * Increase the reference of @gekde by one.
 *
 * Returns: (transfer full): @gekde.
 */
NcGalaxyWLEllipticityKDE *
nc_galaxy_wl_ellipticity_kde_ref (NcGalaxyWLEllipticityKDE *gekde)
{
  return g_object_ref (gekde);
}

/**
 * nc_galaxy_wl_ellipticity_kde_free:
 * @gekde: a #NcGalaxyWLEllipticityKDE
 *
 * Decrease the reference count of @gekde by one.
 *
 */
void
nc_galaxy_wl_ellipticity_kde_free (NcGalaxyWLEllipticityKDE *gekde)
{
  g_object_unref (gekde);
}

/**
 * nc_galaxy_wl_ellipticity_kde_clear:
 * @gekde: a #NcGalaxyWLEllipticityKDE
 *
 * Decrease the reference count of @gekde by one, and sets the pointer *@gekde to
 * NULL.
 *
 */
void
nc_galaxy_wl_ellipticity_kde_clear (NcGalaxyWLEllipticityKDE **gekde)
{
  g_clear_object (gekde);
}

/**
 * nc_galaxy_wl_ellipticity_kde_set_obs:
 * @gekde: a #NcGalaxyWLEllipticityKDE
 * @obs: a #NcmMatrix
 *
 * Sets the observable matrix @obs.
 */
void
nc_galaxy_wl_ellipticity_kde_set_obs (NcGalaxyWLEllipticityKDE *gekde, NcmMatrix *obs)
{
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;

  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);

  ncm_matrix_clear (&self->obs);

  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_ellipticity_kde_peek_obs:
 * @gekde: a #NcGalaxyWLEllipticityKDE
 *
 * Gets the observable matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_ellipticity_kde_peek_obs (NcGalaxyWLEllipticityKDE *gekde)
{
  NcGalaxyWLEllipticityKDEPrivate * const self = gekde->priv;

  return self->obs;
}

