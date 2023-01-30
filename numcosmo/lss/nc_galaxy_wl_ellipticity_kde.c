/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_kde.c
 *
 *  Mon July 27 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_kde.c
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * @short_description: Abstract class describing galaxy weak lensing ellipticity with Kernel Density Estimation
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
  NcmStatsDist1d *kde;
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
nc_galaxy_wl_ellipticity_kde_init (NcGalaxyWLEllipticityKDE *grsg)
{
  NcGalaxyWLEllipticityKDEPrivate * const self = grsg->priv = nc_galaxy_wl_ellipticity_kde_get_instance_private (grsg);
  
  self->obs    = NULL;
  self->kde    = NULL;
  self->e_vec  = NULL;
  self->r      = 0.0;
  self->twolnN = 0.0;
  self->len    = 0;
}

static void
_nc_galaxy_wl_ellipticity_kde_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityKDE *grsg = NC_GALAXY_WL_ELLIPTICITY_KDE (object);
  
  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_KDE (object));
  
  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_wl_ellipticity_kde_set_obs (grsg, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_kde_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityKDE *grsg = NC_GALAXY_WL_ELLIPTICITY_KDE (object);
  
  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_KDE (object));
  
  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_ellipticity_kde_peek_obs (grsg));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_kde_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_kde_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_ellipticity_kde_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_kde_parent_class)->finalize (object);
}

static void _nc_galaxy_wl_ellipticity_kde_m2lnP_initial_prep (NcGalaxyWLDist *gwld,NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
// static void _nc_galaxy_wl_ellipticity_kde_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i);
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
   * FIXME
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
  // wl_dist_class->m2lnP_prep         = &_nc_galaxy_wl_ellipticity_kde_m2lnP_prep;
  wl_dist_class->m2lnP              = &_nc_galaxy_wl_ellipticity_kde_m2lnP;
  wl_dist_class->gen                = &_nc_galaxy_wl_ellipticity_kde_gen;
  wl_dist_class->len                = &_nc_galaxy_wl_ellipticity_kde_len;
}

static void _nc_galaxy_wl_ellipticity_kde_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLEllipticityKDE *gkde               = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = gkde->priv;
  NcmStatsDist1dEPDF *s_kde                    = ncm_stats_dist1d_epdf_new_full (20000, NCM_STATS_DIST1D_EPDF_BW_RoT, 0.1, 1.0e-2);
  NcmStatsDist1d *sd1                          = NCM_STATS_DIST1D (s_kde);
  NcmVector *g_vec                             = ncm_vector_new (self->len);
  NcmMatrix *wl_obs                            = nc_galaxy_wl_ellipticity_kde_peek_obs (NC_GALAXY_WL_ELLIPTICITY_KDE (gwld));
  NcmVector *z_vec                             = nc_galaxy_redshift_spec_peek_z (NC_GALAXY_REDSHIFT_SPEC (gz));
  gdouble min_g_i                              = GSL_POSINF;
  gdouble max_g_i                              = GSL_NEGINF;
  int j = 0;
  gint gal_i;

  ncm_stats_dist1d_free (self->kde);
  ncm_vector_free (self->e_vec);

  ncm_stats_dist1d_set_compute_cdf (sd1, FALSE);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    const gdouble z_i = ncm_vector_get (z_vec, gal_i);
    const gdouble r_i = ncm_matrix_get (wl_obs, gal_i, 0);
    const gdouble g_i = ncm_matrix_get (wl_obs, gal_i, 1);
    const gdouble s_i = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r_i, z_i, z_cluster, z_cluster);

    max_g_i = MAX (max_g_i, g_i);
    min_g_i = MIN (min_g_i, g_i);
    /*if ((g_i > 0.0 && g_i < 0.05 && s_i > 0.0 && s_i < 0.05))*/
    //if (r_i > 0.07)
    {
      ncm_vector_set (g_vec, j, g_i);
      ncm_stats_dist1d_epdf_add_obs (s_kde, s_i);
      j++;
    }
  }

  ncm_stats_dist1d_epdf_set_max (s_kde, max_g_i*1.01);
  if (min_g_i > 0)
    ncm_stats_dist1d_epdf_set_min (s_kde, min_g_i*0.99);
  else
    ncm_stats_dist1d_epdf_set_min (s_kde, min_g_i*1.01);
  ncm_stats_dist1d_prepare (sd1);
  {
    const gdouble h = ncm_stats_dist1d_get_current_h (sd1);
    const gdouble hp = gsl_hypot (h, ncm_matrix_get (wl_obs, 0, 2));

    ncm_stats_dist1d_epdf_set_bw_type (s_kde, NCM_STATS_DIST1D_EPDF_BW_FIXED);
    s_kde->h_fixed = hp;
    ncm_stats_dist1d_prepare (sd1);
  }

  self->kde   = sd1;
  self->e_vec = g_vec;

// ncm_stats_dist1d_free (NCM_STATS_DIST1D (s_kde));
// ncm_vector_free (g_vec);
}

// static void
// _nc_galaxy_wl_ellipticity_kde_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i)
// {
  
// }

static gdouble
_nc_galaxy_wl_ellipticity_kde_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  NcGalaxyWLEllipticityKDE *gkde               = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = gkde->priv;
  gdouble res                                  = 0.0;


  const gdouble e_i = ncm_vector_get (self->e_vec, gal_i);
  const gdouble p = ncm_stats_dist1d_eval_p (self->kde, e_i);
  res += - 2 * log (p);
  
  return res;
}

static gdouble
_nc_galaxy_wl_ellipticity_kde_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLEllipticityKDE *grsg               = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = grsg->priv;
  const gdouble sigma_g                        = ncm_matrix_get (self->obs, 0, 2);
  
  return ncm_rng_gaussian_gen (rng, g_true, sigma_g);
}

static guint
_nc_galaxy_wl_ellipticity_kde_len (NcGalaxyWLDist *gwld)
{
  NcGalaxyWLEllipticityKDE *grsg               = NC_GALAXY_WL_ELLIPTICITY_KDE (gwld);
  NcGalaxyWLEllipticityKDEPrivate * const self = grsg->priv;
  
  return self->len;
}

/**
 * nc_galaxy_wl_ellipticity_kde_new:
 * @pos: a #NcGalaxyWLEllipticityKDEPos
 *
 * Creates a new #NcGalaxyWLEllipticityKDE using
 * @pos as the position type.
 *
 * Returns: (transfer full): a new NcGalaxyWLEllipticityKDE.
 */
NcGalaxyWLEllipticityKDE *
nc_galaxy_wl_ellipticity_kde_new ()
{
  NcGalaxyWLEllipticityKDE *grsg = g_object_new (NC_TYPE_GALAXY_WL_ELLIPTICITY_KDE,
                                                    NULL);
  
  return grsg;
}

/**
 * nc_galaxy_wl_ellipticity_kde_ref:
 * @grsg: a #NcGalaxyWLEllipticityKDE
 *
 * Increase the reference of @grsg by one.
 *
 * Returns: (transfer full): @grsg.
 */
NcGalaxyWLEllipticityKDE *
nc_galaxy_wl_ellipticity_kde_ref (NcGalaxyWLEllipticityKDE *grsg)
{
  return g_object_ref (grsg);
}

/**
 * nc_galaxy_wl_ellipticity_kde_free:
 * @grsg: a #NcGalaxyWLEllipticityKDE
 *
 * Decrease the reference count of @grsg by one.
 *
 */
void
nc_galaxy_wl_ellipticity_kde_free (NcGalaxyWLEllipticityKDE *grsg)
{
  g_object_unref (grsg);
}

/**
 * nc_galaxy_wl_ellipticity_kde_clear:
 * @grsg: a #NcGalaxyWLEllipticityKDE
 *
 * Decrease the reference count of @grsg by one, and sets the pointer *@grsg to
 * NULL.
 *
 */
void
nc_galaxy_wl_ellipticity_kde_clear (NcGalaxyWLEllipticityKDE **grsg)
{
  g_clear_object (grsg);
}

/**
 * nc_galaxy_wl_ellipticity_kde_set_obs:
 * @grsg: a #NcGalaxyWLEllipticityKDE
 * @obs: a #NcmMatrix
 *
 * Sets the observables matrix @obs.
 */
void
nc_galaxy_wl_ellipticity_kde_set_obs (NcGalaxyWLEllipticityKDE *grsg, NcmMatrix *obs)
{
  NcGalaxyWLEllipticityKDEPrivate * const self = grsg->priv;
  
  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);
  
  ncm_matrix_clear (&self->obs);
  
  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_ellipticity_kde_peek_obs:
 * @grsg: a #NcGalaxyWLEllipticityKDE
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_ellipticity_kde_peek_obs (NcGalaxyWLEllipticityKDE *grsg)
{
  NcGalaxyWLEllipticityKDEPrivate * const self = grsg->priv;
  
  return self->obs;
}

