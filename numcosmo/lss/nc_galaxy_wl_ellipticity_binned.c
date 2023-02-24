/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_binned.c
 *
 *  Fri February 24 10:19:05 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_binned.c
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2023 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_galaxy_wl_ellipticity_binned
 * @title: NcGalaxyWLEllipticityBinned
 * @short_description: Class describing galaxy weak lensing ellipticity with binning
 * @stability: Unstable
 *
 *
 * Class defining a galaxy weak lensing ellipticity with binning.
 * probability distribution $P_\mathrm{wl}(g)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_ellipticity_binned.h"
#include "nc_enum_types.h"

#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLEllipticityBinnedPrivate
{
  NcmMatrix *obs;
  NcmVector *bins;
  NcmObjArray *bin_obs;
  gdouble r;
  gdouble twolnN;
  guint len;
};

enum
{
  PROP_0,
  PROP_OBS,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLEllipticityBinned, nc_galaxy_wl_ellipticity_binned, NC_TYPE_GALAXY_WL_DIST);

static void
nc_galaxy_wl_ellipticity_binned_init (NcGalaxyWLEllipticityBinned *gebin)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv = nc_galaxy_wl_ellipticity_binned_get_instance_private (gebin);

  self->obs     = NULL;
  self->bins    = NULL;
  self->bin_obs = NULL;
  self->r       = 0.0;
  self->twolnN  = 0.0;
  self->len     = 0;
}

static void
_nc_galaxy_wl_ellipticity_binned_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityBinned *gebin = NC_GALAXY_WL_ELLIPTICITY_BINNED (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_BINNED (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_wl_ellipticity_binned_set_obs (gebin, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_binned_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityBinned *gebin = NC_GALAXY_WL_ELLIPTICITY_BINNED (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_BINNED (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_ellipticity_binned_peek_obs (gebin));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_binned_dispose (GObject *object)
{
  NcGalaxyWLEllipticityBinned *gebin              = NC_GALAXY_WL_ELLIPTICITY_BINNED (object);
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  ncm_matrix_clear (&self->obs);
  ncm_stats_dist1d_epdf_clear (&self->kde);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_binned_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_ellipticity_binned_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_binned_parent_class)->finalize (object);
}

static void _nc_galaxy_wl_ellipticity_binned_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
static gdouble _nc_galaxy_wl_ellipticity_binned_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
static gdouble _nc_galaxy_wl_ellipticity_binned_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
static guint _nc_galaxy_wl_ellipticity_binned_len (NcGalaxyWLDist *gwld);

static void
nc_galaxy_wl_ellipticity_binned_class_init (NcGalaxyWLEllipticityBinnedClass *klass)
{
  NcGalaxyWLDistClass *wl_dist_class = NC_GALAXY_WL_DIST_CLASS (klass);
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_ellipticity_binned_set_property;
  object_class->get_property = &_nc_galaxy_wl_ellipticity_binned_get_property;
  object_class->dispose      = &_nc_galaxy_wl_ellipticity_binned_dispose;
  object_class->finalize     = &_nc_galaxy_wl_ellipticity_binned_finalize;

  /**
   * NcGalaxyWLEllipticityBinned:obs:
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

  wl_dist_class->m2lnP_initial_prep = &_nc_galaxy_wl_ellipticity_binned_m2lnP_initial_prep;
  wl_dist_class->m2lnP = &_nc_galaxy_wl_ellipticity_binned_m2lnP;
  wl_dist_class->gen   = &_nc_galaxy_wl_ellipticity_binned_gen;
  wl_dist_class->len   = &_nc_galaxy_wl_ellipticity_binned_len;
}

static void
_nc_galaxy_wl_ellipticity_binned_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLEllipticityBinned *gebin              = NC_GALAXY_WL_ELLIPTICITY_BINNED (gwld);
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;
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
_nc_galaxy_wl_ellipticity_binned_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  NcGalaxyWLEllipticityBinned *gebin              = NC_GALAXY_WL_ELLIPTICITY_BINNED (gwld);
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;
  const gdouble e_i                            = ncm_vector_get (self->e_vec, gal_i);
  const gdouble p                              = ncm_stats_dist1d_eval_p (NCM_STATS_DIST1D (self->kde), e_i);
  gdouble res                                  = 0.0;

  res += -2 * log (p);

  return res;
}

static gdouble
_nc_galaxy_wl_ellipticity_binned_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLEllipticityBinned *gebin              = NC_GALAXY_WL_ELLIPTICITY_BINNED (gwld);
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  return ncm_stats_dist1d_gen (NCM_STATS_DIST1D (self->kde), rng);
}

static guint
_nc_galaxy_wl_ellipticity_binned_len (NcGalaxyWLDist *gwld)
{
  NcGalaxyWLEllipticityBinned *gebin              = NC_GALAXY_WL_ELLIPTICITY_BINNED (gwld);
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  return self->len;
}

/**
 * nc_galaxy_wl_ellipticity_binned_new:
 *
 * Creates a new #NcGalaxyWLEllipticityBinned
 *
 * Returns: (transfer full): a new NcGalaxyWLEllipticityBinned.
 */
NcGalaxyWLEllipticityBinned *
nc_galaxy_wl_ellipticity_binned_new ()
{
  NcGalaxyWLEllipticityBinned *gebin = g_object_new (NC_TYPE_GALAXY_WL_ELLIPTICITY_BINNED,
                                                  NULL);

  return gebin;
}

/**
 * nc_galaxy_wl_ellipticity_binned_ref:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 *
 * Increase the reference of @gebin by one.
 *
 * Returns: (transfer full): @gebin.
 */
NcGalaxyWLEllipticityBinned *
nc_galaxy_wl_ellipticity_binned_ref (NcGalaxyWLEllipticityBinned *gebin)
{
  return g_object_ref (gebin);
}

/**
 * nc_galaxy_wl_ellipticity_binned_free:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 *
 * Decrease the reference count of @gebin by one.
 *
 */
void
nc_galaxy_wl_ellipticity_binned_free (NcGalaxyWLEllipticityBinned *gebin)
{
  g_object_unref (gebin);
}

/**
 * nc_galaxy_wl_ellipticity_binned_clear:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 *
 * Decrease the reference count of @gebin by one, and sets the pointer *@gebin to
 * NULL.
 *
 */
void
nc_galaxy_wl_ellipticity_binned_clear (NcGalaxyWLEllipticityBinned **gebin)
{
  g_clear_object (gebin);
}

/**
 * nc_galaxy_wl_ellipticity_binned_set_obs:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 * @obs: a #NcmMatrix
 *
 * Sets the observable matrix @obs.
 */
void
nc_galaxy_wl_ellipticity_binned_set_obs (NcGalaxyWLEllipticityBinned *gebin, NcmMatrix *obs)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);

  ncm_matrix_clear (&self->obs);

  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_ellipticity_binned_set_obs:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 * @obs: a #NcmMatrix
 * @bins: a #NcmVector
 *
 * Sets thje binned data matrix *bin_obs
 */
void
nc_galaxy_wl_ellipticity_binned_set_bin (NcGalaxyWLEllipticityBinned *gebin, NcmVector *bins)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;
  NcmObjArray *bin_obs = ncm_obj_array_sized_new (ncm_vector_len (bins));

  if (ncm_vector_len (bins) == 1)
  {
    /* FIX ME */
  }
  else
  {
    for (bin_i = 0; bin_i < ncm_vector_len (bins) - 1; bin_i++)
    {
      NcmMatrix *bin_data = ncm_matrix_new (0, 3);
      gint bin_min       = ncm_vector_get (bins, bin_i);
      gint bin_max       = ncm_vector_get (bins, bin_i+1);

      for (gal_i = 0;  gal_i < self->len; gal_i ++)
      {
        gdouble r_i   = ncm_matrix_get (self->obs, gal_i, 0);
        gdouble e_i   = ncm_matrix_get (self->obs, gal_i, 1);
        gdouble err_i = ncm_matrix_get (self->obs, gal_i, 2);

        if (r_i >= bin_min && r_i < bin_max)
        {
          NcmMatrix data = ncm_matrix_new (ncm_matrix_nrows (bin_data)+1, 3);

          ncm_matrix_set (data, ncm_matrix_nrows (data)-1, 0, r_i);
          ncm_matrix_set (data, ncm_matrix_nrows (data)-1, 1, 3_i);
          ncm_matrix_set (data, ncm_matrix_nrows (data)-1, 2, err_i);

          if (ncm_matrix_size (bin_data) == 0)
            *bin_data = data;
          else
          {
            for (i = 0; i < ncm_vector_len (data); i++)
            {
              ncm_matrix_set (data, i, 0, ncm_matrix_get (bin_data, i, 0));
              ncm_matrix_set (data, i, 1, ncm_matrix_get (bin_data, i, 1));
              ncm_matrix_set (data, i, 2, ncm_matrix_get (bin_data, i, 2));
            }
            *bin_data = data;
          }
        }
      }
      ncm_obj_array_set (bin_obs, bin_i, bin_data)
    }
  }
  self->bins    = bins
  self->bin_obs = bin_obs
}

/**
 * nc_galaxy_wl_ellipticity_binned_peek_obs:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 *
 * Gets the observable matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_ellipticity_binned_peek_obs (NcGalaxyWLEllipticityBinned *gebin)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  return self->obs;
}

