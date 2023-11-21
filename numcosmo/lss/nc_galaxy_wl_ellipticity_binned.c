/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_binned.c
 *
 *  Fri February 24 10:19:05 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_binned.c
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
  NcmVector *bins;
  NcmObjArray *binobs;
  gdouble r;
  gdouble twolnN;
  guint len;
};

enum
{
  PROP_0,
  PROP_BINS,
  PROP_BINOBS,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLEllipticityBinned, nc_galaxy_wl_ellipticity_binned, NC_TYPE_GALAXY_WL_DIST)

static void
nc_galaxy_wl_ellipticity_binned_init (NcGalaxyWLEllipticityBinned *gebin)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv = nc_galaxy_wl_ellipticity_binned_get_instance_private (gebin);

  self->bins   = NULL;
  self->binobs = ncm_obj_array_new ();
  self->r      = 0.0;
  self->twolnN = 0.0;
  self->len    = 0;
}

static void
_nc_galaxy_wl_ellipticity_binned_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityBinned *gebin = NC_GALAXY_WL_ELLIPTICITY_BINNED (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_BINNED (object));

  switch (prop_id)
  {
    case PROP_BINOBS:
      nc_galaxy_wl_ellipticity_binned_set_binobs (gebin, g_value_get_object (value), g_value_get_object (value));
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
    case PROP_BINS: 
      g_value_set_object (value, nc_galaxy_wl_ellipticity_binned_peek_bins (gebin));
      break;
    case PROP_BINOBS:
      g_value_set_object (value, nc_galaxy_wl_ellipticity_binned_peek_binobs (gebin));
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

  ncm_obj_array_clear (&self->binobs);
  ncm_vector_clear (&self->bins);
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_binned_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_ellipticity_binned_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_binned_parent_class)->finalize (object);
}

// static void _nc_galaxy_wl_ellipticity_binned_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
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
   * NcGalaxyWLEllipticityBinned:binobs:
   *
   * Array with observables matrices for each bin.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_BINOBS,
                                   g_param_spec_boxed ("binobs",
                                                        NULL,
                                                        "Array with observables matrices for each bin",
                                                        NCM_TYPE_OBJ_ARRAY,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  wl_dist_class->m2lnP = &_nc_galaxy_wl_ellipticity_binned_m2lnP;
  wl_dist_class->gen   = &_nc_galaxy_wl_ellipticity_binned_gen;
  wl_dist_class->len   = &_nc_galaxy_wl_ellipticity_binned_len;
}

// static void
// _nc_galaxy_wl_ellipticity_binned_m2lnP_initial_prep (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
// {
// }

static gdouble
_nc_galaxy_wl_ellipticity_binned_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  return 0;
}

static gdouble
_nc_galaxy_wl_ellipticity_binned_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  return 0;
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
 * nc_galaxy_wl_ellipticity_binned_set_binobs:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 * @obs: a #NcmMatrix
 * @bins: a #NcmVector
 *
 * Sets the array with observables matrices for each bin @binobs.
 */
void
nc_galaxy_wl_ellipticity_binned_set_binobs (NcGalaxyWLEllipticityBinned *gebin, NcmMatrix *obs, NcmVector *bins)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;
  gint bin_i;

  for (bin_i = 0; bin_i < ncm_vector_len (bins) - 1; bin_i++)
  {
    NcmMatrix *bin_data = ncm_matrix_new (0, 3);
    gdouble bin_ll         = ncm_vector_get (bins, bin_i);
    gdouble bin_ul         = ncm_vector_get (bins, bin_i+1);
    gint gal_i;
    gint i;

    for (gal_i = 0; gal_i < ncm_matrix_nrows (obs); gal_i++)
    {
      gdouble r_i   = ncm_matrix_get (obs, gal_i, 0);
      gdouble e_i   = ncm_matrix_get (obs, gal_i, 1);
      gdouble err_i = ncm_matrix_get (obs, gal_i, 2);

      if ((r_i >= bin_ll && r_i < bin_ul && bin_i != ncm_vector_len (bins) - 2) || (r_i >= bin_ll && r_i <= bin_ul && bin_i == ncm_vector_len (bins) - 2))
      {
        NcmMatrix *data = ncm_matrix_new (ncm_matrix_nrows (bin_data)+1, 3);

        ncm_matrix_set (data, ncm_matrix_nrows (data) - 1, 0, r_i);
        ncm_matrix_set (data, ncm_matrix_nrows (data) - 1, 1, e_i);
        ncm_matrix_set (data, ncm_matrix_nrows (data) - 1, 2, err_i);

        if (ncm_matrix_nrows (bin_data) != 0)
        {
          for (i = 0; i < ncm_matrix_nrows (bin_data); i++)
          {
            ncm_matrix_set (data, i, 0, ncm_matrix_get (bin_data, i, 0));
            ncm_matrix_set (data, i, 1, ncm_matrix_get (bin_data, i, 1));
            ncm_matrix_set (data, i, 2, ncm_matrix_get (bin_data, i, 2));
          }
        }
        
        ncm_matrix_substitute (&bin_data, data, FALSE);

      }
    }

    ncm_obj_array_add (self->binobs, G_OBJECT (bin_data));
  }

  self->bins = bins;
  self->len  = ncm_matrix_nrows (obs);

}

/**
 * nc_galaxy_wl_ellipticity_binned_peek_binobs:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 *
 * Gets the array with observables matrices for each bin.
 *
 * Returns: (transfer none): the array with observables matrices for each bin.
 */
NcmObjArray *
nc_galaxy_wl_ellipticity_binned_peek_binobs (NcGalaxyWLEllipticityBinned *gebin)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  return self->binobs;
}

/**
 * nc_galaxy_wl_ellipticity_binned_peek_bins:
 * @gebin: a #NcGalaxyWLEllipticityBinned
 *
 * Gets the vector with bins limits.
 *
 * Returns: (transfer none): the vector with bins limits.
 */
NcmVector *
nc_galaxy_wl_ellipticity_binned_peek_bins (NcGalaxyWLEllipticityBinned *gebin)
{
  NcGalaxyWLEllipticityBinnedPrivate * const self = gebin->priv;

  return self->bins;
}
