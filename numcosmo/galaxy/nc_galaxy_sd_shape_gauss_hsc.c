/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape_gauss_hsc.c
 *
 *  Sun Jan 5 12:57:50 2025
 *  Copyright  2025  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_gauss_hsc.c
 * Copyright (C) 2025 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxySDShapeGaussHSC:
 *
 * Class describing a galaxy sample shape distribution with a truncated gaussian p.d.f.
 * convoluted with gaussian noise and accounting for bias. Compatible with Subaru's Hyper
 * Suprime-Cam (HSC) survey data.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape_gauss_hsc.h"
#include "galaxy/nc_galaxy_sd_shape_gauss.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "lss/nc_halo_position.h"
#include "math/ncm_stats_vec.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */


typedef struct _NcGalaxySDShapeGaussHSCPrivate
{
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_hp;
  NcmStatsVec *obs_stats;
} NcGalaxySDShapeGaussHSCPrivate;

struct _NcGalaxySDShapeGaussHSC
{
  NcGalaxySDShape parent_instance;
};

typedef struct _NcGalaxySDShapeGaussHSCData
{
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble std_shape;
  gdouble std_noise;
  gdouble c1;
  gdouble c2;
  gdouble m;
  gdouble sigma;
  gdouble radius;
  gdouble phi;
  NcWLSurfaceMassDensityOptzs optzs;
} NcGalaxySDShapeGaussHSCData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShapeGaussHSC, nc_galaxy_sd_shape_gauss_hsc, NC_TYPE_GALAXY_SD_SHAPE);

static void
nc_galaxy_sd_shape_gauss_hsc_init (NcGalaxySDShapeGaussHSC *gsdshsc)
{
  NcGalaxySDShapeGaussHSCPrivate * const self = nc_galaxy_sd_shape_gauss_hsc_get_instance_private (gsdshsc);

  self->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  self->ctrl_hp    = ncm_model_ctrl_new (NULL);
  self->obs_stats  = ncm_stats_vec_new (6, NCM_STATS_VEC_COV, FALSE);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_gauss_hsc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGaussHSC *gsdshsc = NC_GALAXY_SD_SHAPE_GAUSS_HSC (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (gsdshsc));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_gauss_hsc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGaussHSC *gsdshsc = NC_GALAXY_SD_SHAPE_GAUSS_HSC (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (gsdshsc));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_shape_gauss_hsc_dispose (GObject *object)
{
  NcGalaxySDShapeGaussHSC *gsdshsc            = NC_GALAXY_SD_SHAPE_GAUSS_HSC (object);
  NcGalaxySDShapeGaussHSCPrivate * const self = nc_galaxy_sd_shape_gauss_hsc_get_instance_private (gsdshsc);

  ncm_model_ctrl_clear (&self->ctrl_cosmo);
  ncm_model_ctrl_clear (&self->ctrl_hp);
  ncm_stats_vec_clear (&self->obs_stats);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_hsc_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_shape_gauss_hsc_finalize (GObject *object)
{
  /* NcGalaxySDShapeGaussHSC *gsdshsc          = NC_GALAXY_SD_SHAPE_GAUSS_HSC (object); */
  /* NcGalaxySDShapeGaussHSCPrivate * const self = nc_galaxy_sd_shape_gauss_hsc_get_instance_private (gsdshsc); */

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_hsc_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_gauss_hsc_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng);
static NcGalaxySDShapeIntegrand *_nc_galaxy_sd_shape_gauss_hsc_integ (NcGalaxySDShape *gsds);
static gboolean _nc_galaxy_sd_shape_gauss_hsc_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array);
static void _nc_galaxy_sd_shape_gauss_hsc_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data);
static void _nc_galaxy_sd_shape_gauss_hsc_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho);

static void
nc_galaxy_sd_shape_gauss_hsc_class_init (NcGalaxySDShapeGaussHSCClass *klass)
{
  NcGalaxySDShapeClass *sd_shape_class = NC_GALAXY_SD_SHAPE_CLASS (klass);
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_gauss_hsc_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_gauss_hsc_get_property;
  object_class->dispose     = &_nc_galaxy_sd_shape_gauss_hsc_dispose;
  object_class->finalize    = &_nc_galaxy_sd_shape_gauss_hsc_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian galaxy shape distribution for HSC data", "Gaussian shape for HSC data");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_model_class_check_params_info (model_class);

  sd_shape_class->gen                = &_nc_galaxy_sd_shape_gauss_hsc_gen;
  sd_shape_class->integ              = &_nc_galaxy_sd_shape_gauss_hsc_integ;
  sd_shape_class->prepare_data_array = &_nc_galaxy_sd_shape_gauss_hsc_prepare_data_array;
  sd_shape_class->data_init          = &_nc_galaxy_sd_shape_gauss_hsc_data_init;
  sd_shape_class->direct_estimate    = &_nc_galaxy_sd_shape_gauss_hsc_direct_estimate;
}

static complex double
_gauss_cut_gen (NcmRNG *rng, const gdouble sigma)
{
  gdouble x;
  gdouble y;

  do {
    x = ncm_rng_gaussian_gen (rng, 0.0, sigma);
    y = ncm_rng_gaussian_gen (rng, 0.0, sigma);
  } while (hypot (x, y) > 1.0);

  return x + I * y;
}

static void
_nc_galaxy_sd_shape_gauss_hsc_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  NcGalaxySDShapeGaussHSCData *ldata           = (NcGalaxySDShapeGaussHSCData *) data->ldata;
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  const gdouble ra                             = data->sdpos_data->ra;
  const gdouble dec                            = data->sdpos_data->dec;
  const gdouble z                              = data->sdpos_data->sdz_data->z;
  const gdouble sigma                          = ldata->sigma;
  const gdouble noise1                         = ncm_rng_gaussian_gen (rng, 0.0, ldata->std_noise);
  const gdouble noise2                         = ncm_rng_gaussian_gen (rng, 0.0, ldata->std_noise);
  const gdouble c1                             = ldata->c1;
  const gdouble c2                             = ldata->c2;
  const gdouble m                              = ldata->m;
  complex double e_s                           = _gauss_cut_gen (rng, sigma);
  complex double noise                         = noise1 + I * noise2;
  complex double c                             = c1 + I * c2;
  gdouble theta                                = 0.0;
  gdouble phi                                  = 0.0;
  complex double e_o;
  gdouble radius;

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);

  nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

  if (data->coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
  {
    phi   = M_PI - phi;
    e_s   = conj (e_s);
    noise = conj (noise);
  }

  radius = nc_halo_position_projected_radius (halo_position, cosmo, theta);

  if (z > z_cl)
  {
    const gdouble gt = nc_wl_surface_mass_density_reduced_shear (surface_mass_density,
                                                                 density_profile,
                                                                 cosmo,
                                                                 radius, z, z_cl, z_cl);

    complex double g    = gt * cexp (2.0 * I * phi);
    NcmComplex cplx_E_s = NCM_COMPLEX_INIT (e_s);
    NcmComplex cplx_E_o, cplx_g;

    /* Adding bias */
    g = (1.0 + m) * g + c;

    ncm_complex_set_c (&cplx_g, g);
    nc_galaxy_sd_shape_apply_shear (gsds, &cplx_g, &cplx_E_s, &cplx_E_o);

    e_o = ncm_complex_c (&cplx_E_o);
  }
  else
  {
    e_o = e_s + c;
  }

  e_o += noise;

  data->epsilon_int_1  = creal (e_s);
  data->epsilon_int_2  = cimag (e_s);
  ldata->epsilon_obs_1 = creal (e_o);
  ldata->epsilon_obs_2 = cimag (e_o);
  ldata->radius        = radius;
  ldata->phi           = phi;
}

struct _IntegData
{
  NcGalaxySDShapeGaussHSC *gsdshsc;
  NcGalaxySDShapeData *data;
  NcHICosmo *cosmo;
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcHaloDensityProfile *density_profile;
};

/* LCOV_EXCL_START */
static gpointer
_integ_data_copy (gpointer user_data)
{
  struct _IntegData *new_int_data = g_new0 (struct _IntegData, 1);
  struct _IntegData *int_data     = (struct _IntegData *) user_data;

  new_int_data->gsdshsc              = int_data->gsdshsc;
  new_int_data->data                 = int_data->data;
  new_int_data->cosmo                = nc_hicosmo_ref (int_data->cosmo);
  new_int_data->halo_position        = nc_halo_position_ref (int_data->halo_position);
  new_int_data->surface_mass_density = nc_wl_surface_mass_density_ref (int_data->surface_mass_density);
  new_int_data->density_profile      = nc_halo_density_profile_ref (int_data->density_profile);

  return new_int_data;
}
/* LCOV_EXCL_STOP */

static void
_integ_data_free (gpointer user_data)
{
  struct _IntegData *int_data = (struct _IntegData *) user_data;

  nc_hicosmo_free (int_data->cosmo);
  nc_halo_position_free (int_data->halo_position);
  nc_wl_surface_mass_density_free (int_data->surface_mass_density);
  nc_halo_density_profile_free (int_data->density_profile);

  g_free (int_data);
}

static gdouble
_nc_galaxy_sd_shape_gauss_hsc_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDShapeData *data)
{
  struct _IntegData *int_data        = (struct _IntegData *) callback_data;
  NcGalaxySDShape *gsds              = NC_GALAXY_SD_SHAPE (int_data->gsdshsc);
  NcGalaxySDShapeGaussHSCData *ldata = (NcGalaxySDShapeGaussHSCData *) data->ldata;
  gdouble z_cl                       = nc_halo_position_get_redshift (int_data->halo_position);
  gdouble phi                        = ldata->phi;
  gdouble e1                         = ldata->epsilon_obs_1;
  gdouble e2                         = ldata->epsilon_obs_2;
  gdouble std_noise                  = ldata->std_noise;
  gdouble c1                         = ldata->c1;
  gdouble c2                         = ldata->c2;
  gdouble m                          = ldata->m;
  gdouble sigma                      = ldata->sigma;
  complex double e_o                 = e1 + I * e2;
  NcmComplex cplx_g, cplx_E_o, cplx_E_s;

  if (z > z_cl)
  {
    const gdouble gt = nc_wl_surface_mass_density_reduced_shear_optzs (int_data->surface_mass_density,
                                                                       int_data->density_profile,
                                                                       int_data->cosmo,
                                                                       z, z_cl, &ldata->optzs);

    complex double g = gt * cexp (2.0 * I * phi);

    /* Adding bias */
    g = (1.0 + m) * g + (c1 + I * c2);

    ncm_complex_set_c (&cplx_g, g);
    ncm_complex_set_c (&cplx_E_o, e_o);

    nc_galaxy_sd_shape_apply_shear_inv (gsds, &cplx_g, &cplx_E_o, &cplx_E_s);
  }
  else
  {
    complex double e_s = e_o - (c1 + I * c2);

    ncm_complex_set_zero (&cplx_g);
    ncm_complex_set_c (&cplx_E_o, e_o);
    ncm_complex_set_c (&cplx_E_s, e_s);
  }

  /* TODO: compute the actual convolution */
  {
    const gdouble var_int   = gsl_pow_2 (sigma);
    const gdouble total_var = var_int + gsl_pow_2 (std_noise);
    const gdouble chi2_1    = gsl_pow_2 (ncm_complex_Re (&cplx_E_s)) / total_var;
    const gdouble chi2_2    = gsl_pow_2 (ncm_complex_Im (&cplx_E_s)) / total_var;
    const gdouble lndetjac  = nc_galaxy_sd_shape_lndet_jac (gsds, &cplx_g, &cplx_E_o);

    return exp (-0.5 * (chi2_1 + chi2_2) + lndetjac) / (2.0 * M_PI * total_var);
  }
}

static void
_integ_data_prepare (gpointer user_data, NcmMSet *mset)
{
  struct _IntegData *int_data                  = (struct _IntegData *) user_data;
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));

  g_assert_nonnull (cosmo);
  g_assert_nonnull (halo_position);
  g_assert_nonnull (surface_mass_density);
  g_assert_nonnull (density_profile);

  nc_hicosmo_clear (&int_data->cosmo);
  nc_halo_position_clear (&int_data->halo_position);
  nc_wl_surface_mass_density_clear (&int_data->surface_mass_density);
  nc_halo_density_profile_clear (&int_data->density_profile);

  int_data->cosmo                = nc_hicosmo_ref (cosmo);
  int_data->halo_position        = nc_halo_position_ref (halo_position);
  int_data->surface_mass_density = nc_wl_surface_mass_density_ref (surface_mass_density);
  int_data->density_profile      = nc_halo_density_profile_ref (density_profile);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);
}

static NcGalaxySDShapeIntegrand *
_nc_galaxy_sd_shape_gauss_hsc_integ (NcGalaxySDShape *gsds)
{
  NcGalaxySDShapeGaussHSC *gsdshsc = NC_GALAXY_SD_SHAPE_GAUSS_HSC (gsds);
  struct _IntegData *int_data      = g_new0 (struct _IntegData, 1);
  NcGalaxySDShapeIntegrand *integ  = nc_galaxy_sd_shape_integrand_new (_nc_galaxy_sd_shape_gauss_hsc_integ_f,
                                                                       _integ_data_free,
                                                                       _integ_data_copy,
                                                                       _integ_data_prepare,
                                                                       int_data);

  int_data->gsdshsc = gsdshsc;

  return integ;
}

static gboolean
_nc_galaxy_sd_shape_gauss_hsc_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array)
{
  NcGalaxySDShapeGaussHSC *gsdsghsc            = NC_GALAXY_SD_SHAPE_GAUSS_HSC (gsds);
  NcGalaxySDShapeGaussHSCPrivate * const self  = nc_galaxy_sd_shape_gauss_hsc_get_instance_private (gsdsghsc);
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  gdouble z_cl                                 = nc_halo_position_get_redshift (halo_position);
  gboolean update_radius;
  guint i;

  update_radius = (ncm_model_ctrl_update (self->ctrl_hp, NCM_MODEL (halo_position))) ||
                  (ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo)));

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxySDShapeData *data_i          = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeGaussHSCData *ldata_i = (NcGalaxySDShapeGaussHSCData *) data_i->ldata;

    if ((update_radius) || (ldata_i->radius == 0.0))
    {
      const gdouble ra  = data_i->sdpos_data->ra;
      const gdouble dec = data_i->sdpos_data->dec;
      gdouble theta, phi;

      nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

      if (data_i->coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      ldata_i->radius = nc_halo_position_projected_radius (halo_position, cosmo, theta);
      ldata_i->phi    = phi;
    }

    nc_wl_surface_mass_density_reduced_shear_optzs_prep (surface_mass_density,
                                                         density_profile,
                                                         cosmo,
                                                         ldata_i->radius,
                                                         z_cl,
                                                         z_cl,
                                                         &ldata_i->optzs);
  }

  return TRUE;
}

static void
_nc_galaxy_sd_shape_gauss_hsc_ldata_free (gpointer ldata)
{
  g_free (ldata);
}

static void
_nc_galaxy_sd_shape_gauss_hsc_ldata_read_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDShapeGaussHSCData *ldata = (NcGalaxySDShapeGaussHSCData *) data->ldata;

  nc_galaxy_sd_position_data_read_row (data->sdpos_data, obs, i);

  ldata->epsilon_obs_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_1, i);
  ldata->epsilon_obs_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_2, i);
  ldata->std_shape     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
  ldata->sigma         = nc_galaxy_sd_shape_gauss_sigma_from_std_shape (ldata->std_shape);
  ldata->std_noise     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_NOISE, i);
  ldata->c1            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
  ldata->c2            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
  ldata->m             = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);
}

static void
_nc_galaxy_sd_shape_gauss_hsc_ldata_write_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDShapeGaussHSCData *ldata = (NcGalaxySDShapeGaussHSCData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_1, i, ldata->epsilon_obs_1);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_2, i, ldata->epsilon_obs_2);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i, ldata->std_shape);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_NOISE, i, ldata->std_noise);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i, ldata->c1);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i, ldata->c2);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i, ldata->m);
}

static void
_nc_galaxy_sd_shape_gauss_hsc_ldata_required_columns (NcGalaxySDShapeData *data, GList *columns)
{
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_NOISE));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M));
}

static gdouble
_nc_galaxy_sd_shape_gauss_hsc_ldata_get_radius (NcGalaxySDShapeData *data)
{
  NcGalaxySDShapeGaussHSCData *ldata = (NcGalaxySDShapeGaussHSCData *) data->ldata;

  return ldata->radius;
}

static void
_nc_galaxy_sd_shape_gauss_hsc_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data)
{
  NcGalaxySDShapeGaussHSCData *ldata = g_new0 (NcGalaxySDShapeGaussHSCData, 1);

  data->sdpos_data             = sdpos_data;
  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_shape_gauss_hsc_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_shape_gauss_hsc_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_shape_gauss_hsc_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_shape_gauss_hsc_ldata_required_columns;
  data->ldata_get_radius       = &_nc_galaxy_sd_shape_gauss_hsc_ldata_get_radius;
}

static void
_nc_galaxy_sd_shape_gauss_hsc_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  NcGalaxySDShapeGaussHSC *gsdshsc            = NC_GALAXY_SD_SHAPE_GAUSS_HSC (gsds);
  NcGalaxySDShapeGaussHSCPrivate * const self = nc_galaxy_sd_shape_gauss_hsc_get_instance_private (gsdshsc);
  NcHICosmo *cosmo                            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position               = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxyWLObsEllipConv ellip_conv           = nc_galaxy_sd_shape_get_ellip_conv (gsds);
  guint i;

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  ncm_stats_vec_reset (self->obs_stats, TRUE);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxySDShapeData *data_i          = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeGaussHSCData *ldata_i = (NcGalaxySDShapeGaussHSCData *) data_i->ldata;
    const gdouble ra                     = data_i->sdpos_data->ra;
    const gdouble dec                    = data_i->sdpos_data->dec;
    const gdouble e1                     = ldata_i->epsilon_obs_1;
    const gdouble e2                     = ldata_i->epsilon_obs_2;
    const gdouble std_noise              = ldata_i->std_noise;
    const gdouble m                      = ldata_i->m;
    const gdouble c1                     = ldata_i->c1;
    const gdouble c2                     = ldata_i->c2;
    const gdouble std_shape              = ldata_i->std_shape;
    const gdouble var_tot                = std_shape * std_shape + std_noise * std_noise;
    const gdouble weight                 = 1.0 / var_tot;
    complex double e_o                   = e1 + I * e2;
    complex double hat_g;
    gdouble theta, phi;

    nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

    if (data_i->coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
      phi = M_PI - phi;

    hat_g = e_o * cexp (-2.0 * I * phi);

    ncm_stats_vec_set (self->obs_stats, 0, creal (hat_g));
    ncm_stats_vec_set (self->obs_stats, 1, cimag (hat_g));
    ncm_stats_vec_set (self->obs_stats, 2, std_shape * std_shape);
    ncm_stats_vec_set (self->obs_stats, 3, m);
    ncm_stats_vec_set (self->obs_stats, 4, c1);
    ncm_stats_vec_set (self->obs_stats, 5, c2);

    ncm_stats_vec_update_weight (self->obs_stats, weight);
  }

  {
    const gdouble mean_gt          = ncm_stats_vec_get_mean (self->obs_stats, 0);
    const gdouble mean_gx          = ncm_stats_vec_get_mean (self->obs_stats, 1);
    const gdouble mean_sigma_true2 = ncm_stats_vec_get_mean (self->obs_stats, 2);
    const gdouble mean_m           = ncm_stats_vec_get_mean (self->obs_stats, 3);
    const gdouble mean_c1          = ncm_stats_vec_get_mean (self->obs_stats, 4);
    const gdouble mean_c2          = ncm_stats_vec_get_mean (self->obs_stats, 5);
    const gdouble R                = 1.0 - mean_sigma_true2;

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        *gt      = (mean_gt - mean_c1) / (1.0 + mean_m);
        *gx      = (mean_gx - mean_c2) / (1.0 + mean_m);
        *sigma_t = ncm_stats_vec_get_sd (self->obs_stats, 0) / (1.0 + mean_m);
        *sigma_x = ncm_stats_vec_get_sd (self->obs_stats, 1) / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        *gt      = (0.5 * mean_gt / R - mean_c1) / (1.0 + mean_m);
        *gx      = (0.5 * mean_gx / R - mean_c2) / (1.0 + mean_m);
        *sigma_t = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 0) / R / (1.0 + mean_m);
        *sigma_x = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 1) / R / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }
  }
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxySDShapeGaussHSC
 *
 * Returns: (transfer full): a new NcGalaxySDShapeGaussHSC.
 */
NcGalaxySDShapeGaussHSC *
nc_galaxy_sd_shape_gauss_hsc_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxySDShapeGaussHSC *gsdshsc = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_GAUSS_HSC,
                                                   "ellip-conv", ellip_conv,
                                                   NULL);

  return gsdshsc;
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_ref:
 * @gsdshsc: a #NcGalaxySDShapeGaussHSC
 *
 * Increase the reference of @gsdshsc by one.
 *
 * Returns: (transfer full): @gsdshsc.
 */
NcGalaxySDShapeGaussHSC *
nc_galaxy_sd_shape_gauss_hsc_ref (NcGalaxySDShapeGaussHSC *gsdshsc)
{
  return g_object_ref (gsdshsc);
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_free:
 * @gsdshsc: a #NcGalaxySDShapeGaussHSC
 *
 * Decrease the reference count of @gsdshsc by one.
 *
 */
void
nc_galaxy_sd_shape_gauss_hsc_free (NcGalaxySDShapeGaussHSC *gsdshsc)
{
  g_object_unref (gsdshsc);
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_clear:
 * @gsdshsc: a #NcGalaxySDShapeGaussHSC
 *
 * Decrease the reference count of @gsdshsc by one, and sets the pointer *@gsdshsc to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_gauss_hsc_clear (NcGalaxySDShapeGaussHSC **gsdshsc)
{
  g_clear_object (gsdshsc);
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_gen:
 * @gsdshsc: a #NcGalaxySDShapeGaussHSC
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDShapeData
 * @std_shape: the intrinsic shape dispersion
 * @std_noise: the observational shape dispersion
 * @c1: the first additive bias parameter
 * @c2: the second additive bias parameter
 * @m: the multiplicative bias parameter
 * @coord: the coordinate system #NcGalaxyWLObsCoord
 * @rng: a #NcmRNG
 *
 * Generates a galaxy sample shape.
 *
 */
void
nc_galaxy_sd_shape_gauss_hsc_gen (NcGalaxySDShapeGaussHSC *gsdshsc, NcmMSet *mset, NcGalaxySDShapeData *data, const gdouble std_shape, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcGalaxyWLObsCoord coord, NcmRNG *rng)
{
  NcGalaxySDShapeClass *sd_shape_class = NC_GALAXY_SD_SHAPE_GET_CLASS (gsdshsc);
  NcGalaxySDShapeGaussHSCData *ldata   = (NcGalaxySDShapeGaussHSCData *) data->ldata;

  data->coord      = coord;
  ldata->std_shape = std_shape;
  ldata->std_noise = std_noise;
  ldata->sigma     = nc_galaxy_sd_shape_gauss_sigma_from_std_shape (std_shape);
  ldata->c1        = c1;
  ldata->c2        = c2;
  ldata->m         = m;

  sd_shape_class->gen (NC_GALAXY_SD_SHAPE (gsdshsc), mset, data, rng);
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_data_set:
 * @gsdshsc: a #NcGalaxySDShapeGaussHSC
 * @data: a #NcGalaxySDShapeData
 * @epsilon_obs_1: the observed ellipticity component 1
 * @epsilon_obs_2: the observed ellipticity component 2
 * @std_shape: the intrinsic shape dispersion
 * @std_noise: the observational shape dispersion
 * @c1: the first additive bias parameter
 * @c2: the second additive bias parameter
 * @m: the multiplicative bias parameter
 *
 * Sets the observed ellipticity components and the observational shape dispersion.
 *
 */
void
nc_galaxy_sd_shape_gauss_hsc_data_set (NcGalaxySDShapeGaussHSC *gsdshsc, NcGalaxySDShapeData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_shape, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m)
{
  NcGalaxySDShapeGaussHSCData *ldata = (NcGalaxySDShapeGaussHSCData *) data->ldata;

  ldata->epsilon_obs_1 = epsilon_obs_1;
  ldata->epsilon_obs_2 = epsilon_obs_2;
  ldata->std_shape     = std_shape;
  ldata->sigma         = nc_galaxy_sd_shape_gauss_sigma_from_std_shape (std_shape);
  ldata->std_noise     = std_noise;
  ldata->c1            = c1;
  ldata->c2            = c2;
  ldata->m             = m;
}

/**
 * nc_galaxy_sd_shape_gauss_hsc_data_get:
 * @gsdshsc: a #NcGalaxySDShapeGaussHSC
 * @data: a #NcGalaxySDShapeData
 * @epsilon_obs_1: (out): the observed ellipticity component 1
 * @epsilon_obs_2: (out): the observed ellipticity component 2
 * @std_shape: (out): the intrinsic shape dispersion
 * @std_noise: (out): the observational shape dispersion
 * @c1: (out): the first additive bias parameter
 * @c2: (out): the second additive bias parameter
 * @m: (out): the multiplicative bias parameter
 *
 * Gets the observed ellipticity components and the observational shape dispersion.
 *
 */
void
nc_galaxy_sd_shape_gauss_hsc_data_get (NcGalaxySDShapeGaussHSC *gsdshsc, NcGalaxySDShapeData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_shape, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m)
{
  NcGalaxySDShapeGaussHSCData *ldata = (NcGalaxySDShapeGaussHSCData *) data->ldata;

  *epsilon_obs_1 = ldata->epsilon_obs_1;
  *epsilon_obs_2 = ldata->epsilon_obs_2;
  *std_shape     = ldata->std_shape;
  *std_noise     = ldata->std_noise;
  *c1            = ldata->c1;
  *c2            = ldata->c2;
  *m             = ldata->m;
}

