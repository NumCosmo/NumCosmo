/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape_gauss.c
 *
 *  Mon June 5 14:56:41 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_gauss.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxySDShapeGauss:
 *
 * Class describing a galaxy sample shape distribution with a truncated gaussian p.d.f.
 * convoluted with gaussian noise.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape_gauss.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "lss/nc_halo_position.h"
#include "math/ncm_stats_vec.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */


typedef struct _NcGalaxySDShapeGaussPrivate
{
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_hp;
  NcmModelCtrl *ctrl_shape;
  NcmStatsVec *obs_stats;
  gdouble std_shape;
} NcGalaxySDShapeGaussPrivate;

struct _NcGalaxySDShapeGauss
{
  NcGalaxySDShape parent_instance;
};

typedef struct _NcGalaxySDShapeGaussData
{
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble std_noise;
  gdouble radius;
  gdouble phi;
  NcWLSurfaceMassDensityOptzs optzs;
} NcGalaxySDShapeGaussData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShapeGauss, nc_galaxy_sd_shape_gauss, NC_TYPE_GALAXY_SD_SHAPE);

static void
nc_galaxy_sd_shape_gauss_init (NcGalaxySDShapeGauss *gsdsgauss)
{
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  self->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  self->ctrl_hp    = ncm_model_ctrl_new (NULL);
  self->ctrl_shape = ncm_model_ctrl_new (NULL);
  self->obs_stats  = ncm_stats_vec_new (3, NCM_STATS_VEC_COV, FALSE);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss = NC_GALAXY_SD_SHAPE_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsgauss));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss = NC_GALAXY_SD_SHAPE_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsgauss));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_shape_gauss_dispose (GObject *object)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  ncm_model_ctrl_clear (&self->ctrl_cosmo);
  ncm_model_ctrl_clear (&self->ctrl_hp);
  ncm_model_ctrl_clear (&self->ctrl_shape);
  ncm_stats_vec_clear (&self->obs_stats);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_shape_gauss_finalize (GObject *object)
{
  /* NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object); */
  /* NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss); */

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng);
static NcGalaxySDShapeIntegrand *_nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds);
static gboolean _nc_galaxy_sd_shape_gauss_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array);
static void _nc_galaxy_sd_shape_gauss_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data);
static void _nc_galaxy_sd_shape_gauss_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho);
static void _nc_galaxy_sd_shape_gauss_update_std_shape (NcGalaxySDShape *gsds);

static void
nc_galaxy_sd_shape_gauss_class_init (NcGalaxySDShapeGaussClass *klass)
{
  NcGalaxySDShapeClass *sd_shape_class = NC_GALAXY_SD_SHAPE_CLASS (klass);
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_gauss_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_sd_shape_gauss_dispose;
  object_class->finalize    = &_nc_galaxy_sd_shape_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian galaxy shape distribution", "Gaussian shape");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_SHAPE_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDShapeGauss:sigma:
   *
   * The $\sigma$ parameter of the truncated gaussian distribution.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_SIGMA,
                              "\\sigma",
                              "sigma", 0.0, 1.0, 1.0e-1,
                              NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_SIGMA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  sd_shape_class->gen                = &_nc_galaxy_sd_shape_gauss_gen;
  sd_shape_class->integ              = &_nc_galaxy_sd_shape_gauss_integ;
  sd_shape_class->prepare_data_array = &_nc_galaxy_sd_shape_gauss_prepare_data_array;
  sd_shape_class->data_init          = &_nc_galaxy_sd_shape_gauss_data_init;
  sd_shape_class->direct_estimate    = &_nc_galaxy_sd_shape_gauss_direct_estimate;
}

#define VECTOR  (NCM_MODEL (gsds))
#define SIGMA   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_SHAPE_GAUSS_SIGMA))

static void
_nc_galaxy_sd_shape_gauss_update_std_shape (NcGalaxySDShape *gsds)
{
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (NC_GALAXY_SD_SHAPE_GAUSS (gsds));

  if (ncm_model_ctrl_update (self->ctrl_shape, NCM_MODEL (gsds)))
    self->std_shape = nc_galaxy_sd_shape_gauss_std_shape_from_sigma (SIGMA);
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
_nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  NcGalaxySDShapeGaussData *ldata              = (NcGalaxySDShapeGaussData *) data->ldata;
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  const gdouble ra                             = data->sdpos_data->ra;
  const gdouble dec                            = data->sdpos_data->dec;
  const gdouble z                              = data->sdpos_data->sdz_data->z;
  const gdouble sigma                          = SIGMA;
  const gdouble noise1                         = ncm_rng_gaussian_gen (rng, 0.0, ldata->std_noise);
  const gdouble noise2                         = ncm_rng_gaussian_gen (rng, 0.0, ldata->std_noise);
  complex double e_s                           = _gauss_cut_gen (rng, sigma);
  complex double noise                         = noise1 + I * noise2;
  gdouble theta                                = 0.0;
  gdouble phi                                  = 0.0;
  complex double e_o;
  gdouble e1, e2, e1_int, e2_int;
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

    NcmComplex cplx_g   = NCM_COMPLEX_INIT_REAL (gt);
    NcmComplex cplx_E_s = NCM_COMPLEX_INIT (e_s);
    NcmComplex cplx_E_o;

    nc_galaxy_sd_shape_apply_shear (gsds, &cplx_g, &cplx_E_s, &cplx_E_o);

    e_o = ncm_complex_c (&cplx_E_o);
  }
  else
  {
    e_o = e_s;
  }

  e_o = e_o * cexp (2.0 * I * phi) + noise;
  e_s = e_s * cexp (2.0 * I * phi);

  e1 = creal (e_o);
  e2 = cimag (e_o);

  e1_int = creal (e_s);
  e2_int = cimag (e_s);

  data->epsilon_int_1  = e1_int;
  data->epsilon_int_2  = e2_int;
  ldata->epsilon_obs_1 = e1;
  ldata->epsilon_obs_2 = e2;
  ldata->radius        = radius;
  ldata->phi           = phi;
}

struct _IntegData
{
  NcGalaxySDShapeGauss *gsdsgauss;
  NcGalaxySDShapeData *data;
  NcHICosmo *cosmo;
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcHaloDensityProfile *density_profile;
};

static gpointer
_integ_data_copy (gpointer user_data)
{
  struct _IntegData *new_int_data = g_new0 (struct _IntegData, 1);
  struct _IntegData *int_data     = (struct _IntegData *) user_data;

  new_int_data->gsdsgauss            = int_data->gsdsgauss;
  new_int_data->data                 = int_data->data;
  new_int_data->cosmo                = nc_hicosmo_ref (int_data->cosmo);
  new_int_data->halo_position        = nc_halo_position_ref (int_data->halo_position);
  new_int_data->surface_mass_density = nc_wl_surface_mass_density_ref (int_data->surface_mass_density);
  new_int_data->density_profile      = nc_halo_density_profile_ref (int_data->density_profile);

  return new_int_data;
}

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
_nc_galaxy_sd_shape_gauss_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDShapeData *data)
{
  struct _IntegData *int_data     = (struct _IntegData *) callback_data;
  NcGalaxySDShape *gsds           = NC_GALAXY_SD_SHAPE (int_data->gsdsgauss);
  NcGalaxySDShapeGaussData *ldata = (NcGalaxySDShapeGaussData *) data->ldata;
  gdouble z_cl                    = nc_halo_position_get_redshift (int_data->halo_position);
  gdouble phi                     = ldata->phi;
  gdouble e1                      = ldata->epsilon_obs_1;
  gdouble e2                      = ldata->epsilon_obs_2;
  gdouble sigma                   = SIGMA;
  gdouble std_noise               = ldata->std_noise;
  complex double e_o              = e1 + I * e2;
  NcmComplex cplx_g, cplx_E_o, cplx_E_s;

  if (z > z_cl)
  {
    const gdouble gt = nc_wl_surface_mass_density_reduced_shear_optzs (int_data->surface_mass_density,
                                                                       int_data->density_profile,
                                                                       int_data->cosmo,
                                                                       z, z_cl, &ldata->optzs);

    complex double g = gt * cexp (2.0 * I * phi);

    ncm_complex_set_c (&cplx_g, g);
    ncm_complex_set_c (&cplx_E_o, e_o);

    nc_galaxy_sd_shape_apply_shear_inv (gsds, &cplx_g, &cplx_E_o, &cplx_E_s);
  }
  else
  {
    ncm_complex_set_zero (&cplx_g);
    ncm_complex_set_c (&cplx_E_o, e_o);
    ncm_complex_set_c (&cplx_E_s, e_o);
  }

  /* TODO: compute the actual convolution */
  {
    const gdouble var_int   = gsl_pow_2 (sigma);
    const gdouble total_var = var_int + gsl_pow_2 (std_noise);
    const gdouble chi2_1    = gsl_pow_2 (ncm_complex_Re (&cplx_E_s)) / total_var;
    const gdouble chi2_2    = gsl_pow_2 (ncm_complex_Im (&cplx_E_s)) / total_var;
    const gdouble lndetjac  = nc_galaxy_sd_shape_lndet_jac (gsds, &cplx_g, &cplx_E_o);
    return -0.5 * (chi2_1 + chi2_2) + lndetjac - log (2.0 * M_PI * total_var);
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
_nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds)
{
  NcGalaxySDShapeGauss *gsdsgauss = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  struct _IntegData *int_data     = g_new0 (struct _IntegData, 1);
  NcGalaxySDShapeIntegrand *integ = nc_galaxy_sd_shape_integrand_new (_nc_galaxy_sd_shape_gauss_integ_f,
                                                                      _integ_data_free,
                                                                      _integ_data_copy,
                                                                      _integ_data_prepare,
                                                                      int_data);

  int_data->gsdsgauss = gsdsgauss;

  return integ;
}

static gboolean
_nc_galaxy_sd_shape_gauss_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array)
{
  NcGalaxySDShapeGauss *gsdsgauss              = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self     = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
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
    NcGalaxySDShapeData *data_i       = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeGaussData *ldata_i = (NcGalaxySDShapeGaussData *) data_i->ldata;

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
_nc_galaxy_sd_shape_gauss_ldata_free (gpointer ldata)
{
  g_free (ldata);
}

static void
_nc_galaxy_sd_shape_gauss_ldata_read_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDShapeGaussData *ldata = (NcGalaxySDShapeGaussData *) data->ldata;

  nc_galaxy_sd_position_data_read_row (data->sdpos_data, obs, i);

  ldata->epsilon_obs_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
  ldata->epsilon_obs_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
  ldata->std_noise     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);
}

static void
_nc_galaxy_sd_shape_gauss_ldata_write_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDShapeGaussData *ldata = (NcGalaxySDShapeGaussData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i, ldata->epsilon_obs_1);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i, ldata->epsilon_obs_2);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i, ldata->std_noise);
}

static void
_nc_galaxy_sd_shape_gauss_ldata_required_columns (NcGalaxySDShapeData *data, GList *columns)
{
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE));
}

static gdouble
_nc_galaxy_sd_shape_gauss_get_radius (NcGalaxySDShapeData *data)
{
  NcGalaxySDShapeGaussData *ldata = (NcGalaxySDShapeGaussData *) data->ldata;

  return ldata->radius;
}

static void
_nc_galaxy_sd_shape_gauss_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data)
{
  NcGalaxySDShapeGaussData *ldata = g_new0 (NcGalaxySDShapeGaussData, 1);

  data->sdpos_data             = sdpos_data;
  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_shape_gauss_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_shape_gauss_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_shape_gauss_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_shape_gauss_ldata_required_columns;
  data->ldata_get_radius       = &_nc_galaxy_sd_shape_gauss_get_radius;
}

static void
_nc_galaxy_sd_shape_gauss_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  NcHICosmo *cosmo                         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position            = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxyWLObsEllipConv ellip_conv        = nc_galaxy_sd_shape_get_ellip_conv (gsds);
  guint i;

  _nc_galaxy_sd_shape_gauss_update_std_shape (gsds);
  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  ncm_stats_vec_reset (self->obs_stats, TRUE);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxySDShapeData *data_i       = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeGaussData *ldata_i = (NcGalaxySDShapeGaussData *) data_i->ldata;
    const gdouble ra                  = data_i->sdpos_data->ra;
    const gdouble dec                 = data_i->sdpos_data->dec;
    const gdouble e1                  = ldata_i->epsilon_obs_1;
    const gdouble e2                  = ldata_i->epsilon_obs_2;
    const gdouble std_noise           = ldata_i->std_noise;
    const gdouble std_shape           = self->std_shape;
    const gdouble var_tot             = std_shape * std_shape + std_noise * std_noise;
    const gdouble weight              = 1.0 / var_tot;
    complex double e_o                = e1 + I * e2;
    complex double hat_g;
    gdouble theta, phi;

    nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

    if (data_i->coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
      phi = M_PI - phi;

    hat_g = e_o * cexp (-2.0 * I * phi);

    ncm_stats_vec_set (self->obs_stats, 0, creal (hat_g));
    ncm_stats_vec_set (self->obs_stats, 1, cimag (hat_g));
    ncm_stats_vec_set (self->obs_stats, 2, std_shape * std_shape);

    ncm_stats_vec_update_weight (self->obs_stats, weight);
  }

  {
    const gdouble mean_gt = ncm_stats_vec_get_mean (self->obs_stats, 0);
    const gdouble mean_gx = ncm_stats_vec_get_mean (self->obs_stats, 1);

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        *gt      = mean_gt;
        *gx      = mean_gx;
        *sigma_t = ncm_stats_vec_get_sd (self->obs_stats, 0);
        *sigma_x = ncm_stats_vec_get_sd (self->obs_stats, 1);
        *rho     = ncm_stats_vec_get_cov (self->obs_stats, 0, 1);
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      {
        const gdouble mean_std_shape2 = ncm_stats_vec_get_mean (self->obs_stats, 2);
        const gdouble R               = 1.0 - mean_std_shape2;

        *gt      = 0.5 * mean_gt / R;
        *gx      = 0.5 * mean_gx / R;
        *sigma_t = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 0) / R;
        *sigma_x = 0.5 * ncm_stats_vec_get_sd (self->obs_stats, 1) / R;
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      }
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }
  }
}

/**
 * nc_galaxy_sd_shape_gauss_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxySDShapeGauss
 *
 * Returns: (transfer full): a new NcGalaxySDShapeGauss.
 */
NcGalaxySDShapeGauss *
nc_galaxy_sd_shape_gauss_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxySDShapeGauss *gsdsgauss = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_GAUSS,
                                                  "ellip-conv", ellip_conv,
                                                  NULL);

  return gsdsgauss;
}

/**
 * nc_galaxy_sd_shape_gauss_ref:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Increase the reference of @gsdsgauss by one.
 *
 * Returns: (transfer full): @gsdsgauss.
 */
NcGalaxySDShapeGauss *
nc_galaxy_sd_shape_gauss_ref (NcGalaxySDShapeGauss *gsdsgauss)
{
  return g_object_ref (gsdsgauss);
}

/**
 * nc_galaxy_sd_shape_gauss_free:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Decrease the reference count of @gsdsgauss by one.
 *
 */
void
nc_galaxy_sd_shape_gauss_free (NcGalaxySDShapeGauss *gsdsgauss)
{
  g_object_unref (gsdsgauss);
}

/**
 * nc_galaxy_sd_shape_gauss_clear:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Decrease the reference count of @gsdsgauss by one, and sets the pointer *@gsdsgauss to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_gauss_clear (NcGalaxySDShapeGauss **gsdsgauss)
{
  g_clear_object (gsdsgauss);
}

/**
 * nc_galaxy_sd_shape_gauss_gen:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDShapeData
 * @std_noise: the observational shape dispersion
 * @coord: the coordinate system #NcGalaxyWLObsCoord
 * @rng: a #NcmRNG
 *
 * Generates a galaxy sample shape.
 *
 */
void
nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShapeGauss *gsdsgauss, NcmMSet *mset, NcGalaxySDShapeData *data, const gdouble std_noise, NcGalaxyWLObsCoord coord, NcmRNG *rng)
{
  NcGalaxySDShapeClass *sd_shape_class = NC_GALAXY_SD_SHAPE_GET_CLASS (gsdsgauss);
  NcGalaxySDShapeGaussData *ldata      = (NcGalaxySDShapeGaussData *) data->ldata;

  data->coord      = coord;
  ldata->std_noise = std_noise;

  sd_shape_class->gen (NC_GALAXY_SD_SHAPE (gsdsgauss), mset, data, rng);
}

/**
 * nc_galaxy_sd_shape_gauss_data_set:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 * @data: a #NcGalaxySDShapeData
 * @epsilon_obs_1: the observed ellipticity component 1
 * @epsilon_obs_2: the observed ellipticity component 2
 * @std_noise: the observational shape dispersion
 *
 * Sets the observed ellipticity components and the gaussian noise dispersion.
 *
 */
void
nc_galaxy_sd_shape_gauss_data_set (NcGalaxySDShapeGauss *gsdsgauss, NcGalaxySDShapeData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise)
{
  NcGalaxySDShapeGaussData *ldata = (NcGalaxySDShapeGaussData *) data->ldata;

  ldata->epsilon_obs_1 = epsilon_obs_1;
  ldata->epsilon_obs_2 = epsilon_obs_2;
  ldata->std_noise     = std_noise;
}

/**
 * nc_galaxy_sd_shape_gauss_data_get:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 * @data: a #NcGalaxySDShapeData
 * @epsilon_obs_1: (out): the observed ellipticity component 1
 * @epsilon_obs_2: (out): the observed ellipticity component 2
 * @std_noise: (out): the observational shape dispersion
 *
 * Gets the observed ellipticity components and the gaussian noise dispersion.
 *
 */
void
nc_galaxy_sd_shape_gauss_data_get (NcGalaxySDShapeGauss *gsdsgauss, NcGalaxySDShapeData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_noise)
{
  NcGalaxySDShapeGaussData *ldata = (NcGalaxySDShapeGaussData *) data->ldata;

  *epsilon_obs_1 = ldata->epsilon_obs_1;
  *epsilon_obs_2 = ldata->epsilon_obs_2;
  *std_noise     = ldata->std_noise;
}

/**
 * nc_galaxy_sd_shape_gauss_sigma_from_std_shape:
 * @std_shape: the physical shape standard deviation
 *
 * Computes the internal variance parameter ($\sigma$) of the Gaussian shape model
 * corresponding to a given observed shape dispersion ($s$).
 *
 * In this model, the observed ellipticity components are assumed to be independent and
 * normally distributed, but constrained to lie within a unit disk. As a result, the
 * variance parameter used in the Gaussian likelihood ($\sigma^2$) does not directly
 * equal the physical shape dispersion ($s^2$). This function solves numerically for the
 * internal parameter that reproduces the given true shape dispersion.
 *
 * Returns: the internal shape dispersion $\sigma$ corresponding to the input $s$.
 */
gdouble
nc_galaxy_sd_shape_gauss_sigma_from_std_shape (const gdouble std_shape)
{
  const gdouble std_shape2 = std_shape * std_shape;
  const gdouble beta       = 0.5 / std_shape2;
  register gdouble alpha   = beta;
  gdouble prev_delta[3]    = {0.0, 0.0, 0.0};
  guint iter               = 0;

  while (TRUE)
  {
    gdouble lambda = 1.0;
    gdouble delta;

    if (alpha < 1.0e2)
    {
      const gdouble exp_alpha_m1 = expm1 (alpha);
      const gdouble exp_alpha    = exp_alpha_m1 + 1.0;
      const gdouble falpha       = (alpha - beta) * exp_alpha_m1 + alpha * beta;
      const gdouble dfalpha      = (alpha - beta) * exp_alpha + exp_alpha_m1 + beta;

      delta = falpha / dfalpha;
    }
    else
    {
      const gdouble exp_malpha_m1 = expm1 (-alpha);
      const gdouble exp_malpha    = exp_malpha_m1 + 1.0;
      const gdouble falpha        = -(alpha - beta) * exp_malpha_m1 + alpha * beta * exp_malpha;
      const gdouble dfalpha       = -(alpha - beta) * exp_malpha - exp_malpha_m1 + beta * exp_malpha;

      delta = falpha / dfalpha;
    }

    if ((fabs (delta / prev_delta[0] - 1.0) < GSL_DBL_EPSILON) ||
        (fabs (delta / prev_delta[1] - 1.0) < GSL_DBL_EPSILON) ||
        (fabs (delta / prev_delta[2] - 1.0) < GSL_DBL_EPSILON))
      lambda *= 0.5;

    lambda *= (1.0 / (1.0 + 1.0e-1 * (iter / 20)));
    delta  *= lambda;

    alpha -= delta;

    if (lambda < 1.0e-2)
    {
      g_warning ("nc_galaxy_sd_shape_gauss_new: sigma_from_std_shape failed to converge"); /* LCOV_EXCL_LINE */
      break;                                                                               /* LCOV_EXCL_LINE */
    }

    if (++iter == 10000)
    {
      g_warning ("nc_galaxy_sd_shape_gauss_new: sigma_from_std_shape failed to converge"); /* LCOV_EXCL_LINE */
      break;                                                                               /* LCOV_EXCL_LINE */
    }

    if (fabs (delta) < GSL_DBL_EPSILON * 10.0)
      return sqrt (0.5 / alpha);

    prev_delta[2] = prev_delta[1];
    prev_delta[1] = prev_delta[0];
    prev_delta[0] = delta;
  }

  return NAN;
}

/**
 * nc_galaxy_sd_shape_gauss_std_shape_from_sigma:
 * @sigma: the internal shape dispersion
 *
 * Computes the physical (true) shape standard deviation $s$ corresponding to a given
 * internal Gaussian model dispersion $\sigma$.
 *
 * In the Gaussian likelihood model used for shape measurements, the variance parameter
 * ($\sigma^2$) applies to ellipticity components constrained to lie within
 * a unit disk. This function analytically maps the internal dispersion back to the true
 * shape dispersion, accounting for the truncation of the Gaussian support.
 *
 * Returns: the true observed shape dispersion $s$ corresponding to the input $\sigma$.
 */
gdouble
nc_galaxy_sd_shape_gauss_std_shape_from_sigma (const gdouble sigma)
{
  gdouble s2      = sigma * sigma;
  gdouble exp_arg = 1.0 / (2.0 * s2);
  gdouble denom   = -2.0 * expm1 (exp_arg);

  return sqrt (s2 + 1.0 / denom);
}

