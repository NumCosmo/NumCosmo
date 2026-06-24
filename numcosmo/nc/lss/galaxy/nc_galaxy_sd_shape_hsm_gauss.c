/***************************************************************************
 *            nc_galaxy_sd_shape_hsm_gauss.c
 *
 *  Sun Jan 5 12:57:50 2025
 *  Copyright  2025  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_hsm_gauss.c
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
 * NcGalaxySDShapeHSMGauss:
 *
 * Class describing a galaxy sample shape distribution with a truncated gaussian p.d.f.
 * convoluted with gaussian noise. Uses the HSM shape measurement products with
 * individual shape noise.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/lss/galaxy/nc_galaxy_sd_shape_hsm_gauss.h"
#include "nc/lss/galaxy/nc_galaxy_sd_shape_hsm_gauss_global.h"
#include "nc/lss/galaxy/nc_galaxy_sd_shape.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "ncm/stats/ncm_stats_vec.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */


typedef struct _NcGalaxySDShapeHSMGaussPrivate
{
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_hp;
  NcmStatsVec *obs_stats;
} NcGalaxySDShapeHSMGaussPrivate;

struct _NcGalaxySDShapeHSMGauss
{
  NcGalaxySDShape parent_instance;
};

typedef struct _NcGalaxySDShapeHSMGaussData
{
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble epsilon_obs_t;
  gdouble epsilon_obs_x;
  gdouble std_shape;
  gdouble std_noise;
  gdouble c1;
  gdouble c2;

  /* Sky-frame calibration bias (c1, c2) rotated into the tangential/cross frame
   * by exp(-2i phi); precomputed in prepare so both the integrand and the
   * fixed-node eval work consistently in the tangential frame. */
  gdouble c1_rot;
  gdouble c2_rot;
  gdouble m;
  gdouble sigma;
  gdouble radius;
  gdouble phi;
  NcWLSurfaceMassDensityOptzs optzs;
  NcWLSurfaceMassDensityCritCache *crit_cache_arr;
  NcWLSurfaceMassDensitySigmaCache sigma_cache;
  guint crit_cache_len;
} NcGalaxySDShapeHSMGaussData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShapeHSMGauss, nc_galaxy_sd_shape_hsm_gauss, NC_TYPE_GALAXY_SD_SHAPE);

static void
nc_galaxy_sd_shape_hsm_gauss_init (NcGalaxySDShapeHSMGauss *gsdshsc)
{
  NcGalaxySDShapeHSMGaussPrivate * const self = nc_galaxy_sd_shape_hsm_gauss_get_instance_private (gsdshsc);

  self->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  self->ctrl_hp    = ncm_model_ctrl_new (NULL);
  self->obs_stats  = ncm_stats_vec_new (6, NCM_STATS_VEC_COV, FALSE);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_hsm_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeHSMGauss *gsdshsc = NC_GALAXY_SD_SHAPE_HSM_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS (gsdshsc));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_shape_hsm_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeHSMGauss *gsdshsc = NC_GALAXY_SD_SHAPE_HSM_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS (gsdshsc));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_shape_hsm_gauss_dispose (GObject *object)
{
  NcGalaxySDShapeHSMGauss *gsdshsc            = NC_GALAXY_SD_SHAPE_HSM_GAUSS (object);
  NcGalaxySDShapeHSMGaussPrivate * const self = nc_galaxy_sd_shape_hsm_gauss_get_instance_private (gsdshsc);

  ncm_model_ctrl_clear (&self->ctrl_cosmo);
  ncm_model_ctrl_clear (&self->ctrl_hp);
  ncm_stats_vec_clear (&self->obs_stats);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_hsm_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_shape_hsm_gauss_finalize (GObject *object)
{
  /* NcGalaxySDShapeHSMGauss *gsdshsc          = NC_GALAXY_SD_SHAPE_HSM_GAUSS (object); */
  /* NcGalaxySDShapeHSMGaussPrivate * const self = nc_galaxy_sd_shape_hsm_gauss_get_instance_private (gsdshsc); */

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_hsm_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_hsm_gauss_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng);
static NcGalaxySDShapeIntegrand *_nc_galaxy_sd_shape_hsm_gauss_integ (NcGalaxySDShape *gsds, gboolean use_lnp);
static gboolean _nc_galaxy_sd_shape_hsm_gauss_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gboolean update_radius, gboolean update_optzs);
static void _nc_galaxy_sd_shape_hsm_gauss_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data);
static void _nc_galaxy_sd_shape_hsm_gauss_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho);
static gboolean _nc_galaxy_sd_shape_hsm_gauss_prepare_data_array_at_nodes (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, const GPtrArray *z_nodes_per_galaxy, gboolean update_radius, gboolean update_crit, gboolean update_sigma);
static void _nc_galaxy_sd_shape_hsm_gauss_eval_at_nodes (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, const NcmVector *z_nodes, NcmVector *out);

static void
nc_galaxy_sd_shape_hsm_gauss_class_init (NcGalaxySDShapeHSMGaussClass *klass)
{
  NcGalaxySDShapeClass *sd_shape_class = NC_GALAXY_SD_SHAPE_CLASS (klass);
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class           = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_hsm_gauss_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_hsm_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_sd_shape_hsm_gauss_dispose;
  object_class->finalize    = &_nc_galaxy_sd_shape_hsm_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian galaxy shape distribution for HSC data", "Gaussian shape for HSC data");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_model_class_check_params_info (model_class);

  sd_shape_class->gen                         = &_nc_galaxy_sd_shape_hsm_gauss_gen;
  sd_shape_class->integ                       = &_nc_galaxy_sd_shape_hsm_gauss_integ;
  sd_shape_class->prepare_data_array          = &_nc_galaxy_sd_shape_hsm_gauss_prepare_data_array;
  sd_shape_class->data_init                   = &_nc_galaxy_sd_shape_hsm_gauss_data_init;
  sd_shape_class->direct_estimate             = &_nc_galaxy_sd_shape_hsm_gauss_direct_estimate;
  sd_shape_class->prepare_data_array_at_nodes = &_nc_galaxy_sd_shape_hsm_gauss_prepare_data_array_at_nodes;
  sd_shape_class->eval_at_nodes               = &_nc_galaxy_sd_shape_hsm_gauss_eval_at_nodes;
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
_nc_galaxy_sd_shape_hsm_gauss_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  NcGalaxySDShapeHSMGaussData *ldata           = (NcGalaxySDShapeHSMGaussData *) data->ldata;
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

  /* Express the position angle and the (parity-covariant) intrinsic ellipticity
   * and noise in the data's coordinate frame; both helpers are the identity for
   * the celestial frame. The additive bias is frame-fixed and not conjugated. */
  phi   = nc_wl_ellipticity_frame_position_angle (data->coord, phi);
  e_s   = nc_wl_ellipticity_frame_to_celestial_c (data->coord, e_s);
  noise = nc_wl_ellipticity_frame_to_celestial_c (data->coord, noise);

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

  /* Keep the tangential-frame caches (read by the integrand and fixed-node
   * likelihoods) consistent with the freshly generated observed ellipticity.
   * prepare() only refreshes these when the geometry changes, so a resample -
   * which rewrites epsilon_obs_1/2 in place without a model change - would
   * otherwise leave them stale. */
  {
    const complex double e_o_rotated  = e_o * cexp (-2.0 * I * phi);
    const complex double bias_rotated = c * cexp (-2.0 * I * phi);

    ldata->epsilon_obs_t = creal (e_o_rotated);
    ldata->epsilon_obs_x = cimag (e_o_rotated);
    ldata->c1_rot        = creal (bias_rotated);
    ldata->c2_rot        = cimag (bias_rotated);
  }
}

struct _IntegData
{
  NcGalaxySDShapeHSMGauss *gsdshsc;
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

static inline void
_nc_galaxy_sd_shape_hsm_gauss_pre_integ (gpointer callback_data, const gdouble z, NcGalaxySDShapeData *data,
                                         complex double (*apply_shear_inv) (complex double, complex double),
                                         gdouble (*lndet_jac) (complex double, complex double),
                                         gdouble *total_var, gdouble *chi2_1, gdouble *chi2_2, gdouble *lndetjac)
{
  struct _IntegData *int_data        = (struct _IntegData *) callback_data;
  NcGalaxySDShapeHSMGaussData *ldata = (NcGalaxySDShapeHSMGaussData *) data->ldata;
  gdouble z_cl                       = nc_halo_position_get_redshift (int_data->halo_position);
  gdouble et                         = ldata->epsilon_obs_t;
  gdouble ex                         = ldata->epsilon_obs_x;
  gdouble std_noise                  = ldata->std_noise;
  gdouble c1_rot                     = ldata->c1_rot;
  gdouble c2_rot                     = ldata->c2_rot;
  gdouble m                          = ldata->m;
  gdouble sigma                      = ldata->sigma;
  complex double e_o                 = et + I * ex;
  complex double g, e_s;

  /* Work in the tangential/cross frame (real reduced shear gt, observed
   * ellipticity epsilon_obs_t/x, calibration bias pre-rotated to c1_rot/c2_rot),
   * matching the fixed-node eval_at_nodes path. The likelihood is rotation
   * invariant, so this agrees with a sky-frame computation provided the bias is
   * rotated consistently with the observed ellipticity. */
  if (z > z_cl)
  {
    const gdouble gt = nc_wl_surface_mass_density_reduced_shear_optzs (int_data->surface_mass_density,
                                                                       int_data->density_profile,
                                                                       int_data->cosmo,
                                                                       z, z_cl, &ldata->optzs);

    g = gt;
  }
  else
  {
    const gdouble gt = nc_wl_surface_mass_density_reduced_shear_optzs (int_data->surface_mass_density,
                                                                       int_data->density_profile,
                                                                       int_data->cosmo,
                                                                       z_cl * (1.0 + GSL_DBL_EPSILON), z_cl, &ldata->optzs);
    const gdouble step = exp ((z - z_cl) / 0.001);

    g = step * gt;
  }

  /* Adding bias */
  g = (1.0 + m) * g + (c1_rot + I * c2_rot);

  e_s = apply_shear_inv (g, e_o);

  *total_var = gsl_pow_2 (sigma) + gsl_pow_2 (std_noise);
  *chi2_1    = gsl_pow_2 (creal (e_s)) / *total_var;
  *chi2_2    = gsl_pow_2 (cimag (e_s)) / *total_var;
  *lndetjac  = lndet_jac (g, e_o);
}

static inline gdouble
_nc_galaxy_sd_shape_hsm_gauss_integ_value (gdouble total_var, gdouble chi2_1, gdouble chi2_2, gdouble lndetjac)
{
  return exp (-0.5 * (chi2_1 + chi2_2) + lndetjac) / (2.0 * M_PI * total_var);
}

static inline gdouble
_nc_galaxy_sd_shape_hsm_gauss_ln_integ_value (gdouble total_var, gdouble chi2_1, gdouble chi2_2, gdouble lndetjac)
{
  return -0.5 * (chi2_1 + chi2_2) + lndetjac - log (2.0 * M_PI * total_var);
}

/* Per-convention integrands, selected once in _integ() to keep the convention
 * branch out of the per-evaluation path. */
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_INTEG_F(name, apply_shear_inv, lndet_jac, value)        \
        static gdouble                                                                       \
        name (gpointer callback_data, const gdouble z, NcGalaxySDShapeData * data)           \
        {                                                                                    \
          gdouble total_var, chi2_1, chi2_2, lndetjac;                                       \
          _nc_galaxy_sd_shape_hsm_gauss_pre_integ (callback_data, z, data,                   \
                                                   apply_shear_inv, lndet_jac,               \
                                                   &total_var, &chi2_1, &chi2_2, &lndetjac); \
          return value (total_var, chi2_1, chi2_2, lndetjac);                                \
        }

NC_GALAXY_SD_SHAPE_HSM_GAUSS_INTEG_F (_nc_galaxy_sd_shape_hsm_gauss_integ_f_trace,
                                      nc_wl_ellipticity_apply_shear_inv_trace_c,
                                      nc_wl_ellipticity_lndet_jac_trace_c,
                                      _nc_galaxy_sd_shape_hsm_gauss_integ_value)
NC_GALAXY_SD_SHAPE_HSM_GAUSS_INTEG_F (_nc_galaxy_sd_shape_hsm_gauss_integ_f_trace_det,
                                      nc_wl_ellipticity_apply_shear_inv_trace_det_c,
                                      nc_wl_ellipticity_lndet_jac_trace_det_c,
                                      _nc_galaxy_sd_shape_hsm_gauss_integ_value)
NC_GALAXY_SD_SHAPE_HSM_GAUSS_INTEG_F (_nc_galaxy_sd_shape_hsm_gauss_ln_integ_f_trace,
                                      nc_wl_ellipticity_apply_shear_inv_trace_c,
                                      nc_wl_ellipticity_lndet_jac_trace_c,
                                      _nc_galaxy_sd_shape_hsm_gauss_ln_integ_value)
NC_GALAXY_SD_SHAPE_HSM_GAUSS_INTEG_F (_nc_galaxy_sd_shape_hsm_gauss_ln_integ_f_trace_det,
                                      nc_wl_ellipticity_apply_shear_inv_trace_det_c,
                                      nc_wl_ellipticity_lndet_jac_trace_det_c,
                                      _nc_galaxy_sd_shape_hsm_gauss_ln_integ_value)

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
_nc_galaxy_sd_shape_hsm_gauss_integ (NcGalaxySDShape *gsds, gboolean use_lnp)
{
  NcGalaxySDShapeHSMGauss *gsdshsc        = NC_GALAXY_SD_SHAPE_HSM_GAUSS (gsds);
  struct _IntegData *int_data             = g_new0 (struct _IntegData, 1);
  const NcGalaxyWLObsEllipConv ellip_conv = nc_galaxy_sd_shape_get_ellip_conv (gsds);
  NcGalaxySDShapeIntegrandFunc integ_f;
  NcGalaxySDShapeIntegrand *integ;

  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      integ_f = use_lnp ? _nc_galaxy_sd_shape_hsm_gauss_ln_integ_f_trace : _nc_galaxy_sd_shape_hsm_gauss_integ_f_trace;
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      integ_f = use_lnp ? _nc_galaxy_sd_shape_hsm_gauss_ln_integ_f_trace_det : _nc_galaxy_sd_shape_hsm_gauss_integ_f_trace_det;
      break;
    default:                    /* LCOV_EXCL_LINE */
      g_assert_not_reached ();  /* LCOV_EXCL_LINE */
  }

  integ = nc_galaxy_sd_shape_integrand_new (integ_f,
                                            _integ_data_free,
                                            _integ_data_copy,
                                            _integ_data_prepare,
                                            int_data);

  int_data->gsdshsc = gsdshsc;

  return integ;
}

static gboolean
_nc_galaxy_sd_shape_hsm_gauss_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gboolean update_radius, gboolean update_optzs)
{
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  NcWLSurfaceMassDensityLensCtx lens_ctx;
  gboolean pr_prefactor_ready = FALSE;
  gdouble r_s, rho_s, pr_prefactor = 0.0;
  guint i;

  NCM_UNUSED (gsds);

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);

  if (update_radius)
  {
    pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
    pr_prefactor_ready = TRUE;
  }

  if (update_optzs)
  {
    nc_wl_surface_mass_density_lens_ctx_prep (&lens_ctx, surface_mass_density, cosmo, z_cl);
    nc_halo_density_profile_r_s_rho_s (density_profile, cosmo, z_cl, &r_s, &rho_s);
  }

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxySDShapeData *data_i          = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeHSMGaussData *ldata_i = (NcGalaxySDShapeHSMGaussData *) data_i->ldata;

    if ((update_radius) || (ldata_i->radius == 0.0))
    {
      gdouble e1         = ldata_i->epsilon_obs_1;
      gdouble e2         = ldata_i->epsilon_obs_2;
      const gdouble ra   = data_i->sdpos_data->ra;
      const gdouble dec  = data_i->sdpos_data->dec;
      complex double e_o = e1 + I * e2;
      complex double e_o_rotated;
      complex double bias_rotated;
      gdouble theta, phi;

      nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

      phi = nc_wl_ellipticity_frame_position_angle (data_i->coord, phi);

      e_o_rotated  = e_o * cexp (-2.0 * I * phi);
      bias_rotated = (ldata_i->c1 + I * ldata_i->c2) * cexp (-2.0 * I * phi);

      if (!pr_prefactor_ready)
      {
        pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
        pr_prefactor_ready = TRUE;
      }

      ldata_i->radius        = nc_halo_position_projected_radius_from_prefactor (theta, pr_prefactor);
      ldata_i->phi           = phi;
      ldata_i->epsilon_obs_t = creal (e_o_rotated);
      ldata_i->epsilon_obs_x = cimag (e_o_rotated);
      ldata_i->c1_rot        = creal (bias_rotated);
      ldata_i->c2_rot        = cimag (bias_rotated);
    }

    if (update_optzs)
      nc_wl_surface_mass_density_reduced_shear_optzs_prep_with_lens_ctx (
        density_profile,
        &lens_ctx,
        ldata_i->radius,
        r_s,
        rho_s,
        &ldata_i->optzs
      );
  }

  return TRUE;
}

static gboolean
_nc_galaxy_sd_shape_hsm_gauss_prepare_data_array_at_nodes (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, const GPtrArray *z_nodes_per_galaxy, gboolean update_radius, gboolean update_crit, gboolean update_sigma)
{
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  const gdouble z_cl                           = nc_halo_position_get_redshift (halo_position);
  gboolean lens_ctx_ready                      = FALSE;
  gboolean pr_prefactor_ready                  = FALSE;
  NcWLSurfaceMassDensityLensCtx lens_ctx;
  gdouble pr_prefactor = 0.0;
  guint i;

  NCM_UNUSED (gsds);

  /* This unconditionally (re)computes every per-galaxy node cache. Deciding
   * *when* a recompute is needed is the caller's responsibility (e.g. the
   * NcDataClusterWL that owns the data tracks cosmo/halo/profile changes), so
   * these caches are never invalidated through state shared between consumers. */
  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);

  if (update_radius)
  {
    pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
    pr_prefactor_ready = TRUE;
  }

  if (update_crit)
  {
    nc_wl_surface_mass_density_lens_ctx_prep (&lens_ctx, surface_mass_density, cosmo, z_cl);
    lens_ctx_ready = TRUE;
  }

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxySDShapeData *data_i          = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeHSMGaussData *ldata_i = (NcGalaxySDShapeHSMGaussData *) data_i->ldata;
    NcmVector *z_nodes_i                 = (NcmVector *) g_ptr_array_index (z_nodes_per_galaxy, i);
    const guint n_nodes                  = ncm_vector_len (z_nodes_i);

    if ((update_radius) || (ldata_i->radius == 0.0))
    {
      gdouble e1         = ldata_i->epsilon_obs_1;
      gdouble e2         = ldata_i->epsilon_obs_2;
      const gdouble ra   = data_i->sdpos_data->ra;
      const gdouble dec  = data_i->sdpos_data->dec;
      complex double e_o = e1 + I * e2;
      complex double e_o_rotated;
      complex double bias_rotated;
      gdouble theta, phi;

      nc_halo_position_polar_angles (halo_position, ra, dec, &theta, &phi);

      phi = nc_wl_ellipticity_frame_position_angle (data_i->coord, phi);

      e_o_rotated  = e_o * cexp (-2.0 * I * phi);
      bias_rotated = (ldata_i->c1 + I * ldata_i->c2) * cexp (-2.0 * I * phi);

      if (!pr_prefactor_ready)
      {
        pr_prefactor       = nc_halo_position_projected_radius_prefactor (halo_position, cosmo);
        pr_prefactor_ready = TRUE;
      }

      ldata_i->radius        = nc_halo_position_projected_radius_from_prefactor (theta, pr_prefactor);
      ldata_i->phi           = phi;
      ldata_i->epsilon_obs_t = creal (e_o_rotated);
      ldata_i->epsilon_obs_x = cimag (e_o_rotated);
      ldata_i->c1_rot        = creal (bias_rotated);
      ldata_i->c2_rot        = cimag (bias_rotated);
    }

    if (update_crit || (ldata_i->crit_cache_arr == NULL))
    {
      guint j;

      if (!lens_ctx_ready)
      {
        nc_wl_surface_mass_density_lens_ctx_prep (
          &lens_ctx,
          surface_mass_density,
          cosmo,
          z_cl
        );
        lens_ctx_ready = TRUE;
      }

      if (ldata_i->crit_cache_len != n_nodes)
      {
        g_free (ldata_i->crit_cache_arr);
        ldata_i->crit_cache_arr = g_new0 (NcWLSurfaceMassDensityCritCache, n_nodes);
        ldata_i->crit_cache_len = n_nodes;
      }

      for (j = 0; j < n_nodes; j++)
      {
        const gdouble z_j = ncm_vector_get (z_nodes_i, j);

        nc_wl_surface_mass_density_reduced_shear_crit_cache_prep_with_lens_ctx (
          surface_mass_density,
          cosmo,
          &lens_ctx,
          z_j,
          &ldata_i->crit_cache_arr[j]
        );
      }
    }

    if (update_sigma)
      nc_wl_surface_mass_density_reduced_shear_sigma_cache_prep (
        density_profile,
        cosmo,
        ldata_i->radius,
        z_cl,
        z_cl,
        &ldata_i->sigma_cache
      );
  }

  return TRUE;
}

static inline void
_nc_galaxy_sd_shape_hsm_gauss_eval_at_nodes_conv (NcGalaxySDShapeHSMGaussData *ldata, const NcmVector *z_nodes, NcmVector *out,
                                                  gdouble z_cl, complex double e_o, gdouble total_var,
                                                  complex double (*apply_shear_inv) (complex double, complex double),
                                                  gdouble (*lndet_jac) (complex double, complex double))
{
  const gdouble c1_rot = ldata->c1_rot;
  const gdouble c2_rot = ldata->c2_rot;
  const gdouble m      = ldata->m;
  const guint n_nodes  = ncm_vector_len (z_nodes);
  guint j;

  for (j = 0; j < n_nodes; j++)
  {
    const gdouble z_j = ncm_vector_get (z_nodes, j);
    gdouble gt;
    complex double g, e_s;
    gdouble lndetjac, chi2_1, chi2_2;

    if (z_j > z_cl)
      gt = nc_wl_surface_mass_density_reduced_shear_cache (&ldata->crit_cache_arr[j], &ldata->sigma_cache);
    else
      gt = 0.0;

    g   = (1.0 + m) * gt + (c1_rot + I * c2_rot);
    e_s = apply_shear_inv (g, e_o);

    lndetjac = lndet_jac (g, e_o);
    chi2_1   = gsl_pow_2 (creal (e_s)) / total_var;
    chi2_2   = gsl_pow_2 (cimag (e_s)) / total_var;

    ncm_vector_set (out, j, exp (-0.5 * (chi2_1 + chi2_2) + lndetjac) / (2.0 * M_PI * total_var));
  }
}

static void
_nc_galaxy_sd_shape_hsm_gauss_eval_at_nodes (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, const NcmVector *z_nodes, NcmVector *out)
{
  NcHaloPosition *halo_position           = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxySDShapeHSMGaussData *ldata      = (NcGalaxySDShapeHSMGaussData *) data->ldata;
  const gdouble z_cl                      = nc_halo_position_get_redshift (halo_position);
  const gdouble et                        = ldata->epsilon_obs_t;
  const gdouble ex                        = ldata->epsilon_obs_x;
  const gdouble std_noise                 = ldata->std_noise;
  const gdouble sigma                     = ldata->sigma;
  const complex double e_o                = et + I * ex;
  const gdouble total_var                 = gsl_pow_2 (sigma) + gsl_pow_2 (std_noise);
  const NcGalaxyWLObsEllipConv ellip_conv = nc_galaxy_sd_shape_get_ellip_conv (gsds);

  g_assert_nonnull (ldata->crit_cache_arr);

  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      _nc_galaxy_sd_shape_hsm_gauss_eval_at_nodes_conv (ldata, z_nodes, out, z_cl, e_o, total_var,
                                                        nc_wl_ellipticity_apply_shear_inv_trace_c,
                                                        nc_wl_ellipticity_lndet_jac_trace_c);
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      _nc_galaxy_sd_shape_hsm_gauss_eval_at_nodes_conv (ldata, z_nodes, out, z_cl, e_o, total_var,
                                                        nc_wl_ellipticity_apply_shear_inv_trace_det_c,
                                                        nc_wl_ellipticity_lndet_jac_trace_det_c);
      break;
    default:                    /* LCOV_EXCL_LINE */
      g_assert_not_reached ();  /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_shape_hsm_gauss_ldata_free (gpointer ldata)
{
  NcGalaxySDShapeHSMGaussData *hdata = (NcGalaxySDShapeHSMGaussData *) ldata;

  g_clear_pointer (&hdata->crit_cache_arr, g_free);
  g_free (hdata);
}

static void
_nc_galaxy_sd_shape_hsm_gauss_ldata_read_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDShapeHSMGaussData *ldata = (NcGalaxySDShapeHSMGaussData *) data->ldata;

  nc_galaxy_sd_position_data_read_row (data->sdpos_data, obs, i);

  ldata->epsilon_obs_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_1, i, NULL);
  ldata->epsilon_obs_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_2, i, NULL);
  ldata->std_shape     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_SHAPE, i, NULL);
  ldata->sigma         = nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (ldata->std_shape);
  ldata->std_noise     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_NOISE, i, NULL);
  ldata->c1            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C1, i, NULL);
  ldata->c2            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C2, i, NULL);
  ldata->m             = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_M, i, NULL);
}

static void
_nc_galaxy_sd_shape_hsm_gauss_ldata_write_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDShapeHSMGaussData *ldata = (NcGalaxySDShapeHSMGaussData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_1, i, ldata->epsilon_obs_1, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_2, i, ldata->epsilon_obs_2, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_SHAPE, i, ldata->std_shape, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_NOISE, i, ldata->std_noise, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C1, i, ldata->c1, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C2, i, ldata->c2, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_M, i, ldata->m, NULL);
}

static void
_nc_galaxy_sd_shape_hsm_gauss_ldata_required_columns (NcGalaxySDShapeData *data, GList *columns)
{
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_SHAPE));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_NOISE));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_M));
}

static gdouble
_nc_galaxy_sd_shape_hsm_gauss_ldata_get_radius (NcGalaxySDShapeData *data)
{
  NcGalaxySDShapeHSMGaussData *ldata = (NcGalaxySDShapeHSMGaussData *) data->ldata;

  return ldata->radius;
}

static void
_nc_galaxy_sd_shape_hsm_gauss_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data)
{
  NcGalaxySDShapeHSMGaussData *ldata = g_new0 (NcGalaxySDShapeHSMGaussData, 1);

  data->sdpos_data             = sdpos_data;
  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_shape_hsm_gauss_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_shape_hsm_gauss_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_shape_hsm_gauss_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_shape_hsm_gauss_ldata_required_columns;
  data->ldata_get_radius       = &_nc_galaxy_sd_shape_hsm_gauss_ldata_get_radius;
}

static void
_nc_galaxy_sd_shape_hsm_gauss_direct_estimate (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho)
{
  NcGalaxySDShapeHSMGauss *gsdshsc            = NC_GALAXY_SD_SHAPE_HSM_GAUSS (gsds);
  NcGalaxySDShapeHSMGaussPrivate * const self = nc_galaxy_sd_shape_hsm_gauss_get_instance_private (gsdshsc);
  NcHICosmo *cosmo                            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *halo_position               = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxyWLObsEllipConv ellip_conv           = nc_galaxy_sd_shape_get_ellip_conv (gsds);
  guint i;

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  ncm_stats_vec_reset (self->obs_stats, TRUE);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxySDShapeData *data_i          = g_ptr_array_index (data_array, i);
    NcGalaxySDShapeHSMGaussData *ldata_i = (NcGalaxySDShapeHSMGaussData *) data_i->ldata;
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

    phi = nc_wl_ellipticity_frame_position_angle (data_i->coord, phi);

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
        /* epsilon convention: <epsilon> = g, no responsivity factor R */
        *gt      = (mean_gt - mean_c1) / (1.0 + mean_m);
        *gx      = (mean_gx - mean_c2) / (1.0 + mean_m);
        *sigma_t = ncm_stats_vec_get_sd (self->obs_stats, 0) / (1.0 + mean_m);
        *sigma_x = ncm_stats_vec_get_sd (self->obs_stats, 1) / (1.0 + mean_m);
        *rho     = ncm_stats_vec_get_cor (self->obs_stats, 0, 1);
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        /* distortion convention: <chi> = 2 R g, responsivity applies to both components */
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
 * nc_galaxy_sd_shape_hsm_gauss_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxySDShapeHSMGauss
 *
 * Returns: (transfer full): a new NcGalaxySDShapeHSMGauss.
 */
NcGalaxySDShapeHSMGauss *
nc_galaxy_sd_shape_hsm_gauss_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  NcGalaxySDShapeHSMGauss *gsdshsc = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_HSM_GAUSS,
                                                   "ellip-conv", ellip_conv,
                                                   NULL);

  return gsdshsc;
}

/**
 * nc_galaxy_sd_shape_hsm_gauss_ref:
 * @gsdshsc: a #NcGalaxySDShapeHSMGauss
 *
 * Increase the reference of @gsdshsc by one.
 *
 * Returns: (transfer full): @gsdshsc.
 */
NcGalaxySDShapeHSMGauss *
nc_galaxy_sd_shape_hsm_gauss_ref (NcGalaxySDShapeHSMGauss *gsdshsc)
{
  return g_object_ref (gsdshsc);
}

/**
 * nc_galaxy_sd_shape_hsm_gauss_free:
 * @gsdshsc: a #NcGalaxySDShapeHSMGauss
 *
 * Decrease the reference count of @gsdshsc by one.
 *
 */
void
nc_galaxy_sd_shape_hsm_gauss_free (NcGalaxySDShapeHSMGauss *gsdshsc)
{
  g_object_unref (gsdshsc);
}

/**
 * nc_galaxy_sd_shape_hsm_gauss_clear:
 * @gsdshsc: a #NcGalaxySDShapeHSMGauss
 *
 * Decrease the reference count of @gsdshsc by one, and sets the pointer *@gsdshsc to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_hsm_gauss_clear (NcGalaxySDShapeHSMGauss **gsdshsc)
{
  g_clear_object (gsdshsc);
}

/**
 * nc_galaxy_sd_shape_hsm_gauss_gen:
 * @gsdshsc: a #NcGalaxySDShapeHSMGauss
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDShapeData
 * @std_shape: the intrinsic shape dispersion
 * @std_noise: the observational shape dispersion
 * @c1: the first additive bias parameter
 * @c2: the second additive bias parameter
 * @m: the multiplicative bias parameter
 * @coord: the coordinate system #NcWLEllipticityFrame
 * @rng: a #NcmRNG
 *
 * Generates a galaxy sample shape.
 *
 */
void
nc_galaxy_sd_shape_hsm_gauss_gen (NcGalaxySDShapeHSMGauss *gsdshsc, NcmMSet *mset, NcGalaxySDShapeData *data, const gdouble std_shape, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcWLEllipticityFrame coord, NcmRNG *rng)
{
  NcGalaxySDShapeClass *sd_shape_class = NC_GALAXY_SD_SHAPE_GET_CLASS (gsdshsc);
  NcGalaxySDShapeHSMGaussData *ldata   = (NcGalaxySDShapeHSMGaussData *) data->ldata;

  data->coord      = coord;
  ldata->std_shape = std_shape;
  ldata->std_noise = std_noise;
  ldata->sigma     = nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (std_shape);
  ldata->c1        = c1;
  ldata->c2        = c2;
  ldata->m         = m;

  sd_shape_class->gen (NC_GALAXY_SD_SHAPE (gsdshsc), mset, data, rng);
}

/**
 * nc_galaxy_sd_shape_hsm_gauss_data_set:
 * @gsdshsc: a #NcGalaxySDShapeHSMGauss
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
nc_galaxy_sd_shape_hsm_gauss_data_set (NcGalaxySDShapeHSMGauss *gsdshsc, NcGalaxySDShapeData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_shape, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m)
{
  NcGalaxySDShapeHSMGaussData *ldata = (NcGalaxySDShapeHSMGaussData *) data->ldata;

  ldata->epsilon_obs_1 = epsilon_obs_1;
  ldata->epsilon_obs_2 = epsilon_obs_2;
  ldata->std_shape     = std_shape;
  ldata->sigma         = nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (std_shape);
  ldata->std_noise     = std_noise;
  ldata->c1            = c1;
  ldata->c2            = c2;
  ldata->m             = m;
}

/**
 * nc_galaxy_sd_shape_hsm_gauss_data_get:
 * @gsdshsc: a #NcGalaxySDShapeHSMGauss
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
nc_galaxy_sd_shape_hsm_gauss_data_get (NcGalaxySDShapeHSMGauss *gsdshsc, NcGalaxySDShapeData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_shape, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m)
{
  NcGalaxySDShapeHSMGaussData *ldata = (NcGalaxySDShapeHSMGaussData *) data->ldata;

  *epsilon_obs_1 = ldata->epsilon_obs_1;
  *epsilon_obs_2 = ldata->epsilon_obs_2;
  *std_shape     = ldata->std_shape;
  *std_noise     = ldata->std_noise;
  *c1            = ldata->c1;
  *c2            = ldata->c2;
  *m             = ldata->m;
}

