/***************************************************************************
 *            nc_distance.c
 *
 *  Tue May  8 11:05:35 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:nc_distance
 * @title: Cosmological Distances and Times
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_distance.h"
#include "nc_macros.h"
#include "math/integral.h"
#include "math/ncm_c.h"
#include "math/ncm_spline_cubic_notaknot.h"

typedef struct _ComovingDistanceArgument{
  NcHICosmo *model;
  gint diff;
} ComovingDistanceArgument;

enum
{
  PROP_0,
  PROP_ZF
};

G_DEFINE_TYPE (NcDistance, nc_distance, G_TYPE_OBJECT);

/***************************************************************************
 * Functions that transforms distances
 *
 ****************************************************************************/

/**
 * nc_distance_comoving2luminosity:
 * @k: FIXME
 * @sqrt_Omega_k: FIXME
 * @hubble_dist: FIXME
 * @z: FIXME
 * @cd: FIXME
 *
 * FIXME
 *
 * Returns: an integer.
 */
gdouble
nc_distance_comoving2luminosity (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble cd)
{
  switch (k)
  {
	case 0:
	  return (1.0 + z) * hubble_dist * cd;
	case -1:
	  return (1.0 + z) * hubble_dist * sinh (sqrt_Omega_k * cd) / sqrt_Omega_k;
	case 1:
	  return (1.0 + z) * hubble_dist * fabs(sin  (sqrt_Omega_k * cd)) / sqrt_Omega_k;
	default:
	  g_assert_not_reached();
	  return 0.0;
  }
}

/**
 * nc_distance_luminosity2modulo:
 * @k: FIXME
 * @sqrt_Omega_k: FIXME
 * @hubble_dist: FIXME
 * @z: FIXME
 * @ld: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_luminosity2modulo (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble ld)
{
  return 5 * log10 (ld) + 25;
}

/**
 * nc_distance_comoving2modulo:
 * @k: FIXME
 * @sqrt_Omega_k: FIXME
 * @hubble_dist: FIXME
 * @z: FIXME
 * @cd: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving2modulo (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble cd)
{
  gdouble ld = nc_distance_comoving2luminosity (k, sqrt_Omega_k, hubble_dist, z, cd);
  return nc_distance_luminosity2modulo (k, sqrt_Omega_k, hubble_dist, z, ld);
}

/**
 * nc_distance_luminosity2comoving:
 * @k: FIXME
 * @sqrt_Omega_k: FIXME
 * @hubble_dist: FIXME
 * @z: FIXME
 * @ld: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_luminosity2comoving (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble ld)
{
  switch (k)
  {
	case 0:
	  return ld / ((1.0 + z) * hubble_dist);
	case -1:
	  return asinh(sqrt_Omega_k*ld / (hubble_dist * (1 + z))) / sqrt_Omega_k;
	case 1:
	{
	  gdouble arg = sqrt_Omega_k*ld / (hubble_dist * (1 + z));
	  gdouble n = ceil(arg);
	  arg -= n;
	  return (asin(arg) + n * 2.0 * M_PI)/ sqrt_Omega_k;
	}
	default:
	  g_assert_not_reached();
	  return 0.0;
  }
}

/**
 * nc_distance_modulo2luminosity:
 * @k: FIXME
 * @sqrt_Omega_k: FIXME
 * @hubble_dist: FIXME
 * @z: FIXME
 * @dm: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_modulo2luminosity (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble dm)
{
  return pow (10.0, dm/5.0 - 5.0);
}

/**
 * nc_distance_modulo2comoving:
 * @k: FIXME
 * @sqrt_Omega_k: FIXME
 * @hubble_dist: FIXME
 * @z: FIXME
 * @dm: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_modulo2comoving (gint k, gdouble sqrt_Omega_k, gdouble hubble_dist, gdouble z, gdouble dm)
{
  gdouble dl = nc_distance_modulo2luminosity(k, sqrt_Omega_k, hubble_dist, z, dm);
  return nc_distance_luminosity2comoving (k, sqrt_Omega_k, hubble_dist, z, dl);
}

/***************************************************************************
 * Calculate the curvature scale today
 *
 ****************************************************************************/

/**
 * nc_distance_hubble: (skip)
 * @dist: FIXME
 * @model: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_hubble (NcDistance *dist, NcHICosmo *model)
{
  return ncm_c_c () / (nc_hicosmo_H0 (model) * 1.0e3);
}

/***************************************************************************
 * Calculate the comoving distance integrating numerically the Hubble
 * function
 ****************************************************************************/

static gdouble comoving_distance_integral_argument(gdouble z, gpointer p);
static gdouble dcddz (gdouble y, gdouble x, gpointer userdata);

/**
 * nc_distance_comoving:
 * @dist: pointer to the type defined by #NcDistance
 * @model: pointer to the type defined by #NcHICosmo
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble result, error;
  gsl_function F;

  nc_distance_prepare_if_needed (dist, model);

  if (ncm_model_impl (NCM_MODEL (model)) & NC_HICOSMO_IMPL_cd)
    return nc_hicosmo_cd (model, z);

  if (z <= dist->z_f)
    return ncm_spline_eval (dist->comoving_distance_spline->s, z);

  F.function = &comoving_distance_integral_argument;
  F.params = model;

  if (dist->use_cache)
    nc_integral_cached_0_x (dist->comoving_distance_cache, &F, z, &result, &error);
  else
    nc_integral_locked_a_b (&F, 0.0, z, 0.0, NC_INT_ERROR, &result, &error);

  return result;
}

static gdouble
comoving_distance_integral_argument(gdouble z, gpointer p)
{
  NcHICosmo *model = NC_HICOSMO (p);
  gdouble E2 = nc_hicosmo_E2 (model, z);
  if (GSL_SIGN(E2) == -1.0)
	return GSL_POSINF;
  return 1.0 / sqrt (E2);
}

static gdouble
dcddz (gdouble cd, gdouble z, gpointer userdata)
{
  NcHICosmo *model = NC_HICOSMO (userdata);
  gdouble E2 = nc_hicosmo_E2 (model, z);
  return 1.0 / sqrt (E2);
}

/***************************************************************************
 * Calculate the curvature distance using nc_distance_comoving
 *
 ****************************************************************************/

/**
 * nc_distance_curvature:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_curvature (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble sqrt_Omega_k = sqrt(fabs(nc_hicosmo_Omega_k (model)));
  return sqrt_Omega_k * nc_distance_comoving (dist, model, z);
}

/***************************************************************************
 * Calculate the transverse distance using nc_distance_comoving
 *
 ****************************************************************************/

/**
 * nc_distance_transverse:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_transverse (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble Omega_k = nc_hicosmo_Omega_k (model);
  gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  gdouble comoving_dist = nc_distance_comoving (dist, model, z);
  gint k = fabs (Omega_k) < NC_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);

  if (gsl_isinf(comoving_dist)) return comoving_dist;
  switch (k)
  {
	case 0:
	  return comoving_dist;
	case -1:
	  return sinh (sqrt_Omega_k * comoving_dist) / sqrt_Omega_k;
	case 1:
	  return fabs (sin (sqrt_Omega_k * comoving_dist)) / sqrt_Omega_k; // LOOK
	default:
	  g_assert_not_reached();
	  return 0.0;
  }
}

/***************************************************************************
 * Calculate the transverse distance using nc_distance_comoving
 *
 ****************************************************************************/

/**
 * nc_distance_dtransverse_dz:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_dtransverse_dz (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble Omega_k = nc_hicosmo_Omega_k (model);
  gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  gdouble E = sqrt (nc_hicosmo_E2(model, z));
  gint k = fabs (Omega_k) < NC_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);

  switch (k)
  {
	case 0:
	  return 1.0 / E;
	case -1:
	{
	  gdouble comoving_dist = nc_distance_comoving (dist, model, z);
	  return cosh (sqrt_Omega_k * comoving_dist) / E;
	}
	case 1:
	{
	  gdouble comoving_dist = nc_distance_comoving (dist, model, z);
	  return ncm_c_sign_sin (sqrt_Omega_k * comoving_dist) * cos (sqrt_Omega_k * comoving_dist) / E; // LOOK
	}
	default:
	  g_assert_not_reached();
	  return 0.0;
  }
}

/***************************************************************************
 * Calculate the luminosity distance using nc_distance_comoving
 *
 ****************************************************************************/

/**
 * nc_distance_luminosity:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * The function description goes here.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_luminosity (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble Omega_k = nc_hicosmo_Omega_k (model);
  gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  gdouble comoving_dist = nc_distance_comoving (dist, model, z);
  gdouble hubble_dist = nc_distance_hubble (dist, model);
  gint k = fabs (Omega_k) < NC_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);

  if (gsl_isinf(comoving_dist))
	return comoving_dist;

  return nc_distance_comoving2luminosity (k, sqrt_Omega_k, hubble_dist, z, comoving_dist);
}

/***************************************************************************
 * Calculate the distance modulo using nc_distance_luminosity
 *
 ****************************************************************************/

/**
 * nc_distance_modulo:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_modulo (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble luminosity_distance = nc_distance_luminosity (dist, model, z);
  if (!gsl_finite (luminosity_distance)) // CHECK IT LATER, THINK MORE ABOUT IT
  {
	return luminosity_distance;
  }
  return (5.0 * log10 (luminosity_distance) + 25.0);
}

/***************************************************************************
 * Calculate the angular diameter curvature scale
 *
 ****************************************************************************/

/**
 * nc_distance_angular_diameter_curvature_scale:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @userdata: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_angular_diameter_curvature_scale (NcDistance *dist, NcHICosmo *model, gpointer userdata)
{
  gdouble z_star = nc_distance_decoupling_redshift (dist, model);
  if (gsl_finite (z_star))
	return sqrt(nc_hicosmo_E2 (model, z_star)) *
	nc_distance_transverse (dist, model, z_star) / (1.0 + z_star);
  else
	return GSL_NAN;
}

/***************************************************************************
 * Calculate the shift parameter (sqrt(Omega_m)*D_A(z_\star))
 *
 ****************************************************************************/

/**
 * nc_distance_shift_parameter:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_shift_parameter (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble sqrt_mod_Omega_m = sqrt(fabs(nc_hicosmo_Omega_m (model)));
  gdouble transverse = nc_distance_transverse (dist, model, z);
  return sqrt_mod_Omega_m * transverse;
}

/**
 * nc_distance_shift_parameter_lss:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 *
 * Calculate the shift parameter (sqrt(Omega_m)*D_A(z_\star))
 *
 * Returns: FIXME
 */
gdouble
nc_distance_shift_parameter_lss (NcDistance *dist, NcHICosmo *model)
{
  gdouble sqrt_mod_Omega_m = sqrt(fabs(nc_hicosmo_Omega_m (model)));
  gdouble z_star = nc_distance_decoupling_redshift (dist, model);
  gdouble transverse;
  if (gsl_finite (z_star))
  {
	transverse = nc_distance_transverse (dist, model, z_star);
	return sqrt_mod_Omega_m * transverse;
  }
  else
	return GSL_NAN;
}

/**
 * nc_distance_comoving_a0:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving_a0 (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble sqrt_mod_Omega_k = sqrt(fabs(nc_hicosmo_Omega_k (model)));
  gdouble comoving = nc_distance_comoving (dist, model, z);
  return sqrt_mod_Omega_k * comoving;
}

/***************************************************************************
 * Calculate the d_c (z_lss) / a0
 *
 ****************************************************************************/

/**
 * nc_distance_comoving_a0_lss:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving_a0_lss (NcDistance *dist, NcHICosmo *model)
{
  gdouble sqrt_mod_Omega_k = sqrt(fabs(nc_hicosmo_Omega_k (model)));
  gdouble z_star = nc_distance_decoupling_redshift (dist, model);
  gdouble comoving;
  if (gsl_finite (z_star))
  {
	comoving = nc_distance_comoving (dist, model, z_star);
	return sqrt_mod_Omega_k * comoving;
  }
  else
	return GSL_NAN;
}

/***************************************************************************
 * Calculate the d_c (z_lss)
 *
 ****************************************************************************/

/**
 * nc_distance_comoving_lss:
 * @dist: a #NcDistance
 * @model: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving_lss (NcDistance *dist, NcHICosmo *model)
{
  gdouble z_star = nc_distance_decoupling_redshift (dist, model);
  if (gsl_finite (z_star))
	return nc_distance_comoving (dist, model, z_star);
  else
	return GSL_NAN;
}

/***************************************************************************
 * Decoupling redshift (arXiv:astro-ph/9510117)
 *
 ****************************************************************************/

gdouble
nc_distance_decoupling_redshift (NcDistance *dist, NcHICosmo *model)
{
  if (ncm_model_impl (NCM_MODEL (model)) & NC_HICOSMO_IMPL_z_lss)
	return nc_hicosmo_z_lss (model);
  else
  {
	gdouble h = nc_hicosmo_h (model);
	gdouble h2 = h * h;
	gdouble omega_b_h2 = nc_hicosmo_Omega_b (model) * h2;
	gdouble omega_m_h2 = nc_hicosmo_Omega_m (model) * h2;
	gdouble g1 = 0.0783 * pow (omega_b_h2, -0.238) / (1.0 + 39.5 * pow (omega_b_h2, 0.763));
	gdouble g2 = 0.560 / (1.0 + 21.1 * pow (omega_b_h2, 1.81));
	return 1048.0 * (1.0 + 1.24e-3 * pow (omega_b_h2, -0.738)) * (1.0 + g1 * pow (omega_m_h2, g2));
  }
}

/***************************************************************************
 * Calculate the sound horizon integrating numerically the hubble
 * function
 ****************************************************************************/

static gdouble sound_horizon_integral_argument(gdouble z, gpointer p);

gdouble
nc_distance_sound_horizon (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble result, error;
  gsl_function F;

  if (ncm_model_impl (NCM_MODEL (model)) & NC_HICOSMO_IMPL_cd)
	return nc_hicosmo_cd (model, z);

  F.function = &sound_horizon_integral_argument;
  F.params = model;

  if (dist->use_cache)
	nc_integral_cached_x_inf (dist->sound_horizon_cache, &F, z, &result, &error);
  else
	nc_integral_locked_a_inf (&F, z, NC_INT_ABS_ERROR, NC_INT_ERROR, &result, &error);

  return result;
}

static gdouble
sound_horizon_integral_argument (gdouble z, gpointer p)
{
  NcHICosmo *model = NC_HICOSMO (p);
  gdouble E2 = nc_hicosmo_E2 (model, z);
  gdouble omega_b = nc_hicosmo_Omega_b (model);
  gdouble omega_r = nc_hicosmo_Omega_r (model);
  if (omega_r == 0.0)
	return 0.0;
  if (GSL_SIGN(E2) == -1.0)
	return GSL_POSINF;
  return
	1.0 / sqrt (E2 * (3.0 + 9.0 / 4.0 * (1.0 + ncm_c_neutrino_n_eff () * 0.2271) * omega_b / (omega_r * (1.0 + z))));
}

/***************************************************************************
 * Calculate the sound horizon derivative
 *
 ****************************************************************************/

gdouble
nc_distance_dsound_horizon_dz (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  return -sound_horizon_integral_argument (z, model);
}

/***************************************************************************
 * Calculate the acoustic scale l_A = pi D_T (z_lss) / r (z_lss)
 *
 ****************************************************************************/

gdouble
nc_distance_acoustic_scale (NcDistance *dist, NcHICosmo *model)
{
  gdouble z = nc_distance_decoupling_redshift (dist, model);
  if (gsl_finite (z))
	return M_PI * nc_distance_transverse (dist, model, z) / nc_distance_sound_horizon (dist, model, z);
  else
	return GSL_NAN;
}

/***************************************************************************
 * Drag redshift (arXiv:astro-ph/9510117)
 *
 ****************************************************************************/

gdouble
nc_distance_drag_redshift (NcDistance *dist, NcHICosmo *model)
{
  gdouble h = nc_hicosmo_h (model);
  gdouble h2 = h*h;
  gdouble omega_m_h2 = nc_hicosmo_Omega_m (model) * h2;
  gdouble omega_b_h2 = nc_hicosmo_Omega_b (model) * h2;
  gdouble b1 = 0.313 * pow (omega_m_h2, -0.419) * (1.0 + 0.607 * pow (omega_m_h2, 0.674));
  gdouble b2 = 0.238 * pow (omega_m_h2, 0.223);
  return 1291.0 * pow (omega_m_h2, 0.251) / (1.0 + 0.659 * pow (omega_m_h2, 0.828)) *
	(1.0 + b1 * pow (omega_b_h2, b2));
}

/***************************************************************************
 * Dilation scale D_v(z) -- (arXiv:astro-ph/0501171)
 *
 ****************************************************************************/

gdouble
nc_distance_dilation_scale (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble Dt = nc_distance_transverse (dist, model, z);
  gdouble E = sqrt (nc_hicosmo_E2 (model, z));
  gdouble Dv = cbrt (Dt * Dt * z / E);
  return Dv;
}

/***************************************************************************
 * Bao 'A' scaleD_v(z) sqrt(Omega_m) / z -- (arXiv:astro-ph/0501171)
 *
 ****************************************************************************/

gdouble
nc_distance_bao_A_scale (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble Dv = nc_distance_dilation_scale (dist, model, z);
  gdouble sqrt_Omega_m = sqrt (nc_hicosmo_Omega_m (model));
  return sqrt_Omega_m * Dv / z;
}

/***************************************************************************
 * r(z_d) / D_v(z) -- (arXiv:0705.3323)
 *
 ****************************************************************************/

gdouble
nc_distance_bao_r_Dv (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble zd = nc_distance_drag_redshift (dist, model);
  gdouble r_zd;
  gdouble Dv;

  if (!gsl_finite (zd))
	return GSL_NAN;

  r_zd = nc_distance_sound_horizon (dist, model, zd);
  Dv = nc_distance_dilation_scale (dist, model, z);
  return r_zd / Dv;
}

NcDistance *
nc_distance_new (gdouble z_f)
{
  return g_object_new (NC_TYPE_DISTANCE, "zf", z_f, NULL);
}

/**
 * nc_distance_ref:
 * @dist: a #NcDistance.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcDistance *
nc_distance_ref (NcDistance *dist)
{
  return g_object_ref (dist);
}

/**
 * nc_distance_free:
 * @dist: a #NcDistance.
 *
 * FIXME
 *
 */
void
nc_distance_free (NcDistance *dist)
{
  g_object_unref (dist);
}

/**
 * nc_distance_clear:
 * @dist: a #NcDistance.
 *
 * FIXME
 *
 */
void
nc_distance_clear (NcDistance **dist)
{
  g_clear_object (dist);
}

void
nc_distance_prepare (NcDistance *dist, NcHICosmo *model)
{
  dist->comoving_distance_cache->clear = TRUE;
  dist->time_cache->clear = TRUE;
  dist->lookback_time_cache->clear = TRUE;
  dist->conformal_time_cache->clear = TRUE;

  dist->sound_horizon_cache->clear = TRUE;

  if (dist->comoving_distance_spline == NULL)
  {
    NcmSpline *s = ncm_spline_cubic_notaknot_new ();
    dist->comoving_distance_spline =
      ncm_ode_spline_new (s, dcddz, model, 0.0, 0.0, dist->z_f);
    ncm_spline_free (s);
  }

  ncm_ode_spline_prepare (dist->comoving_distance_spline, model);

  ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (model));

  return;
}

/**
 * nc_distance_prepare_if_needed:
 * @dist: FIXME,
 * @model: FIXME
 *
 * FIXME
 */
static void
nc_distance_init (NcDistance *dist)
{
  dist->use_cache = TRUE;

  dist->comoving_distance_cache = nc_function_cache_new (1, NC_INT_ABS_ERROR, NC_INT_ERROR);

  dist->time_cache = nc_function_cache_new (1, NC_INT_ABS_ERROR, NC_INT_ERROR);
  dist->lookback_time_cache = nc_function_cache_new (1, NC_INT_ABS_ERROR, NC_INT_ERROR);
  dist->conformal_time_cache = nc_function_cache_new (1, NC_INT_ABS_ERROR, NC_INT_ERROR);

  dist->sound_horizon_cache = nc_function_cache_new (1, NC_INT_ABS_ERROR, NC_INT_ERROR);

  dist->comoving_distance_spline = NULL;

  dist->ctrl = ncm_model_ctrl_new (NULL);
}

static void
nc_distance_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_distance_parent_class)->constructed (object);
	{

	}
}

static void
nc_distance_dispose (GObject *object)
{
  NcDistance *dist = NC_DISTANCE (object);

  nc_function_cache_clear (&dist->comoving_distance_cache);
  nc_function_cache_clear (&dist->time_cache);
  nc_function_cache_clear (&dist->lookback_time_cache);
  nc_function_cache_clear (&dist->conformal_time_cache);
  nc_function_cache_clear (&dist->sound_horizon_cache);

  ncm_ode_spline_clear (&dist->comoving_distance_spline);

  ncm_model_ctrl_clear (&dist->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_distance_parent_class)->dispose (object);
}

static void
nc_distance_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_distance_parent_class)->finalize (object);
}

static void
nc_distance_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDistance *dist = NC_DISTANCE (object);
  g_return_if_fail (NC_IS_DISTANCE (object));

  switch (prop_id)
  {
	case PROP_ZF:
	  dist->z_f = g_value_get_double (value);
	  ncm_model_ctrl_force_update (dist->ctrl);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
nc_distance_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDistance *dist = NC_DISTANCE (object);
  g_return_if_fail (NC_IS_DISTANCE (object));

  switch (prop_id)
  {
	case PROP_ZF:
	  g_value_set_double (value, dist->z_f);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
nc_distance_class_init (NcDistanceClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->dispose = nc_distance_dispose;
  object_class->finalize = nc_distance_finalize;
  object_class->set_property = nc_distance_set_property;
  object_class->get_property = nc_distance_get_property;
  object_class->constructed  = nc_distance_constructed;

  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        "",
                                                        "Final cached redshift",
                                                        0.0,
                                                        G_MAXDOUBLE,
                                                        10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/***************************************************************************
 *            cosmic_time.c
 *
 *  Wed Nov 12 17:06:27 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/

static gdouble
_nc_time_integrand (gdouble z, gpointer p)
{
  NcHICosmo *model = NC_HICOSMO (p);
  const gdouble x = 1.0 + z;
  const gdouble E = sqrt (nc_hicosmo_E2 (model, z));
  return 1.0 / (x * E);
}

gdouble
nc_distance_cosmic_time (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble result, error;
  gsl_function F;
  F.function = &_nc_time_integrand;
  F.params = model;

  if (dist->use_cache)
	nc_integral_cached_x_inf (dist->time_cache, &F, z, &result, &error);
  else
	nc_integral_locked_a_inf (&F, z, NC_INT_ABS_ERROR, NC_INT_ERROR, &result, &error);

  return result;
}

gdouble
nc_distance_lookback_time (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble result, error;
  gsl_function F;
  F.function = &_nc_time_integrand;
  F.params = model;

  if (dist->use_cache)
	nc_integral_cached_0_x (dist->lookback_time_cache, &F, z, &result, &error);
  else
	nc_integral_locked_a_b (&F, 0.0, z, 0.0, NC_INT_ERROR, &result, &error);

  return result;
}

gdouble
nc_distance_conformal_lookback_time (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  return nc_distance_comoving (dist, model, z);
}

static gdouble
nc_conformal_time_integrand (gdouble logx, gpointer p)
{
  if (logx > GSL_LOG_DBL_MAX)
	return 0.0;
  else
  {
	NcHICosmo *model = NC_HICOSMO (p);
	const gdouble z = expm1 (logx);
	const gdouble x = 1.0 + z;
	const gdouble E = sqrt (nc_hicosmo_E2 (model, z));

	if (!gsl_finite (E))
	  return 0.0;
	else
	  return x / E;
  }
}

gdouble
nc_distance_conformal_time (NcDistance *dist, NcHICosmo *model, gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &nc_conformal_time_integrand;
  F.params = model;

  if (dist->use_cache)
	nc_integral_cached_x_inf (dist->conformal_time_cache, &F, log1p(z), &result, &error);
  else
	nc_integral_locked_a_inf (&F, log1p(z), NC_INT_ABS_ERROR, NC_INT_ERROR, &result, &error);

  return result;
}

typedef struct _NcDistanceFuncData
{
  NcDistance *dist;
  NcDistanceFunc0 f0;
  NcDistanceFunc1 f1;
} NcDistanceFuncData;

static void
_nc_distance_func0 (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcDistanceFuncData *dist_data = (NcDistanceFuncData *)obj;
  f[0] = dist_data->f0 (dist_data->dist, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)));
}

static void
_nc_distance_func1 (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcDistanceFuncData *dist_data = (NcDistanceFuncData *)obj;
  f[0] = dist_data->f1 (dist_data->dist, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)), x[0]);
}

static void
_nc_distance_free (gpointer obj)
{
  NcDistanceFuncData *dist_data = (NcDistanceFuncData *)obj;
  nc_distance_free (dist_data->dist);
  g_slice_free (NcDistanceFuncData, dist_data);
}

/**
 * nc_distance_func0_new:
 * @dist: FIXME
 * @f0: (scope notified): FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
nc_distance_func0_new (NcDistance *dist, NcDistanceFunc0 f0)
{
  NcDistanceFuncData *dist_data = g_slice_new (NcDistanceFuncData);
  dist_data->dist = nc_distance_ref (dist);
  dist_data->f0 = f0;
  dist_data->f1 = NULL;
  return ncm_mset_func_new (&_nc_distance_func0, 0, 1, dist_data, &_nc_distance_free);
}

/**
 * nc_distance_func1_new:
 * @dist: FIXME
 * @f1: (scope notified): FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
nc_distance_func1_new (NcDistance *dist, NcDistanceFunc1 f1)
{
  NcDistanceFuncData *dist_data = g_slice_new (NcDistanceFuncData);
  dist_data->dist = nc_distance_ref (dist);
  dist_data->f0 = NULL;
  dist_data->f1 = f1;
  return ncm_mset_func_new (&_nc_distance_func1, 1, 1, dist_data, &_nc_distance_free);
}
