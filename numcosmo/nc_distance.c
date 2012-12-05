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
 * @short_description: Calculates the cosmological distances and related quantities.
 *
 * This object implements several distances used in cosmology, here we have
 * the following definitions.
 * 
 * The Hubble scale is simply defined as the inverse of the Hubble function
 * \begin{equation}\label{eq:def:DH}
 * D_H(z) = \frac{c}{H(z)}, \qquad D_{H0} = \frac{c}{H_0}. 
 * \end{equation}
 * where $c$ is the speed of light [ncm_c_c()], $H(z)$ is the Hubble function 
 * and $H_0 \equiv H(0)$ is the Hubble parameter. The comoving distance $D_c$ is 
 * defined as
 * \begin{equation}\label{eq:def:Dc}
 * D_c(z) = \int_0^z\frac{1}{E},
 * \end{equation}
 * where $z$ is the redshift and $E$ is the Hubble function over the Hubble 
 * parameter, i.e., $$E \equiv \frac{H(z)}{H_0}.$$ Note that both quantities are 
 * adimensional, in other words, one can think them as in units of the Hubble 
 * scale today.
 * 
 * The curvature distance is defined as
 * \begin{equation}\label{eq:def:Dk}
 * D_k(z) = \sqrt{\vert\Omega_{k0}\vert}D_c(z),
 * \end{equation}
 * where $\Omega_{k0}$ is the value of the curvature today [nc_hicosmo_Omega_k()].
 * 
 * The transverse distance $D_t$ and its derivative with respect to $z$ are 
 * given by
 * \begin{equation}\label{eq:def:Dt}
 * D_t(z) = \frac{\sinh\left(\sqrt{\Omega_{k0}}D_c(z)\right)}{\sqrt{\Omega_{k0}}}, 
 * \qquad \frac{dD_t}{dz}(z) = \frac{\cosh\left(\sqrt{\Omega_{k0}}D_c(z)\right)}{E(z)}.
 * \end{equation}
 * Using the definition above we have the luminosity distance simply
 * \begin{equation}\label{eq:def:Dl}
 * D_l = (1+z)D_t(z).
 * \end{equation}
 * Then, the distance modulus is defined by the following expression
 * \begin{equation}\label{eq:def:mu}
 * \mu(z) = 5\log_{10}(D_l(z)) + 25,
 * \end{equation}
 * where $\log_{10}$ represents the logarithm in the decimal base. Note that 
 * the distance modulus is usually defined as 
 * $$5\log_{10}(D_{H0}D_l(z)/\text{pc}) - 5,$$
 * where $\text{pc}$ is parsec [ncm_c_pc()], thus, this differs from our 
 * definition by a factor of $5\log_{10}(D_{H0}/\text{Mpc})$, where $\text{Mpc}$ 
 * is mega parsec [ncm_c_Mpc()].
 * 
 * 
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
  NcHICosmo *cosmo;
  gint diff;
} ComovingDistanceArgument;

enum
{
  PROP_0,
  PROP_ZF
};

G_DEFINE_TYPE (NcDistance, nc_distance, G_TYPE_OBJECT);

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

  object_class->dispose = nc_distance_dispose;
  object_class->finalize = nc_distance_finalize;
  object_class->set_property = nc_distance_set_property;
  object_class->get_property = nc_distance_get_property;
  object_class->constructed  = nc_distance_constructed;

  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Final cached redshift",
                                                        0.0,
                                                        G_MAXDOUBLE,
                                                        10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_distance_new:
 * @z_f: final redshift $z_f$.
 *
 * Creates a new #NcDistance object optimized to perform distance calculations
 * to redshift up to $z_f$.
 *
 * Returns: a new #NcDistance.
 */
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

static gdouble dcddz (gdouble y, gdouble x, gpointer userdata);

/**
 * nc_distance_prepare:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 */
void
nc_distance_prepare (NcDistance *dist, NcHICosmo *cosmo)
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
      ncm_ode_spline_new (s, dcddz, cosmo, 0.0, 0.0, dist->z_f);
    ncm_spline_free (s);
  }

  ncm_ode_spline_prepare (dist->comoving_distance_spline, cosmo);

  ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (cosmo));

  return;
}

/**
 * nc_distance_prepare_if_needed:
 * @dist: FIXME,
 * @cosmo: FIXME
 *
 * FIXME
 * 
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

/**
 * nc_distance_hubble:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * Calculate the curvature scale today as defined in Eq $\eqref{eq:def:DH}$ in
 * units of mega parsec $\text{Mpc}$ [ncm_c_Mpc()].
 *
 * Returns: $D_{H0}$.
 */
gdouble
nc_distance_hubble (NcDistance *dist, NcHICosmo *cosmo)
{
  return ncm_c_c () / (nc_hicosmo_H0 (cosmo) * 1.0e3);
}

static gdouble comoving_distance_integral_argument(gdouble z, gpointer p);

/**
 * nc_distance_comoving:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 * 
 * Calculates the comoving distance $D_c$ as defined in Eq. $\eqref{eq:def:Dc}$.
 *
 * Returns: $D_c(z)$.
 */
gdouble
nc_distance_comoving (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble result, error;
  gsl_function F;

  nc_distance_prepare_if_needed (dist, cosmo);

  if (ncm_model_impl (NCM_MODEL (cosmo)) & NC_HICOSMO_IMPL_cd)
    return nc_hicosmo_cd (cosmo, z);

  if (z <= dist->z_f)
    return ncm_spline_eval (dist->comoving_distance_spline->s, z);

  F.function = &comoving_distance_integral_argument;
  F.params = cosmo;

  if (dist->use_cache)
    nc_integral_cached_0_x (dist->comoving_distance_cache, &F, z, &result, &error);
  else
    nc_integral_locked_a_b (&F, 0.0, z, 0.0, NC_INT_ERROR, &result, &error);

  return result;
}

static gdouble
comoving_distance_integral_argument(gdouble z, gpointer p)
{
  NcHICosmo *cosmo = NC_HICOSMO (p);
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
  if (GSL_SIGN(E2) == -1.0)
    return GSL_POSINF;
  return 1.0 / sqrt (E2);
}

static gdouble
dcddz (gdouble cd, gdouble z, gpointer userdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (userdata);
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
  return 1.0 / sqrt (E2);
}

/**
 * nc_distance_curvature:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculate the curvature distance $D_k$ defined in Eq. $\eqref{eq:def:Dk}$.
 *
 * Returns: $D_k(z)$.
 */
gdouble
nc_distance_curvature (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble sqrt_Omega_k = sqrt(fabs(nc_hicosmo_Omega_k (cosmo)));
  return sqrt_Omega_k * nc_distance_comoving (dist, cosmo, z);
}

/**
 * nc_distance_transverse:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculates the transverse distance $D_t$ defined in Eq. $\eqref{eq:def:Dt}$.
 *
 * Returns: $D_T(z)$.
 */
gdouble
nc_distance_transverse (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k = nc_hicosmo_Omega_k (cosmo);
  const gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  const gdouble comoving_dist = nc_distance_comoving (dist, cosmo, z);
  const gint k = fabs (Omega_k) < NC_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);
  gdouble Dt;

  if (gsl_isinf (comoving_dist)) 
    return comoving_dist;

  switch (k)
  {
    case 0:
      Dt = comoving_dist;
      break;
    case -1:
      Dt = sinh (sqrt_Omega_k * comoving_dist) / sqrt_Omega_k;
      break;
    case 1:
      Dt = sin (sqrt_Omega_k * comoving_dist) / sqrt_Omega_k;
      break;
    default:
      g_assert_not_reached();
      Dt = 0.0;
      break;
  }

  g_assert (Dt >= 0);
  
  return Dt;
}

/**
 * nc_distance_dtransverse_dz:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculates the derivative of $D_t$ with respect to $z$ defined in 
 * Eq. $\eqref{eq:def:Dt}$.
 *
 * Returns: $\frac{dD_t(z)}{dz}$.
 */
gdouble
nc_distance_dtransverse_dz (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble Omega_k = nc_hicosmo_Omega_k (cosmo);
  gdouble sqrt_Omega_k = sqrt (fabs (Omega_k));
  gdouble E = sqrt (nc_hicosmo_E2(cosmo, z));
  gint k = fabs (Omega_k) < NC_ZERO_LIMIT ? 0 : (Omega_k > 0.0 ? -1 : 1);

  switch (k)
  {
    case 0:
      return 1.0 / E;
      break;
    case -1:
    {
      gdouble comoving_dist = nc_distance_comoving (dist, cosmo, z);
      return cosh (sqrt_Omega_k * comoving_dist) / E;
      break;
    }
    case 1:
    {
      gdouble comoving_dist = nc_distance_comoving (dist, cosmo, z);
      return ncm_c_sign_sin (sqrt_Omega_k * comoving_dist) * cos (sqrt_Omega_k * comoving_dist) / E; // LOOK
      break;
    }
    default:
      g_assert_not_reached();
      return 0.0;
      break;
  }
}

/**
 * nc_distance_luminosity:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculates the luminosity distance $D_l$  defined in Eq. $\eqref{eq:def:Dl}$.
 *
 * Returns: $D_l(z)$.
 */
gdouble
nc_distance_luminosity (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble Dt = nc_distance_transverse (dist, cosmo, z);
  const gdouble Dl = (1.0 + z) * Dt;

  return Dl;
}

/**
 * nc_distance_modulus:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculates the distance modulus $\mu$ defined in Eq. $\eqref{eq:def:mu}$.
 *
 * Returns: $\mu(z)$.
 */
gdouble
nc_distance_modulus (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble Dl = nc_distance_luminosity (dist, cosmo, z);
  if (!gsl_finite (Dl))
    return Dl;
  return (5.0 * log10 (Dl) + 25.0);
}

/**
 * nc_distance_luminosity_hef:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z_he: the redshift $z_{he}$ in our local frame.
 * @z_cmb: the redshift $z_{CMB}$ in the CMB frame.
 *
 * Calculates the luminosity distance $D_l$ corrected to our local frame.
 *
 * Returns: $D_l(z_{hef},z_{CMB})$.
 */
gdouble 
nc_distance_luminosity_hef (NcDistance *dist, NcHICosmo *cosmo, gdouble z_he, gdouble z_cmb)
{
  const gdouble Dt = nc_distance_transverse (dist, cosmo, z_cmb);
  const gdouble Dl = (1.0 + z_he) * Dt;

  return Dl;
}

/**
 * nc_distance_modulus_hef:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z_he: the redshift $z_{he}$ in our local frame.
 * @z_cmb: the redshift $z_{CMB}$ in the CMB frame.
 *
 * Calculates the distance modulus using the frame corrected luminosity distance
 * [nc_distance_luminosity_hef()].
 *
 * Returns: $\mu(z_{hef},z_{CMB})$.
 */
gdouble 
nc_distance_modulus_hef (NcDistance *dist, NcHICosmo *cosmo, gdouble z_he, gdouble z_cmb)
{
  const gdouble Dl = nc_distance_luminosity_hef (dist, cosmo, z_he, z_cmb);
  if (!gsl_finite (Dl))
    return Dl;
  return (5.0 * log10 (Dl) + 25.0);
}

/**
 * nc_distance_angular_diameter_curvature_scale:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Calculate the angular diameter curvature scale $D_a$ defined in 
 * Eq. \eqref{eq:def:Da} at decoupling redshift $z_\star$ given by 
 * [nc_distance_decoupling_redshift()].
 *
 * Returns: $D_a(z_\star)$.
 */
gdouble
nc_distance_angular_diameter_curvature_scale (NcDistance *dist, NcHICosmo *cosmo)
{
  const gdouble z_star = nc_distance_decoupling_redshift (dist, cosmo);
  if (gsl_finite (z_star))
  {
    return sqrt (nc_hicosmo_E2 (cosmo, z_star)) *
      nc_distance_transverse (dist, cosmo, z_star) / (1.0 + z_star);
  }
  else
    return GSL_NAN;
}

/**
 * nc_distance_shift_parameter:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculate the shift parameter $R$ defined in Eq. $\eqref{eq:def:R}$.
 *
 * Returns: $R(z)$.
 */
gdouble
nc_distance_shift_parameter (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble sqrt_mod_Omega_m = sqrt (fabs (nc_hicosmo_Omega_m (cosmo)));
  const gdouble transverse = nc_distance_transverse (dist, cosmo, z);
  return sqrt_mod_Omega_m * transverse;
}

/**
 * nc_distance_shift_parameter_lss:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Calculate the shift parameter $R$ [nc_distance_shift_parameter()] at the
 * decoupling redshift $z_\star$ [nc_distance_decoupling_redshift()].
 *
 * Returns: $R(z_\star)$.
 */
gdouble
nc_distance_shift_parameter_lss (NcDistance *dist, NcHICosmo *cosmo)
{
  const gdouble sqrt_mod_Omega_m = sqrt(fabs(nc_hicosmo_Omega_m (cosmo)));
  const gdouble z_star = nc_distance_decoupling_redshift (dist, cosmo);
  if (gsl_finite (z_star))
  {
    const gdouble transverse = nc_distance_transverse (dist, cosmo, z_star);
    return sqrt_mod_Omega_m * transverse;
  }
  else
    return GSL_NAN;
}

/**
 * nc_distance_comoving_a0:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving_a0 (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble sqrt_mod_Omega_k = sqrt(fabs(nc_hicosmo_Omega_k (cosmo)));
  const gdouble comoving = nc_distance_comoving (dist, cosmo, z);
  return sqrt_mod_Omega_k * comoving;
}

/**
 * nc_distance_comoving_a0_lss:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * Calculate the d_c (z_lss) / a0
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving_a0_lss (NcDistance *dist, NcHICosmo *cosmo)
{
  const gdouble sqrt_mod_Omega_k = sqrt(fabs(nc_hicosmo_Omega_k (cosmo)));
  const gdouble z_star = nc_distance_decoupling_redshift (dist, cosmo);
  if (gsl_finite (z_star))
  {
    const gdouble comoving = nc_distance_comoving (dist, cosmo, z_star);
    return sqrt_mod_Omega_k * comoving;
  }
  else
    return GSL_NAN;
}

/**
 * nc_distance_comoving_lss:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * Calculate the d_c (z_lss)
 *
 * Returns: FIXME
 */
gdouble
nc_distance_comoving_lss (NcDistance *dist, NcHICosmo *cosmo)
{
  const gdouble z_star = nc_distance_decoupling_redshift (dist, cosmo);
  if (gsl_finite (z_star))
    return nc_distance_comoving (dist, cosmo, z_star);
  else
    return GSL_NAN;
}

/**
 * nc_distance_decoupling_redshift:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * Decoupling redshift (arXiv:astro-ph/9510117).
 *
 * Returns: FIXME
 */
gdouble
nc_distance_decoupling_redshift (NcDistance *dist, NcHICosmo *cosmo)
{
  if (ncm_model_impl (NCM_MODEL (cosmo)) & NC_HICOSMO_IMPL_z_lss)
    return nc_hicosmo_z_lss (cosmo);
  else
  {
    gdouble h = nc_hicosmo_h (cosmo);
    gdouble h2 = h * h;
    gdouble omega_b_h2 = nc_hicosmo_Omega_b (cosmo) * h2;
    gdouble omega_m_h2 = nc_hicosmo_Omega_m (cosmo) * h2;
    gdouble g1 = 0.0783 * pow (omega_b_h2, -0.238) / (1.0 + 39.5 * pow (omega_b_h2, 0.763));
    gdouble g2 = 0.560 / (1.0 + 21.1 * pow (omega_b_h2, 1.81));
    return 1048.0 * (1.0 + 1.24e-3 * pow (omega_b_h2, -0.738)) * (1.0 + g1 * pow (omega_m_h2, g2));
  }
}

static gdouble sound_horizon_integral_argument(gdouble z, gpointer p);

/**
 * nc_distance_sound_horizon:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculate the sound horizon integrating numerically the hubble function.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_sound_horizon (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &sound_horizon_integral_argument;
  F.params = cosmo;

  if (dist->use_cache)
    nc_integral_cached_x_inf (dist->sound_horizon_cache, &F, z, &result, &error);
  else
    nc_integral_locked_a_inf (&F, z, NC_INT_ABS_ERROR, NC_INT_ERROR, &result, &error);

  return result;
}

static gdouble
sound_horizon_integral_argument (gdouble z, gpointer p)
{
  NcHICosmo *cosmo = NC_HICOSMO (p);
  gdouble E2 = nc_hicosmo_E2 (cosmo, z);
  gdouble omega_b = nc_hicosmo_Omega_b (cosmo);
  gdouble omega_r = nc_hicosmo_Omega_r (cosmo);
  if (omega_r == 0.0)
    return 0.0;
  if (GSL_SIGN(E2) == -1.0)
    return GSL_POSINF;
  return
    1.0 / sqrt (E2 * (3.0 + 9.0 / 4.0 * (1.0 + ncm_c_neutrino_n_eff () * 0.2271) * omega_b / (omega_r * (1.0 + z))));
}

/**
 * nc_distance_dsound_horizon_dz:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Calculate the sound horizon derivative.
 *
 * Returns: FIXME
 */
gdouble
nc_distance_dsound_horizon_dz (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  return -sound_horizon_integral_argument (z, cosmo);
}

/**
 * nc_distance_acoustic_scale:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * Calculate the acoustic scale l_A = pi D_T (z_lss) / r (z_lss).
 *
 * Returns: FIXME
 */
gdouble
nc_distance_acoustic_scale (NcDistance *dist, NcHICosmo *cosmo)
{
  gdouble z = nc_distance_decoupling_redshift (dist, cosmo);
  if (gsl_finite (z))
    return M_PI * nc_distance_transverse (dist, cosmo, z) / nc_distance_sound_horizon (dist, cosmo, z);
  else
    return GSL_NAN;
}

/**
 * nc_distance_drag_redshift:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 *
 * Drag redshift (arXiv:astro-ph/9510117).
 *
 * Returns: FIXME
 */
gdouble
nc_distance_drag_redshift (NcDistance *dist, NcHICosmo *cosmo)
{
  gdouble h = nc_hicosmo_h (cosmo);
  gdouble h2 = h*h;
  gdouble omega_m_h2 = nc_hicosmo_Omega_m (cosmo) * h2;
  gdouble omega_b_h2 = nc_hicosmo_Omega_b (cosmo) * h2;
  gdouble b1 = 0.313 * pow (omega_m_h2, -0.419) * (1.0 + 0.607 * pow (omega_m_h2, 0.674));
  gdouble b2 = 0.238 * pow (omega_m_h2, 0.223);
  return 1291.0 * pow (omega_m_h2, 0.251) / (1.0 + 0.659 * pow (omega_m_h2, 0.828)) *
    (1.0 + b1 * pow (omega_b_h2, b2));
}

/**
 * nc_distance_dilation_scale:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Dilation scale D_v(z) -- (arXiv:astro-ph/0501171)
 *
 * Returns: FIXME
 */
gdouble
nc_distance_dilation_scale (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble Dt = nc_distance_transverse (dist, cosmo, z);
  gdouble E = sqrt (nc_hicosmo_E2 (cosmo, z));
  gdouble Dv = cbrt (Dt * Dt * z / E);
  return Dv;
}

/**
 * nc_distance_bao_A_scale:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * Bao 'A' scaleD_v(z) sqrt(Omega_m) / z -- (arXiv:astro-ph/0501171)
 *
 * Returns: FIXME
 */
gdouble
nc_distance_bao_A_scale (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble Dv = nc_distance_dilation_scale (dist, cosmo, z);
  gdouble sqrt_Omega_m = sqrt (nc_hicosmo_Omega_m (cosmo));
  return sqrt_Omega_m * Dv / z;
}

/**
 * nc_distance_bao_r_Dv:
 * @dist: a #NcDistance.
 * @cosmo: a #NcHICosmo.
 * @z: the redshift $z$.
 *
 * r(z_d) / D_v(z) -- (arXiv:0705.3323).
 *
 * Returns: FIXME
 */
gdouble
nc_distance_bao_r_Dv (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble zd = nc_distance_drag_redshift (dist, cosmo);
  gdouble r_zd;
  gdouble Dv;

  if (!gsl_finite (zd))
    return GSL_NAN;

  r_zd = nc_distance_sound_horizon (dist, cosmo, zd);
  Dv = nc_distance_dilation_scale (dist, cosmo, z);
  return r_zd / Dv;
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
  NcHICosmo *cosmo = NC_HICOSMO (p);
  const gdouble x = 1.0 + z;
  const gdouble E = sqrt (nc_hicosmo_E2 (cosmo, z));
  return 1.0 / (x * E);
}

gdouble
nc_distance_cosmic_time (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble result, error;
  gsl_function F;
  F.function = &_nc_time_integrand;
  F.params = cosmo;

  if (dist->use_cache)
    nc_integral_cached_x_inf (dist->time_cache, &F, z, &result, &error);
  else
    nc_integral_locked_a_inf (&F, z, NC_INT_ABS_ERROR, NC_INT_ERROR, &result, &error);

  return result;
}

gdouble
nc_distance_lookback_time (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble result, error;
  gsl_function F;
  F.function = &_nc_time_integrand;
  F.params = cosmo;

  if (dist->use_cache)
    nc_integral_cached_0_x (dist->lookback_time_cache, &F, z, &result, &error);
  else
    nc_integral_locked_a_b (&F, 0.0, z, 0.0, NC_INT_ERROR, &result, &error);

  return result;
}

gdouble
nc_distance_conformal_lookback_time (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  return nc_distance_comoving (dist, cosmo, z);
}

static gdouble
nc_conformal_time_integrand (gdouble logx, gpointer p)
{
  if (logx > GSL_LOG_DBL_MAX)
    return 0.0;
  else
  {
    NcHICosmo *cosmo = NC_HICOSMO (p);
    const gdouble z = expm1 (logx);
    const gdouble x = 1.0 + z;
    const gdouble E = sqrt (nc_hicosmo_E2 (cosmo, z));

    if (!gsl_finite (E))
      return 0.0;
    else
      return x / E;
  }
}

gdouble
nc_distance_conformal_time (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &nc_conformal_time_integrand;
  F.params = cosmo;

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
