/***************************************************************************
 *            nc_distance.c
 *
 *  Tue May  8 11:05:35 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcDistance:
 *
 * Cosmological distance and time related quantities.
 *
 * This object implements several distances used in cosmology, here we have
 * the following definitions.
 *
 * $
 *  \newcommand{\RH}{{R_H}}
 *  \newcommand{\RHc}{{R^\mathrm{c}_H}}
 * $
 *
 * The Hubble radius (or scale) is defined as the inverse of the Hubble
 * function $H(z)$ [nc_hicosmo_H()],
 * \begin{equation}\label{eq:def:RHc}
 * \RH = \frac{c}{H(z)}, \qquad \RH_0 = \frac{c}{H_0},
 * \end{equation}
 * where $c$ is the speed of light [ncm_c_c()], $z$ is the redshift and
 * $H_0 \equiv H(0)$ is the Hubble parameter [nc_hicosmo_H0()]. Similarly,
 * we also define the comoving Hubble radius as
 * \begin{equation}\label{eq:def:DH}
 * \RHc(z) = \frac{c}{aH(z)} = \frac{c(1+z)}{a_0H(z)}, \qquad \RHc_0 = \frac{c}{a_0H_0}
 * \end{equation}
 * where ${}_0$ subscript means that the function is calculated at the
 * present time and the redshift $z$ is defined by the expression
 * $$1 + z = \frac{a_0}{a}.$$
 *
 * The comoving distance $D_c$ is defined as
 * \begin{equation}\label{eq:def:dc}
 * d_c(z) = \RHc_0\int_0^z \frac{dz^\prime}{E (z^\prime)},
 * \end{equation}
 * where $E(z)$ is the normalized Hubble function [nc_hicosmo_E()], i.e.,
 * \begin{equation}\label{eq:def:Ez}
 * E(z) \equiv \frac{H(z)}{H_0}.
 * \end{equation}
 *
 * In this object we will compute the dimensionless version of the distances,
 * for the comoving distance we define
 * \begin{equation}\label{eq:def:Dc}
 * D_c(z) \equiv \frac{d_c(z)}{\RHc_0}.
 * \end{equation}
 * note, however, that $D_c(z)$ coincides with the proper distance today
 * $r(z) \equiv a_0 d_c(z)$ in unit of the Hubble radius, i.e.,
 * $D_c(z) = r(z) / \RH_0$. Therefore, both the comoving distance and
 * the proper distance today can be obtained by multiplying $D_c(z)$
 * by $\RHc_0$ and $\RH_0$ respectively.
 *
 * The transverse comoving distance $D_t$ and its derivative with respect to
 * $z$ are given by
 * \begin{equation}\label{eq:def:Dt}
 * D_t(z) = \frac{\sinh\left[\sqrt{\Omega_{k0}}D_c(z)\right]}{\sqrt{\Omega_{k0}}},
 * \qquad \frac{dD_t}{dz}(z) = \frac{\cosh\left[\sqrt{\Omega_{k0}}D_c(z)\right]}{E(z)},
 * \end{equation}
 * where $\Omega_{k0}$ is the value of the curvature today [nc_hicosmo_Omega_k0()].
 * Using the definition above we have that the luminosity and angular diameter distances
 * are respectively:
 * \begin{equation}\label{eq:def:Dl}
 * D_l = (1 + z)D_t(z), \qquad D_A = D_t (z) / (1 + z),
 * \end{equation}
 * and the distance modulus is given by
 * \begin{equation}\label{eq:def:dmu}
 * \delta\mu(z) = 5\log_{10}(D_l(z)) + 25.
 * \end{equation}
 * Note that the distance modulus is defined as
 * $$\mu(z) = 5\log_{10}[\RH_0D_l(z)/(1\,\text{Mpc})] + 25.$$
 * Thus, this differs from our definition by a factor of
 * $5\log_{10}[\RH_0/(1\,\text{Mpc})]$, i.e.,
 * $$\mu(z) = \delta\mu(z) + 5\log_{10}[\RH_0/(1\,\text{Mpc})],$$
 * where $\text{Mpc}$ is megaparsec [ncm_c_Mpc()].
 *
 * We also implement the following distances between $z_1$ and $z_2$:
 * \begin{align}\label{eq:def:Dt12}
 * D_t(z_1, z_2) &= \frac{\sinh\left\\{\sqrt{\Omega_{k0}}\left[D_c(z_2)-D_c(z_1)\right]\right\\}}{\sqrt{\Omega_{k0}}}, \\\\ \label{eq:def:DA12}
 * D_A(z_1, z_2) &= D_t (z_1, z_2) / (1 + z_2).
 * \end{align}
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_distance.h"
#include "math/ncm_integrate.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_mset_func_list.h"

typedef struct _ComovingDistanceArgument
{
  NcHICosmo *cosmo;
  gint diff;
} ComovingDistanceArgument;

enum
{
  PROP_0,
  PROP_ZF,
  PROP_RECOMB,
  PROP_INV_COMOVING,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDistance, nc_distance, G_TYPE_OBJECT)

static void
nc_distance_init (NcDistance *dist)
{
  dist->use_cache = TRUE;

  dist->comoving_distance_cache = ncm_function_cache_new (1, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR);

  dist->comoving_infinity    = ncm_function_cache_new (1, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR);
  dist->time_cache           = ncm_function_cache_new (1, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR);
  dist->lookback_time_cache  = ncm_function_cache_new (1, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR);
  dist->conformal_time_cache = ncm_function_cache_new (1, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR);

  dist->sound_horizon_cache = ncm_function_cache_new (1, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR);

  dist->comoving_distance_spline = NULL;
  dist->inv_comoving_dist        = NULL;
  dist->cpu_inv_comoving         = FALSE;
  dist->max_comoving             = 0.0;

  dist->recomb = NULL;

  dist->cmethod = NC_DISTANCE_COMOVING_METHOD_LEN;

  dist->ctrl = ncm_model_ctrl_new (NULL);
}

static void
_nc_distance_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDistance *dist = NC_DISTANCE (object);

  g_return_if_fail (NC_IS_DISTANCE (object));

  switch (prop_id)
  {
    case PROP_ZF:
      dist->zf = g_value_get_double (value);
      ncm_model_ctrl_force_update (dist->ctrl);
      break;
    case PROP_RECOMB:
      nc_distance_set_recomb (dist, g_value_get_object (value));
      break;
    case PROP_INV_COMOVING:
      nc_distance_compute_inv_comoving (dist, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_distance_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDistance *dist = NC_DISTANCE (object);

  g_return_if_fail (NC_IS_DISTANCE (object));

  switch (prop_id)
  {
    case PROP_ZF:
      g_value_set_double (value, dist->zf);
      break;
    case PROP_RECOMB:
      g_value_set_object (value, dist->recomb);
      break;
    case PROP_INV_COMOVING:
      g_value_set_boolean (value, dist->cpu_inv_comoving);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble _dcddz (gdouble y, gdouble x, gpointer userdata);

static void
_nc_distance_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_distance_parent_class)->constructed (object);
  {
    NcDistance *dist = NC_DISTANCE (object);
    NcmSpline *s     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

    dist->comoving_distance_spline = ncm_ode_spline_new_full (s, _dcddz, 0.0, 0.0, dist->zf);
    ncm_ode_spline_auto_abstol (dist->comoving_distance_spline, TRUE);
    ncm_spline_free (s);
  }
}

static void
_nc_distance_dispose (GObject *object)
{
  NcDistance *dist = NC_DISTANCE (object);

  ncm_function_cache_clear (&dist->comoving_distance_cache);
  ncm_function_cache_clear (&dist->comoving_infinity);
  ncm_function_cache_clear (&dist->time_cache);
  ncm_function_cache_clear (&dist->lookback_time_cache);
  ncm_function_cache_clear (&dist->conformal_time_cache);
  ncm_function_cache_clear (&dist->sound_horizon_cache);

  ncm_ode_spline_clear (&dist->comoving_distance_spline);
  ncm_spline_clear (&dist->inv_comoving_dist);

  ncm_model_ctrl_clear (&dist->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_distance_parent_class)->dispose (object);
}

static void
_nc_distance_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_distance_parent_class)->finalize (object);
}

static void
nc_distance_class_init (NcDistanceClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_nc_distance_constructed;
  object_class->set_property = &_nc_distance_set_property;
  object_class->get_property = &_nc_distance_get_property;
  object_class->dispose      = &_nc_distance_dispose;
  object_class->finalize     = &_nc_distance_finalize;

  /**
   * NcDistance:zf:
   *
   * The final redhift to compute cosmological distances and related quantities.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Final cached redshift",
                                                        0.0,
                                                        G_MAXDOUBLE,
                                                        10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDistance:recomb:
   *
   * The recombination object, #NcRecomb, used by #NcDistance.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RECOMB,
                                   g_param_spec_object ("recomb",
                                                        NULL,
                                                        "Recombination object",
                                                        NC_TYPE_RECOMB,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INV_COMOVING,
                                   g_param_spec_boolean ("compute-inv-comoving",
                                                         NULL,
                                                         "Whether to compute the inverse comoving function",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_distance_new:
 * @zf: final redshift
 *
 * Creates a new #NcDistance object optimized to perform distance calculations
 * to redshift up to @zf.
 *
 * Returns: a new #NcDistance.
 */
NcDistance *
nc_distance_new (gdouble zf)
{
  return g_object_new (NC_TYPE_DISTANCE, "zf", zf, NULL);
}

/**
 * nc_distance_ref:
 * @dist: a #NcDistance
 *
 * Increases the reference count of @dist atomically.
 *
 * Returns: (transfer full): @dist.
 */
NcDistance *
nc_distance_ref (NcDistance *dist)
{
  return g_object_ref (dist);
}

/**
 * nc_distance_free:
 * @dist: a #NcDistance
 *
 * Atomically decrements the reference count of @dist by one.
 * If the reference count drops to 0, all memory allocated by @dist is released.
 *
 */
void
nc_distance_free (NcDistance *dist)
{
  g_object_unref (dist);
}

/**
 * nc_distance_clear:
 * @dist: a #NcDistance
 *
 * Atomically decrements the reference count of @dist by one.
 * If the reference count drops to 0, all memory allocated by @dist is released.
 * Set pointer to NULL.
 *
 */
void
nc_distance_clear (NcDistance **dist)
{
  g_clear_object (dist);
}

/**
 * nc_distance_require_zf:
 * @dist: a #NcDistance
 * @zf: maximum redshift required
 *
 * Requires the final redshift of at least @zf.
 *
 */
void
nc_distance_require_zf (NcDistance *dist, const gdouble zf)
{
  if (zf > dist->zf)
  {
    ncm_ode_spline_set_xf (dist->comoving_distance_spline, dist->zf);
    dist->zf = zf;

    ncm_model_ctrl_force_update (dist->ctrl);
  }
}

/**
 * nc_distance_set_recomb:
 * @dist: a #NcDistance
 * @recomb: a #NcRecomb
 *
 * This function sets @recomb into @dist.
 *
 */
void
nc_distance_set_recomb (NcDistance *dist, NcRecomb *recomb)
{
  if (dist->recomb != recomb)
  {
    nc_recomb_clear (&dist->recomb);

    if (recomb != NULL)
      dist->recomb = nc_recomb_ref (recomb);

    ncm_model_ctrl_force_update (dist->ctrl);
  }
}

/**
 * nc_distance_compute_inv_comoving:
 * @dist: a #NcDistance
 * @cpu_inv_xi: whether to compute the inverse function $z(\xi)$
 *
 * Enable/Disable the computation of $z(\xi)$
 */
void
nc_distance_compute_inv_comoving (NcDistance *dist, gboolean cpu_inv_xi)
{
  if ((dist->cpu_inv_comoving && !cpu_inv_xi) || (!dist->cpu_inv_comoving && cpu_inv_xi))
  {
    dist->cpu_inv_comoving = cpu_inv_xi;
    ncm_model_ctrl_force_update (dist->ctrl);
  }
}

/**
 * nc_distance_prepare:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * This function prepares the object @dist using @cosmo,
 * such that all the available distances functions can be evaluated,
 * e.g. nc_distance_comoving().
 *
 */
void
nc_distance_prepare (NcDistance *dist, NcHICosmo *cosmo)
{
  ncm_function_cache_empty_cache (dist->comoving_distance_cache);
  ncm_function_cache_empty_cache (dist->comoving_infinity);
  ncm_function_cache_empty_cache (dist->time_cache);
  ncm_function_cache_empty_cache (dist->lookback_time_cache);
  ncm_function_cache_empty_cache (dist->conformal_time_cache);
  ncm_function_cache_empty_cache (dist->sound_horizon_cache);

  if (ncm_model_check_impl_opt (NCM_MODEL (cosmo), NC_HICOSMO_IMPL_Dc))
  {
    dist->cmethod = NC_DISTANCE_COMOVING_METHOD_FROM_MODEL;
  }
  else
  {
    ncm_ode_spline_prepare (dist->comoving_distance_spline, cosmo);
    dist->cmethod = NC_DISTANCE_COMOVING_METHOD_INT_E;
  }

  if (dist->cpu_inv_comoving)
  {
    NcmSpline *s      = ncm_ode_spline_peek_spline (dist->comoving_distance_spline);
    NcmVector *z_vec  = ncm_spline_get_xv (s);
    NcmVector *xi_vec = ncm_spline_get_yv (s);

    ncm_spline_clear (&dist->inv_comoving_dist);
    dist->inv_comoving_dist = ncm_spline_copy_empty (s);

    ncm_spline_set (dist->inv_comoving_dist, xi_vec, z_vec, TRUE);
    dist->max_comoving = ncm_vector_get (xi_vec, ncm_vector_len (xi_vec) - 1);

    ncm_vector_free (z_vec);
    ncm_vector_free (xi_vec);
  }

  if (dist->recomb != NULL)
    nc_recomb_prepare_if_needed (dist->recomb, cosmo);

  ncm_model_ctrl_update (dist->ctrl, NCM_MODEL (cosmo));

  return;
}

/**
 * nc_distance_prepare_if_needed:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * This function prepares the object @dist using @cosmo
 * if it was changed since last preparation.
 *
 */

/**
 * nc_distance_hubble:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Calculate the curvature scale today as defined in Eq $\eqref{eq:def:DH}$ in
 * units of megaparsec (Mpc) [ncm_c_Mpc()].
 *
 * Returns: $D_{H0}$.
 */
gdouble
nc_distance_hubble (NcDistance *dist, NcHICosmo *cosmo)
{
  NCM_UNUSED (dist);

  return ncm_c_c () / (nc_hicosmo_H0 (cosmo) * 1.0e3);
}

static gdouble _comoving_distance_integral_argument (gdouble z, gpointer p);

/**
 * nc_distance_comoving:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Calculate the comoving distance $D_c (z)$ as defined in Eq. $\eqref{eq:def:Dc}$.
 *
 * Returns: $D_c(z)$.
 */
gdouble
nc_distance_comoving (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  switch (dist->cmethod)
  {
    case NC_DISTANCE_COMOVING_METHOD_FROM_MODEL:

      return nc_hicosmo_Dc (cosmo, z);

      break;
    case NC_DISTANCE_COMOVING_METHOD_INT_E:
    {
      if (z <= dist->zf)
      {
        return ncm_spline_eval (ncm_ode_spline_peek_spline (dist->comoving_distance_spline), z);
      }
      else
      {
        gdouble result, error;
        gsl_function F;

        F.function = &_comoving_distance_integral_argument;
        F.params   = cosmo;

        if (dist->use_cache)
          ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, z, &result, &error);
        else
          ncm_integral_locked_a_b (&F, 0.0, z, 0.0, NCM_INTEGRAL_ERROR, &result, &error);

        return result;
      }

      break;
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  return GSL_NAN;
}

static gdouble
_comoving_distance_integral_argument (gdouble z, gpointer p)
{
  NcHICosmo *cosmo = NC_HICOSMO (p);
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);

  if (GSL_SIGN (E2) == -1.0)
    return GSL_POSINF;

  return 1.0 / sqrt (E2);
}

static gdouble
_dcddz (gdouble cd, const gdouble z, gpointer userdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (userdata);
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);

  return 1.0 / sqrt (E2);
}

static gdouble
_nc_distance_sinn (const gdouble r, const gdouble Omega_k0)
{
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);

  switch (k)
  {
    case 0:

      return r;

      break;
    case -1:

      return sinh (sqrt_Omega_k0 * r) / sqrt_Omega_k0;

      break;
    case 1:

      return fabs (sin (sqrt_Omega_k0 * r) / sqrt_Omega_k0);

      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_distance_transverse:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the transverse comoving distance $D_t (z)$ defined in Eq. $\eqref{eq:def:Dt}$.
 *
 * Returns: $D_t(z)$.
 */
gdouble
nc_distance_transverse (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble comoving_dist = nc_distance_comoving (dist, cosmo, z);

  if (gsl_isinf (comoving_dist))
    return comoving_dist;

  return _nc_distance_sinn (comoving_dist, Omega_k0);
}

/**
 * nc_distance_comoving_z1_z2:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z1: redshift $z_1$
 * @z2: redshift $z_2$
 *
 * Compute the comoving distance between $z_1$ and $z_2$, $D_c (z_1, z_2)$.
 *
 * Returns: $D_c(z_1, z_2)$.
 */
gdouble
nc_distance_comoving_z1_z2 (NcDistance *dist, NcHICosmo *cosmo, const gdouble z1, const gdouble z2)
{
  const gdouble dc1 = nc_distance_comoving (dist, cosmo, z1);
  const gdouble dc2 = nc_distance_comoving (dist, cosmo, z2);
  const gdouble dc  = dc2 - dc1;

  if (fabs (2.0 * dc / (dc1 + dc2)) > 1.0e-3)
  {
    return dc;
  }
  else
  {
    gsl_function F = {&_comoving_distance_integral_argument, cosmo};
    gdouble result, error;

    ncm_integral_locked_a_b (&F, z1, z2, 0.0, NCM_INTEGRAL_ERROR, &result, &error);

    return result;
  }
}

/**
 * nc_distance_transverse_z1_z2:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z1: redshift $z_1$
 * @z2: redshift $z_2$
 *
 * Compute the transverse comoving distance between $z_1$ and $z_2$, $D_t (z1, z2)$ defined in Eq. $\eqref{eq:def:Dt12}$.
 *
 * Returns: $D_t(z_1, z_2)$.
 */
gdouble
nc_distance_transverse_z1_z2 (NcDistance *dist, NcHICosmo *cosmo, const gdouble z1, const gdouble z2)
{
  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble dc       = nc_distance_comoving_z1_z2 (dist, cosmo, z1, z2);

  return _nc_distance_sinn (dc, Omega_k0);
}

/**
 * nc_distance_dtransverse_dz:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the derivative of $D_t(z)$ with respect to $z$ defined in
 * Eq. $\eqref{eq:def:Dt}$.
 *
 * Returns: $\frac{dD_t(z)}{dz}$.
 */
gdouble
nc_distance_dtransverse_dz (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  gdouble E             = sqrt (nc_hicosmo_E2 (cosmo, z));
  gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);

  switch (k)
  {
    case 0:

      return 1.0 / E;

      break;
    case -1:
    {
      gdouble comoving_dist = nc_distance_comoving (dist, cosmo, z);

      return cosh (sqrt_Omega_k0 * comoving_dist) / E;

      break;
    }
    case 1:
    {
      gdouble comoving_dist = nc_distance_comoving (dist, cosmo, z);

      return ncm_c_sign_sin (sqrt_Omega_k0 * comoving_dist) * cos (sqrt_Omega_k0 * comoving_dist) / E; /* LOOK */

      break;
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_distance_luminosity:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the luminosity distance $D_l(z)$  defined in Eq. $\eqref{eq:def:Dl}$.
 *
 * Returns: $D_l(z)$.
 */
gdouble
nc_distance_luminosity (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble Dt = nc_distance_transverse (dist, cosmo, z);
  const gdouble Dl = (1.0 + z) * Dt;

  return Dl;
}

/**
 * nc_distance_angular_diameter:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the angular diameter $D_A(z)$, defined in Eq. $\eqref{eq:def:Dl}$.
 *
 * Returns: $D_A(z)$.
 */
gdouble
nc_distance_angular_diameter (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble Dt = nc_distance_transverse (dist, cosmo, z);
  const gdouble DA = Dt / (1.0 + z);

  return DA;
}

/**
 * nc_distance_angular_diameter_z1_z2:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z1: redshift $z_1$
 * @z2: redshift $z_2$
 *
 * Compute the angular diameter $D_A(z_1, z_2)$, defined in Eq. $\eqref{eq:def:DA12}$.
 *
 * Returns: $D_A(z_1, z_2)$.
 */
gdouble
nc_distance_angular_diameter_z1_z2 (NcDistance *dist, NcHICosmo *cosmo, const gdouble z1, const gdouble z2)
{
  const gdouble Dt = nc_distance_transverse_z1_z2 (dist, cosmo, z1, z2);
  const gdouble DA = Dt / (1.0 + z2);

  return DA;
}

/**
 * nc_distance_dmodulus:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the distance modulus $\delta\mu(z)$ defined in Eq. $\eqref{eq:def:dmu}$.
 *
 * Returns: $\delta\mu(z)$.
 */
gdouble
nc_distance_dmodulus (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble Dl = nc_distance_luminosity (dist, cosmo, z);

  if (!gsl_finite (Dl))
    return Dl;

  return (5.0 * log10 (Dl) + 25.0);
}

/**
 * nc_distance_comoving_volume_element:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the comoving volume element per unit solid angle $d\Omega$
 * given @z, namely, $$\frac{\mathrm{d}^2V}{\mathrm{d}z\mathrm{d}\Omega} = \frac{D_t^2(z)}{E(z)} = \frac{(1 + z)^2 D_a^2(z)}{E(z)},$$
 * where $E(z)$ is the normalized Hubble function and $D_t$ is the transverse comoving distance (and $D_a$ is the angular diameter distance).
 *
 * Returns: comoving volume element $\frac{\mathrm{d}^2V}{\mathrm{d}z\mathrm{d}\Omega}$.
 */
gdouble
nc_distance_comoving_volume_element (NcDistance *dist, NcHICosmo *cosmo, gdouble z)
{
  const gdouble E     = sqrt (nc_hicosmo_E2 (cosmo, z));
  gdouble Dt          = nc_distance_transverse (dist, cosmo, z);
  gdouble dV_dzdOmega = gsl_pow_2 (Dt) / E;

  return dV_dzdOmega;
}

/**
 * nc_distance_sigma_critical:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @zs: source redshift $z_\mathrm{source}$
 * @zl: lens redshift $z_\mathrm{lens}$
 *
 * Computes the critical surface density,
 * \begin{equation}\label{eq:def:SigmaC}
 * \Sigma_c = \frac{c^2}{4\pi G} \frac{D_s}{D_l D_{ls}},
 * \end{equation}
 * where $c^2$ is the speed of light squared [ncm_c_c2 ()], $G$ is the gravitational constant in units of $m^3/s^2 M_\odot^{-1}$ [ncm_c_G_mass_solar()],
 * $D_s$ ($D_l$) is the angular diameter distance from the observer to the source (lens), and $D_{ls}$ is the angular diameter distance between
 * the lens and the source.
 *
 * Returns: the critical surface density $\Sigma_c$ in units of $M_\odot / Mpc^2$
 */
gdouble
nc_distance_sigma_critical (NcDistance *dist, NcHICosmo *cosmo, const gdouble zs, const gdouble zl)
{
  if (zs < zl)
  {
    return GSL_POSINF;
  }
  else
  {
    const gdouble a   = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
    const gdouble Ds  = nc_distance_angular_diameter (dist, cosmo, zs);
    const gdouble Dl  = nc_distance_angular_diameter (dist, cosmo, zl);
    const gdouble Dls = nc_distance_angular_diameter_z1_z2 (dist, cosmo, zl, zs);

    const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);

    return a * Ds / (Dl * Dls * RH_Mpc);
  }
}

/**
 * nc_distance_sigma_critical_infinity:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @zl: lens redshift $z_\mathrm{lens}$
 *
 * Computes the critical surface density,
 * \begin{equation}\label{eq:def:SigmaC}
 * \Sigma_c = \frac{c^2}{4\pi G} \frac{D_\infty}{D_l D_{l\infty}},
 * \end{equation}
 * where $c^2$ is the speed of light squared [ncm_c_c2 ()], $G$ is the gravitational constant in units of $m^3/s^2 M_\odot^{-1}$ [ncm_c_G_mass_solar()],
 * $D_\infty$ ($D_l$) is the angular diameter distance from the observer to the source at infinite redshift (lens), and $D_{l\infty}$ is the angular diameter
 * the lens and the source.
 *
 * Returns: the critical surface density $\Sigma_c$ in units of $M_\odot / Mpc^2$
 */
gdouble
nc_distance_sigma_critical_infinity (NcDistance *dist, NcHICosmo *cosmo, const gdouble zl)
{
  /*g_assert_cmpfloat (nc_hicosmo_Omega_k0 (cosmo), >=, 0.0); */
  const gdouble a     = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc (); /* [ M_solar / Mpc ] */
  const gdouble Dinf  = nc_distance_transverse_z_to_infinity (dist, cosmo, 0.0);
  const gdouble Dl    = nc_distance_angular_diameter (dist, cosmo, zl);
  const gdouble Dlinf = nc_distance_transverse_z_to_infinity (dist, cosmo, zl);

  const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);

  return a * Dinf / (Dl * Dlinf * RH_Mpc);
}

/**
 * nc_distance_luminosity_hef:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z_he: redshift $z_{he}$ in our local frame
 * @z_cmb: redshift $z_{CMB}$ in the CMB frame
 *
 * Calculate the luminosity distance $D_l$ corrected to our local frame.
 *
 * Returns: $D_l(z_{hef},z_{CMB})$.
 */
gdouble
nc_distance_luminosity_hef (NcDistance *dist, NcHICosmo *cosmo, const gdouble z_he, const gdouble z_cmb)
{
  const gdouble Dt = nc_distance_transverse (dist, cosmo, z_cmb);
  const gdouble Dl = (1.0 + z_he) * Dt;

  return Dl;
}

/**
 * nc_distance_dmodulus_hef:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z_he: redshift $z_{he}$ in our local frame
 * @z_cmb: redshift $z_{CMB}$ in the CMB frame
 *
 * Calculate the distance modulus [Eq. $\eqref{eq:def:dmu}$] using the frame corrected luminosity
 * distance [nc_distance_luminosity_hef()].
 *
 * Returns: $\delta\mu(z_{hef},z_{CMB})$.
 */
gdouble
nc_distance_dmodulus_hef (NcDistance *dist, NcHICosmo *cosmo, const gdouble z_he, const gdouble z_cmb)
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
 * We define the angular diameter curvature scale $D_a(z_\star)$ as
 * $$D_a(z_\star) = \frac{E(z_\star)}{1 + z_\star} D_t(z_\star),$$
 * where $z_\star$ is the decoupling redshift, given by [nc_distance_decoupling_redshift()],
 * $E(z_\star)$ is the normalized Hubble function [Eq. $\eqref{eq:def:Ez}$] and
 * $D_t(z_\star)$ is the transverse comoving distance [Eq. $\eqref{eq:def:Dt}$] both computed at $z_\star$.
 *
 * Returns: $D_a(z_\star)$.
 */
gdouble
nc_distance_angular_diameter_curvature_scale (NcDistance *dist, NcHICosmo *cosmo)
{
  const gdouble z_star = nc_distance_decoupling_redshift (dist, cosmo);

  if (gsl_finite (z_star))
    return sqrt (nc_hicosmo_E2 (cosmo, z_star)) *
           nc_distance_transverse (dist, cosmo, z_star) / (1.0 + z_star);
  else
    return GSL_NAN;
}

/**
 * nc_distance_shift_parameter:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * The shift parameter $R(z)$ is defined as
 * \begin{equation*}
 * R(z) = \frac{\sqrt{\Omega_{m0} H_0^2}}{c} (1 + z) D_A(z) = \sqrt{\Omega_{m0}} D_t(z),
 * \end{equation*}
 * where $\Omega_{m0}$ is the matter density paremeter [nc_hicosmo_Omega_m0()],
 * $D_A(z) = D_{H_0} D_t(z) / (1 + z)$ is the angular diameter distance and
 * $D_t(z)$ is the tranverse comoving distance [Eq. $\eqref{eq:def:Dt}$].
 *
 * Returns: $R(z)$.
 */
gdouble
nc_distance_shift_parameter (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble sqrt_mod_Omega_m0 = sqrt (fabs (nc_hicosmo_Omega_m0 (cosmo)));
  const gdouble transverse        = nc_distance_transverse (dist, cosmo, z);

  return sqrt_mod_Omega_m0 * transverse;
}

/**
 * nc_distance_shift_parameter_lss:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Compute the shift parameter $R(z)$ [nc_distance_shift_parameter()] at the
 * decoupling redshift $z_\star$ [nc_distance_decoupling_redshift()].
 *
 * Returns: $R(z_\star)$.
 */
gdouble
nc_distance_shift_parameter_lss (NcDistance *dist, NcHICosmo *cosmo)
{
  const gdouble sqrt_mod_Omega_m0 = sqrt (fabs (nc_hicosmo_Omega_m0 (cosmo)));
  const gdouble z_star            = nc_distance_decoupling_redshift (dist, cosmo);

  if (gsl_finite (z_star))
  {
    const gdouble transverse = nc_distance_transverse (dist, cosmo, z_star);

    return sqrt_mod_Omega_m0 * transverse;
  }
  else
  {
    return GSL_NAN;
  }
}

/**
 * nc_distance_comoving_lss:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Compute the comoving distance $D_c(z)$ [Eq. \eqref{eq:def:Dc}] at the
 * decoupling redshift $z_\star$ [nc_distance_decoupling_redshift()].
 *
 * Returns: $D_c(z_\star)$.
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
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * The decoupling redshift $z_\star$ corresponds to the epoch of
 * the last scattering surface of the cosmic microwave background photons.
 *
 * This function computes $z_\star$ using [nc_hicosmo_z_lss()], if @cosmo implements
 * it, or using Hu & Sugiyama fitting formula [Hu (1996)][XHu1996][[arXiv](https://arxiv.org/abs/astro-ph/9510117)],
 * $$ z_\star = 1048 \left(1 + 1.24 \times 10^{-3}  (\Omega_{b0} h^2)^{-0.738}\right) \left(1 + g_1 (\Omega_{m0} h^2)^{g_2}\right),$$
 * where $\Omega_{b0} h^2$ [nc_hicosmo_Omega_b0h2()] and $\Omega_{m0} h^2$ [nc_hicosmo_Omega_m0h2()]
 * are, respectively, the baryonic and matter density parameters times the square
 * of the dimensionless Hubble parameter $h$,  $H_0 = 100 \, h \, \text{km/s} \, \text{Mpc}^{-1}$.
 * The parameters $g_1$ and $g_2$ are given by
 * $$g_1 = \frac{0.0783 (\Omega_{b0} h^2)^{-0.238}}{(1 + 39.5 (\Omega_{b0} h^2)^{0.763})}
 * \; \text{and} \; g_2 = \frac{0.56}{\left(1 + 21.1 (\Omega_{b0} h^2)^{1.81}\right)}.$$
 *
 * Returns: $z_\star$.
 */
gdouble
nc_distance_decoupling_redshift (NcDistance *dist, NcHICosmo *cosmo)
{
  NCM_UNUSED (dist);

  if (ncm_model_check_impl_opt (NCM_MODEL (cosmo), NC_HICOSMO_IMPL_z_lss))
  {
    return nc_hicosmo_z_lss (cosmo);
  }
  else
  {
    gdouble omega_b_h2 = nc_hicosmo_Omega_b0h2 (cosmo);
    gdouble omega_m_h2 = nc_hicosmo_Omega_m0h2 (cosmo);
    gdouble g1         = 0.0783 * pow (omega_b_h2, -0.238) / (1.0 + 39.5 * pow (omega_b_h2, 0.763));
    gdouble g2         = 0.560 / (1.0 + 21.1 * pow (omega_b_h2, 1.81));

    return 1048.0 * (1.0 + 1.24e-3 * pow (omega_b_h2, -0.738)) * (1.0 + g1 * pow (omega_m_h2, g2));
  }
}

static gdouble sound_horizon_integral_argument (gdouble z, gpointer p);

/**
 * nc_distance_sound_horizon:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the sound horizon $r_s$,
 * \begin{equation}
 * \theta_s (z) = \int_{z}^\infty \frac{c_s(z^\prime)}{E(z^\prime)} dz^\prime, \quad r_s (z) = \frac{\mathrm{sinn}\left(\sqrt{\Omega_{k0}}\theta_s\right)}{\sqrt{\Omega_{k0}}},
 * \end{equation}
 * where $c_s = \sqrt{c^{b\gamma 2}_s(z)}$ is the baryon-photon plasma speed of sound
 * (see nc_hicosmo_bgp_cs2() for more details)
 * and $E(z)$ is the normalized Hubble function [nc_hicosmo_E()].
 *
 * Returns: $r_s(z)$.
 */
gdouble
nc_distance_sound_horizon (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (cosmo);
  gdouble result, error;
  gsl_function F;

  g_assert (ncm_model_check_impl_opt (NCM_MODEL (cosmo), NC_HICOSMO_IMPL_bgp_cs2));

  F.function = &sound_horizon_integral_argument;
  F.params   = cosmo;

  if (dist->use_cache)
    ncm_integral_cached_x_inf (dist->sound_horizon_cache, &F, z, &result, &error);
  else
    ncm_integral_locked_a_inf (&F, z, NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR, &result, &error);

  return _nc_distance_sinn (result, Omega_k0);
}

static gdouble
sound_horizon_integral_argument (gdouble z, gpointer p)
{
  NcHICosmo *cosmo = NC_HICOSMO (p);

  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble bgp_cs2 = nc_hicosmo_bgp_cs2 (cosmo, z);

  /*printf ("% 22.15g % 22.15g % 22.15g\n", z, bgp_cs2, E2);*/
  return sqrt (bgp_cs2 / E2);
}

/**
 * nc_distance_dsound_horizon_dz:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Calculate the sound horizon [nc_distance_sound_horizon()] derivative with respect to $z$,
 * $$\frac{d r_s(z)}{dz} = - \frac{c_s(z)}{E(z)}\, ,$$
 * where $c_s = \sqrt{c^{b\gamma 2}_s(z)}$ is the baryon-photon plasma speed of sound
 * (see nc_hicosmo_bgp_cs2() for more details)
 * and $E(z)$ is the normalized Hubble function [nc_hicosmo_E()].
 *
 * Returns: $\frac{d r_s(z)}{dz}$.
 */
gdouble
nc_distance_dsound_horizon_dz (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  NCM_UNUSED (dist);

  return -sound_horizon_integral_argument (z, cosmo);
}

/**
 * nc_distance_acoustic_scale:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Compute the acoustic scale $l_A (z_\star)$ at $z_\star$ [nc_distance_decoupling_redshift()],
 * \begin{equation*}
 * l_A(z_\star) = \pi \frac{D_t (z_\star)}{r_s (z_\star)},
 * \end{equation*}
 * where $D_t(z_\star)$ is the comoving transverse distance [nc_distance_transverse()]
 * and $r_s(z_\star)$ is the sound horizon [nc_distance_sound_horizon()] both computed at $z_\star$.
 *
 * Returns: $l_A(z_\star)$.
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
 * nc_distance_theta100CMB:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Compute the $100\theta_\mathrm{CMB}$ angle at $z_\star$ [nc_distance_decoupling_redshift()],
 * \begin{equation*}
 * 100\theta_\mathrm{CMB} = 100 \times \frac{r_s (z_\star)}{D_t (z_\star)},
 * \end{equation*}
 * where $D_t(z_\star)$ is the comoving transverse distance [nc_distance_transverse()]
 * and $r_s(z_\star)$ is the sound horizon [nc_distance_sound_horizon()] both
 * both computed at $z_\star$.
 *
 * Returns: $100\theta_\mathrm{CMB}$.
 */
gdouble
nc_distance_theta100CMB (NcDistance *dist, NcHICosmo *cosmo)
{
  gdouble z = nc_distance_decoupling_redshift (dist, cosmo);

  if (gsl_finite (z))
    return 100.0 * nc_distance_sound_horizon (dist, cosmo, z) / nc_distance_transverse (dist, cosmo, z);
  else
    return GSL_NAN;
}

/**
 * nc_distance_drag_redshift:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * Drag redshift is the epoch at which baryons were released from photons.
 *
 * If the @dist object constains a NcRecomb object, it calculates the drag
 * redshift through the recombination history. Otherwise, it computes $z_d$
 * using the fitting formula given in [Eisenstein & Hu (1998)][XEisenstein1998] [[arXiv](https://arxiv.org/abs/astro-ph/9709112)],
 * $$z_d = \frac{1291 (\Omega_{m0} h^2)^{0.251}}{(1 + 0.659 (\Omega_{m0} h^2)^{0.828})}
 * \left(1 + b_1 (\Omega_{b0} h^2)^{b_2}\right),$$
 * where $\Omega_{b0} h^2$ [nc_hicosmo_Omega_b0h2()] and $\Omega_{m0} h^2$ [nc_hicosmo_Omega_m0h2()]
 * are, respectively, the baryonic and matter density parameters times the square
 * of the dimensionless Hubble parameter $h$,  $H_0 = 100 \, h \, \text{km/s} \, \text{Mpc}^{-1}$.
 * The parameters $b_1$ and $b_2$ are given by
 * $$b_1 = 0.313 (\Omega_{m0} h^2)^{-0.419} \left(1 + 0.607 (\Omega_{m0} h^2)^{0.674}\right) \;
 * \text{and} \; b_2 = 0.238 (\Omega_{m0} h^2)^{0.223}.$$
 *
 * Returns: $z_d$.
 */
gdouble
nc_distance_drag_redshift (NcDistance *dist, NcHICosmo *cosmo)
{
  if (dist->recomb != NULL)
  {
    return nc_recomb_get_tau_drag_z (dist->recomb, cosmo);
  }
  else
  {
    gdouble omega_m_h2 = nc_hicosmo_Omega_m0h2 (cosmo);
    gdouble omega_b_h2 = nc_hicosmo_Omega_b0h2 (cosmo);
    gdouble b1         = 0.313 * pow (omega_m_h2, -0.419) * (1.0 + 0.607 * pow (omega_m_h2, 0.674));
    gdouble b2         = 0.238 * pow (omega_m_h2, 0.223);

    NCM_UNUSED (dist);

    return 1291.0 * pow (omega_m_h2, 0.251) / (1.0 + 0.659 * pow (omega_m_h2, 0.828)) *
           (1.0 + b1 * pow (omega_b_h2, b2));
  }
}

/**
 * nc_distance_dilation_scale:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * The dilation scale is the spherically averaged distance for perturbations along and
 * orthogonal to the line of sight,
 * $$D_V(z) = \left[D_{H_0}^2 D_t(z)^2 \frac{cz}{H(z)} \right]^{1/3},$$
 * where $D_t(z)$ is the transverse comoving distance [Eq. $\eqref{eq:def:Dt}$], $c$ is the speed of light
 * [ncm_c_c()] and $H(z)$ is the Hubble function [nc_hicosmo_H()].
 * See [Eisenstein et al. (2005)][XEisenstein2005] [[arXiv](https://arxiv.org/abs/astro-ph/0501171)].
 *
 * This function computes the dimensionless dilation scale:
 * $$D_V^\star(z) = \left[D_t(z)^2 \frac{z}{E(z)} \right]^{1/3} = \frac{D_V(z)}{D_{H_0}},$$
 * where $E(z)$ is the normalized Hubble function [nc_hicosmo_E2()].
 *
 * Returns: $D_V^\star(z)$.
 */
gdouble
nc_distance_dilation_scale (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble Dt = nc_distance_transverse (dist, cosmo, z);
  gdouble E  = sqrt (nc_hicosmo_E2 (cosmo, z));
  gdouble Dv = cbrt (Dt * Dt * z / E);

  return Dv;
}

/**
 * nc_distance_bao_A_scale:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: the redshift $z$
 *
 * This function returns the [BAO](https://en.wikipedia.org/wiki/Baryon_acoustic_oscillations) $A$ parameter.
 * It is defined as
 * $$ A \equiv D_V (z) \frac{\sqrt{\Omega_{m0} H_0^2}}{z c},$$
 * where $\Omega_{m0}$ is the matter density parameter [nc_hicosmo_Omega_m0()], $c$ is the speed of light [ncm_c_c()],
 * $H_0$ is the Hubble parameter [nc_hicosmo_H0()] and $D_V(z)$ is the dilation scale [nc_distance_dilation_scale()].
 * See Section 4.5 from [Eisenstein et al. (2005)][XEisenstein2005] [[arXiv](https://arxiv.org/abs/astro-ph/0501171)].
 *
 * Returns: $A(z)$.
 */
gdouble
nc_distance_bao_A_scale (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble Dv            = nc_distance_dilation_scale (dist, cosmo, z);
  gdouble sqrt_Omega_m0 = sqrt (nc_hicosmo_Omega_m0 (cosmo));

  return sqrt_Omega_m0 * Dv / z;
}

/**
 * nc_distance_r_zd:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 * This function computes the sound horizon [nc_distance_sound_horizon ()]
 * at the drag redshift [nc_distance_drag_redshift()].
 *
 * Returns: $r(z_d)$.
 */
gdouble
nc_distance_r_zd (NcDistance *dist, NcHICosmo *cosmo)
{
  gdouble r_zd;

  if (ncm_model_check_impl_opt (NCM_MODEL (cosmo), NC_HICOSMO_IMPL_as_drag))
  {
    r_zd = nc_hicosmo_as_drag (cosmo);
  }
  else
  {
    gdouble zd = nc_distance_drag_redshift (dist, cosmo);

    if (!gsl_finite (zd))
      return GSL_NAN;

    r_zd = nc_distance_sound_horizon (dist, cosmo, zd);
  }

  return r_zd;
}

/**
 * nc_distance_r_zd_Mpc:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 *
 *  Similar as nc_distance_r_zd(), but now in units of $\mathrm{Mpc}$, i.e., $ R_H \times r(z_d)$.
 *
 * Returns: r_zd in $ \left[ \mathrm{Mpc} \right]$ units.
 */
gdouble
nc_distance_r_zd_Mpc (NcDistance *dist, NcHICosmo *cosmo)
{
  return nc_distance_r_zd (dist, cosmo) * nc_hicosmo_RH_Mpc (cosmo);
}

/**
 * nc_distance_bao_r_Dv:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: the redshift $z$
 *
 * This function computes $r(z_d) / D_V(z)$,
 * where $r(z_d)$ is given by nc_distance_r_zd() and $D_V(z)$ by nc_distance_dilation_scale().
 * For more information see [Percival et al. (2007)][XPercival2007] [[arXiv](https://arxiv.org/abs/0705.3323)].
 *
 * Returns: $r(z_d) / D_V(z)$.
 */
gdouble
nc_distance_bao_r_Dv (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble r_zd = nc_distance_r_zd (dist, cosmo);
  gdouble Dv   = nc_distance_dilation_scale (dist, cosmo, z);

  return r_zd / Dv;
}

/**
 * nc_distance_DH_r:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: the redshift $z$
 *
 * Computes the ratio between the Hubble (distance) radius and the sound horizon at the drag epoch,
 * $$\frac{R_H}{r_s(z_d)} = \frac{c}{H(z) r_d}.$$
 *
 * Returns: $R_H/r_d = c / (H(z) r_d)$.
 */
gdouble
nc_distance_DH_r (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble r_zd = nc_distance_r_zd (dist, cosmo);
  gdouble E    = nc_hicosmo_E (cosmo, z);

  return 1.0 / (E * r_zd);
}

/**
 * nc_distance_DA_r:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: the redshift $z$
 *
 * Computes the ratio between the angular-diameter distance and the sound horizon at the drag epoch,
 * $$\frac{D_A(z)}{c \, r_s(z_d)}.$$
 *
 * Returns: $D_A(z) / (c \, r_d)$.
 */
gdouble
nc_distance_DA_r (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble r_zd = nc_distance_r_zd (dist, cosmo);
  gdouble DA   = nc_distance_angular_diameter (dist, cosmo, z);

  return DA / r_zd;
}

/**
 * nc_distance_Dt_r:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: the redshift $z$
 *
 * Computes the ratio between the transverse distance and the sound horizon at the drag epoch,
 * $$\frac{D_t(z)}{c \, r_s(z_d)}.$$
 *
 * Returns: $D_t(z) / (c \, r_d)$.
 */
gdouble
nc_distance_Dt_r (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble r_zd = nc_distance_r_zd (dist, cosmo);
  gdouble Dt   = nc_distance_transverse (dist, cosmo, z);

  return Dt / r_zd;
}

/* Distances from z to infinity */

static gdouble
_nc_distance_comoving_infinity_integrand (gdouble logx, gpointer p)
{
  if (logx > GSL_LOG_DBL_MAX)
  {
    return 0.0;
  }
  else
  {
    NcHICosmo *cosmo = NC_HICOSMO (p);
    const gdouble z  = expm1 (logx);
    const gdouble E  = sqrt (nc_hicosmo_E2 (cosmo, z));

    if (gsl_finite (E))
      return (1.0 + z) / E;
    else
      return 0.0;
  }
}

/**
 * nc_distance_comoving_z_to_infinity:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * Computes the comoving distance from @z to infinity.
 *
 * Returns: $D_c$ from $z$ to $\infty$.
 */
gdouble
nc_distance_comoving_z_to_infinity (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &_nc_distance_comoving_infinity_integrand;
  F.params   = cosmo;

  if (dist->use_cache)
    ncm_integral_cached_x_inf (dist->comoving_infinity, &F, log1p (z), &result, &error);
  else
    ncm_integral_locked_a_inf (&F, log1p (z), NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR, &result, &error);

  return result;
}

/**
 * nc_distance_transverse_z_to_infinity:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Compute the transverse comoving distance $D_t (z)$ (defined in Eq. $\eqref{eq:def:Dt}$)
 * but from @z to infinity.
 *
 * Returns: $D_t(z)$ from $z$ to $\infty$.
 */
gdouble
nc_distance_transverse_z_to_infinity (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble comoving_dist = nc_distance_comoving_z_to_infinity (dist, cosmo, z);

  if (gsl_isinf (comoving_dist))
    return comoving_dist;

  return _nc_distance_sinn (comoving_dist, Omega_k0);
}

/**
 * nc_distance_inv_comoving:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @xi: the comoving distance $\xi$
 *
 * Computes the inverse of $\xi(z)$.
 *
 * Returns: $z(\xi)$
 */
gdouble
nc_distance_inv_comoving (NcDistance *dist, NcHICosmo *cosmo, gdouble xi)
{
  g_assert (dist->cpu_inv_comoving);

  if (xi > dist->max_comoving)
    g_error ("nc_distance_inv_comoving: maximum comoving distance exceeded. "
             "Increase the maximum redshift in the cosmology object.");

  return ncm_spline_eval (dist->inv_comoving_dist, xi);
}

/***************************************************************************
 *            cosmic_time.c
 *
 *  Wed Nov 12 17:06:27 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

static gdouble
_nc_distance_cosmic_time_integrand (gdouble logx, gpointer p)
{
  if (logx > GSL_LOG_DBL_MAX)
  {
    return 0.0;
  }
  else
  {
    NcHICosmo *cosmo = NC_HICOSMO (p);
    const gdouble z  = expm1 (logx);
    const gdouble E  = sqrt (nc_hicosmo_E2 (cosmo, z));

    if (gsl_finite (E))
      return 1.0 / E;
    else
      return 0.0;
  }
}

/**
 * nc_distance_cosmic_time:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the cosmological time, $t(z)$, defined as
 * \begin{equation*}
 * t(z) = \int_z^{\infty} \frac{dx}{(1+x)E(x)}.
 * \end{equation*}
 *
 * Returns: $t(z)$.
 */
gdouble
nc_distance_cosmic_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &_nc_distance_cosmic_time_integrand;
  F.params   = cosmo;

  if (dist->use_cache)
    ncm_integral_cached_x_inf (dist->time_cache, &F, log1p (z), &result, &error);
  else
    ncm_integral_locked_a_inf (&F, log1p (z), NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR, &result, &error);

  return result;
}

/**
 * nc_distance_lookback_time:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This functions computes the look-back time, $t_{lb}(z)$, defined as
 *
 * \begin{equation*}
 * t_{lb}(z) = \int_0^{z} \frac{dx}{(1+x)E(x)}.
 * \end{equation*}
 *
 * Returns: $t_{lb}(z)$.
 */
gdouble
nc_distance_lookback_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &_nc_distance_cosmic_time_integrand;
  F.params   = cosmo;

  if (dist->use_cache)
    ncm_integral_cached_0_x (dist->lookback_time_cache, &F, log1p (z), &result, &error);
  else
    ncm_integral_locked_a_b (&F, 0.0, log1p (z), 0.0, NCM_INTEGRAL_ERROR, &result, &error);

  return result;
}

/**
 * nc_distance_conformal_lookback_time:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the conformal look-back time $\eta_{lb}(z)$.
 * Within the chosen units it becomes the same as the comoving distance [nc_distance_comoving()],
 * given in Eq. $\eqref{eq:def:dc}$.
 *
 * Returns: $\eta_{lb}(z)$.
 */
gdouble
nc_distance_conformal_lookback_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  return nc_distance_comoving (dist, cosmo, z);
}

static gdouble
_nc_distance_conformal_time_integrand (gdouble logx, gpointer p)
{
  if (logx > GSL_LOG_DBL_MAX)
  {
    return 0.0;
  }
  else
  {
    NcHICosmo *cosmo = NC_HICOSMO (p);
    const gdouble z  = expm1 (logx);
    const gdouble x  = 1.0 + z;
    const gdouble E  = sqrt (nc_hicosmo_E2 (cosmo, z));

    if (gsl_finite (E))
      return x / E;
    else
      return 0.0;
  }
}

/**
 * nc_distance_conformal_time:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the cosmological conformal time $\eta(z)$ defined as
 * $$ \eta(z) = \int_z^\infty \frac{dz}{E(z)} \, . $$
 *
 * Returns: $\eta(z)$.
 */
gdouble
nc_distance_conformal_time (NcDistance *dist, NcHICosmo *cosmo, const gdouble z)
{
  gdouble result, error;
  gsl_function F;

  F.function = &_nc_distance_conformal_time_integrand;
  F.params   = cosmo;

  if (dist->use_cache)
    ncm_integral_cached_x_inf (dist->conformal_time_cache, &F, log1p (z), &result, &error);
  else
    ncm_integral_locked_a_inf (&F, log1p (z), NCM_INTEGRAL_ABS_ERROR, NCM_INTEGRAL_ERROR, &result, &error);

  return result;
}

/***************************************************************************
 * Vectorized 'distances'
 ****************************************************************************/

typedef gdouble (*_z_getter) (gconstpointer p, guint i);
typedef void (*_dist_setter) (gpointer p, guint i, gdouble dist);

static gdouble
_z_getter_garray (gconstpointer p, guint i)
{
  return g_array_index (((const GArray *) p), gdouble, i);
}

static void
_dist_setter_garray (gpointer p, guint i, gdouble dist)
{
  g_array_index (((GArray *) p), gdouble, i) = dist;
}

static gdouble
_z_getter_vector_no_stride (gconstpointer p, guint i)
{
  return ncm_vector_fast_get (p, i);
}

static void
_dist_setter_vector_no_stride (gpointer p, guint i, gdouble dist)
{
  ncm_vector_fast_set (NCM_VECTOR (p), i, dist);
}

static gdouble
_z_getter_vector (gconstpointer p, guint i)
{
  return ncm_vector_get (p, i);
}

static void
_dist_setter_vector (gpointer p, guint i, gdouble dist)
{
  ncm_vector_set (NCM_VECTOR (p), i, dist);
}

static void
_nc_distance_comoving_worker (NcDistance *dist, NcHICosmo *cosmo, _z_getter get_z, _dist_setter set_dist, const guint len, gconstpointer z, gpointer res)
{
  guint i;

  switch (dist->cmethod)
  {
    case NC_DISTANCE_COMOVING_METHOD_FROM_MODEL:

      for (i = 0; i < len; i++)
      {
        const gdouble zi = get_z (z, i);
        const gdouble Di = nc_hicosmo_Dc (cosmo, zi);

        set_dist (res, i, Di);
      }

      break;
    case NC_DISTANCE_COMOVING_METHOD_INT_E:
    {
      NcmSpline *s = ncm_ode_spline_peek_spline (dist->comoving_distance_spline);
      gdouble result, error;
      gsl_function F;

      F.function = &_comoving_distance_integral_argument;
      F.params   = cosmo;

      for (i = 0; i < len; i++)
      {
        const gdouble zi = get_z (z, i);

        if (zi <= dist->zf)
          result = ncm_spline_eval (s, zi);
        else if (dist->use_cache)
          ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
        else
          ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);

        set_dist (res, i, result);
      }

      break;
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

static void
_nc_distance_transverse_worker (NcDistance *dist, NcHICosmo *cosmo, _z_getter get_z, _dist_setter set_dist, const guint len, gconstpointer z, gpointer res)
{
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (dist->cmethod)
  {
    case NC_DISTANCE_COMOVING_METHOD_FROM_MODEL:

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);

            set_dist (res, i, Di);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);
            const gdouble Dt = sinh (sqrt_Omega_k0 * Di) / sqrt_Omega_k0;

            set_dist (res, i, Dt);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);
            const gdouble Dt = sin (sqrt_Omega_k0 * Di) / sqrt_Omega_k0;

            set_dist (res, i, Dt);
          }

          break;
      }

      break;

    case NC_DISTANCE_COMOVING_METHOD_INT_E:
    {
      NcmSpline *s = ncm_ode_spline_peek_spline (dist->comoving_distance_spline);
      gdouble result, error;
      gsl_function F;

      F.function = &_comoving_distance_integral_argument;
      F.params   = cosmo;

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            set_dist (res, i, result);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dt;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dt = sinh (sqrt_Omega_k0 * result) / sqrt_Omega_k0;

            set_dist (res, i, Dt);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dt;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dt = sin (sqrt_Omega_k0 * result) / sqrt_Omega_k0;

            set_dist (res, i, Dt);
          }

          break;
      }

      break;
    }

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

static void
_nc_distance_luminosity_worker (NcDistance *dist, NcHICosmo *cosmo, _z_getter get_z, _dist_setter set_dist, const guint len, gconstpointer z, gpointer res)
{
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (dist->cmethod)
  {
    case NC_DISTANCE_COMOVING_METHOD_FROM_MODEL:

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Dl = (1.0 + zi) * nc_hicosmo_Dc (cosmo, zi);

            set_dist (res, i, Dl);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);
            const gdouble Dl = (1.0 + zi) * sinh (sqrt_Omega_k0 * Di) / sqrt_Omega_k0;

            set_dist (res, i, Dl);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);
            const gdouble Dl = (1.0 + zi) * sin (sqrt_Omega_k0 * Di) / sqrt_Omega_k0;

            set_dist (res, i, Dl);
          }

          break;
      }

      break;

    case NC_DISTANCE_COMOVING_METHOD_INT_E:
    {
      NcmSpline *s = ncm_ode_spline_peek_spline (dist->comoving_distance_spline);
      gdouble result, error;
      gsl_function F;

      F.function = &_comoving_distance_integral_argument;
      F.params   = cosmo;

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dl;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dl = (1.0 + zi) * result;

            set_dist (res, i, Dl);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dl;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dl = (1.0 + zi) * sinh (sqrt_Omega_k0 * result) / sqrt_Omega_k0;

            set_dist (res, i, Dl);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dl;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dl = (1.0 + zi) * sin (sqrt_Omega_k0 * result) / sqrt_Omega_k0;

            set_dist (res, i, Dl);
          }

          break;
      }

      break;
    }

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

static void
_nc_distance_angular_diameter_worker (NcDistance *dist, NcHICosmo *cosmo, _z_getter get_z, _dist_setter set_dist, guint len, gconstpointer z, gpointer res)
{
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (dist->cmethod)
  {
    case NC_DISTANCE_COMOVING_METHOD_FROM_MODEL:

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble DA = nc_hicosmo_Dc (cosmo, zi) / (1.0 + zi);

            set_dist (res, i, DA);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);
            const gdouble DA = sinh (sqrt_Omega_k0 * Di) / sqrt_Omega_k0 / (1.0 + zi);

            set_dist (res, i, DA);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            const gdouble Di = nc_hicosmo_Dc (cosmo, zi);
            const gdouble DA = sin (sqrt_Omega_k0 * Di) / sqrt_Omega_k0 / (1.0 + zi);

            set_dist (res, i, DA);
          }

          break;
      }

      break;

    case NC_DISTANCE_COMOVING_METHOD_INT_E:
    {
      NcmSpline *s = ncm_ode_spline_peek_spline (dist->comoving_distance_spline);
      gdouble result, error;
      gsl_function F;

      F.function = &_comoving_distance_integral_argument;
      F.params   = cosmo;

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble DA;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            DA = result / (1.0 + zi);

            set_dist (res, i, DA);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble DA;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            DA = sinh (sqrt_Omega_k0 * result) / sqrt_Omega_k0 / (1.0 + zi);

            set_dist (res, i, DA);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble DA;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            DA = sin (sqrt_Omega_k0 * result) / sqrt_Omega_k0 / (1.0 + zi);

            set_dist (res, i, DA);
          }

          break;
      }

      break;
    }

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

static void
_nc_distance_dmodulus_worker (NcDistance *dist, NcHICosmo *cosmo, _z_getter get_z, _dist_setter set_dist, const guint len, gconstpointer z, gpointer res)
{
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (dist->cmethod)
  {
    case NC_DISTANCE_COMOVING_METHOD_FROM_MODEL:

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi  = get_z (z, i);
            const gdouble Dl  = (1.0 + zi) * nc_hicosmo_Dc (cosmo, zi);
            const gdouble dmu = (5.0 * log10 (Dl) + 25.0);

            set_dist (res, i, dmu);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi  = get_z (z, i);
            const gdouble Di  = nc_hicosmo_Dc (cosmo, zi);
            const gdouble Dl  = (1.0 + zi) * sinh (sqrt_Omega_k0 * Di) / sqrt_Omega_k0;
            const gdouble dmu = (5.0 * log10 (Dl) + 25.0);

            set_dist (res, i, dmu);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi  = get_z (z, i);
            const gdouble Di  = nc_hicosmo_Dc (cosmo, zi);
            const gdouble Dl  = (1.0 + zi) * sin (sqrt_Omega_k0 * Di) / sqrt_Omega_k0;
            const gdouble dmu = (5.0 * log10 (Dl) + 25.0);

            set_dist (res, i, dmu);
          }

          break;
      }

      break;

    case NC_DISTANCE_COMOVING_METHOD_INT_E:
    {
      NcmSpline *s = ncm_ode_spline_peek_spline (dist->comoving_distance_spline);
      gdouble result, error;
      gsl_function F;

      F.function = &_comoving_distance_integral_argument;
      F.params   = cosmo;

      switch (k)
      {
        case 0:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dl, dmu;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dl  = (1.0 + zi) * result;
            dmu = (5.0 * log10 (Dl) + 25.0);

            set_dist (res, i, dmu);
          }

          break;

        case -1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dl, dmu;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dl  = (1.0 + zi) * sinh (sqrt_Omega_k0 * result) / sqrt_Omega_k0;
            dmu = (5.0 * log10 (Dl) + 25.0);

            set_dist (res, i, dmu);
          }

          break;

        case 1:

          for (i = 0; i < len; i++)
          {
            const gdouble zi = get_z (z, i);
            gdouble Dl, dmu;

            if (zi <= dist->zf)
            {
              result = ncm_spline_eval (s, zi);
            }
            else
            {
              if (dist->use_cache)
                ncm_integral_cached_0_x (dist->comoving_distance_cache, &F, zi, &result, &error);
              else
                ncm_integral_locked_a_b (&F, 0.0, zi, 0.0, NCM_INTEGRAL_ERROR, &result, &error);
            }

            Dl  = (1.0 + zi) * sin (sqrt_Omega_k0 * result) / sqrt_Omega_k0;
            dmu = (5.0 * log10 (Dl) + 25.0);

            set_dist (res, i, dmu);
          }

          break;
      }

      break;
    }

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      break; /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_distance_comoving_array:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: (element-type gdouble): a vector of redshifts
 *
 * Compute the comoving distance $D_c(z)$ for each redshift in @z.
 *
 * Returns: (element-type gdouble) (transfer full): a vector with the comoving distances.
 */
GArray *
nc_distance_comoving_array (NcDistance *dist, NcHICosmo *cosmo, const GArray *z)
{
  GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z->len);

  g_array_set_size (res, z->len);
  _nc_distance_comoving_worker (dist, cosmo, _z_getter_garray, _dist_setter_garray, z->len, z, res);

  return res;
}

/**
 * nc_distance_transverse_array:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: (element-type gdouble): a vector of redshifts
 *
 * Compute the transverse comoving distance $D_t(z)$ for each redshift in @z.
 *
 * Returns: (element-type gdouble) (transfer full): a vector with the transverse distances.
 */
GArray *
nc_distance_transverse_array (NcDistance *dist, NcHICosmo *cosmo, const GArray *z)
{
  GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z->len);

  g_array_set_size (res, z->len);

  _nc_distance_transverse_worker (dist, cosmo, _z_getter_garray, _dist_setter_garray, z->len, z, res);

  return res;
}

/**
 * nc_distance_luminosity_array:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: (element-type gdouble): a vector of redshifts
 *
 * Compute the luminosity distance $D_L(z)$ for each redshift in @z.
 *
 * Returns: (element-type gdouble) (transfer full): a vector with the luminosity distances.
 */
GArray *
nc_distance_luminosity_array (NcDistance *dist, NcHICosmo *cosmo, const GArray *z)
{
  GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z->len);

  g_array_set_size (res, z->len);
  _nc_distance_luminosity_worker (dist, cosmo, _z_getter_garray, _dist_setter_garray, z->len, z, res);

  return res;
}

/**
 * nc_distance_angular_diameter_array:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: (element-type gdouble): a vector of redshifts
 *
 * Compute the angular diameter distance $D_A(z)$ for each redshift in @z.
 *
 * Returns: (element-type gdouble) (transfer full): a vector with the angular diameter distances.
 */
GArray *
nc_distance_angular_diameter_array (NcDistance *dist, NcHICosmo *cosmo, const GArray *z)
{
  GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z->len);

  g_array_set_size (res, z->len);
  _nc_distance_angular_diameter_worker (dist, cosmo, _z_getter_garray, _dist_setter_garray, z->len, z, res);

  return res;
}

/**
 * nc_distance_dmodulus_array:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: (element-type gdouble): a vector of redshifts
 *
 * Compute the distance modulus $\delta\mu(z)$ for each redshift in @z.
 *
 * Returns: (element-type gdouble) (transfer full): a vector with the distance modulus.
 */
GArray *
nc_distance_dmodulus_array (NcDistance *dist, NcHICosmo *cosmo, const GArray *z)
{
  GArray *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z->len);

  g_array_set_size (res, z->len);
  _nc_distance_dmodulus_worker (dist, cosmo, _z_getter_garray, _dist_setter_garray, z->len, z, res);

  return res;
}

/**
 * nc_distance_comoving_vector:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: a #NcmVector of redshifts
 * @Dc: a #NcmVector to store the comoving distances
 *
 * Compute the comoving distance $D_c(z)$ for each redshift in @z.
 *
 */
void
nc_distance_comoving_vector (NcDistance *dist, NcHICosmo *cosmo, NcmVector *z, NcmVector *Dc)
{
  const guint len       = ncm_vector_len (z);
  const guint z_stride  = ncm_vector_stride (z);
  const guint Dc_stride = ncm_vector_stride (Dc);

  g_assert_cmpuint (len, ==, ncm_vector_len (Dc));

  if (G_LIKELY ((z_stride == 1) && (Dc_stride == 1)))
    _nc_distance_comoving_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector_no_stride, len, z, Dc);
  else if (z_stride == 1)
    _nc_distance_comoving_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector, len, z, Dc);
  else if (Dc_stride == 1)
    _nc_distance_comoving_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector_no_stride, len, z, Dc);
  else
    _nc_distance_comoving_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector, len, z, Dc);
}

/**
 * nc_distance_transverse_vector:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: a #NcmVector of redshifts
 * @Dt: a #NcmVector to store the transverse distances
 *
 * Compute the transverse distance $D_t(z)$ for each redshift in @z.
 *
 */
void
nc_distance_transverse_vector (NcDistance *dist, NcHICosmo *cosmo, NcmVector *z, NcmVector *Dt)
{
  const guint len       = ncm_vector_len (z);
  const guint z_stride  = ncm_vector_stride (z);
  const guint Dt_stride = ncm_vector_stride (Dt);

  g_assert_cmpuint (len, ==, ncm_vector_len (Dt));

  if (G_LIKELY ((z_stride == 1) && (Dt_stride == 1)))
    _nc_distance_transverse_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector_no_stride, len, z, Dt);
  else if (z_stride == 1)
    _nc_distance_transverse_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector, len, z, Dt);
  else if (Dt_stride == 1)
    _nc_distance_transverse_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector_no_stride, len, z, Dt);
  else
    _nc_distance_transverse_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector, len, z, Dt);
}

/**
 * nc_distance_luminosity_vector:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: a #NcmVector of redshifts
 * @Dl: a #NcmVector to store the luminosity distances
 *
 * Compute the luminosity distance $D_t(z)$ for each redshift in @z.
 *
 */
void
nc_distance_luminosity_vector (NcDistance *dist, NcHICosmo *cosmo, NcmVector *z, NcmVector *Dl)
{
  const guint len       = ncm_vector_len (z);
  const guint z_stride  = ncm_vector_stride (z);
  const guint Dl_stride = ncm_vector_stride (Dl);

  g_assert_cmpuint (len, ==, ncm_vector_len (Dl));

  if (G_LIKELY ((z_stride == 1) && (Dl_stride == 1)))
    _nc_distance_luminosity_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector_no_stride, len, z, Dl);
  else if (z_stride == 1)
    _nc_distance_luminosity_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector, len, z, Dl);
  else if (Dl_stride == 1)
    _nc_distance_luminosity_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector_no_stride, len, z, Dl);
  else
    _nc_distance_luminosity_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector, len, z, Dl);
}

/**
 * nc_distance_angular_diameter_vector:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: a #NcmVector of redshifts
 * @DA: a #NcmVector to store the angular diameter distances
 *
 * Compute the angular diameter distance $D_A(z)$ for each redshift in @z.
 *
 */
void
nc_distance_angular_diameter_vector (NcDistance *dist, NcHICosmo *cosmo, NcmVector *z, NcmVector *DA)
{
  const guint len       = ncm_vector_len (z);
  const guint z_stride  = ncm_vector_stride (z);
  const guint DA_stride = ncm_vector_stride (DA);

  g_assert_cmpuint (len, ==, ncm_vector_len (DA));

  if (G_LIKELY ((z_stride == 1) && (DA_stride == 1)))
    _nc_distance_angular_diameter_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector_no_stride, len, z, DA);
  else if (z_stride == 1)
    _nc_distance_angular_diameter_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector, len, z, DA);
  else if (DA_stride == 1)
    _nc_distance_angular_diameter_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector_no_stride, len, z, DA);
  else
    _nc_distance_angular_diameter_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector, len, z, DA);
}

/**
 * nc_distance_dmodulus_vector:
 * @dist: a #NcDistance
 * @cosmo: a #NcHICosmo
 * @z: a #NcmVector of redshifts
 * @dmod: a #NcmVector to store the distance modulus distances
 *
 * Compute the angular diameter distance $D_A(z)$ for each redshift in @z.
 *
 */
void
nc_distance_dmodulus_vector (NcDistance *dist, NcHICosmo *cosmo, NcmVector *z, NcmVector *dmod)
{
  const guint len         = ncm_vector_len (z);
  const guint z_stride    = ncm_vector_stride (z);
  const guint dmod_stride = ncm_vector_stride (dmod);

  g_assert_cmpuint (len, ==, ncm_vector_len (dmod));

  if (G_LIKELY ((z_stride == 1) && (dmod_stride == 1)))
    _nc_distance_dmodulus_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector_no_stride, len, z, dmod);
  else if (z_stride == 1)
    _nc_distance_dmodulus_worker (dist, cosmo, _z_getter_vector_no_stride, _dist_setter_vector, len, z, dmod);
  else if (dmod_stride == 1)
    _nc_distance_dmodulus_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector_no_stride, len, z, dmod);
  else
    _nc_distance_dmodulus_worker (dist, cosmo, _z_getter_vector, _dist_setter_vector, len, z, dmod);
}

#define _NC_DISTANCE_FUNC0_TO_FLIST(fname)                                                                                   \
        static void _nc_distance_flist_ ## fname (NcmMSetFuncList * flist, NcmMSet * mset, const gdouble * x, gdouble * res) \
        {                                                                                                                    \
          NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));                                            \
          NcDistance *dist = NC_DISTANCE (ncm_mset_func_list_peek_obj (flist));                                              \
          nc_distance_prepare_if_needed (dist, cosmo);                                                                       \
          res[0] = nc_distance_ ## fname (dist, cosmo);                                                                      \
        }

_NC_DISTANCE_FUNC0_TO_FLIST (decoupling_redshift)
_NC_DISTANCE_FUNC0_TO_FLIST (drag_redshift)
_NC_DISTANCE_FUNC0_TO_FLIST (shift_parameter_lss)
_NC_DISTANCE_FUNC0_TO_FLIST (comoving_lss)
_NC_DISTANCE_FUNC0_TO_FLIST (acoustic_scale)
_NC_DISTANCE_FUNC0_TO_FLIST (theta100CMB)
_NC_DISTANCE_FUNC0_TO_FLIST (angular_diameter_curvature_scale)
_NC_DISTANCE_FUNC0_TO_FLIST (r_zd)
_NC_DISTANCE_FUNC0_TO_FLIST (r_zd_Mpc)

#define _NC_DISTANCE_FUNC1_TO_FLIST(fname)                                                                                   \
        static void _nc_distance_flist_ ## fname (NcmMSetFuncList * flist, NcmMSet * mset, const gdouble * x, gdouble * res) \
        {                                                                                                                    \
          NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));                                            \
          NcDistance *dist = NC_DISTANCE (ncm_mset_func_list_peek_obj (flist));                                              \
          nc_distance_prepare_if_needed (dist, cosmo);                                                                       \
          res[0] = nc_distance_ ## fname (dist, cosmo, x[0]);                                                                \
        }

_NC_DISTANCE_FUNC1_TO_FLIST (comoving)
_NC_DISTANCE_FUNC1_TO_FLIST (transverse)
_NC_DISTANCE_FUNC1_TO_FLIST (luminosity)
_NC_DISTANCE_FUNC1_TO_FLIST (angular_diameter)
_NC_DISTANCE_FUNC1_TO_FLIST (dmodulus)
_NC_DISTANCE_FUNC1_TO_FLIST (dilation_scale)
_NC_DISTANCE_FUNC1_TO_FLIST (bao_A_scale)
_NC_DISTANCE_FUNC1_TO_FLIST (bao_r_Dv)
_NC_DISTANCE_FUNC1_TO_FLIST (DH_r)
_NC_DISTANCE_FUNC1_TO_FLIST (DA_r)
_NC_DISTANCE_FUNC1_TO_FLIST (sound_horizon)

void
_nc_distance_register_functions (void)
{
  ncm_mset_func_list_register ("decoupling_redshift",              "z_\\mathrm{dec}",          "NcDistance", "Decoupling redshift",              NC_TYPE_DISTANCE, _nc_distance_flist_decoupling_redshift,              0, 1);

  ncm_mset_func_list_register ("drag_redshift",                    "z_\\mathrm{drag}",         "NcDistance", "Drag redshift",                    NC_TYPE_DISTANCE, _nc_distance_flist_drag_redshift,                    0, 1);
  ncm_mset_func_list_register ("shift_parameter_lss",              "R_\\mathrm{lss}",          "NcDistance", "Shift parameter at lss",           NC_TYPE_DISTANCE, _nc_distance_flist_shift_parameter_lss,              0, 1);
  ncm_mset_func_list_register ("comoving_lss",                     "d_\\mathrm{lss}",          "NcDistance", "Comoving scale of lss",            NC_TYPE_DISTANCE, _nc_distance_flist_comoving_lss,                     0, 1);
  ncm_mset_func_list_register ("acoustic_scale",                   "r_\\mathrm{ac}",           "NcDistance", "Acoustic scale at lss",            NC_TYPE_DISTANCE, _nc_distance_flist_acoustic_scale,                   0, 1);
  ncm_mset_func_list_register ("theta100CMB",                      "100\\theta_\\mathrm{CMB}", "NcDistance", "CMB angular scale times 100",      NC_TYPE_DISTANCE, _nc_distance_flist_theta100CMB,                      0, 1);
  ncm_mset_func_list_register ("angular_diameter_curvature_scale", "z_\\mathrm{dec}",          "NcDistance", "Angular diameter curvature scale", NC_TYPE_DISTANCE, _nc_distance_flist_angular_diameter_curvature_scale, 0, 1);
  ncm_mset_func_list_register ("r_zd",                             "r_\\mathrm{dec}",          "NcDistance", "Sound horizon at drag redshift",   NC_TYPE_DISTANCE, _nc_distance_flist_r_zd,                             0, 1);
  ncm_mset_func_list_register ("r_zd_Mpc",                         "r_\\mathrm{dec}R_H",       "NcDistance", "Sound horizon at drag redshift in Mpc", NC_TYPE_DISTANCE, _nc_distance_flist_r_zd_Mpc,                    0, 1);

  ncm_mset_func_list_register ("comoving",         "d_\\mathrm{c}",               "NcDistance", "Comoving distance",          NC_TYPE_DISTANCE, _nc_distance_flist_comoving,         1, 1);
  ncm_mset_func_list_register ("transverse",       "d_\\mathrm{t}",               "NcDistance", "Transverse distance",        NC_TYPE_DISTANCE, _nc_distance_flist_transverse,       1, 1);
  ncm_mset_func_list_register ("luminosity",       "d_\\mathrm{l}",               "NcDistance", "Luminosity distance",        NC_TYPE_DISTANCE, _nc_distance_flist_luminosity,       1, 1);
  ncm_mset_func_list_register ("angular_diameter", "d_\\mathrm{A}",               "NcDistance", "Angular diameter distance",  NC_TYPE_DISTANCE, _nc_distance_flist_angular_diameter, 1, 1);
  ncm_mset_func_list_register ("dmodulus",         "\\delta\\mu",                 "NcDistance", "delta-Distance modulus",     NC_TYPE_DISTANCE, _nc_distance_flist_dmodulus,         1, 1);
  ncm_mset_func_list_register ("dilation_scale",   "D_\\mathrm{A}",               "NcDistance", "Dilation scale",             NC_TYPE_DISTANCE, _nc_distance_flist_dilation_scale,   1, 1);
  ncm_mset_func_list_register ("bao_A_scale",      "D_\\mathrm{BAO}",             "NcDistance", "BAO A scale",                NC_TYPE_DISTANCE, _nc_distance_flist_bao_A_scale,      1, 1);
  ncm_mset_func_list_register ("bao_r_Dv",         "r_\\mathrm{c}/D_\\mathrm{V}", "NcDistance", "BAO r_Dv",                   NC_TYPE_DISTANCE, _nc_distance_flist_bao_r_Dv,         1, 1);
  ncm_mset_func_list_register ("DH_r",             "H / r_{c\\mathrm{zd}}",       "NcDistance", "BAO H/(c r_zd)",             NC_TYPE_DISTANCE, _nc_distance_flist_DH_r,             1, 1);
  ncm_mset_func_list_register ("DA_r",             "\\delta{}A / r",              "NcDistance", "BAO dA/r",                   NC_TYPE_DISTANCE, _nc_distance_flist_DA_r,             1, 1);
  ncm_mset_func_list_register ("sound_horizon",    "r_\\mathrm{sound}",           "NcDistance", "Sound horizon",              NC_TYPE_DISTANCE, _nc_distance_flist_sound_horizon,    1, 1);
}

