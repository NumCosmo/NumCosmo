/***************************************************************************
 *            nc_multiplicity_func_despali.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * NcMultiplicityFuncDespali:
 *
 * Dark matter halo -- Despali multiplicity function.
 *
 * Reference: arxiv:1507.05627v2
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include <numcosmo/numcosmo.h>
#include "lss/nc_multiplicity_func_despali.h"
#include "math/ncm_spline_cubic_d2.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

struct _NcMultiplicityFuncDespaliPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble (*eval) (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
  gboolean EO;
  gboolean CMF;
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_EO,
  PROP_CMF,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncDespali, nc_multiplicity_func_despali, NC_TYPE_MULTIPLICITY_FUNC)

static gdouble
_nc_multiplicity_func_despali_eval_error (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  g_error ("method eval not correctly initialized by %s.", G_OBJECT_TYPE_NAME (mulf));

  return 0.0;
}

static void
nc_multiplicity_func_despali_init (NcMultiplicityFuncDespali *md)
{
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  self->mdef  = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->eval  = &_nc_multiplicity_func_despali_eval_error;
  self->EO    = FALSE;
  self->CMF   = FALSE;
  self->Delta = 0.0;
}

static void
_nc_multiplicity_func_despali_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncDespali *md = NC_MULTIPLICITY_FUNC_DESPALI (object);

  /* NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md); */

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_DESPALI (object));

  switch (prop_id)
  {
    case PROP_EO:
      nc_multiplicity_func_despali_set_eo (md, g_value_get_boolean (value));
      break;
    case PROP_CMF:
      nc_multiplicity_func_despali_set_cmf (md, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_multiplicity_func_despali_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (object);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_DESPALI (object));

  switch (prop_id)
  {
    case PROP_EO:
      g_value_set_boolean (value, self->EO);
      break;
    case PROP_CMF:
      g_value_set_boolean (value, self->CMF);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_multiplicity_func_despali_dispose (GObject *object)
{
  /* NcMultiplicityFuncDespali *md = NC_MULTIPLICITY_FUNC_DESPALI (object); */
  /* NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md); */

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_despali_parent_class)->dispose (object);
}

static void
_nc_multiplicity_func_despali_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_despali_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_despali_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static void _nc_multiplicity_func_despali_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static NcMultiplicityFuncMassDef _nc_multiplicity_func_despali_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_despali_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_despali_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_despali_class_init (NcMultiplicityFuncDespaliClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass *parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_despali_set_property;
  object_class->get_property = &_nc_multiplicity_func_despali_get_property;
  object_class->dispose      = &_nc_multiplicity_func_despali_dispose;
  object_class->finalize     = &_nc_multiplicity_func_despali_finalize;

  g_object_class_install_property (object_class,
                                   PROP_EO,
                                   g_param_spec_boolean ("E0",
                                                         NULL,
                                                         "Whether the halo finder uses elliptical overdensity",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  g_object_class_install_property (object_class,
                                   PROP_CMF,
                                   g_param_spec_boolean ("CMF",
                                                         NULL,
                                                         "Whether the use of the cluster mass function",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->set_mdef  = &_nc_multiplicity_func_despali_set_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_despali_set_Delta;
  parent_class->get_mdef  = &_nc_multiplicity_func_despali_get_mdef;
  parent_class->get_Delta = &_nc_multiplicity_func_despali_get_Delta;
  parent_class->eval      = &_nc_multiplicity_func_despali_eval;
}

static gdouble
_nc_multiplicity_func_despali_virial_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Despali: MNRAS 456, 2486–2504 (2016) */
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);
  gdouble f_Despali_virial                      = 0;
  const gdouble delta_c                         = nc_multiplicity_func_despali_delta_c (md, cosmo, z);
  const gdouble nu                              = pow (delta_c / sigma, 2.0);

  if (self->EO)
  {
    const gdouble A        = 0.3953;
    const gdouble a        = 0.7057;
    const gdouble p        = 0.2206;
    const gdouble nu_prime = a * nu;

    f_Despali_virial = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
  }
  else
  {
    if (self->CMF)
    {
      const gdouble A        = 0.8199;
      const gdouble a        = 0.3141;
      const gdouble p        = 0.0;
      const gdouble nu_prime = a * nu;

      f_Despali_virial = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
    }
    else
    {
      const gdouble A        = 0.3292;
      const gdouble a        = 0.7665;
      const gdouble p        = 0.2488;
      const gdouble nu_prime = a * nu;

      f_Despali_virial = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
    }
  }

  return f_Despali_virial;
}

static gdouble
_nc_multiplicity_func_despali_mean_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Despali: MNRAS 456, 2486–2504 (2016) */
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);
  gdouble f_Despali_mean                        = 0;
  const gdouble delta_c                         = nc_multiplicity_func_despali_delta_c (md, cosmo, z);
  const gdouble nu                              = pow (delta_c / sigma, 2.0);
  const gdouble delta_vir                       = nc_multiplicity_func_despali_delta_vir (md, cosmo, z);
  const gdouble Delta                           = self->Delta;
  const gdouble E2                              = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m                         = nc_hicosmo_E2Omega_m (cosmo, z) / E2;

  if (self->EO)
  {
    const gdouble A0 = 0.3953;
    const gdouble A1 = -0.1768;

    const gdouble a0 = 0.7057;
    const gdouble a1 = 0.2125;
    const gdouble a2 = 0.3268;

    const gdouble p0 = 0.2206;
    const gdouble p1 = 0.1937;
    const gdouble p2 = -0.04570;

    const gdouble x = log10 (Delta * Omega_m / delta_vir);
    const gdouble A = A0 + A1 * x;
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    f_Despali_mean = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
  }
  else
  {
    const gdouble A0 = 0.3292;
    const gdouble A1 = -0.1362;

    const gdouble a0 = 0.7665;
    const gdouble a1 = 0.2263;
    const gdouble a2 = 0.4332;

    const gdouble p0 = 0.2488;
    const gdouble p1 = 0.2554;
    const gdouble p2 = -0.1151;

    const gdouble x = log10 (Delta * Omega_m / delta_vir);
    const gdouble A = A0 + A1 * x;
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    f_Despali_mean = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
  }

  return f_Despali_mean;
}

static gdouble
_nc_multiplicity_func_despali_crit_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Despali: MNRAS 456, 2.0486–2504 (2016) */
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);
  gdouble f_Despali_crit                        = 0;
  const gdouble delta_c                         = nc_multiplicity_func_despali_delta_c (md, cosmo, z);
  const gdouble nu                              = pow (delta_c / sigma, 2.0);
  const gdouble delta_vir                       = nc_multiplicity_func_despali_delta_vir (md, cosmo, z);
  const gdouble Delta                           = self->Delta;

  if (self->EO)
  {
    const gdouble A0 = 0.3953;
    const gdouble A1 = -0.1768;

    const gdouble a0 = 0.7057;
    const gdouble a1 = 0.2125;
    const gdouble a2 = 0.3268;

    const gdouble p0 = 0.2206;
    const gdouble p1 = 0.1937;
    const gdouble p2 = -0.04570;

    const gdouble x = log10 (Delta / delta_vir);
    const gdouble A = A0 + A1 * x;
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    f_Despali_crit = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
  }
  else
  {
    const gdouble A0 = 0.3292;
    const gdouble A1 = -0.1362;

    const gdouble a0 = 0.7665;
    const gdouble a1 = 0.2263;
    const gdouble a2 = 0.4332;

    const gdouble p0 = 0.2488;
    const gdouble p1 = 0.2554;
    const gdouble p2 = -0.1151;

    const gdouble x = log10 (Delta / delta_vir);
    const gdouble A = A0 + A1 * x;
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    f_Despali_crit = 2.0 * A * (1.0 + pow (nu_prime, -p)) * pow (nu_prime, 0.5) / ncm_c_sqrt_2pi () * exp (-nu_prime * 0.5);
  }

  return f_Despali_crit;
}

static void
_nc_multiplicity_func_despali_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      self->eval = &_nc_multiplicity_func_despali_mean_eval;
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      self->eval = &_nc_multiplicity_func_despali_crit_eval;
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      self->eval = &_nc_multiplicity_func_despali_virial_eval;
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncDespali does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef
_nc_multiplicity_func_despali_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_despali_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Despali: MNRAS 456, 2486–2504 (2016) */
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  return self->eval (mulf, cosmo, sigma, z);
}

/**
 * nc_multiplicity_func_despali_new:
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncDespali.
 */
NcMultiplicityFuncDespali *
nc_multiplicity_func_despali_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_DESPALI,
                       NULL);
}

/**
 * nc_multiplicity_func_despali_new_full:
 * @mdef: a #NcMultiplicityFuncMassDef
 * @Delta: parameter that multiplies the background mass density (mean ou critical)
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncDespali.
 */
NcMultiplicityFuncDespali *
nc_multiplicity_func_despali_new_full (NcMultiplicityFuncMassDef mdef, gdouble Delta)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_DESPALI,
                       "mass-def", mdef,
                       "Delta",    Delta,
                       NULL);
}

/**
 * nc_multiplicity_func_despali_ref:
 * @md: a #NcMultiplicityFuncDespali
 *
 * Increases the reference count of @md by one.
 *
 * Returns: (transfer full): @md
 */
NcMultiplicityFuncDespali *
nc_multiplicity_func_despali_ref (NcMultiplicityFuncDespali *md)
{
  return g_object_ref (md);
}

/**
 * nc_multiplicity_func_despali_free:
 * @md: a #NcMultiplicityFuncDespali
 *
 * Atomically decrements the reference count of @md by one. If the reference count drops to 0,
 * all memory allocated by @md is released.
 *
 */
void
nc_multiplicity_func_despali_free (NcMultiplicityFuncDespali *md)
{
  g_object_unref (md);
}

/**
 * nc_multiplicity_func_despali_clear:
 * @md: a #NcMultiplicityFuncDespali
 *
 * Atomically decrements the reference count of @md by one. If the reference count drops to 0,
 * all memory allocated by @md is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_despali_clear (NcMultiplicityFuncDespali **md)
{
  g_clear_object (md);
}

/**
 * nc_multiplicity_func_despali_set_Delta:
 * @md: a #NcMultiplicityFuncDespali.
 * @Delta: value of #NcMultiplicityFuncDespali:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncDespali:Delta property.
 *
 */
static void
_nc_multiplicity_func_despali_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  self->Delta = Delta;
}

/**
 * _nc_multiplicity_func_despali_get_Delta:
 * @md: a #NcMultiplicityFuncDespali.
 *
 * Returns: the value of #NcMultiplicityFuncDespali:Delta property.
 */
gdouble
_nc_multiplicity_func_despali_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncDespali *md                 = NC_MULTIPLICITY_FUNC_DESPALI (mulf);
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  return self->Delta;
}

/**
 * nc_multiplicity_func_despali_delta_c:
 * @md: a #NcMultiplicityFuncDespali.
 * @cosmo: a #NcHICosmo
 * @z: a @gdouble
 *
 * Calculates the critical density using Kitayama & Suto 1996 interpolation
 * (https://arxiv.org/pdf/astro-ph/9604141)
 *
 */
gdouble
nc_multiplicity_func_despali_delta_c (NcMultiplicityFuncDespali *md, NcHICosmo *cosmo, gdouble z)
{
  /* NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md); */
  /* NcMultiplicityFunc *mulf                      = NC_MULTIPLICITY_FUNC (md); */
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;

  return 3.0 / 20.0 * pow (12.0 * M_PI, 2.0 / 3.0) * (1 + 0.012299 * log10 (Omega_m));
}

/**
 * nc_multiplicity_func_despali_delta_vir:
 * @md: a #NcMultiplicityFuncDespali.
 * @cosmo: a #NcHICosmo
 * @z: a @gdouble
 *
 * Calculates the virial delta using Bryan and Norman 1998 interpolation
 * (https://arxiv.org/pdf/astro-ph/9710107)
 *
 */
gdouble
nc_multiplicity_func_despali_delta_vir (NcMultiplicityFuncDespali *md, NcHICosmo *cosmo, gdouble z)
{
  /* NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md); */
  /* NcMultiplicityFunc *mulf                      = NC_MULTIPLICITY_FUNC (md); */
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
  const gdouble x       = Omega_m - 1.0;

  const gboolean is_flat = ncm_cmp ((nc_hicosmo_Omega_k0 (cosmo)), (0.0), 1.0e5, 0.0) == 0;
  if (is_flat)
  {
    return 18.0 * pow (M_PI, 2.0) + 82.0 * x - 39.0 * x * x;
  }
  else
  {
    g_error ("Interpolation does not work in this regime.");
}
}

/**
 * nc_multiplicity_func_despali_set_eo:
 * @md: a #NcMultiplicityFuncDespali
 * @on: Whether the halo finder uses eliptical overdensidy.
 *
 * Sets array of #Set if halo finder uses eliptical overdensidy.
 *
 */
void
nc_multiplicity_func_despali_set_eo (NcMultiplicityFuncDespali *md, gboolean on)
{
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  self->EO = on;
}

/**
 * nc_multiplicity_func_despali_get_eo:
 * @md: a #NcMultiplicityFuncDespali
 *
 * Gets if the eo option is on.
 *
 * Returns: TRUE or FALSE.
 */
gboolean
nc_multiplicity_func_despali_get_eo (NcMultiplicityFuncDespali *md)
{
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  return self->EO;
}

/**
 * nc_multiplicity_func_despali_set_cmf:
 * @md: a #NcMultiplicityFuncDespali
 * @on: Whether the we use cluster mass function.
 *
 * Sets array of #Set if  uses eliptical  mass function.
 *
 */
void
nc_multiplicity_func_despali_set_cmf (NcMultiplicityFuncDespali *md, gboolean on)
{
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  self->CMF = on;
}

/**
 * nc_multiplicity_func_despali_get_cmf:
 * @md: a #NcMultiplicityFuncDespali
 *
 * Gets if the cmf option is on.
 *
 * Returns: TRUE or FALSE.
 */
gboolean
nc_multiplicity_func_despali_get_cmf (NcMultiplicityFuncDespali *md)
{
  NcMultiplicityFuncDespaliPrivate * const self = nc_multiplicity_func_despali_get_instance_private (md);

  return self->CMF;
}

