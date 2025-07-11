/***************************************************************************
 *            nc_halo_bias_despali.c
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
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
 * NcHaloBiasDespali:
 *
 * Despali halo bias function type.
 *
 * Object implementation to compute the halo bias function given
 * the Despali mass function. A description of the mechanism
 * is given below. Check nc_halo_bias.c for a description
 * of halo biases and nc_halo_bias_despali.c for the Despali
 * mass function.
 *
 * The Despali bias was obtained empirically and is given by
 * \begin{align}
 * b(\nu) &= 1 - A \frac{\nu^a}{\nu^a + \delta_c^a} + B \nu^b + C \nu^c
 * , \end{align}
 * where $b(\nu)$ is the Despali bias, $\delta_c$ is the
 * critical threshold, $\nu = \frac{\delta_c}{\sigma}$ and the free parameters
 * $(A,a,B,b,C,c)$ depend on the value of the overdensity chosen.
 *
 * The user must provide input the values: @NcHaloMassFunction, @delta_c, @B, @b and @c @delta_c  - nc_halo_bias_ps_new_full().
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_despali.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcHaloBiasDespali, nc_halo_bias_despali, NC_TYPE_HALO_BIAS)

enum
{
  PROP_0,
  PROP_EO,
  PROP_CMF,
  PROP_SIZE
};

static gdouble _nc_halo_bias_despali_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
static gdouble _nc_halo_bias_despali_virial_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
static gdouble _nc_halo_bias_despali_mean_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
static gdouble _nc_halo_bias_despali_crit_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_halo_bias_despali_init (NcHaloBiasDespali *biasf_despali)
{
  biasf_despali->eo  = FALSE;
  biasf_despali->cmf = FALSE;
}

static void
_nc_halo_bias_despali_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_despali_parent_class)->finalize (object);
}

static void
_nc_halo_bias_despali_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloBiasDespali *biasf_despali = NC_HALO_BIAS_DESPALI (object);

  g_return_if_fail (NC_IS_HALO_BIAS_DESPALI (object));

  switch (prop_id)
  {
    case PROP_EO:
      nc_halo_bias_despali_set_eo (biasf_despali, g_value_get_boolean (value));
      break;
    case PROP_CMF:
      nc_halo_bias_despali_set_cmf (biasf_despali, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_bias_despali_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasDespali *biasf_despali = NC_HALO_BIAS_DESPALI (object);

  g_return_if_fail (NC_IS_HALO_BIAS_DESPALI (object));

  switch (prop_id)
  {
    case PROP_EO:
      g_value_set_boolean (value, biasf_despali->eo);
      break;
    case PROP_CMF:
      g_value_set_boolean (value, biasf_despali->cmf);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_bias_despali_class_init (NcHaloBiasDespaliClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcHaloBiasClass *parent_class = NC_HALO_BIAS_CLASS (klass);

  object_class->finalize     = _nc_halo_bias_despali_finalize;
  object_class->set_property = _nc_halo_bias_despali_set_property;
  object_class->get_property = _nc_halo_bias_despali_get_property;

  g_object_class_install_property (object_class,
                                   PROP_EO,
                                   g_param_spec_boolean ("eo",
                                                         NULL,
                                                         "Whether the halo finder uses eliptical overdensity",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  g_object_class_install_property (object_class,
                                   PROP_CMF,
                                   g_param_spec_boolean ("cmf",
                                                         NULL,
                                                         "Whether the use of the cluster mass function",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->eval = &_nc_halo_bias_despali_eval;
}

/**
 * nc_halo_bias_despali_new:
 * @mfp: a #NcHaloMassFunction
 *
 * Creates a new #NcHaloBiasDespali object with undefined parameters.
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBiasDespali *
nc_halo_bias_despali_new (NcHaloMassFunction *mfp)
{
  return g_object_new (NC_TYPE_HALO_BIAS_DESPALI,
                       "mass-function", mfp,
                       NULL);
}

/**
 * nc_halo_bias_despali_new_full:
 * @mfp: a #NcHaloMassFunction
 * @eo: Empirical parameter for Despali bias function
 * @cmf: Empirical parameter for Despali bias function
 *
 * Creates a new #NcHaloBiasDespali object with the input parameters.
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBiasDespali *
nc_halo_bias_despali_new_full (NcHaloMassFunction *mfp, gboolean eo, gboolean cmf)
{
  return g_object_new (NC_TYPE_HALO_BIAS_DESPALI,
                       "mass-function", mfp,
                       "eo", eo,
                       "cmf", cmf,
                       NULL);
}

/**
 * nc_halo_bias_despali_ref:
 * @biasf_despali: a #NcHaloBiasDespali
 *
 * Increases the reference count of the object.
 *
 * Returns: (transfer full): The object itself.
 */
NcHaloBiasDespali *
nc_halo_bias_despali_ref (NcHaloBiasDespali *biasf_despali)
{
  return NC_HALO_BIAS_DESPALI (g_object_ref (biasf_despali));
}

/**
 * nc_halo_bias_despali_free:
 * @biasf_despali: a #NcHaloBiasDespali
 *
 * Decreases the reference count of the object. If the reference count
 * reaches zero, the object is destroyed.
 *
 */
void
nc_halo_bias_despali_free (NcHaloBiasDespali *biasf_despali)
{
  g_object_unref (biasf_despali);
}

/**
 * nc_halo_bias_despali_clear:
 * @biasf_despali: a #NcHaloBiasDespali
 *
 * If *@biasf_despali is not %NULL, decreases the reference count of the object and sets
 * *@biasf_despali to %NULL. If the reference count reaches zero, the object is
 * destroyed.
 *
 */
void
nc_halo_bias_despali_clear (NcHaloBiasDespali **biasf_despali)
{
  g_clear_object (biasf_despali);
}

static gdouble
_nc_halo_bias_despali_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  /* NcHaloBiasDespali *biasf_despali = NC_HALO_BIAS_DESPALI (biasf); */
  NcMultiplicityFunc *mulf       = nc_halo_mass_function_peek_multiplicity_function (biasf->mfp);
  NcMultiplicityFuncMassDef mdef = nc_multiplicity_func_get_mdef (mulf);

  gdouble eval = 0.0;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      eval = _nc_halo_bias_despali_mean_eval (biasf, cosmo, sigma, z);
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      eval = _nc_halo_bias_despali_crit_eval (biasf, cosmo, sigma, z);
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      eval = _nc_halo_bias_despali_virial_eval (biasf, cosmo, sigma, z);
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcHaloBiasDespali does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  return eval;
}

static gdouble
_nc_halo_bias_despali_virial_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  NcHaloBiasDespali *biasf_despali = NC_HALO_BIAS_DESPALI (biasf);
  gdouble bias_Despali_virial      = 0;


  const gdouble delta_c = nc_halo_bias_despali_delta_c (biasf_despali, cosmo, z);
  const gdouble nu      = pow (delta_c / sigma, 2.0);

  if (biasf_despali->eo)
  {
    /* const gdouble A        = 0.3953; */
    const gdouble a        = 0.7057;
    const gdouble p        = 0.2206;
    const gdouble nu_prime = a * nu;

    bias_Despali_virial = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
  }
  else
  {
    if (biasf_despali->cmf)
    {
      /* const gdouble A        = 0.8199; */
      const gdouble a        = 0.3141;
      const gdouble p        = 0.0;
      const gdouble nu_prime = a * nu;

      bias_Despali_virial = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
    }
    else
    {
      /* const gdouble A        = 0.3295; */
      const gdouble a        = 0.7689;
      const gdouble p        = 0.2536;
      const gdouble nu_prime = a * nu;

      bias_Despali_virial = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
    }
  }

  return bias_Despali_virial;
}

static gdouble
_nc_halo_bias_despali_mean_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  NcHaloBiasDespali *biasf_despali = NC_HALO_BIAS_DESPALI (biasf);
  NcMultiplicityFunc *mulf         = nc_halo_mass_function_peek_multiplicity_function (biasf->mfp);


  gdouble bias_Despali_mean = 0;
  const gdouble delta_c     = nc_halo_bias_despali_delta_c (biasf_despali, cosmo, z);
  const gdouble nu          = pow (delta_c / sigma, 2.0);
  const gdouble delta_vir   = nc_halo_bias_despali_delta_vir (biasf_despali, cosmo, z);
  const gdouble Delta       = nc_multiplicity_func_get_Delta (mulf);
  const gdouble Omega_m     = nc_hicosmo_E2Omega_m (cosmo, z);

  if (biasf_despali->eo)
  {
    /* const gdouble A0 = 0.3953; */
    /* const gdouble A1 = -0.1768; */

    const gdouble a0 = 0.7057;
    const gdouble a1 = 0.2125;
    const gdouble a2 = 0.3268;

    const gdouble p0 = 0.2206;
    const gdouble p1 = 0.1937;
    const gdouble p2 = -0.04570;

    const gdouble x = log10 (Delta * Omega_m / delta_vir);
    /* const gdouble A = A0 + A1 * x; */
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    bias_Despali_mean = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
  }
  else
  {
    /* const gdouble A0 = 0.3292; */
    /* const gdouble A1 = -0.1362; */

    const gdouble a0 = 0.7665;
    const gdouble a1 = 0.2263;
    const gdouble a2 = 0.4332;

    const gdouble p0 = 0.2488;
    const gdouble p1 = 0.2554;
    const gdouble p2 = -0.1151;

    const gdouble x = log10 (Delta * Omega_m / delta_vir);
    /* const gdouble A = A0 + A1 * x; */
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    bias_Despali_mean = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
  }

  return bias_Despali_mean;
}

static gdouble
_nc_halo_bias_despali_crit_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  NcHaloBiasDespali *biasf_despali = NC_HALO_BIAS_DESPALI (biasf);
  NcMultiplicityFunc *mulf         = nc_halo_mass_function_peek_multiplicity_function (biasf->mfp);

  gdouble bias_Despali_crit = 0;
  const gdouble delta_c     = nc_halo_bias_despali_delta_c (biasf_despali, cosmo, z);
  const gdouble nu          = pow (delta_c / sigma, 2.0);
  const gdouble delta_vir   = nc_halo_bias_despali_delta_vir (biasf_despali, cosmo, z);
  const gdouble Delta       = nc_multiplicity_func_get_Delta (mulf);

  if (biasf_despali->eo)
  {
    /* const gdouble A0 = 0.3953; */
    /* const gdouble A1 = -0.1768; */

    const gdouble a0 = 0.7057;
    const gdouble a1 = 0.2125;
    const gdouble a2 = 0.3268;

    const gdouble p0 = 0.2206;
    const gdouble p1 = 0.1937;
    const gdouble p2 = -0.04570;

    const gdouble x = log10 (Delta / delta_vir);
    /* const gdouble A = A0 + A1 * x; */
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    bias_Despali_crit = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
  }
  else
  {
    /* const gdouble A0 = 0.3292; */
    /* const gdouble A1 = -0.1362; */

    const gdouble a0 = 0.7665;
    const gdouble a1 = 0.2263;
    const gdouble a2 = 0.4332;

    const gdouble p0 = 0.2488;
    const gdouble p1 = 0.2554;
    const gdouble p2 = -0.1151;

    const gdouble x = log10 (Delta / delta_vir);
    /* const gdouble A = A0 + A1 * x; */
    const gdouble a = a0 + a1 * x + a2 * x * x;
    const gdouble p = p0 + p1 * x + p2 * x * x;

    const gdouble nu_prime = a * nu;

    bias_Despali_crit = 1 + nu_prime / (2 * delta_c) + p / (delta_c * (pow (nu_prime, p) + 1)) - 3 / (2 * delta_c);
  }

  return bias_Despali_crit;
}

/* _NC_BIAS_FUNCTION_DESPALI_DATASET_1001_3162_DELTA = {1.686, 0.183, 1.5, 2.4, 200.0}; */

/**
 * nc_halo_bias_despali_delta_c:
 * @biasf_despali: a #NcHaloBiasDespali.
 * @cosmo: a #NcHICosmo
 * @z: a @gdouble
 *
 * Calculates the critical density using Kitayama & Suto 1996 interpolation
 * (https://arxiv.org/pdf/astro-ph/9604141).
 *
 */
gdouble
nc_halo_bias_despali_delta_c (NcHaloBiasDespali *biasf_despali, NcHICosmo *cosmo, gdouble z)
{
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;

  return 3.0 / 20.0 * pow (12.0 * M_PI, 2.0 / 3.0) * (1 + 0.012299 * log10 (Omega_m));
}

/**
 * nc_halo_bias_despali_delta_vir:
 * @biasf_despali: a #NcHaloBiasDespali.
 * @cosmo: a #NcHICosmo
 * @z: a @gdouble
 *
 * Calculates the virial delta using Bryan and Norman 1998 interpolation(https://arxiv.org/pdf/astro-ph/9710107)
 *
 */
gdouble
nc_halo_bias_despali_delta_vir (NcHaloBiasDespali *biasf_despali, NcHICosmo *cosmo, gdouble z)
{
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
  const gdouble x       = Omega_m - 1.0;

  if (nc_hicosmo_Omega_k0 (cosmo) == 0)
    return 18.0 * pow (M_PI, 2.0) + 82.0 * x - 39.0 * x * x;

  else
    g_error ("Interpolation does not work in this regime.");
}

/**
 * nc_halo_bias_despali_set_eo:
 * @biasf_despali: a #NcHaloBiasDespali
 * @on: Whether the halo finder uses eliptical overdensidy.
 *
 * Sets array of #Set if halo finder uses eliptical overdensidy.
 *
 */
void
nc_halo_bias_despali_set_eo (NcHaloBiasDespali *biasf_despali, gboolean on)
{
  biasf_despali->eo = on;
}

/**
 * nc_halo_bias_despali_get_eo:
 * @biasf_despali: a #NcHaloBiasDespali
 *
 * Gets if the eo option is on.
 *
 * Returns: TRUE or FALSE.
 */
gboolean
nc_halo_bias_despali_get_eo (NcHaloBiasDespali *biasf_despali)
{
  return biasf_despali->eo;
}

/**
 * nc_halo_bias_despali_set_cmf:
 * @biasf_despali: a #NcHaloBiasDespali
 * @on: Whether the we use cluster mass function.
 *
 * Sets array of #Set if  uses eliptical  mass function.
 *
 */
void
nc_halo_bias_despali_set_cmf (NcHaloBiasDespali *biasf_despali, gboolean on)
{
  biasf_despali->cmf = on;
}

/**
 * nc_halo_bias_despali_get_cmf:
 * @biasf_despali: a #NcHaloBiasDespali
 *
 * Gets if the cmf option is on.
 *
 * Returns: TRUE or FALSE.
 */
gboolean
nc_halo_bias_despali_get_cmf (NcHaloBiasDespali *biasf_despali)
{
  return biasf_despali->cmf;
}

