/***************************************************************************
 *            nc_halo_bias_castro.c
 *
 *  Thu Aug 14 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 */

/**
 * NcHaloBiasCastro:
 *
 * Linear halo bias calibrated together with the Castro halo mass function.
 *
 * The peak-background split prediction obtained from #NcMultiplicityFuncCastro,
 * $$ b_\mathrm{PBS} = 1 - \frac{1}{\delta_c}\frac{\mathrm{d}\ln f}{\mathrm{d}\ln\nu}, $$
 * multiplied by a calibrated correction that depends on $\Omega_m(z)$, on the slope
 * $s = \mathrm{d}\ln\sigma_R/\mathrm{d}\ln R$ and on the clustering amplitude $S_8$,
 * $$ b = b_\mathrm{PBS}\,A_0\,(1 + a_1\Omega_m(z))\,(1 + b_1 s + b_2 s^2)\,(1 + c_1 S_8). $$
 *
 * The derivative is the total one along the mass direction, which is what the
 * correction was calibrated against. It is not obtained by differentiating an
 * interpolation: the model is closed-form in $\sigma_R$ and $s$, so the two partial
 * derivatives are taken by perturbing those numbers directly through
 * nc_multiplicity_func_castro_eval_full(), and the chain rule closes with
 * $\mathrm{d}s/\mathrm{d}\ln R$ obtained from the transform itself.
 *
 * The mass function this bias is attached to must use a #NcMultiplicityFuncCastro,
 * since the correction coefficients were fitted alongside that multiplicity function.
 *
 * See the <a href="../../theory/castro_hmf.html">Castro Halo Mass Function and Bias</a>
 * theory page, and [Castro et al. (2024)](https://arxiv.org/abs/2409.01877).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/halo/nc_halo_bias_castro.h"
#include "nc/lss/halo/nc_multiplicity_func_castro.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcHaloBiasCastroPrivate
{
  gdouble A0;
  gdouble a1;
  gdouble b1;
  gdouble b2;
  gdouble c1;
} NcHaloBiasCastroPrivate;

enum
{
  PROP_0,
  PROP_A0,
  PROP_a1,
  PROP_b1,
  PROP_b2,
  PROP_c1,
  PROP_SIZE,
};

struct _NcHaloBiasCastro
{
  NcHaloBias parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloBiasCastro, nc_halo_bias_castro, NC_TYPE_HALO_BIAS)

/*
 * Step used for the two partial derivatives. They are taken on the closed-form
 * model, whose arguments are plain numbers rather than interpolated quantities,
 * so a central difference converges as the step squared with no interpolation
 * noise to trade against: 1e-5 sits at ~1e-10, close to the round-off floor.
 */
#define NC_HALO_BIAS_CASTRO_DIFF_STEP (1.0e-5)

static void
nc_halo_bias_castro_init (NcHaloBiasCastro *biasf)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  self->A0 = 0.0;
  self->a1 = 0.0;
  self->b1 = 0.0;
  self->b2 = 0.0;
  self->c1 = 0.0;
}

static void
_nc_halo_bias_castro_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloBiasCastro *biasf = NC_HALO_BIAS_CASTRO (object);

  g_return_if_fail (NC_IS_HALO_BIAS_CASTRO (object));

  switch (prop_id)
  {
    case PROP_A0:
      nc_halo_bias_castro_set_A0 (biasf, g_value_get_double (value));
      break;
    case PROP_a1:
      nc_halo_bias_castro_set_a1 (biasf, g_value_get_double (value));
      break;
    case PROP_b1:
      nc_halo_bias_castro_set_b1 (biasf, g_value_get_double (value));
      break;
    case PROP_b2:
      nc_halo_bias_castro_set_b2 (biasf, g_value_get_double (value));
      break;
    case PROP_c1:
      nc_halo_bias_castro_set_c1 (biasf, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_bias_castro_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasCastro *biasf = NC_HALO_BIAS_CASTRO (object);

  g_return_if_fail (NC_IS_HALO_BIAS_CASTRO (object));

  switch (prop_id)
  {
    case PROP_A0:
      g_value_set_double (value, nc_halo_bias_castro_get_A0 (biasf));
      break;
    case PROP_a1:
      g_value_set_double (value, nc_halo_bias_castro_get_a1 (biasf));
      break;
    case PROP_b1:
      g_value_set_double (value, nc_halo_bias_castro_get_b1 (biasf));
      break;
    case PROP_b2:
      g_value_set_double (value, nc_halo_bias_castro_get_b2 (biasf));
      break;
    case PROP_c1:
      g_value_set_double (value, nc_halo_bias_castro_get_c1 (biasf));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_bias_castro_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_castro_parent_class)->finalize (object);
}

static gdouble _nc_halo_bias_castro_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble lnM, gdouble z);

static void
_nc_halo_bias_castro_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_halo_bias_castro_parent_class)->constructed (object);
  {
    NcHaloBias *bias      = NC_HALO_BIAS (object);
    NcmPowspecFilter *psf = nc_halo_mass_function_peek_psf (nc_halo_bias_peek_mass_function (bias));

    /*
     * The slope-running term of the total derivative needs d2lnvar/dlnR2. Declare
     * it before anything prepares the filter; the requirement only ever raises the
     * order, so it cannot disturb another user of the same filter.
     */
    ncm_powspec_filter_require_nderivs (psf, 2);
  }
}

static void
nc_halo_bias_castro_class_init (NcHaloBiasCastroClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcHaloBiasClass *parent_class = NC_HALO_BIAS_CLASS (klass);

  object_class->constructed  = &_nc_halo_bias_castro_constructed;
  object_class->set_property = &_nc_halo_bias_castro_set_property;
  object_class->get_property = &_nc_halo_bias_castro_get_property;
  object_class->finalize     = &_nc_halo_bias_castro_finalize;

  /**
   * NcHaloBiasCastro:A0:
   *
   * Overall amplitude of the correction to the peak-background split prediction.
   */
  g_object_class_install_property (object_class,
                                   PROP_A0,
                                   g_param_spec_double ("A0",
                                                        NULL,
                                                        "Correction amplitude",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.150,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloBiasCastro:a1:
   *
   * Coefficient of the $\Omega_m(z)$ term of the correction.
   */
  g_object_class_install_property (object_class,
                                   PROP_a1,
                                   g_param_spec_double ("a1",
                                                        NULL,
                                                        "Omega_m(z) coefficient",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0929,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloBiasCastro:b1:
   *
   * Coefficient of the linear slope term of the correction.
   */
  g_object_class_install_property (object_class,
                                   PROP_b1,
                                   g_param_spec_double ("b1",
                                                        NULL,
                                                        "Linear slope coefficient",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.256,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloBiasCastro:b2:
   *
   * Coefficient of the quadratic slope term of the correction.
   */
  g_object_class_install_property (object_class,
                                   PROP_b2,
                                   g_param_spec_double ("b2",
                                                        NULL,
                                                        "Quadratic slope coefficient",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.173,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloBiasCastro:c1:
   *
   * Coefficient of the $S_8$ term of the correction.
   */
  g_object_class_install_property (object_class,
                                   PROP_c1,
                                   g_param_spec_double ("c1",
                                                        NULL,
                                                        "S8 coefficient",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, -0.0372,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->eval = &_nc_halo_bias_castro_eval;
}

/*
 * The multiplicity function must be a Castro one: the correction coefficients were
 * fitted alongside it, so pairing them with another f(sigma) is not the model.
 */
static NcMultiplicityFuncCastro *
_nc_halo_bias_castro_peek_mulf (NcHaloBias *biasf)
{
  NcMultiplicityFunc *mulf = nc_halo_mass_function_peek_multiplicity_function (nc_halo_bias_peek_mass_function (biasf));

  if (!NC_IS_MULTIPLICITY_FUNC_CASTRO (mulf))
    g_error ("NcHaloBiasCastro: the mass function must use a NcMultiplicityFuncCastro, "
             "got %s. The bias correction was calibrated together with that "
             "multiplicity function and is not applicable to another one.",
             G_OBJECT_TYPE_NAME (mulf));

  return NC_MULTIPLICITY_FUNC_CASTRO (mulf);
}

/**
 * nc_halo_bias_castro_S8:
 * @biasf: a #NcHaloBiasCastro
 * @cosmo: a #NcHICosmo
 *
 * Clustering amplitude $S_8 = \sigma_8\sqrt{\Omega_{m,0}/0.3}$, with $\sigma_8$ the
 * present-day variance on $8\,h^{-1}\,\mathrm{Mpc}$. NumCosmo works in $\mathrm{Mpc}$,
 * so the radius is converted with the reduced Hubble parameter.
 *
 * Returns: $S_8$.
 */
gdouble
nc_halo_bias_castro_S8 (NcHaloBiasCastro *biasf, NcHICosmo *cosmo)
{
  NcHaloBias *bias      = NC_HALO_BIAS (biasf);
  NcmPowspecFilter *psf = nc_halo_mass_function_peek_psf (nc_halo_bias_peek_mass_function (bias));
  const gdouble h       = nc_hicosmo_h (cosmo);
  const gdouble lnR8    = log (8.0 / h);
  const gdouble sigma8  = ncm_powspec_filter_eval_sigma_lnr (psf, 0.0, lnR8);

  return sigma8 * sqrt (nc_hicosmo_Omega_m0 (cosmo) / 0.3);
}

/**
 * nc_halo_bias_castro_correction:
 * @biasf: a #NcHaloBiasCastro
 * @cosmo: a #NcHICosmo
 * @dlnsigma_dlnR: the slope $\mathrm{d}\ln\sigma_R/\mathrm{d}\ln R$
 * @z: redshift $z$
 *
 * Calibrated correction multiplying the peak-background split prediction,
 * $A_0(1 + a_1\Omega_m(z))(1 + b_1 s + b_2 s^2)(1 + c_1 S_8)$.
 *
 * Returns: the correction factor.
 */
gdouble
nc_halo_bias_castro_correction (NcHaloBiasCastro *biasf, NcHICosmo *cosmo, gdouble dlnsigma_dlnR, gdouble z)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);
  const gdouble Om_z                   = nc_hicosmo_E2Omega_m (cosmo, z) / nc_hicosmo_E2 (cosmo, z);
  const gdouble S8                     = nc_halo_bias_castro_S8 (biasf, cosmo);

  const gdouble f0 = 1.0 + self->a1 * Om_z;
  const gdouble f1 = 1.0 + self->b1 * dlnsigma_dlnR + self->b2 * gsl_pow_2 (dlnsigma_dlnR);
  const gdouble f2 = 1.0 + self->c1 * S8;

  return self->A0 * f0 * f1 * f2;
}

/**
 * nc_halo_bias_castro_pbs:
 * @biasf: a #NcHaloBiasCastro
 * @cosmo: a #NcHICosmo
 * @sigma: standard deviation of the matter density contrast $\sigma_R$
 * @dlnsigma_dlnR: the slope $\mathrm{d}\ln\sigma_R/\mathrm{d}\ln R$
 * @z: redshift $z$
 *
 * Peak-background split prediction
 * $b_\mathrm{PBS} = 1 - \delta_c^{-1}\,\mathrm{d}\ln f/\mathrm{d}\ln\nu$, with the
 * derivative taken at fixed slope. This is the partial derivative: for a scale-free
 * spectrum it is already the total one, and nc_halo_bias_eval() adds the slope-running
 * term where the spectrum is curved.
 *
 * Returns: $b_\mathrm{PBS}$.
 */
gdouble
nc_halo_bias_castro_pbs (NcHaloBiasCastro *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z)
{
  NcMultiplicityFuncCastro *mc = _nc_halo_bias_castro_peek_mulf (NC_HALO_BIAS (biasf));
  const gdouble delta_c        = nc_multiplicity_func_castro_delta_c (mc, cosmo, z);
  const gdouble h              = NC_HALO_BIAS_CASTRO_DIFF_STEP;

  /* nu = delta_c / sigma, so a step in ln(nu) is minus the same step in ln(sigma). */
  const gdouble lnfp = nc_multiplicity_func_castro_eval_lnf (mc, cosmo, sigma * exp (-h), dlnsigma_dlnR, z);
  const gdouble lnfm = nc_multiplicity_func_castro_eval_lnf (mc, cosmo, sigma * exp (h), dlnsigma_dlnR, z);

  const gdouble dlnf_dlnnu = (lnfp - lnfm) / (2.0 * h);

  return 1.0 - dlnf_dlnnu / delta_c;
}

static gdouble
_nc_halo_bias_castro_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble lnM, gdouble z)
{
  NcHaloBiasCastro *bc         = NC_HALO_BIAS_CASTRO (biasf);
  NcMultiplicityFuncCastro *mc = _nc_halo_bias_castro_peek_mulf (biasf);
  NcHaloMassFunction *mfp      = nc_halo_bias_peek_mass_function (biasf);
  NcmPowspecFilter *psf        = nc_halo_mass_function_peek_psf (mfp);
  const gdouble lnR            = nc_halo_mass_function_lnM_to_lnR (mfp, cosmo, lnM);
  const gdouble s              = 0.5 * ncm_powspec_filter_eval_dlnvar_dlnr (psf, z, lnR);
  const gdouble delta_c        = nc_multiplicity_func_castro_delta_c (mc, cosmo, z);
  const gdouble h              = NC_HALO_BIAS_CASTRO_DIFF_STEP;

  /*
   * Total derivative along the mass direction, which is the definition the
   * correction was calibrated against:
   *
   *   dlnf/dlnnu = dlnf/dlnnu|_s + dlnf/ds * ds/dlnnu,
   *
   * with dlnnu/dlnR = -s and ds/dlnR = 0.5 d2lnvar/dlnR2. Both partials are taken
   * on the closed-form model, so neither differentiates an interpolation.
   */
  const gdouble lnf_sp             = nc_multiplicity_func_castro_eval_lnf (mc, cosmo, sigma * exp (-h), s, z);
  const gdouble lnf_sm             = nc_multiplicity_func_castro_eval_lnf (mc, cosmo, sigma * exp (h), s, z);
  const gdouble dlnf_dlnnu_partial = (lnf_sp - lnf_sm) / (2.0 * h);

  const gdouble lnf_ap  = nc_multiplicity_func_castro_eval_lnf (mc, cosmo, sigma, s + h, z);
  const gdouble lnf_am  = nc_multiplicity_func_castro_eval_lnf (mc, cosmo, sigma, s - h, z);
  const gdouble dlnf_ds = (lnf_ap - lnf_am) / (2.0 * h);

  const gdouble ds_dlnR  = 0.5 * ncm_powspec_filter_eval_dnlnvar_dlnrn (psf, z, lnR, 2);
  const gdouble ds_dlnnu = (s != 0.0) ? -ds_dlnR / s : 0.0;

  const gdouble dlnf_dlnnu = dlnf_dlnnu_partial + dlnf_ds * ds_dlnnu;
  const gdouble b_pbs      = 1.0 - dlnf_dlnnu / delta_c;

  return b_pbs * nc_halo_bias_castro_correction (bc, cosmo, s, z);
}

/**
 * nc_halo_bias_castro_new:
 * @mfp: a #NcHaloMassFunction
 *
 * Creates a new #NcHaloBiasCastro with the calibrated correction coefficients.
 * The mass function @mfp must use a #NcMultiplicityFuncCastro.
 *
 * Returns: A new #NcHaloBiasCastro.
 */
NcHaloBiasCastro *
nc_halo_bias_castro_new (NcHaloMassFunction *mfp)
{
  return g_object_new (NC_TYPE_HALO_BIAS_CASTRO,
                       "mass-function", mfp,
                       NULL);
}

/**
 * nc_halo_bias_castro_ref:
 * @biasf: a #NcHaloBiasCastro
 *
 * Increases the reference count of @biasf by one.
 *
 * Returns: (transfer full): @biasf
 */
NcHaloBiasCastro *
nc_halo_bias_castro_ref (NcHaloBiasCastro *biasf)
{
  return g_object_ref (biasf);
}

/**
 * nc_halo_bias_castro_free:
 * @biasf: a #NcHaloBiasCastro
 *
 * Atomically decrements the reference count of @biasf by one. If the reference count
 * drops to 0, all memory allocated by @biasf is released.
 *
 */
void
nc_halo_bias_castro_free (NcHaloBiasCastro *biasf)
{
  g_object_unref (biasf);
}

/**
 * nc_halo_bias_castro_clear:
 * @biasf: a #NcHaloBiasCastro
 *
 * Atomically decrements the reference count of @biasf by one. If the reference count
 * drops to 0, all memory allocated by @biasf is released. Set the pointer to NULL.
 *
 */
void
nc_halo_bias_castro_clear (NcHaloBiasCastro **biasf)
{
  g_clear_object (biasf);
}

/**
 * nc_halo_bias_castro_set_A0:
 * @biasf: a #NcHaloBiasCastro
 * @A0: the correction amplitude $A_0$
 *
 * Sets the correction amplitude.
 *
 */
void
nc_halo_bias_castro_set_A0 (NcHaloBiasCastro *biasf, gdouble A0)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  self->A0 = A0;
}

/**
 * nc_halo_bias_castro_get_A0:
 * @biasf: a #NcHaloBiasCastro
 *
 * Gets the correction amplitude $A_0$.
 *
 * Returns: the correction amplitude $A_0$.
 */
gdouble
nc_halo_bias_castro_get_A0 (NcHaloBiasCastro *biasf)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  return self->A0;
}

/**
 * nc_halo_bias_castro_set_a1:
 * @biasf: a #NcHaloBiasCastro
 * @a1: the $\Omega_m(z)$ coefficient $a_1$
 *
 * Sets the $\Omega_m(z)$ coefficient.
 *
 */
void
nc_halo_bias_castro_set_a1 (NcHaloBiasCastro *biasf, gdouble a1)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  self->a1 = a1;
}

/**
 * nc_halo_bias_castro_get_a1:
 * @biasf: a #NcHaloBiasCastro
 *
 * Gets the $\Omega_m(z)$ coefficient $a_1$.
 *
 * Returns: the $\Omega_m(z)$ coefficient $a_1$.
 */
gdouble
nc_halo_bias_castro_get_a1 (NcHaloBiasCastro *biasf)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  return self->a1;
}

/**
 * nc_halo_bias_castro_set_b1:
 * @biasf: a #NcHaloBiasCastro
 * @b1: the linear slope coefficient $b_1$
 *
 * Sets the linear slope coefficient.
 *
 */
void
nc_halo_bias_castro_set_b1 (NcHaloBiasCastro *biasf, gdouble b1)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  self->b1 = b1;
}

/**
 * nc_halo_bias_castro_get_b1:
 * @biasf: a #NcHaloBiasCastro
 *
 * Gets the linear slope coefficient $b_1$.
 *
 * Returns: the linear slope coefficient $b_1$.
 */
gdouble
nc_halo_bias_castro_get_b1 (NcHaloBiasCastro *biasf)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  return self->b1;
}

/**
 * nc_halo_bias_castro_set_b2:
 * @biasf: a #NcHaloBiasCastro
 * @b2: the quadratic slope coefficient $b_2$
 *
 * Sets the quadratic slope coefficient.
 *
 */
void
nc_halo_bias_castro_set_b2 (NcHaloBiasCastro *biasf, gdouble b2)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  self->b2 = b2;
}

/**
 * nc_halo_bias_castro_get_b2:
 * @biasf: a #NcHaloBiasCastro
 *
 * Gets the quadratic slope coefficient $b_2$.
 *
 * Returns: the quadratic slope coefficient $b_2$.
 */
gdouble
nc_halo_bias_castro_get_b2 (NcHaloBiasCastro *biasf)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  return self->b2;
}

/**
 * nc_halo_bias_castro_set_c1:
 * @biasf: a #NcHaloBiasCastro
 * @c1: the $S_8$ coefficient $c_1$
 *
 * Sets the $S_8$ coefficient.
 *
 */
void
nc_halo_bias_castro_set_c1 (NcHaloBiasCastro *biasf, gdouble c1)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  self->c1 = c1;
}

/**
 * nc_halo_bias_castro_get_c1:
 * @biasf: a #NcHaloBiasCastro
 *
 * Gets the $S_8$ coefficient $c_1$.
 *
 * Returns: the $S_8$ coefficient $c_1$.
 */
gdouble
nc_halo_bias_castro_get_c1 (NcHaloBiasCastro *biasf)
{
  NcHaloBiasCastroPrivate * const self = nc_halo_bias_castro_get_instance_private (biasf);

  return self->c1;
}

