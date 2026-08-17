/***************************************************************************
 *            nc_multiplicity_func_castro.c
 *
 *  Thu Aug 14 10:00:00 2026
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
 * NcMultiplicityFuncCastro:
 *
 * Dark matter halo -- Castro multiplicity function.
 *
 * Non-universal multiplicity function calibrated on spherical-overdensity virial
 * masses. Unlike the universal fits, its shape parameters depend on the local slope
 * of the variance, $\mathrm{d}\ln\sigma_R/\mathrm{d}\ln R$, and on $\Omega_m(z)$, so
 * the object needs the filtered power spectrum in addition to $\sigma_R$ itself; see
 * nc_multiplicity_func_set_psf(). A #NcHaloMassFunction injects its own filter
 * automatically.
 *
 * Two calibrations are available through the #NcMultiplicityFuncCastro:model
 * property. They share the whole coefficient structure but are distinct fits, so
 * the dynamical dark energy one does not reduce to the other when its dark energy
 * term vanishes. Both amplitudes are exact normalisations, $\int f(\nu)\,
 * \mathrm{d}\nu/\nu = 1$.
 *
 * See the <a href="../../theory/castro_hmf.html">Castro Halo Mass Function and Bias</a>
 * theory page for the full expressions, and
 * [Castro et al. (2023)](https://arxiv.org/abs/2208.02174) and
 * [Castro et al. (2025)](https://arxiv.org/abs/2504.07608).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/halo/nc_multiplicity_func_castro.h"
#include "nc/background/nc_hicosmo_de.h"
#include "numcosmo/nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcMultiplicityFuncCastroParams
{
  gdouble a1;
  gdouble a2;
  gdouble az;
  gdouble p1;
  gdouble p2;
  gdouble q1;
  gdouble q2;
  gdouble qz;
  gdouble alphaa;
} NcMultiplicityFuncCastroParams;

/* Castro et al. (2023), table 2. */
static const NcMultiplicityFuncCastroParams _nc_multiplicity_func_castro_c23[] = {
  /* ROCKSTAR     */
  { 0.7962, 0.1449, -0.0658, -0.5612, -0.4743, 0.3688, -0.2804, 0.0251, 0.0 },
  /* AHF          */ { 0.7937, 0.1119, -0.0693, -0.5689, -0.4522, 0.3652, -0.2628, 0.0376, 0.0 },
  /* SUBFIND      */ { 0.7953, 0.1667, -0.0642, -0.6265, -0.4907, 0.3215, -0.2993, 0.0330, 0.0 },
  /* VELOCIraptor */ { 0.7987, 0.1227, -0.0523, -0.5912, -0.4088, 0.3634, -0.2732, 0.0715, 0.0 },
};

/* Castro et al. (2025). q2 is not free here, it is derived from p2 below. */
static const NcMultiplicityFuncCastroParams _nc_multiplicity_func_castro_c25 = {
  0.8501, 0.234, -0.0640, -0.9880, -0.49, 0.559, 0.0, 0.02701, -0.178
};

typedef struct _NcMultiplicityFuncCastroPrivate
{
  NcMultiplicityFuncMassDef mdef;
  NcMultiplicityFuncCastroModel model;
  NcMultiplicityFuncCastroHaloFinder halo_finder;
  gdouble Delta;

  gdouble (*lnf) (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z);
} NcMultiplicityFuncCastroPrivate;

enum
{
  PROP_0,
  PROP_MODEL,
  PROP_HALO_FINDER,
  PROP_SIZE,
};

struct _NcMultiplicityFuncCastro
{
  NcMultiplicityFunc parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncCastro, nc_multiplicity_func_castro, NC_TYPE_MULTIPLICITY_FUNC)

static void
nc_multiplicity_func_castro_init (NcMultiplicityFuncCastro *mc)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  self->mdef        = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->model       = NC_MULTIPLICITY_FUNC_CASTRO_MODEL_LEN;
  self->halo_finder = NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_LEN;
  self->Delta       = 0.0;
  self->lnf         = NULL;
}

static void
_nc_multiplicity_func_castro_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncCastro *mc = NC_MULTIPLICITY_FUNC_CASTRO (object);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_CASTRO (object));

  switch (prop_id)
  {
    case PROP_MODEL:
      nc_multiplicity_func_castro_set_model (mc, g_value_get_enum (value));
      break;
    case PROP_HALO_FINDER:
      nc_multiplicity_func_castro_set_halo_finder (mc, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_multiplicity_func_castro_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncCastro *mc = NC_MULTIPLICITY_FUNC_CASTRO (object);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_CASTRO (object));

  switch (prop_id)
  {
    case PROP_MODEL:
      g_value_set_enum (value, nc_multiplicity_func_castro_get_model (mc));
      break;
    case PROP_HALO_FINDER:
      g_value_set_enum (value, nc_multiplicity_func_castro_get_halo_finder (mc));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_multiplicity_func_castro_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_castro_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_castro_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static NcMultiplicityFuncMassDef _nc_multiplicity_func_castro_get_mdef (NcMultiplicityFunc *mulf);
static void _nc_multiplicity_func_castro_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static gdouble _nc_multiplicity_func_castro_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_castro_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble lnR, gdouble z);
static gdouble _nc_multiplicity_func_castro_c23_lnf (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z);
static gdouble _nc_multiplicity_func_castro_c25_lnf (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z);

static void
nc_multiplicity_func_castro_class_init (NcMultiplicityFuncCastroClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass *parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_castro_set_property;
  object_class->get_property = &_nc_multiplicity_func_castro_get_property;
  object_class->finalize     = &_nc_multiplicity_func_castro_finalize;

  /**
   * NcMultiplicityFuncCastro:model:
   *
   * Calibration to evaluate, see #NcMultiplicityFuncCastroModel.
   */
  g_object_class_install_property (object_class,
                                   PROP_MODEL,
                                   g_param_spec_enum ("model",
                                                      NULL,
                                                      "Castro calibration",
                                                      NC_TYPE_MULTIPLICITY_FUNC_CASTRO_MODEL,
                                                      NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncCastro:halo-finder:
   *
   * Halo finder whose best-fit parameters are used. Only
   * #NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23 has per-finder calibrations; the
   * Castro et al. (2025) model provides a single parameter set and ignores this
   * property. See #NcMultiplicityFuncCastroHaloFinder.
   */
  g_object_class_install_property (object_class,
                                   PROP_HALO_FINDER,
                                   g_param_spec_enum ("halo-finder",
                                                      NULL,
                                                      "Halo finder calibration",
                                                      NC_TYPE_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER,
                                                      NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_ROCKSTAR,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->set_mdef  = &_nc_multiplicity_func_castro_set_mdef;
  parent_class->get_mdef  = &_nc_multiplicity_func_castro_get_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_castro_set_Delta;
  parent_class->get_Delta = &_nc_multiplicity_func_castro_get_Delta;
  parent_class->eval      = &_nc_multiplicity_func_castro_eval;
}

static void
_nc_multiplicity_func_castro_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      /* nothing to do, Castro is calibrated on virial masses */
      break;
    /* LCOV_EXCL_START */
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      g_error ("NcMultiplicityFuncCastro does not support mean mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncCastro does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncCastro does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
      /* LCOV_EXCL_STOP */
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef
_nc_multiplicity_func_castro_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return self->mdef;
}

static void
_nc_multiplicity_func_castro_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  self->Delta = Delta;
}

static gdouble
_nc_multiplicity_func_castro_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return self->Delta;
}

/*
 * Local slope dln(sigma_R)/dln(R), the source of the non-universality. It comes
 * from the filter rather than from sigma alone, which is why the object needs one.
 */
static gdouble
_nc_multiplicity_func_castro_dlnsigma_dlnR (NcMultiplicityFunc *mulf, gdouble lnR, gdouble z)
{
  NcmPowspecFilter *psf = nc_multiplicity_func_peek_psf (mulf);

  if (psf == NULL)
    g_error ("NcMultiplicityFuncCastro is non-universal and needs the filtered power "
             "spectrum: set it with nc_multiplicity_func_set_psf (), or use the object "
             "inside a NcHaloMassFunction, which injects its own.");

  return 0.5 * ncm_powspec_filter_eval_dlnvar_dlnr (psf, z, lnR);
}

static gdouble
_nc_multiplicity_func_castro_Omega_m (NcHICosmo *cosmo, gdouble z)
{
  return nc_hicosmo_E2Omega_m (cosmo, z) / nc_hicosmo_E2 (cosmo, z);
}

/*
 * log (1 + e^u) and log (e^x + e^y), both without ever forming the exponentials
 * themselves. The multiplicity function underflows to zero for large peak heights,
 * but its logarithm stays perfectly finite there, and it is the logarithm that the
 * peak-background split bias differentiates.
 */
static gdouble
_nc_multiplicity_func_castro_log1p_exp (const gdouble u)
{
  if (u > 0.0)
    return u + log1p (exp (-u));

  return log1p (exp (u));
}

static gdouble
_nc_multiplicity_func_castro_logaddexp (const gdouble x, const gdouble y)
{
  const gdouble m = GSL_MAX (x, y);

  if (!gsl_finite (m))
    return m;

  return m + log (exp (x - m) + exp (y - m));
}

static const NcMultiplicityFuncCastroParams *
_nc_multiplicity_func_castro_peek_params (NcMultiplicityFuncCastroPrivate * const self)
{
  if (self->model == NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C25)
    return &_nc_multiplicity_func_castro_c25;

  g_assert_cmpuint (self->halo_finder, <, NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_LEN);

  return &_nc_multiplicity_func_castro_c23[self->halo_finder];
}

static gdouble
_nc_multiplicity_func_castro_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble lnR, gdouble z)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return exp (self->lnf (mulf, cosmo, sigma,
                         _nc_multiplicity_func_castro_dlnsigma_dlnR (mulf, lnR, z), z));
}

static gdouble
_nc_multiplicity_func_castro_c23_lnf (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);
  const NcMultiplicityFuncCastroParams *bf     = _nc_multiplicity_func_castro_peek_params (self);

  const gdouble s    = dlnsigma_dlnR;
  const gdouble Om_z = _nc_multiplicity_func_castro_Omega_m (cosmo, z);
  const gdouble nu   = nc_multiplicity_func_castro_delta_c (mc, cosmo, z) / sigma;

  const gdouble aR = bf->a1 + bf->a2 * gsl_pow_2 (s + 0.6125);
  const gdouble qR = bf->q1 + bf->q2 * (s + 0.5);
  const gdouble a  = aR * pow (Om_z, bf->az);
  const gdouble p  = bf->p1 + bf->p2 * (s + 0.5);
  const gdouble q  = qR * pow (Om_z, bf->qz);

  const gdouble A_inv = pow (2.0, -0.5 - p + 0.5 * q) / sqrt (M_PI) *
                        (pow (2.0, p) * tgamma (0.5 * q) + tgamma (0.5 * q - p));

  const gdouble anu2   = a * nu * nu;
  const gdouble lnanu2 = log (anu2);

  return -log (A_inv)
         + 0.5 * (M_LN2 + lnanu2 - log (M_PI))
         - 0.5 * anu2
         + _nc_multiplicity_func_castro_log1p_exp (-p * lnanu2)
         + (q - 1.0) * 0.5 * lnanu2;
}

static gdouble
_nc_multiplicity_func_castro_c25_lnf (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z)
{
  NcMultiplicityFuncCastro *mc                 = NC_MULTIPLICITY_FUNC_CASTRO (mulf);
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);
  const NcMultiplicityFuncCastroParams *bf     = _nc_multiplicity_func_castro_peek_params (self);

  if (!NC_IS_HICOSMO_DE (cosmo))
    g_error ("NcMultiplicityFuncCastro: the %s calibration needs the dark energy density "
             "and equation of state at turn-around, which requires a NcHICosmoDE; got %s. "
             "Use NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23 instead.",
             "Castro et al. (2025)", G_OBJECT_TYPE_NAME (cosmo));

  {
    NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);

    const gdouble s    = dlnsigma_dlnR;
    const gdouble Om_z = _nc_multiplicity_func_castro_Omega_m (cosmo, z);
    const gdouble nu   = nc_multiplicity_func_castro_delta_c (mc, cosmo, z) / sigma;

    const gdouble z_ta      = nc_multiplicity_func_castro_z_ta (mc, z);
    const gdouble Om_de_zta = nc_hicosmo_de_E2Omega_de (cosmo_de, z_ta) / nc_hicosmo_E2 (cosmo, z_ta);
    const gdouble w_de_zta  = nc_hicosmo_de_w_de (cosmo_de, z_ta);

    /* Coefficients tied to the free ones by the calibration. */
    const gdouble q2     = 1.0020 * bf->p2 + 0.7363;
    const gdouble r1     = -0.9917 * bf->q1 + 0.5560;
    const gdouble r2     = -1.1159 * q2 - 0.1189;
    const gdouble vp     = 1.0667 * bf->q1 + 1.7577;
    const gdouble alphap = 4.7686 * bf->alphaa + 0.0051;

    const gdouble aR = bf->a1 + bf->a2 * gsl_pow_2 (s + 0.6125);
    const gdouble qR = bf->q1 + q2 * (s + 0.5);
    const gdouble q  = qR * pow (Om_z, bf->qz);
    const gdouble r  = r1 + r2 * (s + 0.5);

    const gdouble qdec = 1.5 * (w_de_zta + 1.0) * Om_de_zta;
    const gdouble a    = aR * pow (Om_z, bf->az) * (1.0 + bf->alphaa * qdec);
    const gdouble p    = (bf->p1 + bf->p2 * (s + 0.5)) * (1.0 + alphap * qdec);

    const gdouble g_qr = tgamma (0.5 * q + 0.5 * r + 1.0);
    const gdouble g_pq = tgamma (-p + 0.5 * q + 1.0);
    const gdouble vp2p = pow (vp, 2.0 * p);
    const gdouble num  = pow (2.0, -p + 0.5 * q + 0.5) *
                         (-pow (2.0, p + 0.5 * r) * q * g_qr
                          + pow (2.0, p + 0.5 * r + 1.0) * p * g_qr
                          - q * vp2p * g_pq
                          - r * vp2p * g_pq);
    const gdouble den   = sqrt (M_PI) * (2.0 * p - q) * (q + r);
    const gdouble A_inv = num / den;

    const gdouble lnvst = log (nu) + 0.5 * log (a);
    const gdouble vst   = exp (lnvst);

    /* log (vst^r + (vp / vst)^{2p}), each term taken in the exponent. */
    const gdouble lnshape = _nc_multiplicity_func_castro_logaddexp (r * lnvst,
                                                                    2.0 * p * (log (vp) - lnvst));

    return M_LN2
           - log (A_inv)
           + lnshape
           + lnvst - 0.5 * log (2.0 * M_PI)
           - 0.5 * vst * vst
           + (q - 1.0) * lnvst;
  }
}

/**
 * nc_multiplicity_func_castro_eval_full:
 * @mc: a #NcMultiplicityFuncCastro
 * @cosmo: a #NcHICosmo
 * @sigma: standard deviation of the matter density contrast $\sigma_R$
 * @dlnsigma_dlnR: the slope $\mathrm{d}\ln\sigma_R/\mathrm{d}\ln R$
 * @z: redshift $z$
 *
 * Evaluates the multiplicity function with the slope supplied directly instead of
 * taken from the filter. The model is a closed-form function of @sigma,
 * @dlnsigma_dlnR and the background at @z, so this is the whole of it: no filter is
 * consulted and no interpolation is involved.
 *
 * This is what makes the two partial derivatives $\partial\ln f/\partial\ln\nu$ and
 * $\partial\ln f/\partial s$ cheap and accurate to obtain, since each can be taken by
 * perturbing a plain number rather than differentiating through an interpolated
 * $\sigma_R$. #NcHaloBiasCastro relies on that.
 *
 * Returns: the value of the multiplicity function
 */
gdouble
nc_multiplicity_func_castro_eval_full (NcMultiplicityFuncCastro *mc, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return exp (self->lnf (NC_MULTIPLICITY_FUNC (mc), cosmo, sigma, dlnsigma_dlnR, z));
}

/**
 * nc_multiplicity_func_castro_eval_lnf:
 * @mc: a #NcMultiplicityFuncCastro
 * @cosmo: a #NcHICosmo
 * @sigma: standard deviation of the matter density contrast $\sigma_R$
 * @dlnsigma_dlnR: the slope $\mathrm{d}\ln\sigma_R/\mathrm{d}\ln R$
 * @z: redshift $z$
 *
 * Logarithm of the multiplicity function, computed without ever forming the
 * function itself.
 *
 * This is the numerically meaningful quantity in the large peak-height regime:
 * $f$ carries a factor $e^{-a\nu^2/2}$ and underflows to zero once $\nu$ grows
 * past a few hundred, while $\ln f$ stays finite and smooth. Anything that
 * differentiates $\ln f$ -- the peak-background split bias in particular -- should
 * use this rather than take the logarithm of nc_multiplicity_func_castro_eval_full().
 *
 * Returns: $\ln f$
 */
gdouble
nc_multiplicity_func_castro_eval_lnf (NcMultiplicityFuncCastro *mc, NcHICosmo *cosmo, gdouble sigma, gdouble dlnsigma_dlnR, gdouble z)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return self->lnf (NC_MULTIPLICITY_FUNC (mc), cosmo, sigma, dlnsigma_dlnR, z);
}

/**
 * nc_multiplicity_func_castro_delta_c:
 * @mc: a #NcMultiplicityFuncCastro
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Linear collapse threshold $\delta_c(\Omega_m(z))$ following Kitayama & Suto (1996),
 * as used by the Castro calibrations.
 *
 * Returns: $\delta_c$.
 */
gdouble
nc_multiplicity_func_castro_delta_c (NcMultiplicityFuncCastro *mc, NcHICosmo *cosmo, gdouble z)
{
  const gdouble Om_z = _nc_multiplicity_func_castro_Omega_m (cosmo, z);

  return 3.0 / 20.0 * pow (12.0 * M_PI, 2.0 / 3.0) * (1.0 + 0.012299 * log10 (Om_z));
}

/**
 * nc_multiplicity_func_castro_z_ta:
 * @mc: a #NcMultiplicityFuncCastro
 * @z: collapse redshift $z$
 *
 * Turn-around redshift for a halo collapsing at @z, in the Einstein-de Sitter
 * approximation used by the calibration, $1 + z_\mathrm{ta} = (1+z)\,2^{2/3}$.
 *
 * Returns: $z_\mathrm{ta}$.
 */
gdouble
nc_multiplicity_func_castro_z_ta (NcMultiplicityFuncCastro *mc, gdouble z)
{
  return (1.0 + z) * pow (2.0, 2.0 / 3.0) - 1.0;
}

/**
 * nc_multiplicity_func_castro_set_model:
 * @mc: a #NcMultiplicityFuncCastro
 * @model: a #NcMultiplicityFuncCastroModel
 *
 * Sets the calibration to evaluate.
 *
 */
void
nc_multiplicity_func_castro_set_model (NcMultiplicityFuncCastro *mc, NcMultiplicityFuncCastroModel model)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  switch (model)
  {
    case NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C23:
      self->lnf = &_nc_multiplicity_func_castro_c23_lnf;
      break;
    case NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C25:
      self->lnf = &_nc_multiplicity_func_castro_c25_lnf;
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  self->model = model;
}

/**
 * nc_multiplicity_func_castro_get_model:
 * @mc: a #NcMultiplicityFuncCastro
 *
 * Gets the calibration in use.
 *
 * Returns: the #NcMultiplicityFuncCastroModel.
 */
NcMultiplicityFuncCastroModel
nc_multiplicity_func_castro_get_model (NcMultiplicityFuncCastro *mc)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return self->model;
}

/**
 * nc_multiplicity_func_castro_set_halo_finder:
 * @mc: a #NcMultiplicityFuncCastro
 * @halo_finder: a #NcMultiplicityFuncCastroHaloFinder
 *
 * Sets the halo finder whose best-fit parameters are used. Ignored by
 * #NC_MULTIPLICITY_FUNC_CASTRO_MODEL_C25, which has a single parameter set.
 *
 */
void
nc_multiplicity_func_castro_set_halo_finder (NcMultiplicityFuncCastro *mc, NcMultiplicityFuncCastroHaloFinder halo_finder)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  g_assert_cmpuint (halo_finder, <, NC_MULTIPLICITY_FUNC_CASTRO_HALO_FINDER_LEN);

  self->halo_finder = halo_finder;
}

/**
 * nc_multiplicity_func_castro_get_halo_finder:
 * @mc: a #NcMultiplicityFuncCastro
 *
 * Gets the halo finder in use.
 *
 * Returns: the #NcMultiplicityFuncCastroHaloFinder.
 */
NcMultiplicityFuncCastroHaloFinder
nc_multiplicity_func_castro_get_halo_finder (NcMultiplicityFuncCastro *mc)
{
  NcMultiplicityFuncCastroPrivate * const self = nc_multiplicity_func_castro_get_instance_private (mc);

  return self->halo_finder;
}

/**
 * nc_multiplicity_func_castro_new:
 *
 * Creates a new #NcMultiplicityFuncCastro with the Castro et al. (2023)
 * calibration and ROCKSTAR parameters.
 *
 * Returns: A new #NcMultiplicityFuncCastro.
 */
NcMultiplicityFuncCastro *
nc_multiplicity_func_castro_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_CASTRO,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL,
                       NULL);
}

/**
 * nc_multiplicity_func_castro_new_full:
 * @model: a #NcMultiplicityFuncCastroModel
 * @halo_finder: a #NcMultiplicityFuncCastroHaloFinder
 *
 * Creates a new #NcMultiplicityFuncCastro using @model and @halo_finder.
 *
 * Returns: A new #NcMultiplicityFuncCastro.
 */
NcMultiplicityFuncCastro *
nc_multiplicity_func_castro_new_full (NcMultiplicityFuncCastroModel model, NcMultiplicityFuncCastroHaloFinder halo_finder)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_CASTRO,
                       "mass-def",    NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL,
                       "model",       model,
                       "halo-finder", halo_finder,
                       NULL);
}

/**
 * nc_multiplicity_func_castro_ref:
 * @mc: a #NcMultiplicityFuncCastro
 *
 * Increases the reference count of @mc by one.
 *
 * Returns: (transfer full): @mc
 */
NcMultiplicityFuncCastro *
nc_multiplicity_func_castro_ref (NcMultiplicityFuncCastro *mc)
{
  return g_object_ref (mc);
}

/**
 * nc_multiplicity_func_castro_free:
 * @mc: a #NcMultiplicityFuncCastro
 *
 * Atomically decrements the reference count of @mc by one. If the reference count
 * drops to 0, all memory allocated by @mc is released.
 *
 */
void
nc_multiplicity_func_castro_free (NcMultiplicityFuncCastro *mc)
{
  g_object_unref (mc);
}

/**
 * nc_multiplicity_func_castro_clear:
 * @mc: a #NcMultiplicityFuncCastro
 *
 * Atomically decrements the reference count of @mc by one. If the reference count
 * drops to 0, all memory allocated by @mc is released. Set the pointer to NULL.
 *
 */
void
nc_multiplicity_func_castro_clear (NcMultiplicityFuncCastro **mc)
{
  g_clear_object (mc);
}

