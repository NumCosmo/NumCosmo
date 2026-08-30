/***************************************************************************
 *            ncm_powspec_analytic.c
 *
 *  Tue August 26 12:00:00 2026
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
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmPowspecAnalytic:
 *
 * Closed-form power spectrum, evaluated directly.
 *
 * \begin{equation}
 *   P(k, z) = A \, k^{n_s} \, T(k)^2 \, B(k) \, D(z)^2 .
 * \end{equation}
 *
 * Every other #NcmPowspec in the library is spline-backed: even
 * #NcPowspecMLTransfer, whose transfer function is a closed form, returns
 * `ncm_spline_eval (Pk, log (k))`. This one has no spline anywhere in
 * ncm_powspec_eval(), so it is exactly reproducible in arbitrary precision and
 * can serve as the known value a quadrature is measured against.
 *
 * It also carries no cosmology. The #NcmModel handed to ncm_powspec_eval() is
 * ignored; $\Omega_m$ is a parameter of the growth factor and nothing else.
 * That is deliberate -- the $\sigma_8$ normalization of a physical power
 * spectrum is a filtered-variance quadrature, the one link in the chain with
 * no closed form.
 *
 * The default parameters are a least-squares fit to NumCosmo's own $P(k)$
 * ($h = 0.6736$, $\Omega_m = 0.3$, Eisenstein-Hu) over
 * $k \in [10^{-5}, 10]\,\mathrm{Mpc}^{-1}$: with the BBKS shape they reproduce
 * it to 12% at worst and 4% rms, with the turnover at
 * $0.0104\,\mathrm{Mpc}^{-1}$ against $0.0109$ and a logarithmic slope of
 * $-2.635$ at $k = 10$ against $-2.622$. It is a look-alike, not a cosmology.
 *
 * Units follow #NcmPowspec: $k$ in $\mathrm{Mpc}^{-1}$ and $P$ in
 * $\mathrm{Mpc}^3$.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "ncm/powspec/tests/ncm_powspec_analytic.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_trig.h>
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmPowspecAnalytic
{
  /*< private >*/
  NcmPowspec parent_instance;

  NcmPowspecAnalyticShape shape;
  NcmPowspecAnalyticGrowth growth;

  gdouble amplitude;
  gdouble n_s;
  gdouble k_eq;
  gdouble Omega_m;
  gdouble a_t;

  gdouble bao_amplitude;
  gdouble bao_rd;
  gdouble bao_sigma;

  /* D(a = 1) of the unnormalized growth, so eval_growth() returns D(0) = 1. */
  gdouble growth_norm;
};

enum
{
  PROP_0,
  PROP_SHAPE,
  PROP_GROWTH,
  PROP_AMPLITUDE,
  PROP_N_S,
  PROP_K_EQ,
  PROP_OMEGA_M,
  PROP_A_T,
  PROP_BAO_AMPLITUDE,
  PROP_BAO_RD,
  PROP_BAO_SIGMA,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmPowspecAnalytic, ncm_powspec_analytic, NCM_TYPE_POWSPEC)

/*
 * The k range the two constructors advertise through NcmPowspec:kmin/kmax.
 *
 * The closed form is valid at every positive k -- there is no table to run off
 * the end of -- so this is only what the object reports to callers that ask.
 * NcmPowspec's own defaults (1e-5 to 1 Mpc^-1) are too narrow for that to be
 * useful here: the xcor kernels reach about 5 Mpc^-1 at l = 200, and a caller
 * reading kmax off this object would silently clip its own integration range.
 */
#define NCM_POWSPEC_ANALYTIC_DEFAULT_KMIN (1.0e-6)
#define NCM_POWSPEC_ANALYTIC_DEFAULT_KMAX (1.0e4)

/* BBKS transfer-function coefficients, Bardeen et al. (1986) eq. (G3). */
#define NCM_POWSPEC_ANALYTIC_BBKS_C1 (2.34)
#define NCM_POWSPEC_ANALYTIC_BBKS_C2 (3.89)
#define NCM_POWSPEC_ANALYTIC_BBKS_C3 (16.1)
#define NCM_POWSPEC_ANALYTIC_BBKS_C4 (5.46)
#define NCM_POWSPEC_ANALYTIC_BBKS_C5 (6.71)

static void
ncm_powspec_analytic_init (NcmPowspecAnalytic *psa)
{
  psa->shape         = NCM_POWSPEC_ANALYTIC_SHAPE_LEN;
  psa->growth        = NCM_POWSPEC_ANALYTIC_GROWTH_LEN;
  psa->amplitude     = 0.0;
  psa->n_s           = 0.0;
  psa->k_eq          = 0.0;
  psa->Omega_m       = 0.0;
  psa->a_t           = 0.0;
  psa->bao_amplitude = 0.0;
  psa->bao_rd        = 0.0;
  psa->bao_sigma     = 0.0;
  psa->growth_norm   = 1.0;
}

static void
_ncm_powspec_analytic_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (object);

  g_return_if_fail (NCM_IS_POWSPEC_ANALYTIC (object));

  switch (prop_id)
  {
    case PROP_SHAPE:
      psa->shape = g_value_get_enum (value);
      break;
    case PROP_GROWTH:
      psa->growth = g_value_get_enum (value);
      break;
    case PROP_AMPLITUDE:
      psa->amplitude = g_value_get_double (value);
      break;
    case PROP_N_S:
      psa->n_s = g_value_get_double (value);
      break;
    case PROP_K_EQ:
      psa->k_eq = g_value_get_double (value);
      break;
    case PROP_OMEGA_M:
      psa->Omega_m = g_value_get_double (value);
      break;
    case PROP_A_T:
      psa->a_t = g_value_get_double (value);
      break;
    case PROP_BAO_AMPLITUDE:
      psa->bao_amplitude = g_value_get_double (value);
      break;
    case PROP_BAO_RD:
      psa->bao_rd = g_value_get_double (value);
      break;
    case PROP_BAO_SIGMA:
      psa->bao_sigma = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_analytic_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (object);

  g_return_if_fail (NCM_IS_POWSPEC_ANALYTIC (object));

  switch (prop_id)
  {
    case PROP_SHAPE:
      g_value_set_enum (value, psa->shape);
      break;
    case PROP_GROWTH:
      g_value_set_enum (value, psa->growth);
      break;
    case PROP_AMPLITUDE:
      g_value_set_double (value, psa->amplitude);
      break;
    case PROP_N_S:
      g_value_set_double (value, psa->n_s);
      break;
    case PROP_K_EQ:
      g_value_set_double (value, psa->k_eq);
      break;
    case PROP_OMEGA_M:
      g_value_set_double (value, psa->Omega_m);
      break;
    case PROP_A_T:
      g_value_set_double (value, psa->a_t);
      break;
    case PROP_BAO_AMPLITUDE:
      g_value_set_double (value, psa->bao_amplitude);
      break;
    case PROP_BAO_RD:
      g_value_set_double (value, psa->bao_rd);
      break;
    case PROP_BAO_SIGMA:
      g_value_set_double (value, psa->bao_sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Unnormalized flat matter-plus-Lambda linear growth,
 *
 *   D(a) = a 2F1(1/3, 1; 11/6; -(1 - Om)/Om a^3) .
 *
 * The argument is <= 0 and reaches -9 already at Omega_m = 0.1, well outside
 * the |x| < 1 that gsl_sf_hyperg_2F1() accepts, so it is reached through
 * Pfaff's transformation
 *
 *   2F1(a, b; c; z) = (1 - z)^-a 2F1(a, c - b; c; z/(z - 1)) ,
 *
 * whose argument z/(z-1) lies in [0, 1) for every z <= 0. The transformed
 * series converges at the endpoint as well, since c - a - b = 2/3 > 0.
 */
static gdouble
_ncm_powspec_analytic_growth_lcdm_raw (NcmPowspecAnalytic *psa, const gdouble a)
{
  const gdouble x = -(1.0 - psa->Omega_m) / psa->Omega_m * gsl_pow_3 (a);
  const gdouble w = x / (x - 1.0);

  return a * pow (1.0 - x, -1.0 / 3.0) *
         gsl_sf_hyperg_2F1 (1.0 / 3.0, 5.0 / 6.0, 11.0 / 6.0, w);
}

/*
 * d/da of the above. The derivative of the hypergeometric raises every
 * parameter by one,
 *
 *   d/dx 2F1(a, b; c; x) = (ab/c) 2F1(a+1, b+1; c+1; x) ,
 *
 * so with x = -L a^3, D'(a) = F(x) - 3 L a^3 F'(x), and the raised function is
 * reached by the same Pfaff transformation (c - a - b = 2/3 > 0 there too).
 */
static gdouble
_ncm_powspec_analytic_growth_lcdm_raw_deriv (NcmPowspecAnalytic *psa, const gdouble a)
{
  const gdouble L = (1.0 - psa->Omega_m) / psa->Omega_m;
  const gdouble x = -L *gsl_pow_3 (a);

  const gdouble w = x / (x - 1.0);
  const gdouble F = pow (1.0 - x, -1.0 / 3.0) *
                    gsl_sf_hyperg_2F1 (1.0 / 3.0, 5.0 / 6.0, 11.0 / 6.0, w);
  const gdouble Fp = (2.0 / 11.0) * pow (1.0 - x, -4.0 / 3.0) *
                     gsl_sf_hyperg_2F1 (4.0 / 3.0, 5.0 / 6.0, 17.0 / 6.0, w);

  return F - 3.0 * L * gsl_pow_3 (a) * Fp;
}

static gdouble
_ncm_powspec_analytic_growth_rational_raw (NcmPowspecAnalytic *psa, const gdouble a)
{
  return a * pow (1.0 + gsl_pow_3 (a / psa->a_t), -1.0 / 3.0);
}

static gdouble
_ncm_powspec_analytic_growth_rational_raw_deriv (NcmPowspecAnalytic *psa, const gdouble a)
{
  const gdouble u = gsl_pow_3 (a / psa->a_t);

  return pow (1.0 + u, -4.0 / 3.0);
}

static void
_ncm_powspec_analytic_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_powspec_analytic_parent_class)->constructed (object);
  {
    NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (object);

    /* No validation here: every shape parameter is range-checked by its own
     * GParamSpec (k-eq, a-t and Omega-m are bounded away from zero, the two
     * enums by their GType), so an out-of-range value never reaches this far.
     * All that is left is the growth normalization, fixed once so that
     * eval_growth() returns D(0) = 1. */
    switch (psa->growth)
    {
      case NCM_POWSPEC_ANALYTIC_GROWTH_NONE:
        psa->growth_norm = 1.0;
        break;
      case NCM_POWSPEC_ANALYTIC_GROWTH_LCDM:
        psa->growth_norm = _ncm_powspec_analytic_growth_lcdm_raw (psa, 1.0);
        break;
      case NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL:
        psa->growth_norm = _ncm_powspec_analytic_growth_rational_raw (psa, 1.0);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }
  }
}

static void
_ncm_powspec_analytic_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_analytic_parent_class)->finalize (object);
}

static void _ncm_powspec_analytic_prepare (NcmPowspec *powspec, NcmModel *model);
static gdouble _ncm_powspec_analytic_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _ncm_powspec_analytic_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
static gdouble _ncm_powspec_analytic_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static gdouble _ncm_powspec_analytic_deriv_k (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _ncm_powspec_analytic_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

static void
ncm_powspec_analytic_class_init (NcmPowspecAnalyticClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_ncm_powspec_analytic_set_property;
  object_class->get_property = &_ncm_powspec_analytic_get_property;
  object_class->constructed  = &_ncm_powspec_analytic_constructed;
  object_class->finalize     = &_ncm_powspec_analytic_finalize;

  /**
   * NcmPowspecAnalytic:shape:
   *
   * Transfer-function shape, a #NcmPowspecAnalyticShape.
   */
  g_object_class_install_property (object_class,
                                   PROP_SHAPE,
                                   g_param_spec_enum ("shape",
                                                      NULL,
                                                      "Transfer function shape",
                                                      NCM_TYPE_POWSPEC_ANALYTIC_SHAPE,
                                                      NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:growth:
   *
   * Growth-factor shape, a #NcmPowspecAnalyticGrowth.
   */
  g_object_class_install_property (object_class,
                                   PROP_GROWTH,
                                   g_param_spec_enum ("growth",
                                                      NULL,
                                                      "Growth factor shape",
                                                      NCM_TYPE_POWSPEC_ANALYTIC_GROWTH,
                                                      NCM_POWSPEC_ANALYTIC_GROWTH_LCDM,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:amplitude:
   *
   * Overall amplitude $A$, in $\mathrm{Mpc}^{3+n_s}$. The default is fitted to
   * NumCosmo's own $P(k)$; see #NcmPowspecAnalytic.
   */
  g_object_class_install_property (object_class,
                                   PROP_AMPLITUDE,
                                   g_param_spec_double ("amplitude",
                                                        NULL,
                                                        "Overall amplitude",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 2.08e7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:n-s:
   *
   * Primordial spectral index $n_s$.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_S,
                                   g_param_spec_double ("n-s",
                                                        NULL,
                                                        "Spectral index",
                                                        -3.0, 3.0, 0.9875,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:k-eq:
   *
   * Transfer-function scale $k_\mathrm{eq}$, in $\mathrm{Mpc}^{-1}$. With the
   * BBKS shape the turnover of $P$ falls near $0.1 k_\mathrm{eq}$, not at
   * $k_\mathrm{eq}$ itself.
   */
  g_object_class_install_property (object_class,
                                   PROP_K_EQ,
                                   g_param_spec_double ("k-eq",
                                                        NULL,
                                                        "Transfer function scale in 1/Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 0.10594,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:Omega-m:
   *
   * Matter fraction used by %NCM_POWSPEC_ANALYTIC_GROWTH_LCDM. It is a
   * parameter of the growth factor and of nothing else -- this object carries
   * no cosmology.
   */
  g_object_class_install_property (object_class,
                                   PROP_OMEGA_M,
                                   g_param_spec_double ("Omega-m",
                                                        NULL,
                                                        "Matter fraction, growth factor only",
                                                        G_MINDOUBLE, 1.0, 0.3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:a-t:
   *
   * Transition scale factor of %NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL.
   */
  g_object_class_install_property (object_class,
                                   PROP_A_T,
                                   g_param_spec_double ("a-t",
                                                        NULL,
                                                        "Growth transition scale factor",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:bao-amplitude:
   *
   * Amplitude of the optional oscillatory factor
   * $B(k) = 1 + a\,\mathrm{sinc}(k r_d) e^{-(k\sigma_d)^2}$. Zero, the default,
   * turns it off. It is a third oscillation scale, not a physical BAO feature.
   */
  g_object_class_install_property (object_class,
                                   PROP_BAO_AMPLITUDE,
                                   g_param_spec_double ("bao-amplitude",
                                                        NULL,
                                                        "Amplitude of the oscillatory factor",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:bao-rd:
   *
   * Oscillation scale $r_d$ of the optional oscillatory factor, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_BAO_RD,
                                   g_param_spec_double ("bao-rd",
                                                        NULL,
                                                        "Oscillation scale in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 147.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspecAnalytic:bao-sigma:
   *
   * Damping scale $\sigma_d$ of the optional oscillatory factor, in Mpc.
   */
  g_object_class_install_property (object_class,
                                   PROP_BAO_SIGMA,
                                   g_param_spec_double ("bao-sigma",
                                                        NULL,
                                                        "Oscillation damping scale in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  powspec_class->prepare    = &_ncm_powspec_analytic_prepare;
  powspec_class->eval       = &_ncm_powspec_analytic_eval;
  powspec_class->eval_vec   = &_ncm_powspec_analytic_eval_vec;
  powspec_class->deriv_z    = &_ncm_powspec_analytic_deriv_z;
  powspec_class->deriv_k    = &_ncm_powspec_analytic_deriv_k;
  powspec_class->get_nknots = &_ncm_powspec_analytic_get_nknots;
}

/*
 * Nothing to precompute: every evaluation is a closed form and the growth
 * normalization is fixed at construction. The model is ignored.
 */
static void
_ncm_powspec_analytic_prepare (NcmPowspec *powspec, NcmModel *model)
{
}

/* T(k). */
static gdouble
_ncm_powspec_analytic_transfer (NcmPowspecAnalytic *psa, const gdouble k)
{
  switch (psa->shape)
  {
    case NCM_POWSPEC_ANALYTIC_SHAPE_POWER_LAW:

      return 1.0;

    case NCM_POWSPEC_ANALYTIC_SHAPE_BBKS:
    {
      const gdouble q  = k / psa->k_eq;
      const gdouble c1 = NCM_POWSPEC_ANALYTIC_BBKS_C1 * q;
      const gdouble M  = 1.0
                         + NCM_POWSPEC_ANALYTIC_BBKS_C2 * q
                         + gsl_pow_2 (NCM_POWSPEC_ANALYTIC_BBKS_C3 * q)
                         + gsl_pow_3 (NCM_POWSPEC_ANALYTIC_BBKS_C4 * q)
                         + gsl_pow_4 (NCM_POWSPEC_ANALYTIC_BBKS_C5 * q);

      /* log1p (x) / x is accurate to the last bit down to the smallest
       * normal, but x = 0 is a genuine 0/0 whose limit is one. */
      const gdouble L = (c1 == 0.0) ? 1.0 : log1p (c1) / c1;

      return L * pow (M, -0.25);
    }

    case NCM_POWSPEC_ANALYTIC_SHAPE_RATIONAL:

      return 1.0 / (1.0 + gsl_pow_2 (k / psa->k_eq));

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/* d ln T / d k. */
static gdouble
_ncm_powspec_analytic_dlnT_dk (NcmPowspecAnalytic *psa, const gdouble k)
{
  switch (psa->shape)
  {
    case NCM_POWSPEC_ANALYTIC_SHAPE_POWER_LAW:

      return 0.0;

    case NCM_POWSPEC_ANALYTIC_SHAPE_BBKS:
    {
      const gdouble q  = k / psa->k_eq;
      const gdouble c1 = NCM_POWSPEC_ANALYTIC_BBKS_C1 * q;
      const gdouble M  = 1.0
                         + NCM_POWSPEC_ANALYTIC_BBKS_C2 * q
                         + gsl_pow_2 (NCM_POWSPEC_ANALYTIC_BBKS_C3 * q)
                         + gsl_pow_3 (NCM_POWSPEC_ANALYTIC_BBKS_C4 * q)
                         + gsl_pow_4 (NCM_POWSPEC_ANALYTIC_BBKS_C5 * q);
      const gdouble Mp = NCM_POWSPEC_ANALYTIC_BBKS_C2
                         + 2.0 * gsl_pow_2 (NCM_POWSPEC_ANALYTIC_BBKS_C3) * q
                         + 3.0 * gsl_pow_3 (NCM_POWSPEC_ANALYTIC_BBKS_C4) * gsl_pow_2 (q)
                         + 4.0 * gsl_pow_4 (NCM_POWSPEC_ANALYTIC_BBKS_C5) * gsl_pow_3 (q);
      /* d/dq of ln [log1p (c1 q) / (c1 q)], with the 0/0 limit -c1/2 at q = 0. */
      const gdouble dlnL_dq = (q == 0.0) ?
                              -0.5 * NCM_POWSPEC_ANALYTIC_BBKS_C1 :
                              NCM_POWSPEC_ANALYTIC_BBKS_C1 / ((1.0 + c1) * log1p (c1)) - 1.0 / q;

      return (dlnL_dq - 0.25 * Mp / M) / psa->k_eq;
    }

    case NCM_POWSPEC_ANALYTIC_SHAPE_RATIONAL:
    {
      const gdouble u = gsl_pow_2 (k / psa->k_eq);

      return -2.0 * k / (gsl_pow_2 (psa->k_eq) * (1.0 + u));
    }

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/* B(k), one when the oscillatory factor is off. */
static gdouble
_ncm_powspec_analytic_bao (NcmPowspecAnalytic *psa, const gdouble k)
{
  const gdouble x = k * psa->bao_rd;

  if (psa->bao_amplitude == 0.0)
    return 1.0;

  return 1.0 + psa->bao_amplitude * gsl_sf_sinc (x / M_PI) *
         exp (-gsl_pow_2 (k * psa->bao_sigma));
}

/* dB/dk. */
static gdouble
_ncm_powspec_analytic_dbao_dk (NcmPowspecAnalytic *psa, const gdouble k)
{
  const gdouble r = psa->bao_rd;
  const gdouble x = k * r;
  const gdouble d = exp (-gsl_pow_2 (k * psa->bao_sigma));
  gdouble sinc, dsinc_dx;

  if (psa->bao_amplitude == 0.0)
    return 0.0;

  sinc = gsl_sf_sinc (x / M_PI);
  /* d/dx [sin x / x] = (x cos x - sin x)/x^2, with limit 0 at x = 0. */
  dsinc_dx = (x == 0.0) ? 0.0 : (cos (x) - sinc) / x;

  return psa->bao_amplitude * d *
         (dsinc_dx * r - sinc * 2.0 * k * gsl_pow_2 (psa->bao_sigma));
}

static gdouble
_ncm_powspec_analytic_growth (NcmPowspecAnalytic *psa, const gdouble z)
{
  const gdouble a = 1.0 / (1.0 + z);

  /* D is normalized to D(0) = 1, so z = 0 is exactly one and needs no special
   * function. Worth the branch: every Tier A evaluation is at z = 0 by
   * construction (the analytic xcor kernels fold all evolution into the window
   * and call ncm_powspec_eval() with z = 0.0), and the LCDM branch below costs
   * a gsl_sf_hyperg_2F1 per call -- measured at 29% of an xcor solve driven by
   * this object, all of it spent returning 1.0. */
  if (z == 0.0)
    return 1.0;

  switch (psa->growth)
  {
    case NCM_POWSPEC_ANALYTIC_GROWTH_NONE:

      return 1.0;

    case NCM_POWSPEC_ANALYTIC_GROWTH_LCDM:

      return _ncm_powspec_analytic_growth_lcdm_raw (psa, a) / psa->growth_norm;

    case NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL:

      return _ncm_powspec_analytic_growth_rational_raw (psa, a) / psa->growth_norm;

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/* dD/dz = (dD/da) da/dz, with da/dz = -a^2. */
static gdouble
_ncm_powspec_analytic_dgrowth_dz (NcmPowspecAnalytic *psa, const gdouble z)
{
  const gdouble a = 1.0 / (1.0 + z);
  gdouble dD_da;

  switch (psa->growth)
  {
    case NCM_POWSPEC_ANALYTIC_GROWTH_NONE:

      return 0.0;

    case NCM_POWSPEC_ANALYTIC_GROWTH_LCDM:
      dD_da = _ncm_powspec_analytic_growth_lcdm_raw_deriv (psa, a);
      break;

    case NCM_POWSPEC_ANALYTIC_GROWTH_RATIONAL:
      dD_da = _ncm_powspec_analytic_growth_rational_raw_deriv (psa, a);
      break;

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }

  return -dD_da *gsl_pow_2 (a) / psa->growth_norm;
}

static gdouble
_ncm_powspec_analytic_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (powspec);
  const gdouble T         = _ncm_powspec_analytic_transfer (psa, k);
  const gdouble D         = _ncm_powspec_analytic_growth (psa, z);

  return psa->amplitude * pow (k, psa->n_s) * gsl_pow_2 (T) *
         _ncm_powspec_analytic_bao (psa, k) * gsl_pow_2 (D);
}

static void
_ncm_powspec_analytic_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (powspec);
  const gdouble D2        = gsl_pow_2 (_ncm_powspec_analytic_growth (psa, z));
  const guint len         = ncm_vector_len (k);
  guint i;

  g_assert_cmpuint (len, ==, ncm_vector_len (Pk));

  for (i = 0; i < len; i++)
  {
    const gdouble ki = ncm_vector_get (k, i);
    const gdouble T  = _ncm_powspec_analytic_transfer (psa, ki);

    ncm_vector_set (Pk, i, psa->amplitude * pow (ki, psa->n_s) * gsl_pow_2 (T) *
                    _ncm_powspec_analytic_bao (psa, ki) * D2);
  }
}

/* dP/dz = 2 P (dD/dz) / D. */
static gdouble
_ncm_powspec_analytic_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (powspec);
  const gdouble D         = _ncm_powspec_analytic_growth (psa, z);
  const gdouble P         = _ncm_powspec_analytic_eval (powspec, model, z, k);

  return 2.0 * P * _ncm_powspec_analytic_dgrowth_dz (psa, z) / D;
}

/* dP/dk = P [n_s/k + 2 dlnT/dk + (dB/dk)/B]. */
static gdouble
_ncm_powspec_analytic_deriv_k (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecAnalytic *psa = NCM_POWSPEC_ANALYTIC (powspec);
  const gdouble P         = _ncm_powspec_analytic_eval (powspec, model, z, k);
  const gdouble B         = _ncm_powspec_analytic_bao (psa, k);

  return P * (psa->n_s / k
              + 2.0 * _ncm_powspec_analytic_dlnT_dk (psa, k)
              + _ncm_powspec_analytic_dbao_dk (psa, k) / B);
}

/*
 * Only a hint, for callers that build a spline from this object. Nothing here
 * is sampled, so the numbers describe the shape's structure rather than a
 * stored table.
 */
static void
_ncm_powspec_analytic_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  *Nz = 20;
  *Nk = 1000;
}

/**
 * ncm_powspec_analytic_new:
 * @shape: a #NcmPowspecAnalyticShape
 * @growth: a #NcmPowspecAnalyticGrowth
 *
 * Creates a new #NcmPowspecAnalytic with the default amplitude, spectral
 * index, scale and matter fraction, which together imitate NumCosmo's own
 * $P(k)$ (see #NcmPowspecAnalytic).
 *
 * Returns: (transfer full): a new #NcmPowspecAnalytic.
 */
NcmPowspecAnalytic *
ncm_powspec_analytic_new (NcmPowspecAnalyticShape shape, NcmPowspecAnalyticGrowth growth)
{
  return g_object_new (NCM_TYPE_POWSPEC_ANALYTIC,
                       "shape", shape,
                       "growth", growth,
                       "kmin", NCM_POWSPEC_ANALYTIC_DEFAULT_KMIN,
                       "kmax", NCM_POWSPEC_ANALYTIC_DEFAULT_KMAX,
                       NULL);
}

/**
 * ncm_powspec_analytic_new_full:
 * @shape: a #NcmPowspecAnalyticShape
 * @growth: a #NcmPowspecAnalyticGrowth
 * @amplitude: the overall amplitude $A$
 * @n_s: the spectral index $n_s$
 * @k_eq: the transfer-function scale $k_\mathrm{eq}$, in $\mathrm{Mpc}^{-1}$
 * @Omega_m: the matter fraction used by the growth factor
 *
 * Creates a new #NcmPowspecAnalytic with every shape parameter given.
 *
 * Returns: (transfer full): a new #NcmPowspecAnalytic.
 */
NcmPowspecAnalytic *
ncm_powspec_analytic_new_full (NcmPowspecAnalyticShape shape, NcmPowspecAnalyticGrowth growth, gdouble amplitude, gdouble n_s, gdouble k_eq, gdouble Omega_m)
{
  return g_object_new (NCM_TYPE_POWSPEC_ANALYTIC,
                       "shape", shape,
                       "growth", growth,
                       "amplitude", amplitude,
                       "n-s", n_s,
                       "k-eq", k_eq,
                       "Omega-m", Omega_m,
                       "kmin", NCM_POWSPEC_ANALYTIC_DEFAULT_KMIN,
                       "kmax", NCM_POWSPEC_ANALYTIC_DEFAULT_KMAX,
                       NULL);
}

/**
 * ncm_powspec_analytic_ref:
 * @psa: a #NcmPowspecAnalytic
 *
 * Increases the reference count of @psa by one.
 *
 * Returns: (transfer full): @psa.
 */
NcmPowspecAnalytic *
ncm_powspec_analytic_ref (NcmPowspecAnalytic *psa)
{
  return g_object_ref (psa);
}

/**
 * ncm_powspec_analytic_free:
 * @psa: a #NcmPowspecAnalytic
 *
 * Decreases the reference count of @psa by one.
 *
 */
void
ncm_powspec_analytic_free (NcmPowspecAnalytic *psa)
{
  g_object_unref (psa);
}

/**
 * ncm_powspec_analytic_clear:
 * @psa: a #NcmPowspecAnalytic
 *
 * If *@psa is not %NULL, decreases its reference count by one and sets *@psa
 * to %NULL.
 *
 */
void
ncm_powspec_analytic_clear (NcmPowspecAnalytic **psa)
{
  g_clear_object (psa);
}

/**
 * ncm_powspec_analytic_get_shape:
 * @psa: a #NcmPowspecAnalytic
 *
 * Returns: the transfer-function shape.
 */
NcmPowspecAnalyticShape
ncm_powspec_analytic_get_shape (NcmPowspecAnalytic *psa)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), NCM_POWSPEC_ANALYTIC_SHAPE_LEN);

  return psa->shape;
}

/**
 * ncm_powspec_analytic_get_growth:
 * @psa: a #NcmPowspecAnalytic
 *
 * Returns: the growth-factor shape.
 */
NcmPowspecAnalyticGrowth
ncm_powspec_analytic_get_growth (NcmPowspecAnalytic *psa)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), NCM_POWSPEC_ANALYTIC_GROWTH_LEN);

  return psa->growth;
}

/**
 * ncm_powspec_analytic_get_amplitude:
 * @psa: a #NcmPowspecAnalytic
 *
 * Returns: the overall amplitude $A$.
 */
gdouble
ncm_powspec_analytic_get_amplitude (NcmPowspecAnalytic *psa)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return psa->amplitude;
}

/**
 * ncm_powspec_analytic_get_n_s:
 * @psa: a #NcmPowspecAnalytic
 *
 * Returns: the spectral index $n_s$.
 */
gdouble
ncm_powspec_analytic_get_n_s (NcmPowspecAnalytic *psa)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return psa->n_s;
}

/**
 * ncm_powspec_analytic_get_k_eq:
 * @psa: a #NcmPowspecAnalytic
 *
 * Returns: the transfer-function scale $k_\mathrm{eq}$, in $\mathrm{Mpc}^{-1}$.
 */
gdouble
ncm_powspec_analytic_get_k_eq (NcmPowspecAnalytic *psa)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return psa->k_eq;
}

/**
 * ncm_powspec_analytic_get_Omega_m:
 * @psa: a #NcmPowspecAnalytic
 *
 * Returns: the matter fraction used by the growth factor.
 */
gdouble
ncm_powspec_analytic_get_Omega_m (NcmPowspecAnalytic *psa)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return psa->Omega_m;
}

/**
 * ncm_powspec_analytic_eval_transfer:
 * @psa: a #NcmPowspecAnalytic
 * @k: the mode $k$, in $\mathrm{Mpc}^{-1}$
 *
 * Evaluates the transfer function $T(k)$ alone.
 *
 * Returns: $T(k)$.
 */
gdouble
ncm_powspec_analytic_eval_transfer (NcmPowspecAnalytic *psa, const gdouble k)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return _ncm_powspec_analytic_transfer (psa, k);
}

/**
 * ncm_powspec_analytic_eval_bao:
 * @psa: a #NcmPowspecAnalytic
 * @k: the mode $k$, in $\mathrm{Mpc}^{-1}$
 *
 * Evaluates the optional oscillatory factor $B(k)$ alone. It is one whenever
 * #NcmPowspecAnalytic:bao-amplitude is zero.
 *
 * Returns: $B(k)$.
 */
gdouble
ncm_powspec_analytic_eval_bao (NcmPowspecAnalytic *psa, const gdouble k)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return _ncm_powspec_analytic_bao (psa, k);
}

/**
 * ncm_powspec_analytic_eval_growth:
 * @psa: a #NcmPowspecAnalytic
 * @z: the redshift $z$
 *
 * Evaluates the growth factor $D(z)$ alone, normalized so that $D(0) = 1$.
 *
 * Returns: $D(z)$.
 */
gdouble
ncm_powspec_analytic_eval_growth (NcmPowspecAnalytic *psa, const gdouble z)
{
  g_return_val_if_fail (NCM_IS_POWSPEC_ANALYTIC (psa), 0.0);

  return _ncm_powspec_analytic_growth (psa, z);
}

