/***************************************************************************
 *            ncm_sf_sbessel.c
 *
 *  Wed Mar 10 17:15:25 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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
 * NcmSFSBessel:
 *
 * Double precision spherical bessel implementation.
 *
 * Implementation of double precision spherical Bessel functions. This module leverages
 * the multiple precision spherical Bessel functions implementation for precise
 * computations. It involves converting the arguments to multiple precision, performing
 * the calculations, and then converting the results back to double precision, ensuring
 * accuracy in the computation of spherical Bessel functions with the convenience of
 * double precision output.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sf_sbessel.h"
#include "math/ncm_mpsf_sbessel.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <mpfr.h>
#endif /* NUMCOSMO_GIR_SCAN */

/**
 * NcmSFSBesselArray:
 *
 * Spherical Bessel function array evaluator with automatic cutoff.
 *
 * This object efficiently evaluates spherical Bessel functions j_l(x) for multiple l values
 * using the Steed/Barnett algorithm. It includes automatic cutoff logic to prevent numerical
 * instability for high l values where j_l(x) becomes negligibly small.
 */

struct _NcmSFSBesselArray
{
  /*< private >*/
  GObject parent_instance;
  guint lmax;
  gdouble threshold;
  gdouble *ell_cutoff_x;
  guint ell_cutoff_size;
};

enum
{
  PROP_0,
  PROP_LMAX,
  PROP_THRESHOLD,
};

G_DEFINE_TYPE (NcmSFSBesselArray, ncm_sf_sbessel_array, G_TYPE_OBJECT)

static void
ncm_sf_sbessel_array_init (NcmSFSBesselArray *sba)
{
  sba->lmax            = 0;
  sba->threshold       = 0.0;
  sba->ell_cutoff_x    = NULL;
  sba->ell_cutoff_size = 0;
}

static void _ncm_sf_sbessel_array_build_cutoff_cache (NcmSFSBesselArray *sba);

static void
_ncm_sf_sbessel_array_constructed (GObject *object)
{
  NcmSFSBesselArray *sba = NCM_SF_SBESSEL_ARRAY (object);

  /* Chain up */
  G_OBJECT_CLASS (ncm_sf_sbessel_array_parent_class)->constructed (object);

  /* Build cutoff cache */
  _ncm_sf_sbessel_array_build_cutoff_cache (sba);
}

static void
_ncm_sf_sbessel_array_dispose (GObject *object)
{
  NcmSFSBesselArray *sba = NCM_SF_SBESSEL_ARRAY (object);

  g_clear_pointer (&sba->ell_cutoff_x, g_free);

  /* Chain up */
  G_OBJECT_CLASS (ncm_sf_sbessel_array_parent_class)->dispose (object);
}

static void
_ncm_sf_sbessel_array_finalize (GObject *object)
{
  /* Chain up */
  G_OBJECT_CLASS (ncm_sf_sbessel_array_parent_class)->finalize (object);
}

static void
_ncm_sf_sbessel_array_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSFSBesselArray *sba = NCM_SF_SBESSEL_ARRAY (object);

  switch (prop_id)
  {
    case PROP_LMAX:
      sba->lmax = g_value_get_uint (value);
      break;
    case PROP_THRESHOLD:
      sba->threshold = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sf_sbessel_array_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSFSBesselArray *sba = NCM_SF_SBESSEL_ARRAY (object);

  switch (prop_id)
  {
    case PROP_LMAX:
      g_value_set_uint (value, sba->lmax);
      break;
    case PROP_THRESHOLD:
      g_value_set_double (value, sba->threshold);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_sf_sbessel_array_class_init (NcmSFSBesselArrayClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_sf_sbessel_array_constructed;
  object_class->set_property = &_ncm_sf_sbessel_array_set_property;
  object_class->get_property = &_ncm_sf_sbessel_array_get_property;
  object_class->dispose      = &_ncm_sf_sbessel_array_dispose;
  object_class->finalize     = &_ncm_sf_sbessel_array_finalize;

  /**
   * NcmSFSBesselArray:lmax:
   *
   * Maximum l value for which the array can compute spherical Bessel functions. This
   * value is relevant to automatic cutoff logic.
   */
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum l value",
                                                      0, G_MAXUINT, 10000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSFSBesselArray:threshold:
   *
   * Threshold value below which j_l(x) is considered negligible and set to zero.
   */
  g_object_class_install_property (object_class,
                                   PROP_THRESHOLD,
                                   g_param_spec_double ("threshold",
                                                        NULL,
                                                        "Threshold value",
                                                        0.0, 1.0, 1.0e-100,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
_ncm_sf_sbessel_array_build_cutoff_cache (NcmSFSBesselArray *sba)
{
  const gdouble log_threshold = log (sba->threshold);
  const guint max_ell         = sba->lmax;
  guint ell;

  /* Allocate cache array indexed by ell */
  g_clear_pointer (&sba->ell_cutoff_x, g_free);

  sba->ell_cutoff_size = max_ell + 1;
  sba->ell_cutoff_x    = g_new (gdouble, sba->ell_cutoff_size);

  /* Compute x value where j_ell(x) = threshold for each ell */
  sba->ell_cutoff_x[0] = 0.0; /* ell=0 case */

  for (ell = 1; ell <= max_ell; ell++)
  {
    const gdouble ell_d       = (gdouble) ell;
    const gdouble ln_j_l_peak = log (fabs (0.447 / pow (ell_d, 2.0 / 3.0))); /* Approximate peak value of j_ell */
    const gdouble x_ell       = exp (((1.0 + ell_d) * M_LN2 + lgamma (1.5 + ell_d) + log_threshold - 0.5 * M_PI - ln_j_l_peak) / ell_d);

    sba->ell_cutoff_x[ell] = x_ell;
  }
}

/* Fast ell cutoff lookup using binary search on cached values */
static gint
_ncm_sf_sbessel_array_get_ell_cutoff (NcmSFSBesselArray *sba, gdouble x)
{
  guint lo = 0;
  guint hi = sba->ell_cutoff_size - 1;

  /* If x is beyond our cache, return max ell */
  if (x >= sba->ell_cutoff_x[hi])
    return hi;

  /* Binary search for position where x fits in array */
  /* Find largest index where ell_cutoff_x[index] <= x */
  while (lo < hi)
  {
    guint mid = (lo + hi + 1) / 2;

    if (sba->ell_cutoff_x[mid] <= x)
      lo = mid;
    else
      hi = mid - 1;
  }

  /* Return index + 1 as the ell threshold */
  return lo + 1;
}

/**
 * ncm_sf_sbessel_array_new:
 *
 * Creates a new #NcmSFSBesselArray for computing spherical Bessel functions up to l =
 * 10000 with automatic cutoff at 1e-100.
 *
 * Returns: (transfer full): a new #NcmSFSBesselArray
 */
NcmSFSBesselArray *
ncm_sf_sbessel_array_new ()
{
  NcmSFSBesselArray *sba = g_object_new (NCM_TYPE_SF_SBESSEL_ARRAY,
                                         NULL);

  return sba;
}

/**
 * ncm_sf_sbessel_array_new_full:
 * @lmax: maximum l value
 * @threshold: threshold value below which j_l(x) is negligible
 *
 * Creates a new #NcmSFSBesselArray for computing spherical Bessel functions
 * up to l = @lmax with automatic cutoff at @threshold.
 *
 * Returns: (transfer full): a new #NcmSFSBesselArray
 */
NcmSFSBesselArray *
ncm_sf_sbessel_array_new_full (guint lmax, gdouble threshold)
{
  NcmSFSBesselArray *sba = g_object_new (NCM_TYPE_SF_SBESSEL_ARRAY,
                                         "lmax", lmax,
                                         "threshold", threshold,
                                         NULL);

  return sba;
}

/**
 * ncm_sf_sbessel_array_ref:
 * @sba: a #NcmSFSBesselArray
 *
 * Increases the reference count of @sba by one.
 *
 * Returns: (transfer full): @sba
 */
NcmSFSBesselArray *
ncm_sf_sbessel_array_ref (NcmSFSBesselArray *sba)
{
  return g_object_ref (sba);
}

/**
 * ncm_sf_sbessel_array_free:
 * @sba: a #NcmSFSBesselArray
 *
 * Decreases the reference count of @sba by one.
 */
void
ncm_sf_sbessel_array_free (NcmSFSBesselArray *sba)
{
  g_object_unref (sba);
}

/**
 * ncm_sf_sbessel_array_clear:
 * @sba: a #NcmSFSBesselArray
 *
 * If @sba is different from NULL, decreases the reference count of
 * @sba by one and sets @sba to NULL.
 */
void
ncm_sf_sbessel_array_clear (NcmSFSBesselArray **sba)
{
  g_clear_object (sba);
}

static inline gdouble
_ncm_sf_sbessel_asymptotic_3 (gint ell, gdouble x)
{
  const gdouble x_inv  = 1.0 / x;
  const gdouble x_inv2 = x_inv * x_inv;
  const gdouble x_inv3 = x_inv2 * x_inv;

  const gdouble llp1 = (gdouble) ell * (ell + 1);

  const gdouble a1 = 0.5 * llp1;
  const gdouble a2 = 0.125 * llp1 * (llp1 - 2.0);
  const gdouble a3 = (1.0 / 48.0) * llp1 * (llp1 - 2.0) * (llp1 - 6.0);

  const gdouble phi = x - 0.5 * G_PI * ell;
  const gdouble s   = sin (phi);
  const gdouble c   = cos (phi);

  return x_inv * (s
                  - a1 * x_inv  * c
                  - a2 * x_inv2 * s
                  + a3 * x_inv3 * c);
}

static void
_ncm_sf_sbessel_array_eval_asymptotic (gint lmax, gdouble x, gdouble *jl_x)
{
  const gdouble x_inv = 1.0 / x;
  gint L;

  /* seed top two values */
  jl_x[lmax] = ncm_sf_sbessel (lmax, x);

  if (lmax > 0)
    jl_x[lmax - 1] = ncm_sf_sbessel (lmax - 1, x);

  /* downward recurrence */
  if (lmax > 1)
  {
    gdouble PL  = lmax * x_inv;
    gdouble XP2 = jl_x[lmax];
    gdouble FP;
    gint LL = lmax;

#pragma GCC unroll 4

    for ( ; LL > 1; --LL)
    {
      jl_x[LL - 1] = PL * jl_x[LL] + XP2;
      FP           = PL * jl_x[LL - 1] - jl_x[LL];
      XP2          = FP;
      PL          -= x_inv;
    }
  }

  /* normalization using j0 */
  {
    const gdouble norm = fabs (jl_x[0]);

    if (norm > 0.0)
    {
      const gdouble W = x_inv / norm;

      for (L = 0; L <= lmax; ++L)
        jl_x[L] *= W;
    }
  }
}

/**
 * ncm_sf_sbessel_array_eval:
 * @sba: a #NcmSFSBesselArray
 * @ell: maximum l value to compute (must be <= lmax)
 * @x: argument value
 * @jl_x: (out) (array) (transfer none): output array of size (ell+1) for j_l(x) values
 *
 * Computes spherical Bessel functions j_l(x) for l = 0 to min(ell, lmax, cutoff(x))
 * using the Steed/Barnett algorithm with automatic cutoff for numerical stability.
 * Values beyond the cutoff are set to zero.
 */
void
ncm_sf_sbessel_array_eval (NcmSFSBesselArray *sba, guint ell, gdouble x, gdouble *jl_x)
{
  gint lmax;

  g_return_if_fail (NCM_IS_SF_SBESSEL_ARRAY (sba));
  g_return_if_fail (ell <= sba->lmax);
  g_return_if_fail (jl_x != NULL);

  /* Initialize output array to zero */
  memset (jl_x, 0, sizeof (gdouble) * (ell + 1));

  /* Apply cutoff */
  lmax = GSL_MIN ((gint) ell, _ncm_sf_sbessel_array_get_ell_cutoff (sba, x));

  if (G_UNLIKELY (x == 0.0))
  {
    jl_x[0] = 1.0;

    return;
  }
  else if (x < 2.0 * GSL_ROOT4_DBL_EPSILON)
  {
    /* first two terms of Taylor series */
    gdouble inv_fact = 1.0; /* 1/(1 3 5 ... (2l+1)) */
    gdouble x_l      = 1.0; /* x^l */
    gint l;

    for (l = 0; l <= lmax; l++)
    {
      jl_x[l]   = x_l * inv_fact;
      jl_x[l]  *= 1.0 - 0.5 * x * x / (2.0 * l + 3.0);
      inv_fact /= 2.0 * l + 3.0;
      x_l      *= x;
    }

    return;
  }
  else if ((x > (gdouble) (lmax + 1.0)) && FALSE)
  {
    /* Use asymptotic expansion for large x */
    _ncm_sf_sbessel_array_eval_asymptotic (lmax, x, jl_x);
  }
  else
  {
    /* Steed/Barnett algorithm [Comp. Phys. Comm. 21, 297 (1981)] */
    gdouble x_inv = 1.0 / x;
    gdouble W     = 2.0 * x_inv;
    gdouble F     = 1.0;
    gdouble FP    = (lmax + 1.0) * x_inv;
    gdouble B     = 2.0 * FP + x_inv;
    gdouble end   = B + 20000.0 * W;
    gdouble D     = 1.0 / B;
    gdouble del   = -D;

    FP += del;

    /* continued fraction */
    do {
      B   += W;
      D    = 1.0 / (B - D);
      del *= (B * D - 1.);
      FP  += del;

      if (D < 0.0)
        F = -F;

      if (B > end)
      {
        g_warning ("ncm_sf_sbessel_array_eval: continued fraction did not converge for lmax=%d, x=%g", lmax, x);

        return;
      }
    } while (fabs (del) >= fabs (FP) * GSL_DBL_EPSILON);

    FP *= F;

    if (lmax > 0)
    {
      /* downward recursion */
      gdouble XP2 = FP;
      gdouble PL  = lmax * x_inv;
      gint L      = lmax;
      gint LP;

      jl_x[lmax] = F;
#pragma GCC unroll 4

      for (LP = 1; LP <= lmax; LP++)
      {
        jl_x[L - 1] = PL * jl_x[L] + XP2;
        FP          = PL * jl_x[L - 1] - jl_x[L];
        XP2         = FP;
        PL         -= x_inv;
        --L;
      }

      F = jl_x[0];
    }

    /* normalization */
    W       = x_inv / hypot (FP, F);
    jl_x[0] = W * F;

    if (lmax > 0)
    {
      gint L;

      for (L = 1; L <= lmax; L++)
      {
        jl_x[L] *= W;
      }
    }

    return;
  }
}

/**
 * ncm_sf_sbessel_array_eval_ell_cutoff:
 * @sba: a #NcmSFSBesselArray
 * @x: argument value
 *
 * Determines the maximum l value for which j_l(x) is above the threshold.
 * For l values above this cutoff, j_l(x) is negligibly small and set to zero.
 *
 * Returns: the cutoff l value
 */
guint
ncm_sf_sbessel_array_eval_ell_cutoff (NcmSFSBesselArray *sba, gdouble x)
{
  g_return_val_if_fail (NCM_IS_SF_SBESSEL_ARRAY (sba), 0);

  return _ncm_sf_sbessel_array_get_ell_cutoff (sba, x);
}

/**
 * ncm_sf_sbessel_array_get_lmax:
 * @sba: a #NcmSFSBesselArray
 *
 * Gets the maximum l value for this array.
 *
 * Returns: the lmax value
 */
guint
ncm_sf_sbessel_array_get_lmax (NcmSFSBesselArray *sba)
{
  g_return_val_if_fail (NCM_IS_SF_SBESSEL_ARRAY (sba), 0);

  return sba->lmax;
}

/**
 * ncm_sf_sbessel_array_get_threshold:
 * @sba: a #NcmSFSBesselArray
 *
 * Gets the threshold value for this array.
 *
 * Returns: the threshold value
 */
gdouble
ncm_sf_sbessel_array_get_threshold (NcmSFSBesselArray *sba)
{
  g_return_val_if_fail (NCM_IS_SF_SBESSEL_ARRAY (sba), 0.0);

  return sba->threshold;
}

/**
 * ncm_sf_sbessel:
 * @l: Spherical Bessel order $\ell$
 * @x: Spherical Bessel argument $x$
 *
 * Computes Spherical Bessel function $j_\ell(x)$.
 *
 * Returns: the value of $j_\ell(x)$.
 */
gdouble
ncm_sf_sbessel (gulong l, gdouble x)
{
  MPFR_DECL_INIT (res, 53); /* Should it be 53? FIXME */

  gdouble res_d;

  ncm_mpsf_sbessel_d (l, x, res, GMP_RNDN);
  res_d = mpfr_get_d (res, GMP_RNDN);

  return res_d;
}

static void
_taylor_jl (const glong l, const gdouble x, const gdouble x2, const gdouble x3, const gdouble jl, const gdouble jlp1, gdouble *deriv)
{
  const gdouble llm1    = l * (l - 1.0);
  const gdouble llm1lm2 = llm1 * (l - 2.0);

  deriv[0] = jl;
  deriv[1] = (l * jl - x * jlp1) / x;
  deriv[2] = (((llm1 - x2) * jl + 2.0 * x * jlp1) / x2) / (1.0 * 2.0);
  deriv[3] = (((llm1lm2 - (l - 2.0) * x2) * jl - x * (l * (l + 1.0) + 6.0 - x2) * jlp1) / x3) / (1.0 * 2.0 * 3.0);
}

/**
 * ncm_sf_sbessel_taylor:
 * @l: Spherical Bessel order $\ell$
 * @x: Spherical Bessel argument $x$
 * @djl: (out) (array fixed-size=4): Output power series coefficients
 *
 * Computes Spherical Bessel function power series
 * coefficients up to order three, i.e.,
 * $$\left(j_\ell(x),\; j'_\ell(x), \frac{j''_\ell(x)}{2!}, \frac{j'''_\ell(x)}{3!}\right).$$
 */
void
ncm_sf_sbessel_taylor (gulong l, gdouble x, gdouble *djl)
{
  const gdouble jl   = ncm_sf_sbessel (l, x);
  const gdouble jlp1 = ncm_sf_sbessel (l + 1, x);
  const gdouble x2   = x * x;
  const gdouble x3   = x2 * x;

  _taylor_jl (l, x, x2, x3, jl, jlp1, djl);

  return;
}

static gdouble
_ncm_sf_sbessel_spline_calc (gdouble x, gpointer data)
{
  gulong *l = (gulong *) data;

  return ncm_sf_sbessel (*l, x);
}

/**
 * ncm_sf_sbessel_spline:
 * @l: Spherical Bessel order $\ell$.
 * @xi: Spherical Bessel interval lower-bound $x_i$.
 * @xf: Spherical Bessel interval lower-bound $x_f$.
 * @reltol: Interpolation error tolerance.
 *
 * Computes a spline approximation of the Spherical Bessel
 * $j_\ell$ in the interval $[x_i, x_f]$.
 *
 * Returns: (transfer full): A #NcmSpline with the Spherical Bessel approximation.
 */
NcmSpline *
ncm_sf_sbessel_spline (gulong l, gdouble xi, gdouble xf, gdouble reltol)
{
  NcmSpline *s = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  gsl_function F;

  F.function = &_ncm_sf_sbessel_spline_calc;
  F.params   = &l;

  ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, xi, xf, 0, reltol);

  return s;
}

