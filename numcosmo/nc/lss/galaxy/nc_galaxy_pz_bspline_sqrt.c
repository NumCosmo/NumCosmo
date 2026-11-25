/***************************************************************************
 *            nc_galaxy_pz_bspline_sqrt.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_pz_bspline_sqrt.c
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
 * NcGalaxyPzBSplineSqrt:
 *
 * The redshift density of one galaxy as the square of a B-spline.
 *
 * The density is
 * $$
 * p(z) = s(z)^2,
 * $$
 * where $s$ is the #NcmBSpline of order %NC_GALAXY_PZ_BSPLINE_SQRT_ORDER in
 * #NcGalaxyPzBSplineSqrt:bspline. It is non-negative and $C^4$, and its logarithm is
 * $2\ln\vert s(z)\vert$. The range is that of $s$. The normalization is the exact
 * integral of $s^2$ (6-point Gauss-Legendre on each polynomial piece), computed at
 * construction, and sampling inverts the same integral. The Chebyshev moments of
 * nc_galaxy_pz_cheb_moments() are exact as well: $p\,T_k$ is a polynomial of degree
 * $10 + k$ on each piece, integrated by Gauss-Legendre with enough points.
 *
 * nc_galaxy_pz_bspline_sqrt_new_from_table() builds the object from a table
 * $(z_j, p_j)$: $s$ is then the interpolating B-spline through $(z_j, \sqrt{p_j})$.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_pz_bspline_sqrt.h"
#include "ncm/spline/ncm_spline_bspline.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyPzBSplineSqrt
{
  /*< private >*/
  NcGalaxyPz parent_instance;
  NcmBSpline *s;
  const gdouble *brk; /* breakpoints of s: n_brk of its knots, from the index order - 1 */
  guint n_brk;
  GArray *cum; /* integral of s^2 from brk[0] to each breakpoint; the last is the norm */
};

G_DEFINE_TYPE (NcGalaxyPzBSplineSqrt, nc_galaxy_pz_bspline_sqrt, NC_TYPE_GALAXY_PZ)

enum
{
  PROP_0,
  PROP_BSPLINE,
};

/* Points of the Gauss-Legendre rule: 6 points are exact up to degree 11, s^2 has degree 10 */
#define NC_GALAXY_PZ_BSPLINE_SQRT_GL_N (6)

/* The Gauss-Legendre rule, allocated once in class_init */
static gsl_integration_glfixed_table *_nc_galaxy_pz_bspline_sqrt_gl = NULL;

static void
nc_galaxy_pz_bspline_sqrt_init (NcGalaxyPzBSplineSqrt *gpzs)
{
  gpzs->s     = NULL;
  gpzs->brk   = NULL;
  gpzs->n_brk = 0;
  gpzs->cum   = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void
_nc_galaxy_pz_bspline_sqrt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (object);

  g_return_if_fail (NC_IS_GALAXY_PZ_BSPLINE_SQRT (object));

  switch (prop_id)
  {
    case PROP_BSPLINE:
      gpzs->s = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_pz_bspline_sqrt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (object);

  g_return_if_fail (NC_IS_GALAXY_PZ_BSPLINE_SQRT (object));

  switch (prop_id)
  {
    case PROP_BSPLINE:
      g_value_set_object (value, gpzs->s);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _nc_galaxy_pz_bspline_sqrt_s2_integ (NcGalaxyPzBSplineSqrt *gpzs, const gdouble a, const gdouble b);

/*
 * Stores the breakpoints of s and computes the integral of s^2 from the first breakpoint
 * to each of them, the last being the normalization. The breakpoints are the knots of s
 * from the index order - 1 to the index equal to the number of control points. Aborts
 * when the B-spline is missing, when its order is not NC_GALAXY_PZ_BSPLINE_SQRT_ORDER, or
 * when the normalization is not positive.
 */
static void
_nc_galaxy_pz_bspline_sqrt_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_pz_bspline_sqrt_parent_class)->constructed (object);
  {
    NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (object);
    gdouble total               = 0.0;
    guint j;

    if (gpzs->s == NULL)
      g_error ("NcGalaxyPzBSplineSqrt: the B-spline is required.");

    if (ncm_bspline_get_order (gpzs->s) != NC_GALAXY_PZ_BSPLINE_SQRT_ORDER)
      g_error ("NcGalaxyPzBSplineSqrt: the B-spline must have order %d, got %u.",
               NC_GALAXY_PZ_BSPLINE_SQRT_ORDER, ncm_bspline_get_order (gpzs->s));

    gpzs->brk   = ncm_vector_data (ncm_bspline_peek_knots (gpzs->s)) + NC_GALAXY_PZ_BSPLINE_SQRT_ORDER - 1;
    gpzs->n_brk = ncm_vector_len (ncm_bspline_peek_control_points (gpzs->s)) - NC_GALAXY_PZ_BSPLINE_SQRT_ORDER + 2;

    g_array_set_size (gpzs->cum, gpzs->n_brk);
    g_array_index (gpzs->cum, gdouble, 0) = total;

    for (j = 0; j + 1 < gpzs->n_brk; j++)
    {
      total                                    += _nc_galaxy_pz_bspline_sqrt_s2_integ (gpzs, gpzs->brk[j], gpzs->brk[j + 1]);
      g_array_index (gpzs->cum, gdouble, j + 1) = total;
    }

    if (!(total > 0.0))
      g_error ("NcGalaxyPzBSplineSqrt: the normalization is %g; it must be positive.", total);
  }
}

static void
_nc_galaxy_pz_bspline_sqrt_dispose (GObject *object)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (object);

  ncm_bspline_clear (&gpzs->s);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_pz_bspline_sqrt_parent_class)->dispose (object);
}

static void
_nc_galaxy_pz_bspline_sqrt_finalize (GObject *object)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (object);

  g_clear_pointer (&gpzs->cum, g_array_unref);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_pz_bspline_sqrt_parent_class)->finalize (object);
}

static void _nc_galaxy_pz_bspline_sqrt_get_range (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max);
static gdouble _nc_galaxy_pz_bspline_sqrt_eval_p (NcGalaxyPz *gpz, const gdouble z);
static gdouble _nc_galaxy_pz_bspline_sqrt_eval_ln_p (NcGalaxyPz *gpz, const gdouble z);
static gdouble _nc_galaxy_pz_bspline_sqrt_get_norm (NcGalaxyPz *gpz);
static gdouble _nc_galaxy_pz_bspline_sqrt_gen (NcGalaxyPz *gpz, NcmRNG *rng);
static void _nc_galaxy_pz_bspline_sqrt_cheb_moments (NcGalaxyPz *gpz, const gdouble a, const gdouble b, const guint n, GArray *moments);

static void
nc_galaxy_pz_bspline_sqrt_class_init (NcGalaxyPzBSplineSqrtClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcGalaxyPzClass *gpz_class = NC_GALAXY_PZ_CLASS (klass);

  object_class->set_property = &_nc_galaxy_pz_bspline_sqrt_set_property;
  object_class->get_property = &_nc_galaxy_pz_bspline_sqrt_get_property;
  object_class->constructed  = &_nc_galaxy_pz_bspline_sqrt_constructed;
  object_class->dispose      = &_nc_galaxy_pz_bspline_sqrt_dispose;
  object_class->finalize     = &_nc_galaxy_pz_bspline_sqrt_finalize;

  /**
   * NcGalaxyPzBSplineSqrt:bspline:
   *
   * The B-spline $s$, of order %NC_GALAXY_PZ_BSPLINE_SQRT_ORDER.
   */
  g_object_class_install_property (object_class,
                                   PROP_BSPLINE,
                                   g_param_spec_object ("bspline",
                                                        NULL,
                                                        "B-spline whose square is the density",
                                                        NCM_TYPE_BSPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gpz_class->get_range    = &_nc_galaxy_pz_bspline_sqrt_get_range;
  gpz_class->eval_p       = &_nc_galaxy_pz_bspline_sqrt_eval_p;
  gpz_class->eval_ln_p    = &_nc_galaxy_pz_bspline_sqrt_eval_ln_p;
  gpz_class->get_norm     = &_nc_galaxy_pz_bspline_sqrt_get_norm;
  gpz_class->gen          = &_nc_galaxy_pz_bspline_sqrt_gen;
  gpz_class->cheb_moments = &_nc_galaxy_pz_bspline_sqrt_cheb_moments;

  _nc_galaxy_pz_bspline_sqrt_gl = gsl_integration_glfixed_table_alloc (NC_GALAXY_PZ_BSPLINE_SQRT_GL_N);
}

/*
 * Integral of s^2 over [a, b] inside one polynomial piece of s, where s^2 has degree 10
 * and the 6-point Gauss-Legendre rule is exact.
 */
static gdouble
_nc_galaxy_pz_bspline_sqrt_s2_integ (NcGalaxyPzBSplineSqrt *gpzs, const gdouble a, const gdouble b)
{
  const gsl_integration_glfixed_table *gl = _nc_galaxy_pz_bspline_sqrt_gl;
  const gdouble c                         = 0.5 * (a + b);
  const gdouble h                         = 0.5 * (b - a);
  gdouble res                             = 0.0;
  guint j;

  for (j = 0; j < NC_GALAXY_PZ_BSPLINE_SQRT_GL_N / 2; j++)
  {
    const gdouble s_m = ncm_bspline_eval (gpzs->s, c - h * gl->x[j]);
    const gdouble s_p = ncm_bspline_eval (gpzs->s, c + h * gl->x[j]);

    res += gl->w[j] * (s_m * s_m + s_p * s_p);
  }

  return h * res;
}

static void
_nc_galaxy_pz_bspline_sqrt_get_range (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (gpz);

  ncm_bspline_get_range (gpzs->s, z_min, z_max);
}

static gdouble
_nc_galaxy_pz_bspline_sqrt_eval_p (NcGalaxyPz *gpz, const gdouble z)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (gpz);
  const gdouble s_z           = ncm_bspline_eval (gpzs->s, z);

  return s_z * s_z;
}

static gdouble
_nc_galaxy_pz_bspline_sqrt_eval_ln_p (NcGalaxyPz *gpz, const gdouble z)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (gpz);

  return 2.0 * log (fabs (ncm_bspline_eval (gpzs->s, z)));
}

static gdouble
_nc_galaxy_pz_bspline_sqrt_get_norm (NcGalaxyPz *gpz)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (gpz);

  return g_array_index (gpzs->cum, gdouble, gpzs->n_brk - 1);
}

/*
 * Draws z with density s^2 / norm: picks the polynomial piece by its share of the
 * normalization, then solves int_{brk_j}^{z} s^2 = target in it by Newton steps kept
 * inside a bisection bracket.
 */
static gdouble
_nc_galaxy_pz_bspline_sqrt_gen (NcGalaxyPz *gpz, NcmRNG *rng)
{
  NcGalaxyPzBSplineSqrt *gpzs = NC_GALAXY_PZ_BSPLINE_SQRT (gpz);
  const gdouble total         = g_array_index (gpzs->cum, gdouble, gpzs->n_brk - 1);
  const gdouble target        = ncm_rng_uniform_gen (rng, 0.0, total);
  guint lo                    = 0;
  guint hi                    = gpzs->n_brk - 1;
  gdouble a, b, z, rem;
  guint iter;

  while (hi - lo > 1)
  {
    const guint mid = (lo + hi) / 2;

    if (target < g_array_index (gpzs->cum, gdouble, mid))
      hi = mid;
    else
      lo = mid;
  }

  a   = gpzs->brk[lo];
  b   = gpzs->brk[lo + 1];
  rem = target - g_array_index (gpzs->cum, gdouble, lo);

  {
    const gdouble x0 = a;
    gdouble z_lo     = a;
    gdouble z_hi     = b;

    z = 0.5 * (a + b);

    for (iter = 0; iter < 200; iter++)
    {
      const gdouble F  = _nc_galaxy_pz_bspline_sqrt_s2_integ (gpzs, x0, z) - rem;
      const gdouble sz = ncm_bspline_eval (gpzs->s, z);
      const gdouble dF = sz * sz;
      gdouble z_new;

      if (F > 0.0)
        z_hi = z;
      else
        z_lo = z;

      z_new = (dF > 0.0) ? z - F / dF : 0.5 * (z_lo + z_hi);

      if (!((z_new > z_lo) && (z_new < z_hi)))
        z_new = 0.5 * (z_lo + z_hi);

      if (fabs (z_new - z) <= 4.0 * GSL_DBL_EPSILON * fabs (z))
      {
        z = z_new;
        break;
      }

      z = z_new;
    }
  }

  return z;
}

/*
 * The moments of s^2 on [a, b] (nc_galaxy_pz_cheb_moments() checked the interval and sized
 * moments to n + 1): on each polynomial piece of s that meets [a, b], the Gauss-Legendre rule
 * of q points with 2 q - 1 >= 10 + n integrates s^2 T_k exactly for every k <= n.
 */
static void
_nc_galaxy_pz_bspline_sqrt_cheb_moments (NcGalaxyPz *gpz, const gdouble a, const gdouble b, const guint n, GArray *moments)
{
  NcGalaxyPzBSplineSqrt *gpzs        = NC_GALAXY_PZ_BSPLINE_SQRT (gpz);
  const guint q                      = (n + 12) / 2;
  gsl_integration_glfixed_table *glt = gsl_integration_glfixed_table_alloc (q);
  gdouble *m                         = &g_array_index (moments, gdouble, 0);
  guint lo                           = 0;
  guint hi                           = gpzs->n_brk - 1;
  guint j, k, i;

  for (k = 0; k <= n; k++)
    m[k] = 0.0;

  /* The piece [brk[lo], brk[lo + 1]] that holds a */
  while (hi - lo > 1)
  {
    const guint mid = (lo + hi) / 2;

    if (a < gpzs->brk[mid])
      hi = mid;
    else
      lo = mid;
  }

  for (j = lo; (j + 1 < gpzs->n_brk) && (gpzs->brk[j] < b); j++)
  {
    const gdouble z_lo = GSL_MAX (a, gpzs->brk[j]);
    const gdouble z_hi = GSL_MIN (b, gpzs->brk[j + 1]);

    if (!(z_hi > z_lo))
      continue;

    for (i = 0; i < q; i++)
    {
      gdouble zq, wq, t, T0, T1, s_z, pw;

      gsl_integration_glfixed_point (z_lo, z_hi, i, &zq, &wq, glt);

      s_z   = ncm_bspline_eval (gpzs->s, zq);
      pw    = wq * s_z * s_z;
      t     = (2.0 * zq - (a + b)) / (b - a);
      T0    = 1.0;
      T1    = t;
      m[0] += pw;

      if (n == 0)
        continue;

      m[1] += pw * t;

      for (k = 2; k <= n; k++)
      {
        const gdouble T2 = 2.0 * t * T1 - T0;

        m[k] += pw * T2;
        T0    = T1;
        T1    = T2;
      }
    }
  }

  gsl_integration_glfixed_table_free (glt);
}

/**
 * nc_galaxy_pz_bspline_sqrt_new:
 * @s: a #NcmBSpline of order %NC_GALAXY_PZ_BSPLINE_SQRT_ORDER
 *
 * Creates the density $s(z)^2$, referencing @s. Aborts when the order of @s is not
 * %NC_GALAXY_PZ_BSPLINE_SQRT_ORDER or when the integral of $s^2$ over its range is not
 * positive.
 *
 * Returns: (transfer full): a new #NcGalaxyPzBSplineSqrt.
 */
NcGalaxyPzBSplineSqrt *
nc_galaxy_pz_bspline_sqrt_new (NcmBSpline *s)
{
  return g_object_new (NC_TYPE_GALAXY_PZ_BSPLINE_SQRT,
                       "bspline", s,
                       NULL);
}

/**
 * nc_galaxy_pz_bspline_sqrt_new_from_table:
 * @zv: the redshifts $z_j$, increasing
 * @pv: the values $p_j$ at @zv
 *
 * Creates the density $s(z)^2$ where $s$ is the interpolating B-spline of order
 * %NC_GALAXY_PZ_BSPLINE_SQRT_ORDER through $(z_j, \sqrt{p_j})$. The vectors are read and
 * not kept. Aborts on a $p_j$ that is negative or not a number, on fewer than
 * %NC_GALAXY_PZ_BSPLINE_SQRT_ORDER points, or when every $p_j$ is zero.
 *
 * Returns: (transfer full): a new #NcGalaxyPzBSplineSqrt.
 */
NcGalaxyPzBSplineSqrt *
nc_galaxy_pz_bspline_sqrt_new_from_table (NcmVector *zv, NcmVector *pv)
{
  const guint n = ncm_vector_len (pv);
  NcGalaxyPzBSplineSqrt *gpzs;
  NcmSplineBSpline *sbs;
  NcmBSpline *s;
  NcmVector *sv;
  guint j;

  if (n < NC_GALAXY_PZ_BSPLINE_SQRT_ORDER)
    g_error ("nc_galaxy_pz_bspline_sqrt_new_from_table: the p(z) table has %u points, at least %d are required.",
             n, NC_GALAXY_PZ_BSPLINE_SQRT_ORDER);

  sv = ncm_vector_new (n);

  for (j = 0; j < n; j++)
  {
    const gdouble p_j = ncm_vector_get (pv, j);

    if (!(p_j >= 0.0))
      g_error ("nc_galaxy_pz_bspline_sqrt_new_from_table: the p(z) table value %u is %g at z = %g; values must be non-negative.",
               j, p_j, ncm_vector_get (zv, j));

    ncm_vector_fast_set (sv, j, sqrt (p_j));
  }

  sbs  = ncm_spline_bspline_new_full (NC_GALAXY_PZ_BSPLINE_SQRT_ORDER, zv, sv, TRUE);
  s    = ncm_spline_bspline_dup_bspline (sbs);
  gpzs = nc_galaxy_pz_bspline_sqrt_new (s);

  ncm_bspline_free (s);
  ncm_spline_free (NCM_SPLINE (sbs));
  ncm_vector_free (sv);

  return gpzs;
}

/**
 * nc_galaxy_pz_bspline_sqrt_peek_bspline:
 * @gpzs: a #NcGalaxyPzBSplineSqrt
 *
 * Gets #NcGalaxyPzBSplineSqrt:bspline.
 *
 * Returns: (transfer none): the B-spline $s$.
 */
NcmBSpline *
nc_galaxy_pz_bspline_sqrt_peek_bspline (NcGalaxyPzBSplineSqrt *gpzs)
{
  return gpzs->s;
}

