/***************************************************************************
 *            ncm_spectral_panels_product.c
 *
 *  Sat Oct 10 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * NcmSpectralPanelsProduct:
 *
 * Exact integrals of products of two panel sets.
 *
 * ncm_spectral_panels_product_integrate() computes, for each row $r$,
 * $$I_r = \int k^p P_r(k / s_{1r}) Q_r(k / s_{2r}) \, dk$$
 * over the range both panel sets cover, $p$ being #NcmSpectralPanelsProduct:p. On each
 * cell of the common refinement of the two sets of edges both factors are polynomials, so
 * the integral is the bilinear form of their coefficients with the tables
 * $\int_{-1}^{1} T_i T_j T_m \, dt$, $m \le p$, and carries no quadrature error.
 *
 * The rows whose two scales are equal share one refinement in $x = k / s_r$. Otherwise each
 * row is integrated alone, in the variable of the side whose scale differs from 1, or of the
 * first side when both do.
 *
 * The error estimate propagates $d(PQ) = |P| \, dQ + |Q| \, dP$, with $dP$ on a cell the
 * sum of the absolute coefficients of its panel above the order of the level below, and
 * $|P|$, $|Q|$ sampled at the five Gauss-Legendre nodes of every cell.
 *
 * The tables are built at construction and calls take their scratch space from a pool, so
 * one object may serve calls in several threads.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_panels_product.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/core/ncm_memory_pool.h"

#include <math.h>
#include <string.h>

/* The tables are built once per object up to this number of coefficients */
#define NCM_SPECTRAL_PANELS_PRODUCT_TABLE_N (65)

/* Five-node Gauss-Legendre rule on [-1, 1], for the error terms */
#define NCM_SPECTRAL_PANELS_PRODUCT_GL5_N (5)

static const gdouble _ncm_spectral_panels_product_gl5_x[NCM_SPECTRAL_PANELS_PRODUCT_GL5_N] = {
  -0.9061798459386640, -0.5384693101056831, 0.0,
  0.5384693101056831, 0.9061798459386640
};

static const gdouble _ncm_spectral_panels_product_gl5_w[NCM_SPECTRAL_PANELS_PRODUCT_GL5_N] = {
  0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
  0.4786286704993665, 0.2369268850561891
};

enum
{
  PROP_0,
  PROP_P,
};

struct _NcmSpectralPanelsProduct
{
  /*< private >*/
  GObject parent_instance;
  guint p;
  GArray *tables;           /* gdouble, the p + 1 tables A_m[i][j] of row stride TABLE_N */
  NcmMemoryPool *work_pool; /* NcmSpectralPanelsProductWork, one per concurrent call */
};

/* The scratch space of one call */
typedef struct _NcmSpectralPanelsProductWork
{
  NcmSpectral *spectral;
  GArray *edges;   /* gdouble, the common refinement of a group */
  GArray *w;       /* gdouble, x^p in the Chebyshev polynomials of a cell, p + 1 */
  GArray *w_next;  /* gdouble, scratch of the expansion of x^p */
  GArray *G;       /* gdouble, the bilinear form of a cell, n1 n2 */
  GArray *A_local; /* gdouble, tables for counts above TABLE_N */
  guint nmax_local;
  GArray *sigma; /* gdouble, per row, k = sigma x */
  GArray *s1;    /* gdouble, per row, the scales of the two sides */
  GArray *s2;
  GArray *sum;   /* gdouble, per row, the integral */
  GArray *res;   /* gdouble, per row, the error estimate */
  GArray *rows1; /* guint, per cell, the panel of each side holding it */
  GArray *rows2;
  NcmMatrix *cell1;
  NcmMatrix *cell2;
} NcmSpectralPanelsProductWork;

/* What a group of rows needs */
typedef struct _NcmSpectralPanelsProductGroup
{
  NcmSpectralPanels *P;
  NcmSpectralPanels *Q;
  gboolean isauto;
  gboolean scaled1; /* whether any scale of a side differs from 1 */
  gboolean scaled2;
  guint r0;
  guint n;
  gdouble c1; /* the variable of each side, x_j = c_j x */
  gdouble c2;
  gboolean want_err;
} NcmSpectralPanelsProductGroup;

G_DEFINE_TYPE (NcmSpectralPanelsProduct, ncm_spectral_panels_product, G_TYPE_OBJECT)

static gpointer _ncm_spectral_panels_product_work_new (gpointer userdata);
static void _ncm_spectral_panels_product_work_free (gpointer p);
static void _ncm_spectral_panels_product_tables_fill (gdouble *A, const guint n_tables, const guint nmax);

static void
ncm_spectral_panels_product_init (NcmSpectralPanelsProduct *prod)
{
  prod->p         = 0;
  prod->tables    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  prod->work_pool = ncm_memory_pool_new (_ncm_spectral_panels_product_work_new, NULL,
                                         _ncm_spectral_panels_product_work_free);
}

static void
_ncm_spectral_panels_product_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectralPanelsProduct *prod = NCM_SPECTRAL_PANELS_PRODUCT (object);

  g_return_if_fail (NCM_IS_SPECTRAL_PANELS_PRODUCT (object));

  switch (prop_id)
  {
    case PROP_P:
      prod->p = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_panels_product_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectralPanelsProduct *prod = NCM_SPECTRAL_PANELS_PRODUCT (object);

  g_return_if_fail (NCM_IS_SPECTRAL_PANELS_PRODUCT (object));

  switch (prop_id)
  {
    case PROP_P:
      g_value_set_uint (value, prod->p);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spectral_panels_product_constructed (GObject *object)
{
  NcmSpectralPanelsProduct *prod = NCM_SPECTRAL_PANELS_PRODUCT (object);
  const guint nmax               = NCM_SPECTRAL_PANELS_PRODUCT_TABLE_N;

  /* Chain up : start */
  G_OBJECT_CLASS (ncm_spectral_panels_product_parent_class)->constructed (object);

  g_array_set_size (prod->tables, (prod->p + 1) * nmax * nmax);
  _ncm_spectral_panels_product_tables_fill (&g_array_index (prod->tables, gdouble, 0), prod->p + 1, nmax);
}

static void
_ncm_spectral_panels_product_finalize (GObject *object)
{
  NcmSpectralPanelsProduct *prod = NCM_SPECTRAL_PANELS_PRODUCT (object);

  g_array_unref (prod->tables);
  ncm_memory_pool_free (prod->work_pool, TRUE);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spectral_panels_product_parent_class)->finalize (object);
}

static void
ncm_spectral_panels_product_class_init (NcmSpectralPanelsProductClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spectral_panels_product_set_property;
  object_class->get_property = &_ncm_spectral_panels_product_get_property;
  object_class->constructed  = &_ncm_spectral_panels_product_constructed;
  object_class->finalize     = &_ncm_spectral_panels_product_finalize;

  /**
   * NcmSpectralPanelsProduct:p:
   *
   * The power of $k$ in the integrand.
   */
  g_object_class_install_property (object_class,
                                   PROP_P,
                                   g_param_spec_uint ("p",
                                                      NULL,
                                                      "Power of k",
                                                      0, 64, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* Half the integral of T_m over [-1, 1]: 1 / (1 - m^2) for even m, 0 for odd m */
static inline gdouble
_ncm_spectral_panels_product_half_int (const guint m)
{
  return (m % 2 == 0) ? 1.0 / (1.0 - (gdouble) m * (gdouble) m) : 0.0;
}

/*
 * The tables A_m[i][j] = int T_i T_j T_m dt over [-1, 1], m below n_tables and i, j below
 * nmax, row stride nmax. From T_i T_j = (T_{i+j} + T_{|i-j|}) / 2 applied twice,
 * A_m[i][j] = (h(s + m) + h(|s - m|) + h(d + m) + h(|d - m|)) / 2 with s = i + j,
 * d = |i - j| and h half the integral of T.
 */
static void
_ncm_spectral_panels_product_tables_fill (gdouble *A, const guint n_tables, const guint nmax)
{
  guint i, j, m;

  for (m = 0; m < n_tables; m++)
  {
    for (i = 0; i < nmax; i++)
    {
      for (j = 0; j < nmax; j++)
      {
        const guint s = i + j;
        const guint d = (i > j) ? i - j : j - i;

        A[(m * nmax + i) * nmax + j] = 0.5 * (_ncm_spectral_panels_product_half_int (s + m) +
                                              _ncm_spectral_panels_product_half_int ((s > m) ? s - m : m - s) +
                                              _ncm_spectral_panels_product_half_int (d + m) +
                                              _ncm_spectral_panels_product_half_int ((d > m) ? d - m : m - d));
      }
    }
  }
}

static gpointer
_ncm_spectral_panels_product_work_new (gpointer userdata)
{
  NcmSpectralPanelsProductWork *work = g_new0 (NcmSpectralPanelsProductWork, 1);

  work->spectral   = ncm_spectral_new ();
  work->edges      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->w          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->w_next     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->G          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->A_local    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->nmax_local = 0;
  work->sigma      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->s1         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->s2         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->sum        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->res        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  work->rows1      = g_array_new (FALSE, FALSE, sizeof (guint));
  work->rows2      = g_array_new (FALSE, FALSE, sizeof (guint));
  work->cell1      = NULL;
  work->cell2      = NULL;

  return work;
}

static void
_ncm_spectral_panels_product_work_free (gpointer p)
{
  NcmSpectralPanelsProductWork *work = p;

  ncm_spectral_free (work->spectral);
  g_array_unref (work->edges);
  g_array_unref (work->w);
  g_array_unref (work->w_next);
  g_array_unref (work->G);
  g_array_unref (work->A_local);
  g_array_unref (work->sigma);
  g_array_unref (work->s1);
  g_array_unref (work->s2);
  g_array_unref (work->sum);
  g_array_unref (work->res);
  g_array_unref (work->rows1);
  g_array_unref (work->rows2);
  ncm_matrix_clear (&work->cell1);
  ncm_matrix_clear (&work->cell2);

  g_free (work);
}

/**
 * ncm_spectral_panels_product_new:
 * @p: power of $k$ in the integrand
 *
 * Creates a product integrator for $\int k^p P Q \, dk$.
 *
 * Returns: (transfer full): a new #NcmSpectralPanelsProduct.
 */
NcmSpectralPanelsProduct *
ncm_spectral_panels_product_new (guint p)
{
  return g_object_new (NCM_TYPE_SPECTRAL_PANELS_PRODUCT,
                       "p", p,
                       NULL);
}

/**
 * ncm_spectral_panels_product_ref:
 * @prod: a #NcmSpectralPanelsProduct
 *
 * Increases the reference count of @prod by one.
 *
 * Returns: (transfer full): @prod.
 */
NcmSpectralPanelsProduct *
ncm_spectral_panels_product_ref (NcmSpectralPanelsProduct *prod)
{
  return g_object_ref (prod);
}

/**
 * ncm_spectral_panels_product_free:
 * @prod: a #NcmSpectralPanelsProduct
 *
 * Decreases the reference count of @prod by one.
 */
void
ncm_spectral_panels_product_free (NcmSpectralPanelsProduct *prod)
{
  g_object_unref (prod);
}

/**
 * ncm_spectral_panels_product_clear:
 * @prod: a #NcmSpectralPanelsProduct
 *
 * Decreases the reference count of *@prod by one and sets *@prod to %NULL.
 */
void
ncm_spectral_panels_product_clear (NcmSpectralPanelsProduct **prod)
{
  g_clear_object (prod);
}

/**
 * ncm_spectral_panels_product_get_p:
 * @prod: a #NcmSpectralPanelsProduct
 *
 * Returns: #NcmSpectralPanelsProduct:p.
 */
guint
ncm_spectral_panels_product_get_p (NcmSpectralPanelsProduct *prod)
{
  return prod->p;
}

static void _ncm_spectral_panels_product_scales (NcmVector *s, guint fdim, GArray *out, gboolean *scaled);
static void _ncm_spectral_panels_product_group (NcmSpectralPanelsProduct *prod, NcmSpectralPanelsProductWork *work, const NcmSpectralPanelsProductGroup *group);

/**
 * ncm_spectral_panels_product_integrate:
 * @prod: a #NcmSpectralPanelsProduct
 * @P: the first panel set
 * @s1: (nullable): the scales of the rows of @P, %NULL for all 1
 * @Q: the second panel set
 * @s2: (nullable): the scales of the rows of @Q, %NULL for all 1
 * @result: the integral of each row
 * @err: (nullable): the error estimate of each row
 *
 * Computes $I_r = \int k^p P_r(k / s_{1r}) Q_r(k / s_{2r}) \, dk$ for every row over the
 * range both panel sets cover, zero when they do not overlap. @P and @Q have the same
 * #NcmSpectralPanels:fdim, the length of @result, @err, @s1 and @s2. When @Q is @P with
 * the same scales the product is computed as a square. A cell of the common refinement
 * outside a single panel of either side is an error.
 */
void
ncm_spectral_panels_product_integrate (NcmSpectralPanelsProduct *prod, NcmSpectralPanels *P, NcmVector *s1,
                                       NcmSpectralPanels *Q, NcmVector *s2, NcmVector *result, NcmVector *err)
{
  const guint fdim                   = ncm_spectral_panels_get_fdim (P);
  NcmSpectralPanelsProductWork **ptr = ncm_memory_pool_get (prod->work_pool);
  NcmSpectralPanelsProductWork *work = *ptr;
  NcmSpectralPanelsProductGroup group;
  gboolean equal = TRUE;
  guint r;

  if (ncm_spectral_panels_get_fdim (Q) != fdim)
    g_error ("ncm_spectral_panels_product_integrate: the panel sets have %u and %u rows.",
             fdim, ncm_spectral_panels_get_fdim (Q));

  if ((ncm_vector_len (result) != fdim) || ((err != NULL) && (ncm_vector_len (err) != fdim)) ||
      ((s1 != NULL) && (ncm_vector_len (s1) != fdim)) || ((s2 != NULL) && (ncm_vector_len (s2) != fdim)))
    g_error ("ncm_spectral_panels_product_integrate: the vectors must have fdim = %u entries.", fdim);

  group.P        = P;
  group.Q        = Q;
  group.want_err = (err != NULL);

  _ncm_spectral_panels_product_scales (s1, fdim, work->s1, &group.scaled1);
  _ncm_spectral_panels_product_scales (s2, fdim, work->s2, &group.scaled2);

  g_array_set_size (work->sigma, fdim);
  g_array_set_size (work->sum, fdim);
  g_array_set_size (work->res, fdim);
  memset (work->sum->data, 0, sizeof (gdouble) * fdim);
  memset (work->res->data, 0, sizeof (gdouble) * fdim);

  for (r = 0; r < fdim; r++)
    equal = equal && (g_array_index (work->s1, gdouble, r) == g_array_index (work->s2, gdouble, r));

  group.isauto = (P == Q) && equal;

  if (equal)
  {
    /* One refinement in x = k / s for every row */
    memcpy (work->sigma->data, work->s1->data, sizeof (gdouble) * fdim);
    group.r0 = 0;
    group.n  = fdim;
    group.c1 = 1.0;
    group.c2 = 1.0;
    _ncm_spectral_panels_product_group (prod, work, &group);
  }
  else
  {
    /* One refinement per row, in the variable of the side whose scale differs from 1, or
     * of the first side when both do */
    const gboolean first = group.scaled1;

    memcpy (work->sigma->data, (first ? work->s1 : work->s2)->data, sizeof (gdouble) * fdim);
    group.n = 1;

    for (r = 0; r < fdim; r++)
    {
      const gdouble sigma = g_array_index (work->sigma, gdouble, r);

      group.r0 = r;
      group.c1 = first ? 1.0 : sigma;
      group.c2 = first ? (group.scaled2 ? sigma / g_array_index (work->s2, gdouble, r) : sigma) : 1.0;
      _ncm_spectral_panels_product_group (prod, work, &group);
    }
  }

  for (r = 0; r < fdim; r++)
    ncm_vector_set (result, r, g_array_index (work->sum, gdouble, r));

  if (err != NULL)
    for (r = 0; r < fdim; r++)
      ncm_vector_set (err, r, g_array_index (work->res, gdouble, r));

  ncm_memory_pool_return (ptr);
}

/* The scales of one side, 1 when s is NULL, and whether any differs from 1 */
static void
_ncm_spectral_panels_product_scales (NcmVector *s, guint fdim, GArray *out, gboolean *scaled)
{
  guint r;

  g_array_set_size (out, fdim);
  *scaled = FALSE;

  for (r = 0; r < fdim; r++)
  {
    const gdouble sr = (s != NULL) ? ncm_vector_get (s, r) : 1.0;

    g_array_index (out, gdouble, r) = sr;
    *scaled                         = *scaled || (sr != 1.0);
  }
}

static void _ncm_spectral_panels_product_edges (const NcmSpectralPanelsProductGroup *group, gdouble x_min, gdouble x_max, GArray *edges);
static gboolean _ncm_spectral_panels_product_restrict (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, NcmMatrix **cell);
static void _ncm_spectral_panels_product_weights (const guint p, const gdouble mid, const gdouble half, GArray *w, GArray *w_next);
static const gdouble *_ncm_spectral_panels_product_tables (NcmSpectralPanelsProduct *prod, NcmSpectralPanelsProductWork *work, const guint n1, const guint n2, guint *nmax);
static gdouble _ncm_spectral_panels_product_pow (const gdouble s, const guint n);
static void _ncm_spectral_panels_product_err (NcmSpectralPanelsProduct *prod, NcmSpectralPanelsProductWork *work, const NcmSpectralPanelsProductGroup *group);

/*
 * Integrates the rows [r0, r0 + n) on the common refinement of the two edge sets in x,
 * with k = sigma_r x and side j in x_j = c_j x. On a cell both sides are polynomials in the
 * cell's variable t, x = mid + half t, and x^p = sum_m w_m T_m(t), so
 * int x^p P Q dx = half a^T G b with G = sum_m w_m A_m, one G for every row of the cell;
 * k^p dk = sigma^(p + 1) x^p dx.
 */
static void
_ncm_spectral_panels_product_group (NcmSpectralPanelsProduct *prod, NcmSpectralPanelsProductWork *work,
                                    const NcmSpectralPanelsProductGroup *group)
{
  const guint p = prod->p;
  gdouble lo1, hi1, lo2, hi2, x_min, x_max;
  guint ie, r;

  if ((ncm_spectral_panels_get_n_panels (group->P) == 0) || (ncm_spectral_panels_get_n_panels (group->Q) == 0))
    return;

  ncm_spectral_panels_get_range (group->P, &lo1, &hi1);
  ncm_spectral_panels_get_range (group->Q, &lo2, &hi2);

  x_min = MAX (lo1 / group->c1, lo2 / group->c2);
  x_max = MIN (hi1 / group->c1, hi2 / group->c2);

  if (x_min >= x_max)
    return;

  _ncm_spectral_panels_product_edges (group, x_min, x_max, work->edges);

  for (ie = 0; ie + 1 < work->edges->len; ie++)
  {
    const gdouble a    = g_array_index (work->edges, gdouble, ie);
    const gdouble b    = g_array_index (work->edges, gdouble, ie + 1);
    const gdouble mid  = 0.5 * (a + b);
    const gdouble half = 0.5 * (b - a);
    NcmMatrix *cm1, *cm2;

    if (!_ncm_spectral_panels_product_restrict (group->P, work->spectral, a * group->c1, b * group->c1, &work->cell1))
      g_error ("ncm_spectral_panels_product_integrate: cell [%.17g, %.17g] is not inside a single panel "
               "of the first panel set.", a, b);

    cm1 = work->cell1;

    if (group->isauto)
    {
      cm2 = cm1;
    }
    else
    {
      if (!_ncm_spectral_panels_product_restrict (group->Q, work->spectral, a * group->c2, b * group->c2, &work->cell2))
        g_error ("ncm_spectral_panels_product_integrate: cell [%.17g, %.17g] is not inside a single panel "
                 "of the second panel set.", a, b);

      cm2 = work->cell2;
    }

    {
      const guint n1     = ncm_matrix_ncols (cm1);
      const guint n2     = ncm_matrix_ncols (cm2);
      const guint tda1   = ncm_matrix_tda (cm1);
      const guint tda2   = ncm_matrix_tda (cm2);
      const gdouble *dc1 = ncm_matrix_data (cm1);
      const gdouble *dc2 = ncm_matrix_data (cm2);
      const gdouble *w;
      const gdouble *A;
      gdouble *G;
      guint nmax, i, j, m;

      _ncm_spectral_panels_product_weights (p, mid, half, work->w, work->w_next);
      A = _ncm_spectral_panels_product_tables (prod, work, n1, n2, &nmax);
      w = &g_array_index (work->w, gdouble, 0);

      g_array_set_size (work->G, n1 * n2);
      G = &g_array_index (work->G, gdouble, 0);

      for (i = 0; i < n1; i++)
      {
        const guint row = i * nmax;

        for (j = 0; j < n2; j++)
        {
          gdouble g = w[0] * A[row + j];

          for (m = 1; m <= p; m++)
            g += w[m] * A[m * nmax * nmax + row + j];

          G[i * n2 + j] = g;
        }
      }

      for (r = group->r0; r < group->r0 + group->n; r++)
      {
        const gdouble s_p1 = _ncm_spectral_panels_product_pow (g_array_index (work->sigma, gdouble, r), p + 1);
        const gdouble *ca  = &dc1[r * tda1];
        const gdouble *cb  = &dc2[r * tda2];
        gdouble acc        = 0.0;

        for (i = 0; i < n1; i++)
        {
          const gdouble *Gi = &G[i * n2];
          gdouble rr        = 0.0;

          for (j = 0; j < n2; j++)
            rr += Gi[j] * cb[j];

          acc += ca[i] * rr;
        }

        g_array_index (work->sum, gdouble, r) += s_p1 * half * acc;
      }
    }
  }

  if (group->want_err)
    _ncm_spectral_panels_product_err (prod, work, group);
}

static gint
_ncm_spectral_panels_product_cmp_gdouble (gconstpointer pa, gconstpointer pb)
{
  const gdouble a = *(const gdouble *) pa;
  const gdouble b = *(const gdouble *) pb;

  return (a > b) - (a < b);
}

/* Appends the interior edges of sp, in x = x_sp / c, that fall strictly inside (x_min, x_max) */
static void
_ncm_spectral_panels_product_append_edges (NcmSpectralPanels *sp, const gdouble c, gdouble x_min, gdouble x_max, GArray *edges)
{
  const guint n_panels = ncm_spectral_panels_get_n_panels (sp);
  guint i;

  for (i = 0; i < n_panels; i++)
  {
    gdouble a, b;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);
    b /= c;

    if ((b > x_min) && (b < x_max))
      g_array_append_val (edges, b);
  }
}

/* The common refinement of the two edge sets on [x_min, x_max], ascending, without duplicates */
static void
_ncm_spectral_panels_product_edges (const NcmSpectralPanelsProductGroup *group, gdouble x_min, gdouble x_max, GArray *edges)
{
  guint w = 1;
  guint i;

  g_array_set_size (edges, 0);
  g_array_append_val (edges, x_min);

  _ncm_spectral_panels_product_append_edges (group->P, group->c1, x_min, x_max, edges);

  if (!group->isauto)
    _ncm_spectral_panels_product_append_edges (group->Q, group->c2, x_min, x_max, edges);

  g_array_sort (edges, _ncm_spectral_panels_product_cmp_gdouble);
  g_array_append_val (edges, x_max);

  for (i = 1; i < edges->len; i++)
    if (g_array_index (edges, gdouble, i) > g_array_index (edges, gdouble, w - 1))
      g_array_index (edges, gdouble, w++) = g_array_index (edges, gdouble, i);

  g_array_set_size (edges, w);
}

/*
 * The coefficients of every row of sp on [a, b], which lies inside one panel. A cell
 * reaches a panel edge through the conversion between the two sides' variables, so an
 * end can land a rounding outside the panel; within that slack it is clamped to the edge.
 */
static gboolean
_ncm_spectral_panels_product_restrict (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, NcmMatrix **cell)
{
  const guint fdim = ncm_spectral_panels_get_fdim (sp);
  const guint i    = ncm_spectral_panels_find_panel (sp, 0.5 * (a + b));
  gdouble pa, pb;
  guint N;
  const gdouble *c    = ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);
  const gdouble slack = 1.0e-12 * (pb - pa);

  if ((a < pa - slack) || (b > pb + slack))
    return FALSE;

  a = MAX (a, pa);
  b = MIN (b, pb);

  if ((*cell != NULL) && ((ncm_matrix_nrows (*cell) != fdim) || (ncm_matrix_ncols (*cell) != N)))
    ncm_matrix_clear (cell);

  if (*cell == NULL)
    *cell = ncm_matrix_new (fdim, N);

  if ((a == pa) && (b == pb))
  {
    guint r;

    for (r = 0; r < fdim; r++)
      memcpy (ncm_matrix_ptr (*cell, r, 0), c + (gsize) r * N, sizeof (gdouble) * N);
  }
  else
  {
    ncm_spectral_chebyshev_rebase_rows_data (spectral, c, fdim, N, N, pa, pb, a, b, *cell);
  }

  return TRUE;
}

/*
 * x^p = sum_m w_m T_m(t) with x = mid + half t. Up to p = 2 the closed forms; above it each
 * multiplication by mid + half T_1 uses T_1 T_0 = T_1 and T_1 T_k = (T_{k-1} + T_{k+1}) / 2.
 */
static void
_ncm_spectral_panels_product_weights (const guint p, const gdouble mid, const gdouble half, GArray *w, GArray *w_next)
{
  gdouble *c;
  guint deg, k;

  g_array_set_size (w, p + 1);
  c = &g_array_index (w, gdouble, 0);

  switch (p)
  {
    case 0:
      c[0] = 1.0;

      return;

    case 1:
      c[0] = mid;
      c[1] = half;

      return;

    default:
      c[0] = mid * mid + 0.5 * half * half;
      c[1] = 2.0 * mid * half;
      c[2] = 0.5 * half * half;
      break;
  }

  g_array_set_size (w_next, p + 1);

  for (deg = 2; deg < p; deg++)
  {
    gdouble *cn = &g_array_index (w_next, gdouble, 0);

    for (k = 0; k <= deg + 1; k++)
    {
      const gdouble c_km1 = (k >= 1) ? c[k - 1] : 0.0;
      const gdouble c_k   = (k <= deg) ? c[k] : 0.0;
      const gdouble c_kp1 = (k + 1 <= deg) ? c[k + 1] : 0.0;
      const gdouble t1c   = (k == 0) ? 0.5 * c_kp1 : ((k == 1) ? c_km1 + 0.5 * c_kp1 : 0.5 * (c_km1 + c_kp1));

      cn[k] = mid * c_k + half * t1c;
    }

    memcpy (c, cn, sizeof (gdouble) * (deg + 2));
  }
}

/* The tables for counts n1 and n2: the object's when both fit, else ones built in the scratch space */
static const gdouble *
_ncm_spectral_panels_product_tables (NcmSpectralPanelsProduct *prod, NcmSpectralPanelsProductWork *work,
                                     const guint n1, const guint n2, guint *nmax)
{
  const guint n = MAX (n1, n2);

  if (n <= NCM_SPECTRAL_PANELS_PRODUCT_TABLE_N)
  {
    *nmax = NCM_SPECTRAL_PANELS_PRODUCT_TABLE_N;

    return &g_array_index (prod->tables, gdouble, 0);
  }

  if (n > work->nmax_local)
  {
    g_array_set_size (work->A_local, (prod->p + 1) * n * n);
    _ncm_spectral_panels_product_tables_fill (&g_array_index (work->A_local, gdouble, 0), prod->p + 1, n);
    work->nmax_local = n;
  }

  *nmax = work->nmax_local;

  return &g_array_index (work->A_local, gdouble, 0);
}

/* s^n, n >= 1, by repeated multiplication */
static gdouble
_ncm_spectral_panels_product_pow (const gdouble s, const guint n)
{
  gdouble v = s;
  guint i;

  for (i = 1; i < n; i++)
    v *= s;

  return v;
}

static void _ncm_spectral_panels_product_cell_rows (NcmSpectralPanels *sp, GArray *edges, const gdouble c, GArray *rows);
static gdouble _ncm_spectral_panels_product_tail (NcmSpectralPanels *sp, guint panel, guint r);
static gdouble _ncm_spectral_panels_product_eval_k (NcmSpectralPanels *sp, guint r, const gdouble k, const gdouble s, const gboolean scaled);

/*
 * Adds the error terms of the rows [r0, r0 + n) on the cells of the group: |P| dQ + |Q| dP
 * at the five Gauss-Legendre nodes of every cell, dP the tail of the panel of P holding the
 * cell. On a square the two terms are one, 2 |P| dP.
 */
static void
_ncm_spectral_panels_product_err (NcmSpectralPanelsProduct *prod, NcmSpectralPanelsProductWork *work,
                                  const NcmSpectralPanelsProductGroup *group)
{
  const guint p = prod->p;
  guint ie, ig, r, m;

  _ncm_spectral_panels_product_cell_rows (group->P, work->edges, group->c1, work->rows1);

  if (!group->isauto)
    _ncm_spectral_panels_product_cell_rows (group->Q, work->edges, group->c2, work->rows2);

  for (ie = 0; ie + 1 < work->edges->len; ie++)
  {
    const gdouble cell_lo = g_array_index (work->edges, gdouble, ie);
    const gdouble cell_hi = g_array_index (work->edges, gdouble, ie + 1);
    const gdouble mid     = 0.5 * (cell_lo + cell_hi);
    const gdouble half    = 0.5 * (cell_hi - cell_lo);
    const guint panel1    = g_array_index (work->rows1, guint, ie);
    const guint panel2    = group->isauto ? panel1 : g_array_index (work->rows2, guint, ie);

    for (r = group->r0; r < group->r0 + group->n; r++)
    {
      const gdouble sigma = g_array_index (work->sigma, gdouble, r);
      const gdouble s_p1  = _ncm_spectral_panels_product_pow (sigma, p + 1);
      const gdouble dW1   = _ncm_spectral_panels_product_tail (group->P, panel1, r);
      const gdouble dW2   = group->isauto ? dW1 : _ncm_spectral_panels_product_tail (group->Q, panel2, r);

      for (ig = 0; ig < NCM_SPECTRAL_PANELS_PRODUCT_GL5_N; ig++)
      {
        const gdouble x     = mid + half * _ncm_spectral_panels_product_gl5_x[ig];
        const gdouble k     = sigma * x;
        const gdouble absW1 = fabs (_ncm_spectral_panels_product_eval_k (group->P, r, k, g_array_index (work->s1, gdouble, r), group->scaled1));
        gdouble w           = half * _ncm_spectral_panels_product_gl5_w[ig];

        for (m = 0; m < p; m++)
          w *= x;

        if (group->isauto)
        {
          g_array_index (work->res, gdouble, r) += 2.0 * s_p1 * w * absW1 * dW1;
        }
        else
        {
          const gdouble absW2 = fabs (_ncm_spectral_panels_product_eval_k (group->Q, r, k, g_array_index (work->s2, gdouble, r), group->scaled2));

          g_array_index (work->res, gdouble, r) += s_p1 * w * (absW1 * dW2 + absW2 * dW1);
        }
      }
    }
  }
}

/* For each cell, the panel of sp holding it, by one marching index over both sorted sets */
static void
_ncm_spectral_panels_product_cell_rows (NcmSpectralPanels *sp, GArray *edges, const gdouble c, GArray *rows)
{
  const guint n_panels = ncm_spectral_panels_get_n_panels (sp);
  guint ie, j = 0;

  g_array_set_size (rows, 0);

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble cell_lo = c * g_array_index (edges, gdouble, ie);

    while (j + 1 < n_panels)
    {
      gdouble a, b;
      guint N;

      ncm_spectral_panels_peek_panel (sp, j, &a, &b, &N);

      if (b > cell_lo)
        break;

      j++;
    }

    g_array_append_val (rows, j);
  }
}

/* The sum of the absolute coefficients of row r of a panel above the order of the level below */
static gdouble
_ncm_spectral_panels_product_tail (NcmSpectralPanels *sp, guint panel, guint r)
{
  gdouble a, b, tail = 0.0;
  guint N, j;
  const gdouble *c   = ncm_spectral_panels_peek_panel (sp, panel, &a, &b, &N);
  const gdouble *row = c + (gsize) r * N;

  for (j = (N + 1) / 2; j < N; j++)
    tail += fabs (row[j]);

  return tail;
}

/* Row r of sp at k, in its variable k / s when the side is scaled; zero outside its range */
static gdouble
_ncm_spectral_panels_product_eval_k (NcmSpectralPanels *sp, guint r, const gdouble k, const gdouble s, const gboolean scaled)
{
  const gdouble x = scaled ? k / s : k;
  gdouble lo, hi;

  ncm_spectral_panels_get_range (sp, &lo, &hi);

  if ((x < lo) || (x > hi))
    return 0.0;

  return ncm_spectral_panels_eval_row (sp, r, x);
}

