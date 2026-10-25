/***************************************************************************
 *            ncm_stats_rows.c
 *
 *  Fri Oct 24 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_rows.c
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
 * NcmStatsRows:
 *
 * Analyses of a sample of rows, each a #NcmVector, given with ncm_stats_rows_set_rows():
 * the burn-in estimates ncm_stats_rows_max_ess_time(), ncm_stats_rows_heidel_diag() and
 * ncm_stats_rows_estimate_const_break(), and the robust covariances
 * ncm_stats_rows_cov_robust_diag() and ncm_stats_rows_cov_robust_ogk(). The work buffers
 * sized by the number of columns are allocated once. Those sized by the number of rows are
 * allocated by the first analysis that needs them and freed, with the reference to the
 * sample, by ncm_stats_rows_reset(), so the intended use is ncm_stats_rows_set_rows(), the
 * analyses, then ncm_stats_rows_reset().
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_rows.h"
#include "ncm/stats/ncm_stats_vec.h"
#include "ncm/stats/ncm_stats_acorr.h"
#include "ncm/algebra/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsRows
{
  GObject parent_instance;
  guint nrows;
  guint ncols;
  GPtrArray *sample;
  /* Work buffers sized by ncols. */
  NcmVector *esss;
  NcmVector *Ivals;
  NcmVector *spec0;
  NcmVector *cumsum;
  NcmVector *sigma_x;
  NcmVector *sigma_z;
  NcmVector *c;
  NcmMatrix *E;
  NcmMatrix *c_cov;
  NcmStatsAcorr *acorr;
  NcmStatsAcorr *acorr1;
  NcmStatsVec *chunk;
  GArray *ar_order;

  /* Work buffers sized by nrows, allocated by the first analysis that needs them and freed
   * by ncm_stats_rows_reset(); NULL when not allocated. */
  NcmVector *data;
  NcmVector *work;
  GArray *work_int;
  NcmMatrix *y;
  NcmMatrix *z;
  NcmVector *ones;
  NcmVector *yfit;
};

enum
{
  PROP_0,
  PROP_NCOLS,
};

G_DEFINE_TYPE (NcmStatsRows, ncm_stats_rows, G_TYPE_OBJECT)

/* Frees the buffers sized by nrows. */
static void
_ncm_stats_rows_free_row_buffers (NcmStatsRows *rows)
{
  ncm_vector_clear (&rows->data);
  ncm_vector_clear (&rows->work);
  g_clear_pointer (&rows->work_int, g_array_unref);
  ncm_matrix_clear (&rows->y);
  ncm_matrix_clear (&rows->z);
  ncm_vector_clear (&rows->ones);
  ncm_vector_clear (&rows->yfit);
}

static void
ncm_stats_rows_init (NcmStatsRows *rows)
{
  rows->nrows    = 0;
  rows->ncols    = 0;
  rows->sample   = NULL;
  rows->esss     = NULL;
  rows->Ivals    = NULL;
  rows->spec0    = NULL;
  rows->cumsum   = NULL;
  rows->sigma_x  = NULL;
  rows->sigma_z  = NULL;
  rows->c        = NULL;
  rows->E        = NULL;
  rows->c_cov    = NULL;
  rows->acorr    = NULL;
  rows->acorr1   = NULL;
  rows->chunk    = NULL;
  rows->ar_order = g_array_new (FALSE, FALSE, sizeof (guint));
  rows->data     = NULL;
  rows->work     = NULL;
  rows->work_int = NULL;
  rows->y        = NULL;
  rows->z        = NULL;
  rows->ones     = NULL;
  rows->yfit     = NULL;
}

static void
_ncm_stats_rows_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsRows *rows = NCM_STATS_ROWS (object);

  g_return_if_fail (NCM_IS_STATS_ROWS (object));

  switch (prop_id)
  {
    case PROP_NCOLS:
      rows->ncols = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_rows_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsRows *rows = NCM_STATS_ROWS (object);

  g_return_if_fail (NCM_IS_STATS_ROWS (object));

  switch (prop_id)
  {
    case PROP_NCOLS:
      g_value_set_uint (value, rows->ncols);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* The buffers sized by ncols; ncols is construct-only. */
static void
_ncm_stats_rows_constructed (GObject *object)
{
  G_OBJECT_CLASS (ncm_stats_rows_parent_class)->constructed (object);
  {
    NcmStatsRows *rows = NCM_STATS_ROWS (object);
    const guint ncols  = rows->ncols;

    g_assert_cmpuint (ncols, >, 0);

    rows->esss    = ncm_vector_new (ncols);
    rows->Ivals   = ncm_vector_new (ncols);
    rows->spec0   = ncm_vector_new (ncols);
    rows->cumsum  = ncm_vector_new (ncols);
    rows->sigma_x = ncm_vector_new (ncols);
    rows->sigma_z = ncm_vector_new (ncols);
    rows->c       = ncm_vector_new (1);
    rows->E       = ncm_matrix_new (ncols, ncols);
    rows->c_cov   = ncm_matrix_new (1, 1);
    rows->acorr   = ncm_stats_acorr_new_full (ncols, NCM_STATS_ACORR_DEFAULT_MAX_LAG, NCM_STATS_ACORR_DEFAULT_MAX_LEVELS, NCM_STATS_ACORR_METHOD_AR);
    rows->acorr1  = ncm_stats_acorr_new_full (1, NCM_STATS_ACORR_DEFAULT_MAX_LAG, NCM_STATS_ACORR_DEFAULT_MAX_LEVELS, NCM_STATS_ACORR_METHOD_AR);
    rows->chunk   = ncm_stats_vec_new (ncols, NCM_STATS_VEC_VAR, FALSE);

    g_array_set_size (rows->ar_order, ncols);
  }
}

static void
_ncm_stats_rows_dispose (GObject *object)
{
  NcmStatsRows *rows = NCM_STATS_ROWS (object);

  g_clear_pointer (&rows->sample, g_ptr_array_unref);
  ncm_vector_clear (&rows->esss);
  ncm_vector_clear (&rows->Ivals);
  ncm_vector_clear (&rows->spec0);
  ncm_vector_clear (&rows->cumsum);
  ncm_vector_clear (&rows->sigma_x);
  ncm_vector_clear (&rows->sigma_z);
  ncm_vector_clear (&rows->c);
  ncm_matrix_clear (&rows->E);
  ncm_matrix_clear (&rows->c_cov);
  ncm_stats_acorr_clear (&rows->acorr);
  ncm_stats_acorr_clear (&rows->acorr1);
  ncm_stats_vec_clear (&rows->chunk);

  _ncm_stats_rows_free_row_buffers (rows);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_rows_parent_class)->dispose (object);
}

static void
_ncm_stats_rows_finalize (GObject *object)
{
  NcmStatsRows *rows = NCM_STATS_ROWS (object);

  g_clear_pointer (&rows->ar_order, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_rows_parent_class)->finalize (object);
}

static void
ncm_stats_rows_class_init (NcmStatsRowsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_stats_rows_set_property;
  object_class->get_property = &_ncm_stats_rows_get_property;
  object_class->constructed  = &_ncm_stats_rows_constructed;
  object_class->dispose      = &_ncm_stats_rows_dispose;
  object_class->finalize     = &_ncm_stats_rows_finalize;

  /**
   * NcmStatsRows:ncols:
   *
   * The number of columns of the sample, the length of each row.
   */
  g_object_class_install_property (object_class,
                                   PROP_NCOLS,
                                   g_param_spec_uint ("ncols",
                                                      NULL,
                                                      "Number of columns",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_stats_rows_new:
 * @ncols: number of columns of the sample
 *
 * Creates a #NcmStatsRows for samples of @ncols columns. The sample is given with
 * ncm_stats_rows_set_rows().
 *
 * Returns: (transfer full): a new #NcmStatsRows.
 */
NcmStatsRows *
ncm_stats_rows_new (guint ncols)
{
  NcmStatsRows *rows = g_object_new (NCM_TYPE_STATS_ROWS,
                                     "ncols", ncols,
                                     NULL);

  return rows;
}

/**
 * ncm_stats_rows_ref:
 * @rows: a #NcmStatsRows
 *
 * Increases the reference count of @rows by one.
 *
 * Returns: (transfer full): @rows.
 */
NcmStatsRows *
ncm_stats_rows_ref (NcmStatsRows *rows)
{
  return g_object_ref (rows);
}

/**
 * ncm_stats_rows_free:
 * @rows: a #NcmStatsRows
 *
 * Decreases the reference count of @rows by one.
 */
void
ncm_stats_rows_free (NcmStatsRows *rows)
{
  g_object_unref (rows);
}

/**
 * ncm_stats_rows_clear:
 * @rows: a #NcmStatsRows
 *
 * Decreases the reference count of *@rows by one and sets *@rows to %NULL.
 */
void
ncm_stats_rows_clear (NcmStatsRows **rows)
{
  g_clear_object (rows);
}

/**
 * ncm_stats_rows_set_rows:
 * @rows: a #NcmStatsRows
 * @sample: (element-type NcmVector): the rows, each a #NcmVector of #NcmStatsRows:ncols elements
 *
 * Sets the sample, which @rows references until ncm_stats_rows_reset() or the next call, and
 * the number of rows to its length; see ncm_stats_rows_set_nrows() to analyze fewer rows.
 */
void
ncm_stats_rows_set_rows (NcmStatsRows *rows, GPtrArray *sample)
{
  g_clear_pointer (&rows->sample, g_ptr_array_unref);
  rows->sample = g_ptr_array_ref (sample);

  ncm_stats_rows_set_nrows (rows, sample->len);
}

/**
 * ncm_stats_rows_set_nrows:
 * @rows: a #NcmStatsRows
 * @nrows: number of rows to analyze, the first @nrows of the sample
 *
 * Sets the number of rows the analyses use, at most the length of the sample.
 */
void
ncm_stats_rows_set_nrows (NcmStatsRows *rows, guint nrows)
{
  g_assert (rows->sample != NULL);
  g_assert_cmpuint (nrows, <=, rows->sample->len);

  rows->nrows = nrows;
}

/**
 * ncm_stats_rows_reset:
 * @rows: a #NcmStatsRows
 *
 * Releases the reference to the sample, frees the work buffers sized by the number of rows,
 * and sets the number of rows to zero. The buffers sized by the number of columns are kept.
 */
void
ncm_stats_rows_reset (NcmStatsRows *rows)
{
  g_clear_pointer (&rows->sample, g_ptr_array_unref);
  rows->nrows = 0;

  _ncm_stats_rows_free_row_buffers (rows);
}

/**
 * ncm_stats_rows_get_nrows:
 * @rows: a #NcmStatsRows
 *
 * Returns: the number of rows the analyses use.
 */
guint
ncm_stats_rows_get_nrows (NcmStatsRows *rows)
{
  return rows->nrows;
}

/**
 * ncm_stats_rows_get_ncols:
 * @rows: a #NcmStatsRows
 *
 * Returns: the number of columns of the sample.
 */
guint
ncm_stats_rows_get_ncols (NcmStatsRows *rows)
{
  return rows->ncols;
}

/**
 * ncm_stats_rows_peek_row:
 * @rows: a #NcmStatsRows
 * @i: row index, below ncm_stats_rows_get_nrows()
 *
 * Returns: (transfer none): row @i of the sample.
 */
NcmVector *
ncm_stats_rows_peek_row (NcmStatsRows *rows, guint i)
{
  g_assert_cmpuint (i, <, rows->nrows);

  return g_ptr_array_index (rows->sample, i);
}

#define _ROW(rows, i) ((NcmVector *) g_ptr_array_index ((rows)->sample, (i)))

/* Allocates *v with at least n elements when it is missing or shorter. */
static void
_ncm_stats_rows_alloc_vector (NcmVector **v, guint n)
{
  if ((*v == NULL) || (ncm_vector_len (*v) < n))
  {
    ncm_vector_clear (v);
    *v = ncm_vector_new (n);
  }
}

/* Allocates *m as an nrows x ncols matrix when it is missing or of another shape. */
static void
_ncm_stats_rows_alloc_matrix (NcmMatrix **m, guint nrows, guint ncols)
{
  if ((*m == NULL) || (ncm_matrix_nrows (*m) != nrows) || (ncm_matrix_ncols (*m) != ncols))
  {
    ncm_matrix_clear (m);
    *m = ncm_matrix_new (nrows, ncols);
  }
}

/**
 * ncm_stats_rows_max_ess_time:
 * @rows: a #NcmStatsRows
 * @ntests: number of tests
 * @bindex: (out): time index of the best ESS's
 * @wp: (out): worst parameter index
 * @wp_order: (out): worst parameter AR fit order
 * @wp_ess: (out): worst parameter ESS
 *
 * Finds the starting row that maximizes the smallest effective sample size (ESS) over the
 * columns, computed from that row to the last with an AR fit, testing @ntests starting
 * rows (10 when @ntests is zero). Requires at least 10 rows.
 *
 * Returns: (transfer full): a #NcmVector containing the best ess.
 */
NcmVector *
ncm_stats_rows_max_ess_time (NcmStatsRows *rows, const guint ntests, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_ess)
{
  const gint size   = rows->nrows;
  const gint block  = (ntests == 0) ? ((size - 1) / 10 + 1) : ((size - 1) / (gint) ntests + 1);
  NcmVector *esss   = ncm_vector_new (rows->ncols);
  gdouble max_t_ess = 0.0;
  gint i, j = 0;

  g_assert_cmpuint (rows->nrows, >=, 10);

  ncm_stats_acorr_reset (rows->acorr);

  bindex[0] = -1;

  for (i = size - 1; i >= 0; i--)
  {
    ncm_stats_acorr_update (rows->acorr, _ROW (rows, i));

    if ((i == 0) || ((i % block == 0) && (j >= 99)))
    {
      gdouble min_ess = GSL_POSINF;
      guint cur_size  = size - i;
      guint lwp_order = 0;
      guint k, lwp    = 0;

      for (k = 0; k < rows->ncols; k++)
      {
        const gdouble ess   = ncm_stats_acorr_get_ess (rows->acorr, k);
        const gdouble c_ess = GSL_MIN (cur_size, ess);

        lwp_order = ncm_stats_acorr_get_ar_order (rows->acorr, k);

        ncm_vector_set (rows->esss, k, ess);

        if (c_ess < min_ess)
        {
          min_ess = c_ess;
          lwp     = k;
        }
      }

      if (min_ess >= max_t_ess)
      {
        max_t_ess   = min_ess;
        bindex[0]   = i;
        wp_order[0] = lwp_order;
        wp[0]       = lwp;

        ncm_vector_memcpy (esss, rows->esss);
      }
    }

    j++;
  }

  wp_ess[0] = ncm_vector_get (esss, wp[0]);

  return esss;
}

static gdouble
_ncm_stats_rows_heidel_diag_pcramer (const gdouble q)
{
  const gdouble pi_32 = sqrt (gsl_pow_3 (M_PI));
  const guint maxiter = 100;
  const gdouble ffac  = 1.0 / (pi_32 * sqrt (q));
  gdouble p           = 0.0;
  guint i             = 0;

  g_assert_cmpfloat (q, >=, 0.0);

  for (i = 0; i < maxiter; i++)
  {
    gint sig          = 0;
    const gdouble lnf = lgamma_r (i + 0.5, &sig) - lgamma_r (i + 1.0, &sig);
    const gdouble z   = sqrt (4.0 * i + 1.0);
    const gdouble u   = gsl_pow_2 (4.0 * i + 1.0) / (16.0 * q);
    const gdouble ti  = z * exp (-u + lnf + gsl_sf_bessel_lnKnu (0.25, u));

    p += ti;

    if (fabs (ti / p) < GSL_DBL_EPSILON)
      break;
  }

  return GSL_MIN (GSL_MAX (p * ffac, 0.0), 1.0);
}

/**
 * ncm_stats_rows_heidel_diag:
 * @rows: a #NcmStatsRows
 * @ntests: number of tests
 * @pvalue: required p-value
 * @bindex: (out): index of the best p-values
 * @wp: (out): worst parameter index
 * @wp_order: (out): worst parameter AR fit order
 * @wp_pvalue: (out): value of the worst parameter
 *
 * Applies the Heidelberger--Welch convergence diagnostic with @ntests
 * sequential Schruben tests. Uses 10 tests when @ntests is zero and a
 * p-value of $0.05$ when @pvalue is zero.
 *
 * Each test applies the Cramér-von Mises statistic of the Brownian bridge of the rows from
 * a starting index to the last, with the spectral density at zero estimated by an AR fit
 * to the second half of the rows. The values reported are the Cramér-von Mises
 * cumulative distribution at the statistic, that is one minus the p-values. Sets @bindex
 * to the smallest starting index at which every value is at most $1 -$ @pvalue, and to -1
 * when no index qualifies. The returned vector contains the values at @bindex, or for the
 * full sample when no index qualifies; @wp, @wp_order and @wp_pvalue identify the
 * column with the largest value there, its AR order and its value. Requires at least 10
 * rows.
 *
 * See:
 *
 * - [Heidelberger (1981)](https://doi.org/10.1145/358598.358630)
 * - [Schruben (1982)](https://doi.org/10.1287/opre.30.3.569)
 * - [Heidelberger (1983)](https://doi.org/10.1287/opre.31.6.1109)
 *
 * Returns: (transfer full): the Cramér-von Mises cumulative distribution values.
 */
NcmVector *
ncm_stats_rows_heidel_diag (NcmStatsRows *rows, const guint ntests, const gdouble pvalue, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_pvalue)
{
  const guint nrows    = rows->nrows;
  const guint ncols    = rows->ncols;
  const gint half_size = nrows / 2;
  const gint block     = (ntests == 0) ? ((half_size - 1) / 10 + 1) : ((half_size - 1) / (gint) ntests + 1);
  const gdouble onepv  = (pvalue == 0.0) ? 0.95 : (1.0 - pvalue);
  NcmVector *pvals     = ncm_vector_new (ncols);
  gint i;

  g_assert_cmpuint (nrows, >=, 10);
  g_assert_cmpfloat (pvalue, <, 1.0);

  ncm_stats_vec_reset (rows->chunk, TRUE);
  ncm_stats_acorr_reset (rows->acorr);

  for (i = nrows - 1; i >= half_size; i--)
  {
    ncm_stats_vec_append (rows->chunk, _ROW (rows, i), FALSE);
    ncm_stats_acorr_update (rows->acorr, _ROW (rows, i));
  }

  for (i = 0; i < (gint) ncols; i++)
  {
    ncm_vector_set (rows->spec0, i, ncm_stats_acorr_get_spec0 (rows->acorr, i));
    g_array_index (rows->ar_order, guint, i) = ncm_stats_acorr_get_ar_order (rows->acorr, i);
  }

  bindex[0] = -1;
  wp[0]     = 0;

  for (i = half_size - 1; i >= 0; i--)
  {
    ncm_stats_vec_append (rows->chunk, _ROW (rows, i), FALSE);

    if ((i % block) == 0)
    {
      const guint nitens = nrows - i;
      gint j;
      guint p;

      ncm_vector_set_zero (rows->cumsum);
      ncm_vector_set_zero (rows->Ivals);

      for (j = nrows - 1; j >= i; j--)
      {
        NcmVector *row  = _ROW (rows, j);
        const gdouble n = (nrows - j);

        for (p = 0; p < ncols; p++)
        {
          const gdouble mean_p_n = ncm_stats_vec_get_mean (rows->chunk, p) * n;
          const gdouble cumsum_p = ncm_vector_get (rows->cumsum, p) + ncm_vector_get (row, p);

          ncm_vector_set (rows->cumsum, p, cumsum_p);
          ncm_vector_addto (rows->Ivals, p, gsl_pow_2 (cumsum_p - mean_p_n));
        }
      }

      {
        gdouble max_pval = 0.0;
        gint lwp         = 0;

        for (p = 0; p < ncols; p++)
        {
          const gdouble Ival_p = ncm_vector_get (rows->Ivals, p) / (gsl_pow_2 (nitens) * ncm_vector_get (rows->spec0, p));
          const gdouble pval_p = _ncm_stats_rows_heidel_diag_pcramer (Ival_p);

          ncm_vector_set (rows->Ivals, p, pval_p);

          if (pval_p > max_pval)
          {
            max_pval = pval_p;
            lwp      = p;
          }
        }

        if (max_pval <= onepv)
        {
          bindex[0] = i;
          wp[0]     = lwp;

          ncm_vector_memcpy (pvals, rows->Ivals);
        }
      }
    }
  }

  if (bindex[0] == -1)
  {
    ncm_vector_memcpy (pvals, rows->Ivals);
    wp[0] = ncm_vector_get_max_index (pvals);
  }

  wp_pvalue[0] = ncm_vector_get (pvals, wp[0]);
  wp_order[0]  = g_array_index (rows->ar_order, guint, wp[0]);

  return pvals;
}

/**
 * ncm_stats_rows_visual_heidel_diag:
 * @rows: a #NcmStatsRows
 * @p: column index
 * @fi: first index
 * @mean: (out): mean
 * @var: (out): test's variance
 *
 * Computes, for column @p and the rows from the last down to @fi, the cumulative sums in
 * that order, their mean and the variance used by ncm_stats_rows_heidel_diag().
 *
 * Returns: (transfer full): the cumulative sums.
 */
NcmVector *
ncm_stats_rows_visual_heidel_diag (NcmStatsRows *rows, const guint p, const guint fi, gdouble *mean, gdouble *var)
{
  const guint nitens  = rows->nrows - fi;
  gdouble cumsum      = 0.0;
  NcmVector *cumsum_v = ncm_vector_new (nitens);
  gint i, j = 0;

  g_assert_cmpuint (rows->nrows, >=, 10);
  g_assert_cmpuint (fi, <, rows->nrows);

  ncm_stats_acorr_reset (rows->acorr1);

  for (i = rows->nrows - 1; i >= (gint) fi; i--)
  {
    const gdouble p_val = ncm_vector_get (_ROW (rows, i), p);

    cumsum += p_val;
    ncm_vector_set (cumsum_v, j, cumsum);
    j++;

    ncm_stats_acorr_update_var (rows->acorr1, 0, p_val);
  }

  mean[0] = ncm_stats_acorr_get_mean (rows->acorr1, 0);
  var[0]  = ncm_stats_acorr_get_spec0 (rows->acorr1, 0) * nitens;

  return cumsum_v;
}

/* First row at or after pad within the robust cutoff of the robust mean of column p over the
 * rows from pad to the last, counted from pad. The GSL robust workspace is sized to the
 * number of rows fitted, so it is allocated for each call. */
static guint
_ncm_stats_rows_estimate_const_break_int (NcmStatsRows *rows, guint p, guint pad)
{
  g_assert_cmpuint (pad, <, rows->nrows);
  {
    const guint n                     = rows->nrows - pad;
    const gsl_multifit_robust_type *T = gsl_multifit_robust_default;
    gsl_multifit_robust_workspace *w  = gsl_multifit_robust_alloc (T, n, 1);
    gsl_multifit_robust_stats stats;
    gsl_matrix_view X;
    gsl_vector_view y;
    gint status;
    gdouble t0, cutoff;
    guint i;

    for (i = 0; i < n; i++)
    {
      ncm_vector_set (rows->yfit, i, ncm_vector_get (_ROW (rows, i + pad), p));
      ncm_vector_set (rows->ones, i, 1.0);
    }

    X = gsl_matrix_view_array (ncm_vector_data (rows->ones), n, 1);
    y = gsl_vector_view_array (ncm_vector_data (rows->yfit), n);

    gsl_multifit_robust_maxiter (100000, w);
    status = gsl_multifit_robust (&X.matrix, &y.vector, ncm_vector_gsl (rows->c), ncm_matrix_gsl (rows->c_cov), w);

    if ((status != GSL_SUCCESS) && (status != GSL_EMAXITER))
      g_error ("_ncm_stats_rows_estimate_const_break_int: error %d computing gsl_multifit_robust\n", status);

    stats = gsl_multifit_robust_statistics (w);

    t0 = ncm_vector_get (rows->c, 0);

    cutoff = ceil (sqrt (gsl_cdf_chisq_Qinv (1.0 / n, 1.0)));

    for (i = 0; i < n; i++)
    {
      if (fabs ((ncm_vector_get (rows->yfit, i) - t0) / stats.sigma_rob) < cutoff)
        break;
    }

    gsl_multifit_robust_free (w);

    return i;
  }
}

/**
 * ncm_stats_rows_estimate_const_break:
 * @rows: a #NcmStatsRows
 * @p: column index
 *
 * Estimates the mean $\mu$ and standard deviation $\sigma$ of column @p
 * with robust regression and returns the first index $t_0$ within
 * $\alpha\sigma$ of $\mu$, where $\alpha$ is $\sqrt{x}$ rounded up, with $x$ the value
 * exceeded with probability $1/N$ by a $\chi^2_1$ variable, and $N$ the size of the sample.
 * The robust regression is repeated on the rows after each cut until no row is cut.
 *
 * Returns: $t_0$
 */
gdouble
ncm_stats_rows_estimate_const_break (NcmStatsRows *rows, guint p)
{
  const guint nrows = rows->nrows;
  guint t0          = 0;
  guint t1          = 0;

  _ncm_stats_rows_alloc_vector (&rows->yfit, nrows);
  _ncm_stats_rows_alloc_vector (&rows->ones, nrows);

  do {
    t1  = _ncm_stats_rows_estimate_const_break_int (rows, p, t0);
    t0 += t1;

    if (t0 >= nrows)
    {
      t0 = nrows - 1;
      break;
    }
  } while (t1 > 0);

  return t0;
}

/* Sizes the scratch of _ncm_stats_rows_qn() and of the eigenproblem of the OGK estimate. */
static void
_ncm_stats_rows_size_qn (NcmStatsRows *rows)
{
  const guint n_int = MAX (rows->nrows * 5, 2 * rows->ncols);

  _ncm_stats_rows_alloc_vector (&rows->data, MAX (rows->nrows, rows->ncols));
  _ncm_stats_rows_alloc_vector (&rows->work, rows->nrows * 3);

  if (rows->work_int == NULL)
    rows->work_int = g_array_new (FALSE, FALSE, sizeof (gint));

  if (rows->work_int->len < n_int)
    g_array_set_size (rows->work_int, n_int);
}

/* Qn scale of the column held in data, sorted in place. */
static gdouble
_ncm_stats_rows_qn (NcmStatsRows *rows)
{
  gsl_sort (ncm_vector_data (rows->data), 1, rows->nrows);

  return gsl_stats_Qn_from_sorted_data (ncm_vector_data (rows->data),
                                        1, rows->nrows,
                                        ncm_vector_data (rows->work),
                                        &g_array_index (rows->work_int, gint, 0));
}

/**
 * ncm_stats_rows_cov_robust_diag:
 * @rows: a #NcmStatsRows
 *
 * Estimates the variance of each column as the square of the Qn scale estimator of
 * Rousseeuw and Croux. Requires at least 4 rows.
 *
 * Returns: (transfer full): a diagonal matrix of the variances.
 */
NcmMatrix *
ncm_stats_rows_cov_robust_diag (NcmStatsRows *rows)
{
  const guint nrows = rows->nrows;
  const guint ncols = rows->ncols;
  NcmMatrix *cov    = ncm_matrix_new (ncols, ncols);
  guint i;

  if (nrows < 4)
    g_error ("ncm_stats_rows_cov_robust_diag: too few points to estimate the covariance [%u].", nrows);

  _ncm_stats_rows_size_qn (rows);
  ncm_matrix_set_zero (cov);

  for (i = 0; i < ncols; i++)
  {
    guint a;

    for (a = 0; a < nrows; a++)
      ncm_vector_fast_set (rows->data, a, ncm_vector_get (_ROW (rows, a), i));

    ncm_matrix_set (cov, i, i, gsl_pow_2 (_ncm_stats_rows_qn (rows)));
  }

  return cov;
}

/**
 * ncm_stats_rows_cov_robust_ogk:
 * @rows: a #NcmStatsRows
 *
 * Estimates the covariance by the orthogonalized Gnanadesikan-Kettenring (OGK) method of
 * Maronna and Zamar (2002), with the Qn scale estimator. Requires at least 4 rows.
 *
 * Returns: (transfer full): the covariance matrix.
 */
NcmMatrix *
ncm_stats_rows_cov_robust_ogk (NcmStatsRows *rows)
{
  const guint nrows = rows->nrows;
  const guint ncols = rows->ncols;
  NcmMatrix *cov    = ncm_matrix_new (ncols, ncols);
  NcmMatrix *E      = rows->E;
  NcmMatrix *y, *z;
  guint a, i;

  if (nrows < 4)
    g_error ("ncm_stats_rows_cov_robust_ogk: too few points to estimate the covariance [%u].", nrows);

  _ncm_stats_rows_size_qn (rows);
  _ncm_stats_rows_alloc_matrix (&rows->y, nrows, ncols);
  _ncm_stats_rows_alloc_matrix (&rows->z, ncols, nrows);
  y = rows->y;
  z = rows->z;

  /* y holds every row; each column is scaled after its scale is read from it unscaled. */
  for (a = 0; a < nrows; a++)
    ncm_matrix_set_row (y, a, _ROW (rows, a));

  for (i = 0; i < ncols; i++)
  {
    gdouble sigma_i;

    for (a = 0; a < nrows; a++)
      ncm_vector_fast_set (rows->data, a, ncm_matrix_get (y, a, i));

    sigma_i = _ncm_stats_rows_qn (rows);

    ncm_vector_set (rows->sigma_x, i, sigma_i);
    ncm_matrix_mul_col (y, i, 1.0 / sigma_i);
  }

  ncm_matrix_set_identity (cov);

  for (i = 0; i < ncols; i++)
  {
    guint j;

    for (j = i + 1; j < ncols; j++)
    {
      gdouble s_ipj, s_imj;

      for (a = 0; a < nrows; a++)
        ncm_vector_fast_set (rows->data, a, ncm_matrix_get (y, a, i) + ncm_matrix_get (y, a, j));

      s_ipj = _ncm_stats_rows_qn (rows);

      for (a = 0; a < nrows; a++)
        ncm_vector_fast_set (rows->data, a, ncm_matrix_get (y, a, i) - ncm_matrix_get (y, a, j));

      s_imj = _ncm_stats_rows_qn (rows);

      ncm_matrix_set (cov, i, j, 0.25 * (s_ipj * s_ipj - s_imj * s_imj));
    }
  }

  {
    NcmLapackWS *lapack_work = ncm_lapack_ws_new ();
    gint neval;
    gint info;

    info = ncm_lapack_dsyevr ('V', 'A', 'U', ncols,
                              ncm_matrix_data (cov), ncm_matrix_tda (cov),
                              0.0, 0.0, 0.0, 0.0,
                              0.0, &neval, ncm_vector_data (rows->data),
                              ncm_matrix_data (E), ncm_matrix_tda (E),
                              &g_array_index (rows->work_int, gint, 0),
                              lapack_work);

    NCM_LAPACK_CHECK_INFO ("dsyevr", info);

    ncm_lapack_ws_free (lapack_work);
  }

  ncm_matrix_dgemm (z, 'N', 'T', 1.0, E, y, 0.0);

  for (i = 0; i < ncols; i++)
  {
    for (a = 0; a < nrows; a++)
      ncm_vector_fast_set (rows->data, a, ncm_matrix_get (z, i, a));

    ncm_vector_set (rows->sigma_z, i, _ncm_stats_rows_qn (rows));
  }

  for (i = 0; i < ncols; i++)
  {
    const gdouble sigma_z_i = ncm_vector_get (rows->sigma_z, i);
    guint j;

    for (j = 0; j < ncols; j++)
    {
      const gdouble sigma_x_j = ncm_vector_get (rows->sigma_x, j);
      const gdouble E_ij      = ncm_matrix_get (E, i, j);

      ncm_matrix_set (E, i, j, sigma_z_i * E_ij * sigma_x_j);
    }
  }

  /* dsyrk writes one triangle; the caller gets a full covariance matrix. */
  ncm_matrix_dsyrk (cov, 'U', 'T', 1.0, E, 0.0);
  ncm_matrix_copy_triangle (cov, 'U');

  return cov;
}

