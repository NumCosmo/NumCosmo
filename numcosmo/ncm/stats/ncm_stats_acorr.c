/***************************************************************************
 *            ncm_stats_acorr.c
 *
 *  Tue Sep 16 09:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_acorr.c
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
 * NcmStatsAcorr:
 *
 * Integrated autocorrelation time $\tau$ of one or more scalar series, updated one
 * sample at a time.
 *
 * For each series the object accumulates the lagged autocovariances $C_k$ up to
 * #NcmStatsAcorr:max-lag, and repeats that accumulation on successive levels of block
 * averages, level $j$ seeing the means of $2^j$ consecutive samples. Each update costs
 * $O(\mathrm{max\text{-}lag})$ per series and the stored state is
 * $O(\mathrm{max\text{-}lag}\log n)$: the samples are not retained. The accumulated $C_k$
 * are exact, in the sense that feeding a series one sample at a time and computing its
 * autocovariances in one pass over the whole series give the same numbers to rounding.
 *
 * Three estimators turn $C_k$ into $\tau$ — an auto-regressive spectral estimate with the
 * order chosen by AICc, Geyer's initial monotone positive sequence, and Sokal's
 * self-consistent window — selected by #NcmStatsAcorr:method. They are pure functions of
 * $C_k$, and %NCM_STATS_ACORR_METHOD_MAX reports the largest and flags a disagreement. The
 * level whose blocks resolve the correlation within the available lags is selected
 * automatically, which is what removes the ceiling a fixed maximum lag would otherwise put
 * on $\tau$. Every estimate carries a #NcmStatsAcorrDiag saying whether it is to be
 * trusted.
 *
 * Error modes: a series with no variance never moved, so its correlation time is unbounded
 * and $\tau$ is reported at its cap, the number of samples, with
 * %NCM_STATS_ACORR_DIAG_ZERO_VARIANCE, which makes the effective sample size one; an estimator whose truncated sum is not positive
 * gives $\tau = 1$; $\tau$ is capped at the number of samples, since a longer correlation
 * is not measurable from the series. None of these abort.
 *
 * See <a href="../../theory/autocorrelation.html">Autocorrelation Time and Effective
 * Sample Size</a> for the definitions, the identity the accumulator updates, the
 * level-selection rule and the meaning of each condition.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_acorr.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmStatsAcorrLevel
{
  guint64 n;
  GArray *R;
  GArray *head;
  GArray *ring;
  gdouble S;
  gdouble pending_sum;
  guint pending_n;
} NcmStatsAcorrLevel;

typedef struct _NcmStatsAcorrVar
{
  gdouble offset;
  gboolean started;
  GPtrArray *levels;
  NcmVector *acov;
  gboolean dirty;
  gdouble mean;
  gdouble var;
  gdouble tau;
  gdouble tau_ar;
  gdouble tau_geyer;
  gdouble drift_z;
  gdouble var_ratio;
  guint level;
  guint window;
  guint ar_order;
  NcmStatsAcorrDiag diag;
} NcmStatsAcorrVar;

typedef struct _NcmStatsAcorrPrivate
{
  guint len;
  guint max_lag;
  guint max_levels;
  NcmStatsAcorrMethod method;
  NcmStatsAcorrARCrit ar_crit;
  gdouble reliability_factor;
  gdouble drift_threshold;
  NcmStatsAcorrVar *vars;
} NcmStatsAcorrPrivate;

enum
{
  PROP_0,
  PROP_LEN,
  PROP_MAX_LAG,
  PROP_MAX_LEVELS,
  PROP_METHOD,
  PROP_AR_CRITERION,
  PROP_RELIABILITY_FACTOR,
  PROP_DRIFT_THRESHOLD,
};

struct _NcmStatsAcorr
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsAcorr, ncm_stats_acorr, G_TYPE_OBJECT);

/* Minimum number of block values a level must hold to be considered for selection. */
#define _NCM_STATS_ACORR_MIN_BLOCKS (16)
/* A level resolves the correlation when its tau fits this fraction of the lag budget. */
#define _NCM_STATS_ACORR_LAG_MARGIN (5.0)
/* Recenter the accumulator when the mean leaves this many standard deviations. */
#define _NCM_STATS_ACORR_RECENTER_SD (8.0)

/*
 * Level accumulator.
 *
 * ring[i mod max_lag] holds the i-th value of this level, so the value at lag k behind
 * the last one is ring[(n - k) mod max_lag]; head[i] holds the first max_lag values.
 * R[k] is the uncentered lagged sum of the values with the variable offset already
 * subtracted, which is what keeps the centered autocovariance well conditioned.
 */

static NcmStatsAcorrLevel *
_ncm_stats_acorr_level_new (guint max_lag)
{
  NcmStatsAcorrLevel *lev = g_new0 (NcmStatsAcorrLevel, 1);

  lev->R    = g_array_new (FALSE, TRUE, sizeof (gdouble));
  lev->head = g_array_new (FALSE, TRUE, sizeof (gdouble));
  lev->ring = g_array_new (FALSE, TRUE, sizeof (gdouble));

  g_array_set_size (lev->R,    max_lag + 1);
  g_array_set_size (lev->head, max_lag);
  g_array_set_size (lev->ring, max_lag);

  return lev;
}

static void
_ncm_stats_acorr_level_free (gpointer ptr)
{
  NcmStatsAcorrLevel *lev = ptr;

  g_clear_pointer (&lev->R,    g_array_unref);
  g_clear_pointer (&lev->head, g_array_unref);
  g_clear_pointer (&lev->ring, g_array_unref);
  g_free (lev);
}

static void
_ncm_stats_acorr_level_update (NcmStatsAcorrLevel *lev, guint max_lag, const gdouble y)
{
  gdouble * const R    = &g_array_index (lev->R,    gdouble, 0);
  gdouble * const head = &g_array_index (lev->head, gdouble, 0);
  gdouble * const ring = &g_array_index (lev->ring, gdouble, 0);
  const guint64 n      = lev->n;
  const guint kmax     = (n < max_lag) ? (guint) n : max_lag;
  const guint pos      = (guint) (n % max_lag);
  const guint kwrap    = (kmax < pos) ? kmax : pos;
  guint k;

  /* The value at lag k sits at ring[pos - k], wrapping once at k = pos: two contiguous
   * runs rather than a division per lag. */
  for (k = 1; k <= kwrap; k++)
    R[k] += y * ring[pos - k];

  for (k = kwrap + 1; k <= kmax; k++)
    R[k] += y * ring[max_lag + pos - k];

  R[0]   += y * y;
  lev->S += y;

  if (n < max_lag)
    head[n] = y;

  ring[n % max_lag] = y;
  lev->n            = n + 1;
}

/*
 * Suffix sums T_k of the last k values and prefix sums H_k of the first k values, both
 * needed to center the uncentered lagged sums exactly.
 */
static void
_ncm_stats_acorr_level_edge_sums (NcmStatsAcorrLevel *lev, guint max_lag, guint nlag, gdouble *T, gdouble *H)
{
  const gdouble * const head = &g_array_index (lev->head, gdouble, 0);
  const gdouble * const ring = &g_array_index (lev->ring, gdouble, 0);
  const guint pos            = (guint) (lev->n % max_lag);
  const guint kwrap          = (nlag < pos) ? nlag : pos;
  guint k;

  T[0] = 0.0;
  H[0] = 0.0;

  for (k = 1; k <= kwrap; k++)
    T[k] = T[k - 1] + ring[pos - k];

  for (k = kwrap + 1; k <= nlag; k++)
    T[k] = T[k - 1] + ring[max_lag + pos - k];

  for (k = 1; k <= nlag; k++)
    H[k] = H[k - 1] + head[k - 1];
}

/*
 * Centered autocovariances of this level, C_k for k = 0 .. nlag, normalized by the
 * number of values (the biased estimator, which keeps the sequence positive definite).
 */
static void
_ncm_stats_acorr_level_acov (NcmStatsAcorrLevel *lev, guint max_lag, NcmVector *acov)
{
  const gdouble * const R = &g_array_index (lev->R, gdouble, 0);
  const guint64 n         = lev->n;
  const guint nlag        = ncm_vector_len (acov) - 1;
  const gdouble m         = lev->S / (1.0 * n);
  GArray *edges           = g_array_new (FALSE, TRUE, sizeof (gdouble));
  gdouble *T, *H;
  guint k;

  g_assert_cmpuint (nlag, <, n);

  g_array_set_size (edges, 2 * (nlag + 1));
  T = &g_array_index (edges, gdouble, 0);
  H = &g_array_index (edges, gdouble, nlag + 1);

  _ncm_stats_acorr_level_edge_sums (lev, max_lag, nlag, T, H);

  for (k = 0; k <= nlag; k++)
  {
    const gdouble Ck = R[k] - m * (2.0 * lev->S - T[k] - H[k]) + (n - k) * m * m;

    ncm_vector_set (acov, k, Ck / (1.0 * n));
  }

  g_array_unref (edges);
}

/*
 * Move the offset by d. Exact: the same algebra that centers R_k re-centers it.
 */
static void
_ncm_stats_acorr_level_recenter (NcmStatsAcorrLevel *lev, guint max_lag, const gdouble d)
{
  gdouble * const R    = &g_array_index (lev->R,    gdouble, 0);
  gdouble * const head = &g_array_index (lev->head, gdouble, 0);
  gdouble * const ring = &g_array_index (lev->ring, gdouble, 0);
  const guint64 n      = lev->n;
  const guint nlag     = (n <= max_lag) ? (guint) (n - 1) : max_lag;
  const guint nstored  = (n < max_lag) ? (guint) n : max_lag;
  GArray *edges        = NULL;
  gdouble *T, *H;
  guint k;

  /* A level is created and fed in the same call, so it always holds at least one value. */
  g_assert_cmpuint (n, >, 0);

  edges = g_array_new (FALSE, TRUE, sizeof (gdouble));
  g_array_set_size (edges, 2 * (nlag + 1));
  T = &g_array_index (edges, gdouble, 0);
  H = &g_array_index (edges, gdouble, nlag + 1);

  _ncm_stats_acorr_level_edge_sums (lev, max_lag, nlag, T, H);

  for (k = 0; k <= nlag; k++)
    R[k] += -d * (2.0 * lev->S - T[k] - H[k]) + (n - k) * d * d;

  lev->S -= n * d;

  for (k = 0; k < nstored; k++)
  {
    head[k] -= d;
    ring[k] -= d;
  }

  lev->pending_sum -= lev->pending_n * d;

  g_array_unref (edges);
}

static void
_ncm_stats_acorr_var_init (NcmStatsAcorrVar *var, guint max_lag)
{
  var->offset  = 0.0;
  var->started = FALSE;
  var->levels  = g_ptr_array_new_with_free_func (_ncm_stats_acorr_level_free);
  var->acov    = ncm_vector_new (max_lag + 1);
  var->dirty   = TRUE;

  g_ptr_array_add (var->levels, _ncm_stats_acorr_level_new (max_lag));
}

static void
_ncm_stats_acorr_var_clear (NcmStatsAcorrVar *var)
{
  g_clear_pointer (&var->levels, g_ptr_array_unref);
  ncm_vector_clear (&var->acov);
}

static void
_ncm_stats_acorr_var_reset (NcmStatsAcorrVar *var, guint max_lag)
{
  g_ptr_array_set_size (var->levels, 0);
  g_ptr_array_add (var->levels, _ncm_stats_acorr_level_new (max_lag));

  var->offset  = 0.0;
  var->started = FALSE;
  var->dirty   = TRUE;
}

static void
_ncm_stats_acorr_var_add (NcmStatsAcorrVar *var, guint max_lag, guint max_levels, guint level, const gdouble y)
{
  NcmStatsAcorrLevel *lev;

  if (level >= max_levels)
    return;

  while (var->levels->len <= level)
    g_ptr_array_add (var->levels, _ncm_stats_acorr_level_new (max_lag));

  lev = g_ptr_array_index (var->levels, level);

  _ncm_stats_acorr_level_update (lev, max_lag, y);

  lev->pending_sum += y;
  lev->pending_n++;

  if (lev->pending_n == 2)
  {
    const gdouble block = 0.5 * lev->pending_sum;

    lev->pending_sum = 0.0;
    lev->pending_n   = 0;

    _ncm_stats_acorr_var_add (var, max_lag, max_levels, level + 1, block);
  }
}

static void
_ncm_stats_acorr_var_update (NcmStatsAcorrVar *var, guint max_lag, guint max_levels, const gdouble x)
{
  NcmStatsAcorrLevel *lev0;

  if (!var->started)
  {
    var->offset  = x;
    var->started = TRUE;
  }

  _ncm_stats_acorr_var_add (var, max_lag, max_levels, 0, x - var->offset);

  lev0 = g_ptr_array_index (var->levels, 0);

  if (lev0->n >= _NCM_STATS_ACORR_MIN_BLOCKS)
  {
    const gdouble m  = lev0->S / (1.0 * lev0->n);
    const gdouble R0 = g_array_index (lev0->R, gdouble, 0);
    const gdouble v  = R0 / (1.0 * lev0->n) - m * m;

    if (fabs (m) > _NCM_STATS_ACORR_RECENTER_SD * sqrt (fabs (v)))
    {
      guint j;

      for (j = 0; j < var->levels->len; j++)
        _ncm_stats_acorr_level_recenter (g_ptr_array_index (var->levels, j), max_lag, m);

      var->offset += m;
    }
  }

  var->dirty = TRUE;
}

static void
ncm_stats_acorr_init (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  self->len                = 0;
  self->max_lag            = 0;
  self->max_levels         = 0;
  self->method             = NCM_STATS_ACORR_METHOD_MAX;
  self->ar_crit            = NCM_STATS_ACORR_DEFAULT_AR_CRIT;
  self->reliability_factor = NCM_STATS_ACORR_DEFAULT_RELIABILITY_FACTOR;
  self->drift_threshold    = NCM_STATS_ACORR_DEFAULT_DRIFT_THRESHOLD;
  self->vars               = NULL;
}

static void
_ncm_stats_acorr_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_stats_acorr_parent_class)->constructed (object);
  {
    NcmStatsAcorr *acorr              = NCM_STATS_ACORR (object);
    NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
    guint p;

    g_assert_cmpuint (self->len, >, 0);
    g_assert_cmpuint (self->max_lag, >, 3);
    g_assert_cmpuint (self->max_levels, >, 0);

    self->vars = g_new0 (NcmStatsAcorrVar, self->len);

    for (p = 0; p < self->len; p++)
      _ncm_stats_acorr_var_init (&self->vars[p], self->max_lag);
  }
}

static void
_ncm_stats_acorr_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsAcorr *acorr              = NCM_STATS_ACORR (object);
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_return_if_fail (NCM_IS_STATS_ACORR (object));

  switch (prop_id)
  {
    case PROP_LEN:
      self->len = g_value_get_uint (value);
      break;
    case PROP_MAX_LAG:
      self->max_lag = g_value_get_uint (value);
      break;
    case PROP_MAX_LEVELS:
      self->max_levels = g_value_get_uint (value);
      break;
    case PROP_METHOD:
      ncm_stats_acorr_set_method (acorr, g_value_get_enum (value));
      break;
    case PROP_AR_CRITERION:
      ncm_stats_acorr_set_ar_criterion (acorr, g_value_get_enum (value));
      break;
    case PROP_RELIABILITY_FACTOR:
      ncm_stats_acorr_set_reliability_factor (acorr, g_value_get_double (value));
      break;
    case PROP_DRIFT_THRESHOLD:
      ncm_stats_acorr_set_drift_threshold (acorr, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_acorr_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsAcorr *acorr              = NCM_STATS_ACORR (object);
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_return_if_fail (NCM_IS_STATS_ACORR (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, self->len);
      break;
    case PROP_MAX_LAG:
      g_value_set_uint (value, self->max_lag);
      break;
    case PROP_MAX_LEVELS:
      g_value_set_uint (value, self->max_levels);
      break;
    case PROP_METHOD:
      g_value_set_enum (value, self->method);
      break;
    case PROP_AR_CRITERION:
      g_value_set_enum (value, self->ar_crit);
      break;
    case PROP_RELIABILITY_FACTOR:
      g_value_set_double (value, self->reliability_factor);
      break;
    case PROP_DRIFT_THRESHOLD:
      g_value_set_double (value, self->drift_threshold);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_acorr_finalize (GObject *object)
{
  NcmStatsAcorr *acorr              = NCM_STATS_ACORR (object);
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  if (self->vars != NULL)
  {
    guint p;

    for (p = 0; p < self->len; p++)
      _ncm_stats_acorr_var_clear (&self->vars[p]);

    g_clear_pointer (&self->vars, g_free);
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_acorr_parent_class)->finalize (object);
}

static void
ncm_stats_acorr_class_init (NcmStatsAcorrClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_stats_acorr_constructed;
  object_class->set_property = &_ncm_stats_acorr_set_property;
  object_class->get_property = &_ncm_stats_acorr_get_property;
  object_class->finalize     = &_ncm_stats_acorr_finalize;

  /**
   * NcmStatsAcorr:len:
   *
   * Number of series tracked.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      NULL,
                                                      "Number of series",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcmStatsAcorr:max-lag:
   *
   * Number of lags accumulated at each level. It bounds the correlation each level can
   * resolve, not the reported $\tau$: a correlation longer than that is resolved at a
   * coarser level.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_LAG,
                                   g_param_spec_uint ("max-lag",
                                                      NULL,
                                                      "Number of lags per level",
                                                      4, G_MAXUINT, NCM_STATS_ACORR_DEFAULT_MAX_LAG,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcmStatsAcorr:max-levels:
   *
   * Maximum number of block-averaging levels. Level $j$ averages $2^j$ consecutive
   * samples, so the longest resolvable $\tau$ is of order
   * $2^{\mathrm{max\text{-}levels}-1}\,\mathrm{max\text{-}lag}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_LEVELS,
                                   g_param_spec_uint ("max-levels",
                                                      NULL,
                                                      "Maximum number of block-averaging levels",
                                                      1, 64, NCM_STATS_ACORR_DEFAULT_MAX_LEVELS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcmStatsAcorr:method:
   *
   * Estimator used to turn the autocovariances into $\tau$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Estimator",
                                                      NCM_TYPE_STATS_ACORR_METHOD, NCM_STATS_ACORR_METHOD_MAX,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcmStatsAcorr:ar-criterion:
   *
   * Rule that picks the order of the auto-regressive fit.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_AR_CRITERION,
                                   g_param_spec_enum ("ar-criterion",
                                                      NULL,
                                                      "Auto-regressive order criterion",
                                                      NCM_TYPE_STATS_ACORR_AR_CRIT, NCM_STATS_ACORR_DEFAULT_AR_CRIT,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcmStatsAcorr:reliability-factor:
   *
   * A series shorter than this many autocorrelation times gets
   * %NCM_STATS_ACORR_DIAG_SHORT_CHAIN.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELIABILITY_FACTOR,
                                   g_param_spec_double ("reliability-factor",
                                                        NULL,
                                                        "Required number of autocorrelation times",
                                                        1.0, G_MAXDOUBLE, NCM_STATS_ACORR_DEFAULT_RELIABILITY_FACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcmStatsAcorr:drift-threshold:
   *
   * Number of standard errors between the means of the first and of the second half of
   * the series above which %NCM_STATS_ACORR_DIAG_DRIFT is set.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DRIFT_THRESHOLD,
                                   g_param_spec_double ("drift-threshold",
                                                        NULL,
                                                        "Drift z-score threshold",
                                                        0.0, G_MAXDOUBLE, NCM_STATS_ACORR_DEFAULT_DRIFT_THRESHOLD,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

/*
 * One Levinson-Durbin pass up to @pmax. @crit picks the order among those tried. Returns
 * the innovation variance of the selected order and sets the order and the sum of its
 * coefficients; @phi_out and @pacf_out, when given, receive the coefficients of the
 * selected order and the reflection coefficients of every order tried.
 */
static gdouble
_ncm_stats_acorr_ar_fit_pass (NcmVector *acov, const gdouble n, const guint pmax,
                              NcmStatsAcorrARCrit crit, guint *order, gdouble *phisum,
                              GArray *phi_out, GArray *pacf_out)
{
  GArray *phi_arr  = g_array_new (FALSE, TRUE, sizeof (gdouble));
  GArray *best_phi = g_array_new (FALSE, TRUE, sizeof (gdouble));
  const gdouble C0 = ncm_vector_get (acov, 0);
  gdouble *phi, *phi_prev;
  gdouble v, best_v, best_crit, best_sum;
  guint m, best_order;

  g_array_set_size (phi_arr,  2 * (pmax + 1));
  g_array_set_size (best_phi, pmax + 1);
  phi      = &g_array_index (phi_arr, gdouble, 0);
  phi_prev = &g_array_index (phi_arr, gdouble, pmax + 1);

  if (pacf_out != NULL)
    g_array_set_size (pacf_out, 0);

  v          = C0;
  best_v     = C0;
  best_order = 0;
  best_sum   = 0.0;

  switch (crit)
  {
    case NCM_STATS_ACORR_AR_CRIT_NONE:
      best_crit = GSL_POSINF;
      break;
    case NCM_STATS_ACORR_AR_CRIT_FPE:
      best_crit = C0;
      break;
    case NCM_STATS_ACORR_AR_CRIT_AIC:
      best_crit = n * log (C0) + 2.0;
      break;
    case NCM_STATS_ACORR_AR_CRIT_AICC:
      best_crit = n * log (C0) + 2.0 * n / (n - 2.0);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  for (m = 1; m <= pmax; m++)
  {
    gdouble acc = ncm_vector_get (acov, m);
    gdouble k, crit_m, sum;
    guint j;

    for (j = 1; j < m; j++)
      acc -= phi_prev[j] * ncm_vector_get (acov, m - j);

    k = acc / v;

    if (!gsl_finite (k) || (fabs (k) >= 1.0))
      break;

    if (pacf_out != NULL)
      g_array_append_val (pacf_out, k);

    for (j = 1; j < m; j++)
      phi[j] = phi_prev[j] - k * phi_prev[m - j];

    phi[m] = k;
    v     *= 1.0 - k * k;

    if (v <= 0.0) /* LCOV_EXCL_BR_LINE: |k| < 1 keeps v positive; underflow guard only */
      break;  /* LCOV_EXCL_LINE */

    sum = 0.0;

    for (j = 1; j <= m; j++)
      sum += phi[j];

    switch (crit)
    {
      case NCM_STATS_ACORR_AR_CRIT_NONE:
        crit_m = -1.0 * m;
        break;
      case NCM_STATS_ACORR_AR_CRIT_FPE:
        crit_m = v * (n + m) / (n - m);
        break;
      case NCM_STATS_ACORR_AR_CRIT_AIC:
        crit_m = n * log (v) + 2.0 * (m + 1.0);
        break;
      case NCM_STATS_ACORR_AR_CRIT_AICC:
        crit_m = n * log (v) + 2.0 * n * (m + 1.0) / (n - m - 2.0);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    if (crit_m < best_crit)
    {
      best_crit  = crit_m;
      best_v     = v;
      best_order = m;
      best_sum   = sum;

      if (phi_out != NULL)
        memcpy (&g_array_index (best_phi, gdouble, 0), phi, sizeof (gdouble) * (pmax + 1));
    }

    memcpy (phi_prev, phi, sizeof (gdouble) * (pmax + 1));
  }

  if (phi_out != NULL)
  {
    guint j;

    g_array_set_size (phi_out, 0);

    for (j = 1; j <= best_order; j++)
      g_array_append_val (phi_out, g_array_index (best_phi, gdouble, j));
  }

  g_array_unref (phi_arr);
  g_array_unref (best_phi);

  order[0]  = best_order;
  phisum[0] = best_sum;

  return best_v;
}

/*
 * The order is searched up to floor (10 log10 n); whenever the criterion selects that
 * bound the search is repeated with the bound doubled, since a criterion that picks the
 * largest order offered has not been given enough of them. NONE always selects the bound,
 * so it is not escalated.
 */
static gdouble
_ncm_stats_acorr_ar_fit_full (NcmVector *acov, const gdouble n, NcmStatsAcorrARCrit crit,
                              guint *order, gdouble *phisum, GArray *phi_out, GArray *pacf_out)
{
  const guint nlag  = ncm_vector_len (acov) - 1;
  const gdouble lgn = floor (10.0 * log10 (n));
  const guint pcap  = (guint) GSL_MAX (GSL_MIN (1.0 * nlag, n - 4.0), 1.0);
  guint pmax        = (guint) GSL_MAX (GSL_MIN (lgn, 1.0 * pcap), 1.0);
  gdouble best_v;

  while (TRUE)
  {
    best_v = _ncm_stats_acorr_ar_fit_pass (acov, n, pmax, crit, order, phisum, phi_out, pacf_out);

    if ((crit == NCM_STATS_ACORR_AR_CRIT_NONE) || (order[0] < pmax) || (pmax >= pcap))
      break;

    pmax = (guint) GSL_MIN (2.0 * pmax, 1.0 * pcap);
  }

  return best_v;
}

/**
 * ncm_stats_acorr_tau_ar:
 * @acov: autocovariances $C_0 \dots C_L$
 * @nitens: number of samples the autocovariances were computed from
 * @crit: a #NcmStatsAcorrARCrit
 * @ar_order: (out) (optional): order selected
 *
 * Integrated autocorrelation time from an auto-regressive fit of @acov by the
 * Levinson-Durbin recursion, the order chosen by @crit, as $\tau = S(0) / C_0$ with $S(0)$
 * the spectral density of the fitted model at zero frequency. See
 * <a href="../../theory/autocorrelation.html">Autocorrelation Time and Effective Sample
 * Size</a>.
 *
 * Returns: $\tau$, capped at @nitens.
 */
gdouble
ncm_stats_acorr_tau_ar (NcmVector *acov, guint64 nitens, NcmStatsAcorrARCrit crit, guint *ar_order)
{
  const gdouble C0 = ncm_vector_get (acov, 0);
  const gdouble n  = 1.0 * nitens;
  gdouble best_v, best_sum = 0.0;
  guint best_order = 0;

  g_assert_cmpuint (crit, <, NCM_STATS_ACORR_AR_CRIT_LEN);

  if (ar_order != NULL)
    ar_order[0] = 0;

  if ((C0 <= 0.0) || (nitens < 8))
    return 1.0;

  best_v = _ncm_stats_acorr_ar_fit_full (acov, n, crit, &best_order, &best_sum, NULL, NULL);

  if (ar_order != NULL)
    ar_order[0] = best_order;

  {
    const gdouble ivar  = best_v * (n - 1.0) / (n - (best_order + 1.0));
    const gdouble denom = gsl_pow_2 (1.0 - best_sum);
    const gdouble tau   = ivar / (denom * C0);

    if (!gsl_finite (tau)) /* LCOV_EXCL_BR_LINE: the fitted model at the stationarity boundary */
      return n;  /* LCOV_EXCL_LINE */

    if (tau <= 0.0) /* LCOV_EXCL_BR_LINE: ivar, denom and C0 are positive at this point */
      return 1.0;  /* LCOV_EXCL_LINE */

    return GSL_MIN (tau, n);
  }
}

/**
 * ncm_stats_acorr_ar_fit:
 * @acov: autocovariances $C_0 \dots C_L$
 * @nitens: number of samples the autocovariances were computed from
 * @crit: a #NcmStatsAcorrARCrit
 * @phi: (out) (transfer full) (optional): coefficients of the selected order
 * @pacf: (out) (transfer full) (optional): reflection coefficients of every order tried
 * @ivar: (out) (optional): innovation variance of the selected order
 * @order: (out) (optional): order selected
 *
 * Fits an auto-regressive model to @acov and returns what the fit is made of, rather than
 * the $\tau$ that ncm_stats_acorr_tau_ar() builds from it. @phi holds
 * $\phi_{p,1} \dots \phi_{p,p}$ and @pacf the reflection coefficients $\kappa_m$, whose
 * decay with $m$ is what an order-selection rule is reading.
 *
 * Returns: TRUE when the selected order is not zero.
 */
gboolean
ncm_stats_acorr_ar_fit (NcmVector *acov, guint64 nitens, NcmStatsAcorrARCrit crit,
                        NcmVector **phi, NcmVector **pacf, gdouble *ivar, guint *order)
{
  const gdouble C0 = ncm_vector_get (acov, 0);
  const gdouble n  = 1.0 * nitens;
  GArray *phi_a    = g_array_new (FALSE, TRUE, sizeof (gdouble));
  GArray *pacf_a   = g_array_new (FALSE, TRUE, sizeof (gdouble));
  gdouble best_v = C0, best_sum = 0.0;
  guint best_order = 0;

  g_assert_cmpuint (crit, <, NCM_STATS_ACORR_AR_CRIT_LEN);

  if ((C0 > 0.0) && (nitens >= 8))
    best_v = _ncm_stats_acorr_ar_fit_full (acov, n, crit, &best_order, &best_sum, phi_a, pacf_a);

  if (phi != NULL)
    *phi = (phi_a->len > 0) ? ncm_vector_new_array (phi_a) : NULL;

  if (pacf != NULL)
    *pacf = (pacf_a->len > 0) ? ncm_vector_new_array (pacf_a) : NULL;

  if (ivar != NULL)
    ivar[0] = (best_order > 0) ? best_v * (n - 1.0) / (n - (best_order + 1.0)) : best_v;

  if (order != NULL)
    order[0] = best_order;

  g_array_unref (phi_a);
  g_array_unref (pacf_a);

  return (best_order > 0);
}

/**
 * ncm_stats_acorr_tau_geyer:
 * @acov: autocovariances $C_0 \dots C_L$
 * @window: (out) (optional): number of lags summed
 *
 * Integrated autocorrelation time of @acov by Geyer's initial monotone positive
 * sequence: the lag pairs $\Gamma_k = C_{2k} + C_{2k+1}$ are summed while they stay
 * positive, after being made non-increasing.
 *
 * The estimate is conservative for a reversible chain: its expectation is not below the
 * true $\tau$, which is what makes it the companion of the auto-regressive estimate in
 * %NCM_STATS_ACORR_METHOD_MAX.
 *
 * Returns: $\tau$.
 */
gdouble
ncm_stats_acorr_tau_geyer (NcmVector *acov, guint *window)
{
  const guint nlag = ncm_vector_len (acov) - 1;
  const gdouble C0 = ncm_vector_get (acov, 0);
  gdouble sum      = 0.0;
  gdouble prev     = GSL_POSINF;
  guint k;

  if (window != NULL)
    window[0] = 0;

  if (C0 <= 0.0)
    return 1.0;

  for (k = 0; 2 * k + 1 <= nlag; k++)
  {
    gdouble G = ncm_vector_get (acov, 2 * k) + ncm_vector_get (acov, 2 * k + 1);

    if (G <= 0.0)
      break;

    G    = GSL_MIN (G, prev);
    prev = G;

    sum += G;

    if (window != NULL)
      window[0] = 2 * k + 2;
  }

  {
    const gdouble tau = (2.0 * sum - C0) / C0;

    if (!gsl_finite (tau) || (tau <= 0.0))
      return 1.0;

    return tau;
  }
}

/**
 * ncm_stats_acorr_tau_sokal:
 * @acov: autocovariances $C_0 \dots C_L$
 * @c: window factor, %NCM_STATS_ACORR_SOKAL_C is the usual choice
 * @window: (out) (optional): number of lags summed
 *
 * Integrated autocorrelation time of @acov summed up to the smallest window $M$
 * satisfying $M \geq c\,\tau(M)$. If no window satisfies it the whole sequence is
 * summed.
 *
 * The rule assumes a non-negative autocorrelation function. On a sequence whose first
 * lag is negative enough ($\rho_1 \leq -0.4$) the window closes at $M = 1$ and the
 * estimate is $1 + 2\rho_1$, which can fall to zero; milder negative lags widen the window
 * instead. Either way the variance reduction such a sequence carries is not measured;
 * %NCM_STATS_ACORR_METHOD_GEYER and %NCM_STATS_ACORR_METHOD_AR are exact there.
 *
 * Returns: $\tau$.
 */
gdouble
ncm_stats_acorr_tau_sokal (NcmVector *acov, gdouble c, guint *window)
{
  const guint nlag = ncm_vector_len (acov) - 1;
  const gdouble C0 = ncm_vector_get (acov, 0);
  gdouble sum      = 0.0;
  gdouble tau      = 1.0;
  guint M;

  if (window != NULL)
    window[0] = 0;

  if (C0 <= 0.0)
    return 1.0;

  for (M = 1; M <= nlag; M++)
  {
    sum += ncm_vector_get (acov, M) / C0;
    tau  = 1.0 + 2.0 * sum;

    if (window != NULL)
      window[0] = M;

    if (M >= c * tau)
      break;
  }

  if (!gsl_finite (tau) || (tau <= 0.0))
    return 1.0;

  return tau;
}

/**
 * ncm_stats_acorr_acov_fft:
 * @series: a #NcmVector holding the series
 * @max_lag: highest lag required, 0 for all of them
 *
 * Centered autocovariances of @series computed by zero-padded Fourier transform,
 * normalized by the number of samples. This is the $O(n\log n)$ path used when a whole
 * series is available at once; the accumulator computes the same numbers one sample at
 * a time.
 *
 * Returns: (transfer full): a #NcmVector with $C_0 \dots C_L$.
 */
NcmVector *
ncm_stats_acorr_acov_fft (NcmVector *series, guint max_lag)
{
  const guint n            = ncm_vector_len (series);
  const guint nlag         = ((max_lag == 0) || (max_lag > n - 1)) ? n - 1 : max_lag;
  const guint effsize      = ncm_util_fact_size (2 * n);
  guint fftw_default_flags = ncm_cfg_get_fftw_default_flag ();
  NcmVector *acov          = ncm_vector_new (nlag + 1);
  fftw_complex *fft        = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * (effsize / 2 + 1));
  gdouble *data            = (gdouble *) fftw_malloc (sizeof (gdouble) * effsize);
  gdouble mean             = 0.0;
  fftw_plan r2c, c2r;
  guint i;

  g_assert_cmpuint (n, >, 1);

  for (i = 0; i < n; i++)
    mean += ncm_vector_get (series, i);

  mean = mean / (1.0 * n);

  for (i = 0; i < n; i++)
    data[i] = ncm_vector_get (series, i) - mean;

  memset (&data[n], 0, sizeof (gdouble) * (effsize - n));

  ncm_cfg_load_fftw_wisdom ("ncm_stats_acorr_%u", effsize);
  ncm_cfg_lock_plan_fftw ();
  r2c = fftw_plan_dft_r2c_1d (effsize, data, fft, fftw_default_flags | FFTW_DESTROY_INPUT);
  c2r = fftw_plan_dft_c2r_1d (effsize, fft, data, fftw_default_flags | FFTW_DESTROY_INPUT);
  ncm_cfg_unlock_plan_fftw ();
  ncm_cfg_save_fftw_wisdom ("ncm_stats_acorr_%u", effsize);

  for (i = 0; i < n; i++)
    data[i] = ncm_vector_get (series, i) - mean;

  memset (&data[n], 0, sizeof (gdouble) * (effsize - n));

  fftw_execute (r2c);

  for (i = 0; i < effsize / 2 + 1; i++)
    fft[i] = fft[i] * conj (fft[i]);

  fftw_execute (c2r);

  for (i = 0; i <= nlag; i++)
    ncm_vector_set (acov, i, data[i] / (1.0 * effsize * n));

  ncm_cfg_lock_plan_fftw ();
  fftw_destroy_plan (r2c);
  fftw_destroy_plan (c2r);
  ncm_cfg_unlock_plan_fftw ();

  fftw_free (fft);
  fftw_free (data);

  return acov;
}

/*
 * Half-to-half drift z-score from the block means of the coarsest level holding at
 * least _NCM_STATS_ACORR_MIN_BLOCKS values, scaled by the standard error the selected
 * long-run variance implies.
 */
static void
_ncm_stats_acorr_half_diag (NcmStatsAcorrVar *var, guint max_lag, const gdouble spec0, const guint64 n,
                            gdouble *drift_z, gdouble *var_ratio)
{
  guint j;

  drift_z[0]   = 0.0;
  var_ratio[0] = 1.0;

  for (j = var->levels->len; j > 0; j--)
  {
    NcmStatsAcorrLevel *lev = g_ptr_array_index (var->levels, j - 1);

    if (lev->n >= _NCM_STATS_ACORR_MIN_BLOCKS)
    {
      const gdouble *head = &g_array_index (lev->head, gdouble, 0);
      const guint nb      = (lev->n < max_lag) ? (guint) lev->n : max_lag;
      const guint half    = nb / 2;
      gdouble m1 = 0.0, m2 = 0.0, v1 = 0.0, v2 = 0.0;
      guint i;

      for (i = 0; i < half; i++)
        m1 += head[i];

      for (i = half; i < 2 * half; i++)
        m2 += head[i];

      m1 = m1 / (1.0 * half);
      m2 = m2 / (1.0 * half);

      for (i = 0; i < half; i++)
        v1 += gsl_pow_2 (head[i] - m1);

      for (i = half; i < 2 * half; i++)
        v2 += gsl_pow_2 (head[i] - m2);

      if ((v1 > 0.0) && (v2 > 0.0))
        var_ratio[0] = v2 / v1;

      if (spec0 > 0.0)
        drift_z[0] = (m2 - m1) / sqrt (4.0 * spec0 / (1.0 * n));

      return;
    }
  }
}

static gdouble
_ncm_stats_acorr_tau_level (NcmStatsAcorrVar *var, NcmStatsAcorrPrivate *self, NcmStatsAcorrLevel *lev, gdouble *C0_level)
{
  const guint nlag = (lev->n - 1 < self->max_lag) ? (guint) (lev->n - 1) : self->max_lag;
  NcmVector *acov  = ncm_vector_get_subvector (var->acov, 0, nlag + 1);
  gdouble tau;

  _ncm_stats_acorr_level_acov (lev, self->max_lag, acov);

  C0_level[0] = ncm_vector_get (acov, 0);

  var->tau_ar    = ncm_stats_acorr_tau_ar (acov, lev->n, self->ar_crit, &var->ar_order);
  var->tau_geyer = ncm_stats_acorr_tau_geyer (acov, &var->window);

  switch (self->method)
  {
    case NCM_STATS_ACORR_METHOD_AR:
      tau = var->tau_ar;
      break;
    case NCM_STATS_ACORR_METHOD_GEYER:
      tau = var->tau_geyer;
      break;
    case NCM_STATS_ACORR_METHOD_SOKAL:
      tau = ncm_stats_acorr_tau_sokal (acov, NCM_STATS_ACORR_SOKAL_C, &var->window);
      break;
    case NCM_STATS_ACORR_METHOD_MAX:
      tau = GSL_MAX (var->tau_ar, var->tau_geyer);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  ncm_vector_free (acov);

  return tau;
}

static void
_ncm_stats_acorr_prepare (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  NcmStatsAcorrVar *var             = &self->vars[p];
  NcmStatsAcorrLevel *lev0          = g_ptr_array_index (var->levels, 0);
  const guint64 n                   = lev0->n;

  if (!var->dirty)
    return;

  var->dirty     = FALSE;
  var->diag      = NCM_STATS_ACORR_DIAG_OK;
  var->level     = 0;
  var->window    = 0;
  var->ar_order  = 0;
  var->drift_z   = 0.0;
  var->var_ratio = 1.0;
  var->tau       = 1.0;

  if (n == 0)
  {
    var->mean = 0.0;
    var->var  = 0.0;
    var->diag = NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_ZERO_VARIANCE;

    return;
  }

  {
    const gdouble m = lev0->S / (1.0 * n);

    var->mean = var->offset + m;
  }

  if (n < 2)
  {
    var->var  = 0.0;
    var->diag = NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_ZERO_VARIANCE;

    return;
  }

  {
    NcmVector *acov0 = ncm_vector_get_subvector (var->acov, 0, 1);

    _ncm_stats_acorr_level_acov (lev0, self->max_lag, acov0);
    var->var = ncm_vector_get (acov0, 0);
    ncm_vector_free (acov0);
  }

  if (var->var <= 0.0)
  {
    /* A series that never moved is perfectly correlated, so its autocorrelation time is
     * unbounded, not one: report the cap, which makes the effective sample size one and
     * every error derived from it the full sample standard deviation. The flag marks the
     * value as a bound rather than a measurement. */
    var->tau  = 1.0 * n;
    var->diag = NCM_STATS_ACORR_DIAG_ZERO_VARIANCE;

    if (n < self->reliability_factor)
      var->diag |= NCM_STATS_ACORR_DIAG_SHORT_CHAIN;

    return;
  }

  {
    const gdouble margin = self->max_lag / _NCM_STATS_ACORR_LAG_MARGIN;
    gboolean resolved    = FALSE;
    guint j;

    for (j = 0; j < var->levels->len; j++)
    {
      NcmStatsAcorrLevel *lev = g_ptr_array_index (var->levels, j);
      gdouble tau_j, C0_j = 0.0;

      if ((lev->n < _NCM_STATS_ACORR_MIN_BLOCKS) && (j > 0))
        break;

      tau_j      = _ncm_stats_acorr_tau_level (var, self, lev, &C0_j);
      var->level = j;
      var->tau   = tau_j * ldexp (1.0, j) * C0_j / var->var;

      if (tau_j <= margin)
      {
        resolved = TRUE;
        break;
      }
    }

    if (!resolved)
      var->diag |= NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED;
  }

  var->tau = GSL_MAX (var->tau, 1.0 / (1.0 * n));
  var->tau = GSL_MIN (var->tau, 1.0 * n);

  if ((self->method == NCM_STATS_ACORR_METHOD_MAX) &&
      (GSL_MAX (var->tau_ar, var->tau_geyer) > 2.0 * GSL_MIN (var->tau_ar, var->tau_geyer)))
    var->diag |= NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT;

  if (n < self->reliability_factor * var->tau)
    var->diag |= NCM_STATS_ACORR_DIAG_SHORT_CHAIN;

  _ncm_stats_acorr_half_diag (var, self->max_lag, var->tau * var->var, n, &var->drift_z, &var->var_ratio);

  if (fabs (var->drift_z) > self->drift_threshold)
    var->diag |= NCM_STATS_ACORR_DIAG_DRIFT;

  /* A series whose two halves differ in scale by orders of magnitude is not one series.
   * The ratio is used rather than a z-score because the block values it is built from are
   * not independent, and because what this has to catch is orders of magnitude. */
  if ((var->var_ratio > NCM_STATS_ACORR_VARIANCE_SHIFT_FACTOR) ||
      (var->var_ratio < 1.0 / NCM_STATS_ACORR_VARIANCE_SHIFT_FACTOR))
    var->diag |= NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT;
}

/**
 * ncm_stats_acorr_new:
 * @len: number of series
 *
 * Creates a new #NcmStatsAcorr tracking @len series with the default lag budget and
 * estimator.
 *
 * Returns: (transfer full): a new #NcmStatsAcorr.
 */
NcmStatsAcorr *
ncm_stats_acorr_new (guint len)
{
  return ncm_stats_acorr_new_full (len,
                                   NCM_STATS_ACORR_DEFAULT_MAX_LAG,
                                   NCM_STATS_ACORR_DEFAULT_MAX_LEVELS,
                                   NCM_STATS_ACORR_METHOD_MAX);
}

/**
 * ncm_stats_acorr_new_full:
 * @len: number of series
 * @max_lag: lags accumulated per level
 * @max_levels: maximum number of block-averaging levels
 * @method: a #NcmStatsAcorrMethod
 *
 * Creates a new #NcmStatsAcorr with every parameter given.
 *
 * Returns: (transfer full): a new #NcmStatsAcorr.
 */
NcmStatsAcorr *
ncm_stats_acorr_new_full (guint len, guint max_lag, guint max_levels, NcmStatsAcorrMethod method)
{
  NcmStatsAcorr *acorr = g_object_new (NCM_TYPE_STATS_ACORR,
                                       "len", len,
                                       "max-lag", max_lag,
                                       "max-levels", max_levels,
                                       "method", method,
                                       NULL);

  return acorr;
}

/**
 * ncm_stats_acorr_ref:
 * @acorr: a #NcmStatsAcorr
 *
 * Increases the reference count of @acorr by one.
 *
 * Returns: (transfer full): @acorr.
 */
NcmStatsAcorr *
ncm_stats_acorr_ref (NcmStatsAcorr *acorr)
{
  return g_object_ref (acorr);
}

/**
 * ncm_stats_acorr_free:
 * @acorr: a #NcmStatsAcorr
 *
 * Decreases the reference count of @acorr by one.
 *
 */
void
ncm_stats_acorr_free (NcmStatsAcorr *acorr)
{
  g_object_unref (acorr);
}

/**
 * ncm_stats_acorr_clear:
 * @acorr: a #NcmStatsAcorr
 *
 * If *@acorr is different from NULL, decreases the reference count of *@acorr by one
 * and sets *@acorr to NULL.
 *
 */
void
ncm_stats_acorr_clear (NcmStatsAcorr **acorr)
{
  g_clear_object (acorr);
}

/**
 * ncm_stats_acorr_len:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: the number of series tracked.
 */
guint
ncm_stats_acorr_len (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->len;
}

/**
 * ncm_stats_acorr_set_method:
 * @acorr: a #NcmStatsAcorr
 * @method: a #NcmStatsAcorrMethod
 *
 * Sets the estimator. The accumulated autocovariances are unaffected, so the estimate
 * can be changed at any time.
 *
 */
void
ncm_stats_acorr_set_method (NcmStatsAcorr *acorr, NcmStatsAcorrMethod method)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  guint p;

  g_assert_cmpuint (method, <, NCM_STATS_ACORR_METHOD_LEN);

  self->method = method;

  for (p = 0; p < self->len; p++)
    self->vars[p].dirty = TRUE;
}

/**
 * ncm_stats_acorr_get_method:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: the estimator in use.
 */
NcmStatsAcorrMethod
ncm_stats_acorr_get_method (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->method;
}

/**
 * ncm_stats_acorr_get_max_lag:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: the number of lags accumulated per level.
 */
guint
ncm_stats_acorr_get_max_lag (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->max_lag;
}

/**
 * ncm_stats_acorr_get_max_levels:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: the maximum number of levels.
 */
guint
ncm_stats_acorr_get_max_levels (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->max_levels;
}

/**
 * ncm_stats_acorr_set_ar_criterion:
 * @acorr: a #NcmStatsAcorr
 * @crit: a #NcmStatsAcorrARCrit
 *
 * Sets the rule that picks the order of the auto-regressive fit. The accumulated
 * autocovariances are unaffected, so it can be changed at any time.
 *
 */
void
ncm_stats_acorr_set_ar_criterion (NcmStatsAcorr *acorr, NcmStatsAcorrARCrit crit)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  guint p;

  g_assert_cmpuint (crit, <, NCM_STATS_ACORR_AR_CRIT_LEN);

  self->ar_crit = crit;

  for (p = 0; p < self->len; p++)
    self->vars[p].dirty = TRUE;
}

/**
 * ncm_stats_acorr_get_ar_criterion:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: the rule that picks the auto-regressive order.
 */
NcmStatsAcorrARCrit
ncm_stats_acorr_get_ar_criterion (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->ar_crit;
}

/**
 * ncm_stats_acorr_set_reliability_factor:
 * @acorr: a #NcmStatsAcorr
 * @factor: required number of autocorrelation times
 *
 * Sets #NcmStatsAcorr:reliability-factor.
 *
 */
void
ncm_stats_acorr_set_reliability_factor (NcmStatsAcorr *acorr, gdouble factor)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  guint p;

  g_assert_cmpfloat (factor, >=, 1.0);

  self->reliability_factor = factor;

  for (p = 0; p < self->len; p++)
    self->vars[p].dirty = TRUE;
}

/**
 * ncm_stats_acorr_get_reliability_factor:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: #NcmStatsAcorr:reliability-factor.
 */
gdouble
ncm_stats_acorr_get_reliability_factor (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->reliability_factor;
}

/**
 * ncm_stats_acorr_set_drift_threshold:
 * @acorr: a #NcmStatsAcorr
 * @threshold: drift z-score threshold
 *
 * Sets #NcmStatsAcorr:drift-threshold.
 *
 */
void
ncm_stats_acorr_set_drift_threshold (NcmStatsAcorr *acorr, gdouble threshold)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  guint p;

  g_assert_cmpfloat (threshold, >=, 0.0);

  self->drift_threshold = threshold;

  for (p = 0; p < self->len; p++)
    self->vars[p].dirty = TRUE;
}

/**
 * ncm_stats_acorr_get_drift_threshold:
 * @acorr: a #NcmStatsAcorr
 *
 * Returns: #NcmStatsAcorr:drift-threshold.
 */
gdouble
ncm_stats_acorr_get_drift_threshold (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  return self->drift_threshold;
}

/**
 * ncm_stats_acorr_reset:
 * @acorr: a #NcmStatsAcorr
 *
 * Discards every accumulated sample.
 *
 */
void
ncm_stats_acorr_reset (NcmStatsAcorr *acorr)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  guint p;

  for (p = 0; p < self->len; p++)
    _ncm_stats_acorr_var_reset (&self->vars[p], self->max_lag);
}

/**
 * ncm_stats_acorr_update:
 * @acorr: a #NcmStatsAcorr
 * @x: a #NcmVector with one sample of every series
 *
 * Adds one sample of each series. Costs $O(\mathrm{max\text{-}lag})$ per series.
 *
 */
void
ncm_stats_acorr_update (NcmStatsAcorr *acorr, NcmVector *x)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  guint p;

  g_assert_cmpuint (ncm_vector_len (x), >=, self->len);

  for (p = 0; p < self->len; p++)
    _ncm_stats_acorr_var_update (&self->vars[p], self->max_lag, self->max_levels, ncm_vector_get (x, p));
}

/**
 * ncm_stats_acorr_update_var:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @x_p: the sample
 *
 * Adds one sample to series @p alone.
 *
 */
void
ncm_stats_acorr_update_var (NcmStatsAcorr *acorr, guint p, gdouble x_p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);

  _ncm_stats_acorr_var_update (&self->vars[p], self->max_lag, self->max_levels, x_p);
}

/**
 * ncm_stats_acorr_set_series:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @series: a #NcmVector holding the whole series
 *
 * Discards whatever series @p held and accumulates @series in order.
 *
 */
void
ncm_stats_acorr_set_series (NcmStatsAcorr *acorr, guint p, NcmVector *series)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  const guint n                     = ncm_vector_len (series);
  guint i;

  g_assert_cmpuint (p, <, self->len);

  _ncm_stats_acorr_var_reset (&self->vars[p], self->max_lag);

  for (i = 0; i < n; i++)
    _ncm_stats_acorr_var_update (&self->vars[p], self->max_lag, self->max_levels, ncm_vector_get (series, i));
}

/**
 * ncm_stats_acorr_set_series_matrix:
 * @acorr: a #NcmStatsAcorr
 * @series: a #NcmMatrix whose rows are samples and columns are series
 *
 * Discards everything accumulated and feeds every row of @series in order.
 *
 */
void
ncm_stats_acorr_set_series_matrix (NcmStatsAcorr *acorr, NcmMatrix *series)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  const guint nrows                 = ncm_matrix_nrows (series);
  guint i, p;

  g_assert_cmpuint (ncm_matrix_ncols (series), >=, self->len);

  ncm_stats_acorr_reset (acorr);

  for (i = 0; i < nrows; i++)
    for (p = 0; p < self->len; p++)
      _ncm_stats_acorr_var_update (&self->vars[p], self->max_lag, self->max_levels, ncm_matrix_get (series, i, p));
}

/**
 * ncm_stats_acorr_nitens:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the number of samples accumulated in series @p.
 */
guint64
ncm_stats_acorr_nitens (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  NcmStatsAcorrLevel *lev0;

  g_assert_cmpuint (p, <, self->len);
  lev0 = g_ptr_array_index (self->vars[p].levels, 0);

  return lev0->n;
}

/**
 * ncm_stats_acorr_get_mean:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the mean of series @p.
 */
gdouble
ncm_stats_acorr_get_mean (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].mean;
}

/**
 * ncm_stats_acorr_get_var:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the variance of series @p, normalized by the number of samples.
 */
gdouble
ncm_stats_acorr_get_var (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].var;
}

/**
 * ncm_stats_acorr_get_tau:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Integrated autocorrelation time of series @p, in samples, by the estimator set in
 * #NcmStatsAcorr:method. Check ncm_stats_acorr_get_diag() before using it.
 *
 * Returns: $\tau$.
 */
gdouble
ncm_stats_acorr_get_tau (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].tau;
}

/**
 * ncm_stats_acorr_get_tau_method:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @method: a #NcmStatsAcorrMethod
 *
 * Integrated autocorrelation time of series @p by @method, at the level selected for
 * #NcmStatsAcorr:method. Lets one estimator be compared with another without changing
 * the object's configuration.
 *
 * Returns: $\tau$.
 */
gdouble
ncm_stats_acorr_get_tau_method (NcmStatsAcorr *acorr, guint p, NcmStatsAcorrMethod method)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  NcmStatsAcorrVar *var             = &self->vars[p];
  NcmStatsAcorrLevel *lev;
  NcmVector *acov;
  guint nlag;
  gdouble tau_j, C0_j;

  g_assert_cmpuint (p, <, self->len);
  g_assert_cmpuint (method, <, NCM_STATS_ACORR_METHOD_LEN);
  _ncm_stats_acorr_prepare (acorr, p);

  if (var->var <= 0.0)
    return 1.0;

  lev  = g_ptr_array_index (var->levels, var->level);
  nlag = (lev->n - 1 < self->max_lag) ? (guint) (lev->n - 1) : self->max_lag;
  acov = ncm_vector_get_subvector (var->acov, 0, nlag + 1);

  _ncm_stats_acorr_level_acov (lev, self->max_lag, acov);
  C0_j = ncm_vector_get (acov, 0);

  switch (method)
  {
    case NCM_STATS_ACORR_METHOD_AR:
      tau_j = ncm_stats_acorr_tau_ar (acov, lev->n, self->ar_crit, NULL);
      break;
    case NCM_STATS_ACORR_METHOD_GEYER:
      tau_j = ncm_stats_acorr_tau_geyer (acov, NULL);
      break;
    case NCM_STATS_ACORR_METHOD_SOKAL:
      tau_j = ncm_stats_acorr_tau_sokal (acov, NCM_STATS_ACORR_SOKAL_C, NULL);
      break;
    case NCM_STATS_ACORR_METHOD_MAX:
      tau_j = GSL_MAX (ncm_stats_acorr_tau_ar (acov, lev->n, self->ar_crit, NULL),
                       ncm_stats_acorr_tau_geyer (acov, NULL));
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  ncm_vector_free (acov);

  return GSL_MIN (tau_j * ldexp (1.0, var->level) * C0_j / var->var, 1.0 * ncm_stats_acorr_nitens (acorr, p));
}

/**
 * ncm_stats_acorr_get_ess:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Effective sample size $n/\tau$ of series @p.
 *
 * Returns: the effective sample size.
 */
gdouble
ncm_stats_acorr_get_ess (NcmStatsAcorr *acorr, guint p)
{
  const gdouble tau = ncm_stats_acorr_get_tau (acorr, p);
  const guint64 n   = ncm_stats_acorr_nitens (acorr, p);

  if (n == 0)
    return 0.0;

  return (1.0 * n) / tau;
}

/**
 * ncm_stats_acorr_get_spec0:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Long-run variance $S(0) = \tau\,\mathrm{Var}(x)$ of series @p: the spectral density
 * at zero frequency, which is what the variance of the mean is built from.
 *
 * Returns: $S(0)$.
 */
gdouble
ncm_stats_acorr_get_spec0 (NcmStatsAcorr *acorr, guint p)
{
  return ncm_stats_acorr_get_tau (acorr, p) * ncm_stats_acorr_get_var (acorr, p);
}

/**
 * ncm_stats_acorr_get_var_mean:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Variance of the sample mean of series @p, $S(0)/n$. No assumption is made about the
 * samples beyond the series itself being the one averaged.
 *
 * Returns: the variance of the mean.
 */
gdouble
ncm_stats_acorr_get_var_mean (NcmStatsAcorr *acorr, guint p)
{
  const guint64 n = ncm_stats_acorr_nitens (acorr, p);

  if (n == 0)
    return 0.0;

  return ncm_stats_acorr_get_spec0 (acorr, p) / (1.0 * n);
}

/**
 * ncm_stats_acorr_get_sd_mean:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the standard error of the mean of series @p.
 */
gdouble
ncm_stats_acorr_get_sd_mean (NcmStatsAcorr *acorr, guint p)
{
  return sqrt (ncm_stats_acorr_get_var_mean (acorr, p));
}

/**
 * ncm_stats_acorr_get_diag:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Conditions attached to the estimate of series @p.
 *
 * Returns: the #NcmStatsAcorrDiag flags.
 */
NcmStatsAcorrDiag
ncm_stats_acorr_get_diag (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].diag;
}

/**
 * ncm_stats_acorr_diag_to_string:
 * @diag: a #NcmStatsAcorrDiag
 *
 * Returns: (transfer full): a comma separated list of the conditions set in @diag, or
 *    "ok" when none is.
 */
gchar *
ncm_stats_acorr_diag_to_string (NcmStatsAcorrDiag diag)
{
  GString *str = g_string_new ("");

  if (diag & NCM_STATS_ACORR_DIAG_SHORT_CHAIN)
    g_string_append (str, "short-chain,");

  if (diag & NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED)
    g_string_append (str, "window-truncated,");

  if (diag & NCM_STATS_ACORR_DIAG_DRIFT)
    g_string_append (str, "drift,");

  if (diag & NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT)
    g_string_append (str, "method-disagreement,");

  if (diag & NCM_STATS_ACORR_DIAG_ZERO_VARIANCE)
    g_string_append (str, "zero-variance,");

  if (diag & NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT)
    g_string_append (str, "variance-shift,");

  if (str->len == 0)
    g_string_append (str, "ok");
  else
    g_string_truncate (str, str->len - 1);

  return g_string_free (str, FALSE);
}

/**
 * ncm_stats_acorr_get_level:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the block-averaging level the estimate of series @p was taken from.
 */
guint
ncm_stats_acorr_get_level (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].level;
}

/**
 * ncm_stats_acorr_get_window:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the number of lags summed by the windowed estimator of series @p.
 */
guint
ncm_stats_acorr_get_window (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].window;
}

/**
 * ncm_stats_acorr_get_ar_order:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the auto-regressive order AICc selected for series @p.
 */
guint
ncm_stats_acorr_get_ar_order (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].ar_order;
}

/**
 * ncm_stats_acorr_get_drift_z:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Difference between the mean of the first and of the second half of series @p, in
 * units of the standard error of that difference.
 *
 * Returns: the drift z-score.
 */
gdouble
ncm_stats_acorr_get_drift_z (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].drift_z;
}

/**
 * ncm_stats_acorr_get_var_ratio:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Variance of the second half of series @p divided by the variance of the first half,
 * both taken over the block values of the coarsest level holding enough of them. One
 * when there are too few.
 *
 * Returns: the ratio of the two half variances.
 */
gdouble
ncm_stats_acorr_get_var_ratio (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  return self->vars[p].var_ratio;
}

/**
 * ncm_stats_acorr_nlevels:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 *
 * Returns: the number of block-averaging levels series @p has filled.
 */
guint
ncm_stats_acorr_nlevels (NcmStatsAcorr *acorr, guint p)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);

  g_assert_cmpuint (p, <, self->len);

  return self->vars[p].levels->len;
}

/**
 * ncm_stats_acorr_level_nitens:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @level: level index
 *
 * Returns: the number of block values level @level of series @p holds.
 */
guint64
ncm_stats_acorr_level_nitens (NcmStatsAcorr *acorr, guint p, guint level)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  NcmStatsAcorrLevel *lev;

  g_assert_cmpuint (p, <, self->len);
  g_assert_cmpuint (level, <, self->vars[p].levels->len);

  lev = g_ptr_array_index (self->vars[p].levels, level);

  return lev->n;
}

/**
 * ncm_stats_acorr_get_ar_fit:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @phi: (out) (transfer full) (optional): coefficients of the selected order
 * @pacf: (out) (transfer full) (optional): reflection coefficients of every order tried
 * @ivar: (out) (optional): innovation variance of the selected order
 * @order: (out) (optional): order selected
 *
 * Auto-regressive fit of series @p at the level its estimate was taken from, returning
 * what the fit is made of rather than the $\tau$ built from it. See
 * ncm_stats_acorr_ar_fit().
 *
 * Returns: TRUE when the selected order is not zero.
 */
gboolean
ncm_stats_acorr_get_ar_fit (NcmStatsAcorr *acorr, guint p, NcmVector **phi, NcmVector **pacf, gdouble *ivar, guint *order)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  NcmStatsAcorrLevel *lev;
  NcmVector *acov;
  gboolean ret;

  g_assert_cmpuint (p, <, self->len);
  _ncm_stats_acorr_prepare (acorr, p);

  lev  = g_ptr_array_index (self->vars[p].levels, self->vars[p].level);
  acov = ncm_stats_acorr_get_acov (acorr, p, self->vars[p].level);
  ret  = ncm_stats_acorr_ar_fit (acov, lev->n, self->ar_crit, phi, pacf, ivar, order);

  ncm_vector_free (acov);

  return ret;
}

/**
 * ncm_stats_acorr_get_acov:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @level: level index
 *
 * Autocovariances of series @p at level @level, $C_0 \dots C_L$ with
 * $L = \min(\mathrm{max\text{-}lag}, n_\mathrm{level}-1)$, normalized by the number of
 * block values of that level.
 *
 * Returns: (transfer full): the autocovariances.
 */
NcmVector *
ncm_stats_acorr_get_acov (NcmStatsAcorr *acorr, guint p, guint level)
{
  NcmStatsAcorrPrivate * const self = ncm_stats_acorr_get_instance_private (acorr);
  NcmStatsAcorrLevel *lev;
  NcmVector *acov;
  guint nlag;

  g_assert_cmpuint (p, <, self->len);
  g_assert_cmpuint (level, <, self->vars[p].levels->len);

  lev = g_ptr_array_index (self->vars[p].levels, level);

  g_assert_cmpuint (lev->n, >, 1);

  nlag = (lev->n - 1 < self->max_lag) ? (guint) (lev->n - 1) : self->max_lag;
  acov = ncm_vector_new (nlag + 1);

  _ncm_stats_acorr_level_acov (lev, self->max_lag, acov);

  return acov;
}

/**
 * ncm_stats_acorr_get_acf:
 * @acorr: a #NcmStatsAcorr
 * @p: series index
 * @level: level index
 *
 * Autocorrelation function of series @p at level @level, $C_k/C_0$.
 *
 * Returns: (transfer full): the autocorrelation function.
 */
NcmVector *
ncm_stats_acorr_get_acf (NcmStatsAcorr *acorr, guint p, guint level)
{
  NcmVector *acf   = ncm_stats_acorr_get_acov (acorr, p, level);
  const gdouble C0 = ncm_vector_get (acf, 0);

  if (C0 > 0.0)
    ncm_vector_scale (acf, 1.0 / C0);

  return acf;
}

