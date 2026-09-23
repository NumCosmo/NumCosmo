/***************************************************************************
 *            ncm_stats_acorr.h
 *
 *  Tue Sep 16 09:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_acorr.h
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

#ifndef _NCM_STATS_ACORR_H_
#define _NCM_STATS_ACORR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_ACORR (ncm_stats_acorr_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsAcorr, ncm_stats_acorr, NCM, STATS_ACORR, GObject)

/**
 * NcmStatsAcorrMethod:
 * @NCM_STATS_ACORR_METHOD_AR: auto-regressive spectral estimate
 * @NCM_STATS_ACORR_METHOD_GEYER: Geyer initial monotone positive sequence
 * @NCM_STATS_ACORR_METHOD_SOKAL: self-consistent window $M \geq c\tau(M)$
 * @NCM_STATS_ACORR_METHOD_MAX: largest of AR and GEYER
 *
 * Estimator used to turn the accumulated autocovariances into an integrated
 * autocorrelation time.
 *
 */
typedef enum _NcmStatsAcorrMethod /*< enum,underscore_name=NCM_STATS_ACORR_METHOD,prefix=NCM_STATS_ACORR_METHOD >*/
{
  NCM_STATS_ACORR_METHOD_AR = 0,
  NCM_STATS_ACORR_METHOD_GEYER,
  NCM_STATS_ACORR_METHOD_SOKAL,
  NCM_STATS_ACORR_METHOD_MAX,
  /* < private > */
  NCM_STATS_ACORR_METHOD_LEN, /*< skip >*/
} NcmStatsAcorrMethod;

/**
 * NcmStatsAcorrARCrit:
 * @NCM_STATS_ACORR_AR_CRIT_NONE: no selection, the largest order tried is used
 * @NCM_STATS_ACORR_AR_CRIT_FPE: final prediction error
 * @NCM_STATS_ACORR_AR_CRIT_AIC: Akaike information criterion
 * @NCM_STATS_ACORR_AR_CRIT_AICC: Akaike information criterion corrected for small samples
 *
 * Rule that picks the order of the auto-regressive fit.
 *
 */
typedef enum _NcmStatsAcorrARCrit /*< enum,underscore_name=NCM_STATS_ACORR_AR_CRIT,prefix=NCM_STATS_ACORR_AR_CRIT >*/
{
  NCM_STATS_ACORR_AR_CRIT_NONE = 0,
  NCM_STATS_ACORR_AR_CRIT_FPE,
  NCM_STATS_ACORR_AR_CRIT_AIC,
  NCM_STATS_ACORR_AR_CRIT_AICC,
  /* < private > */
  NCM_STATS_ACORR_AR_CRIT_LEN, /*< skip >*/
} NcmStatsAcorrARCrit;

/**
 * NcmStatsAcorrDiag:
 * @NCM_STATS_ACORR_DIAG_OK: no condition detected
 * @NCM_STATS_ACORR_DIAG_SHORT_CHAIN: fewer than #NcmStatsAcorr:reliability-factor
 *    autocorrelation times in the series
 * @NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED: no level resolves the correlation within
 *    #NcmStatsAcorr:max-lag lags, the estimate is a lower bound
 * @NCM_STATS_ACORR_DIAG_DRIFT: the means of the first and of the second half of the
 *    series differ by more than #NcmStatsAcorr:drift-threshold standard errors
 * @NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT: the two estimators of
 *    %NCM_STATS_ACORR_METHOD_MAX differ by more than a factor of two
 * @NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT: the two halves of the series differ in variance by
 *    more than a factor of a hundred, which is what an unremoved burn-in looks like
 * @NCM_STATS_ACORR_DIAG_ZERO_VARIANCE: the series never moved, so $\tau$ is unbounded and
 *    is reported at its cap, the number of samples
 *
 * Conditions attached to an estimate. Every one of them means the estimate is not to be
 * read as a converged autocorrelation time.
 *
 */
typedef enum _NcmStatsAcorrDiag /*< flags,underscore_name=NCM_STATS_ACORR_DIAG,prefix=NCM_STATS_ACORR_DIAG >*/
{
  NCM_STATS_ACORR_DIAG_OK                  = 0,
  NCM_STATS_ACORR_DIAG_SHORT_CHAIN         = 1 << 0,
  NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED    = 1 << 1,
  NCM_STATS_ACORR_DIAG_DRIFT               = 1 << 2,
  NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT = 1 << 3,
  NCM_STATS_ACORR_DIAG_ZERO_VARIANCE       = 1 << 4,
  NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT      = 1 << 5,
} NcmStatsAcorrDiag;

/* Constructors and references */
NcmStatsAcorr *ncm_stats_acorr_new (guint len);
NcmStatsAcorr *ncm_stats_acorr_new_full (guint len, guint max_lag, guint max_levels, NcmStatsAcorrMethod method);
NcmStatsAcorr *ncm_stats_acorr_ref (NcmStatsAcorr *acorr);
void ncm_stats_acorr_free (NcmStatsAcorr *acorr);
void ncm_stats_acorr_clear (NcmStatsAcorr **acorr);

/* Configuration */
guint ncm_stats_acorr_len (NcmStatsAcorr *acorr);
void ncm_stats_acorr_set_method (NcmStatsAcorr *acorr, NcmStatsAcorrMethod method);
NcmStatsAcorrMethod ncm_stats_acorr_get_method (NcmStatsAcorr *acorr);
guint ncm_stats_acorr_get_max_lag (NcmStatsAcorr *acorr);
guint ncm_stats_acorr_get_max_levels (NcmStatsAcorr *acorr);
void ncm_stats_acorr_set_reliability_factor (NcmStatsAcorr *acorr, gdouble factor);
gdouble ncm_stats_acorr_get_reliability_factor (NcmStatsAcorr *acorr);
void ncm_stats_acorr_set_drift_threshold (NcmStatsAcorr *acorr, gdouble threshold);
gdouble ncm_stats_acorr_get_drift_threshold (NcmStatsAcorr *acorr);
void ncm_stats_acorr_set_ar_criterion (NcmStatsAcorr *acorr, NcmStatsAcorrARCrit crit);
NcmStatsAcorrARCrit ncm_stats_acorr_get_ar_criterion (NcmStatsAcorr *acorr);

/* Feeding */
void ncm_stats_acorr_reset (NcmStatsAcorr *acorr);
void ncm_stats_acorr_update (NcmStatsAcorr *acorr, NcmVector *x);
void ncm_stats_acorr_update_var (NcmStatsAcorr *acorr, guint p, gdouble x_p);
void ncm_stats_acorr_set_series (NcmStatsAcorr *acorr, guint p, NcmVector *series);
void ncm_stats_acorr_set_series_matrix (NcmStatsAcorr *acorr, NcmMatrix *series);

/* Results */
guint64 ncm_stats_acorr_nitens (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_mean (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_var (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_tau (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_tau_method (NcmStatsAcorr *acorr, guint p, NcmStatsAcorrMethod method);
gdouble ncm_stats_acorr_get_ess (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_spec0 (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_var_mean (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_sd_mean (NcmStatsAcorr *acorr, guint p);

/* Diagnostics */
NcmStatsAcorrDiag ncm_stats_acorr_get_diag (NcmStatsAcorr *acorr, guint p);
gchar *ncm_stats_acorr_diag_to_string (NcmStatsAcorrDiag diag);
guint ncm_stats_acorr_get_level (NcmStatsAcorr *acorr, guint p);
guint ncm_stats_acorr_get_window (NcmStatsAcorr *acorr, guint p);
guint ncm_stats_acorr_get_ar_order (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_drift_z (NcmStatsAcorr *acorr, guint p);
gdouble ncm_stats_acorr_get_var_ratio (NcmStatsAcorr *acorr, guint p);
gboolean ncm_stats_acorr_get_ar_fit (NcmStatsAcorr *acorr, guint p, NcmVector **phi, NcmVector **pacf, gdouble *ivar, guint *order);

/* Autocovariance access */
NcmVector *ncm_stats_acorr_get_acov (NcmStatsAcorr *acorr, guint p, guint level);
NcmVector *ncm_stats_acorr_get_acf (NcmStatsAcorr *acorr, guint p, guint level);
guint ncm_stats_acorr_nlevels (NcmStatsAcorr *acorr, guint p);
guint64 ncm_stats_acorr_level_nitens (NcmStatsAcorr *acorr, guint p, guint level);

/* Estimators on a given autocovariance sequence */
gdouble ncm_stats_acorr_tau_ar (NcmVector *acov, guint64 nitens, NcmStatsAcorrARCrit crit, guint *ar_order);
gboolean ncm_stats_acorr_ar_fit (NcmVector *acov, guint64 nitens, NcmStatsAcorrARCrit crit, NcmVector **phi, NcmVector **pacf, gdouble *ivar, guint *order);
gdouble ncm_stats_acorr_tau_geyer (NcmVector *acov, guint *window);
gdouble ncm_stats_acorr_tau_sokal (NcmVector *acov, gdouble c, guint *window);
NcmVector *ncm_stats_acorr_acov_fft (NcmVector *series, guint max_lag);

#define NCM_STATS_ACORR_DEFAULT_MAX_LAG (512)
#define NCM_STATS_ACORR_DEFAULT_MAX_LEVELS (24)
#define NCM_STATS_ACORR_DEFAULT_RELIABILITY_FACTOR (50.0)
#define NCM_STATS_ACORR_DEFAULT_DRIFT_THRESHOLD (3.0)
#define NCM_STATS_ACORR_SOKAL_C (5.0)
#define NCM_STATS_ACORR_DEFAULT_AR_CRIT (NCM_STATS_ACORR_AR_CRIT_AICC)
#define NCM_STATS_ACORR_VARIANCE_SHIFT_FACTOR (100.0)

G_END_DECLS

#endif /* _NCM_STATS_ACORR_H_ */

