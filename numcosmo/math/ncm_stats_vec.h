/***************************************************************************
 *            ncm_stats_vec.h
 *
 *  Fri August 02 13:41:15 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_vec.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_STATS_VEC_H_
#define _NCM_STATS_VEC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */
#ifndef NUMCOSMO_GIR_SCAN
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_STATS_VEC             (ncm_stats_vec_get_type ())
#define NCM_STATS_VEC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_VEC, NcmStatsVec))
#define NCM_STATS_VEC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_VEC, NcmStatsVecClass))
#define NCM_IS_STATS_VEC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_VEC))
#define NCM_IS_STATS_VEC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_VEC))
#define NCM_STATS_VEC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_VEC, NcmStatsVecClass))

typedef struct _NcmStatsVecClass NcmStatsVecClass;
typedef struct _NcmStatsVec NcmStatsVec;

struct _NcmStatsVecClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmStatsVecType:
 * @NCM_STATS_VEC_MEAN: Calculates mean only.
 * @NCM_STATS_VEC_VAR: Calculates mean and variance.
 * @NCM_STATS_VEC_COV: Calculates mean, variance and covariance.
 * 
 * FIXME
 * 
 */
typedef enum 
{
  NCM_STATS_VEC_MEAN = 0,
  NCM_STATS_VEC_VAR,
  NCM_STATS_VEC_COV,       
  /* < private > */
  NCM_STATS_VEC_TYPES_LEN, /*< skip >*/
} NcmStatsVecType;

typedef void (*NcmStatsVecUpdateFunc) (NcmStatsVec *svec, const gdouble w, NcmVector *x);

struct _NcmStatsVec
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsVecType t;
  NcmStatsVecUpdateFunc update;
  NcmStatsVec *tmp;
  guint len;
  gboolean save_x;
  gdouble weight;
  gdouble weight2;
  gdouble bias_wt;
  guint nitens;
  NcmVector *x;
  NcmVector *mean;
  NcmVector *var;
  NcmMatrix *cov;
  NcmMatrix *real_cov;
  GPtrArray *saved_x;
  GPtrArray *q_array;
#ifdef NUMCOSMO_HAVE_FFTW3
  guint fft_size;
  guint fft_plan_size;
  gdouble *param_data;
  fftw_complex *param_fft;
  fftw_plan param_r2c;
  fftw_plan param_c2r;
#endif /* NUMCOSMO_HAVE_FFTW3 */
};

/**
 * NcmStatsVecARType:
 * @NCM_STATS_VEC_AR_NONE: Calculates using the required order.
 * @NCM_STATS_VEC_AR_FPE: Uses the FPE criterium to choose the ar order.
 * @NCM_STATS_VEC_AR_AIC: Uses the AIC criterium to choose the ar order.
 * @NCM_STATS_VEC_AR_AICC: Uses the AICc criterium to choose the ar order.
 * 
 * FIXME
 * 
 */
typedef enum 
{
  NCM_STATS_VEC_AR_NONE = 0,
  NCM_STATS_VEC_AR_FPE,
  NCM_STATS_VEC_AR_AIC,
  NCM_STATS_VEC_AR_AICC, 
  /* < private > */
  NCM_STATS_VEC_AR_LEN,  /*< skip >*/
} NcmStatsVecARType;

GType ncm_stats_vec_get_type (void) G_GNUC_CONST;

NcmStatsVec *ncm_stats_vec_new (guint len, NcmStatsVecType t, gboolean save_x);
NcmStatsVec *ncm_stats_vec_ref (NcmStatsVec *svec);
void ncm_stats_vec_free (NcmStatsVec *svec);
void ncm_stats_vec_clear (NcmStatsVec **svec);

void ncm_stats_vec_reset (NcmStatsVec *svec, gboolean rm_saved);
void ncm_stats_vec_update_weight (NcmStatsVec *svec, gdouble w);

void ncm_stats_vec_append_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup);
void ncm_stats_vec_prepend_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup);

void ncm_stats_vec_append (NcmStatsVec *svec, NcmVector *x, gboolean dup);
void ncm_stats_vec_prepend (NcmStatsVec *svec, NcmVector *x, gboolean dup);
void ncm_stats_vec_append_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup);
void ncm_stats_vec_prepend_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup);

void ncm_stats_vec_enable_quantile (NcmStatsVec *svec, gdouble p);
void ncm_stats_vec_disable_quantile (NcmStatsVec *svec);
gdouble ncm_stats_vec_get_quantile (NcmStatsVec *svec, guint i);
gdouble ncm_stats_vec_get_quantile_spread (NcmStatsVec *svec, guint i);

NcmVector *ncm_stats_vec_get_autocorr (NcmStatsVec *svec, guint p);
NcmVector *ncm_stats_vec_get_subsample_autocorr (NcmStatsVec *svec, guint p, guint subsample);
gdouble ncm_stats_vec_get_autocorr_tau (NcmStatsVec *svec, guint p, const guint max_lag);
gdouble ncm_stats_vec_get_subsample_autocorr_tau (NcmStatsVec *svec, guint p, guint subsample, const guint max_lag);

gboolean ncm_stats_vec_fit_ar_model (NcmStatsVec *svec, guint p, const guint order, NcmStatsVecARType ar_crit, NcmVector **rho, NcmVector **pacf, gdouble *ivar, guint *c_order);
gdouble ncm_stats_vec_ar_ess (NcmStatsVec *svec, guint p, NcmStatsVecARType ar_crit, gdouble *spec0, guint *c_order);
gdouble ncm_stats_vec_estimate_const_break (NcmStatsVec *svec, guint p);

NcmVector *ncm_stats_vec_max_ess_time (NcmStatsVec *svec, const guint ntests, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_ess);
NcmVector *ncm_stats_vec_heidel_diag (NcmStatsVec *svec, const guint ntests, const gdouble pvalue, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_pvalue);
NcmVector *ncm_stats_vec_visual_heidel_diag (NcmStatsVec *svec, const guint p, const guint fi, gdouble *mean, gdouble *var);

GPtrArray *ncm_stats_vec_dup_saved_x (NcmStatsVec *svec);

NCM_INLINE NcmVector *ncm_stats_vec_peek_x (NcmStatsVec *svec);
NCM_INLINE void ncm_stats_vec_set (NcmStatsVec *svec, guint i, gdouble x_i);
NCM_INLINE gdouble ncm_stats_vec_get (NcmStatsVec *svec, guint i);
NCM_INLINE void ncm_stats_vec_update (NcmStatsVec *svec);
NCM_INLINE guint ncm_stats_vec_len (NcmStatsVec *svec);
NCM_INLINE gdouble ncm_stats_vec_get_mean (NcmStatsVec *svec, guint i);
NCM_INLINE gdouble ncm_stats_vec_get_var (NcmStatsVec *svec, guint i);
NCM_INLINE gdouble ncm_stats_vec_get_sd (NcmStatsVec *svec, guint i);
NCM_INLINE gdouble ncm_stats_vec_get_cov (NcmStatsVec *svec, guint i, guint j);
NCM_INLINE gdouble ncm_stats_vec_get_cor (NcmStatsVec *svec, guint i, guint j);
NCM_INLINE gdouble ncm_stats_vec_get_weight (NcmStatsVec *svec);
NCM_INLINE void ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *mean, guint offset);
NCM_INLINE NcmVector *ncm_stats_vec_peek_mean (NcmStatsVec *svec);
NCM_INLINE void ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m, guint offset);
NCM_INLINE NcmMatrix *ncm_stats_vec_peek_cov_matrix (NcmStatsVec *svec, guint offset);
NCM_INLINE guint ncm_stats_vec_nrows (NcmStatsVec *svec);
NCM_INLINE guint ncm_stats_vec_nitens (NcmStatsVec *svec);
NCM_INLINE NcmVector *ncm_stats_vec_peek_row (NcmStatsVec *svec, guint i);
NCM_INLINE gdouble ncm_stats_vec_get_param_at (NcmStatsVec *svec, guint i, guint p);

#define NCM_STATS_VEC_HEIDEL_PVAL_COR(pvalue,n) (1.0 - pow (1.0 - (pvalue), 1.0 / ((gdouble)(n))))

G_END_DECLS

#endif /* _NCM_STATS_VEC_H_ */

#ifndef _NCM_STATS_VEC_INLINE_H_
#define _NCM_STATS_VEC_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcmVector *
ncm_stats_vec_peek_x (NcmStatsVec *svec)
{
  return svec->x;
}

NCM_INLINE void 
ncm_stats_vec_set (NcmStatsVec *svec, guint i, gdouble x_i)
{
  ncm_vector_fast_set (svec->x, i, x_i);
}

NCM_INLINE gdouble 
ncm_stats_vec_get (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->x, i);
}

NCM_INLINE void 
ncm_stats_vec_update (NcmStatsVec *svec)
{
  ncm_stats_vec_update_weight (svec, 1.0);
}

NCM_INLINE guint 
ncm_stats_vec_len (NcmStatsVec *svec)
{
  return svec->len;
}

NCM_INLINE gdouble 
ncm_stats_vec_get_mean (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->mean, i);
}

NCM_INLINE gdouble 
ncm_stats_vec_get_var (NcmStatsVec *svec, guint i)
{
  g_assert (svec->t == NCM_STATS_VEC_VAR || svec->t == NCM_STATS_VEC_COV);
  return ncm_vector_fast_get (svec->var, i) * svec->bias_wt;
}

NCM_INLINE gdouble 
ncm_stats_vec_get_sd (NcmStatsVec *svec, guint i)
{
  return sqrt (ncm_stats_vec_get_var (svec, i));
}

NCM_INLINE gdouble 
ncm_stats_vec_get_cov (NcmStatsVec *svec, guint i, guint j)
{
  g_assert (svec->t == NCM_STATS_VEC_COV);
  if (i == j)
    return ncm_stats_vec_get_var (svec, i);
  else
    return ncm_matrix_get (svec->cov, i, j) * svec->bias_wt;
}

NCM_INLINE gdouble 
ncm_stats_vec_get_cor (NcmStatsVec *svec, guint i, guint j)
{
  if (i == j)
    return 1.0;
  else
    return ncm_stats_vec_get_cov (svec, i, j) / (ncm_stats_vec_get_sd (svec, i) * ncm_stats_vec_get_sd (svec, j));
}

NCM_INLINE gdouble 
ncm_stats_vec_get_weight (NcmStatsVec *svec)
{
  return svec->weight;
}

NCM_INLINE void 
ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *x, guint offset)
{
  g_assert (x != NULL);
  g_assert_cmpint (offset, <, svec->len);
  ncm_vector_memcpy2 (x, svec->mean, 0, offset, svec->len - offset);
}

NCM_INLINE NcmVector *
ncm_stats_vec_peek_mean (NcmStatsVec *svec)
{
  return svec->mean;
}

NCM_INLINE void 
ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m, guint offset)
{
  guint i;
  g_assert (m != NULL);
  g_assert_cmpint (offset, <, svec->len);

  if (offset > 0)
  {
    NcmMatrix *m_src = ncm_matrix_get_submatrix (svec->cov, offset, offset, svec->len - offset, svec->len - offset);
    ncm_matrix_memcpy (m, m_src);
    ncm_matrix_free (m_src);
  }
  else
    ncm_matrix_memcpy (m, svec->cov);
  
  for (i = 0; i < svec->len - offset; i++)
    ncm_matrix_set (m, i, i, ncm_vector_fast_get (svec->var, i + offset));

  ncm_matrix_scale (m, svec->bias_wt);
}

NCM_INLINE NcmMatrix *
ncm_stats_vec_peek_cov_matrix (NcmStatsVec *svec, guint offset)
{
  gint effsize = svec->len - offset;
  g_assert_cmpint (effsize, >, 0);
  if (svec->real_cov != NULL)
  {
    if ((gint)ncm_matrix_nrows (svec->real_cov) != effsize)
    {
      ncm_matrix_free (svec->real_cov);
      svec->real_cov = ncm_matrix_new (effsize, effsize);
    }
  }
  else
    svec->real_cov = ncm_matrix_new (effsize, effsize);

  ncm_stats_vec_get_cov_matrix (svec, svec->real_cov, offset);

  return svec->real_cov;
}

NCM_INLINE guint 
ncm_stats_vec_nrows (NcmStatsVec *svec)
{
  g_assert (svec->save_x);
  return svec->saved_x->len;
}

NCM_INLINE guint 
ncm_stats_vec_nitens (NcmStatsVec *svec)
{
  return svec->nitens;
}

NCM_INLINE NcmVector *
ncm_stats_vec_peek_row (NcmStatsVec *svec, guint i)
{
  g_assert (svec->save_x);
  g_assert (i < svec->saved_x->len);

  return g_ptr_array_index (svec->saved_x, i);
}

NCM_INLINE gdouble 
ncm_stats_vec_get_param_at (NcmStatsVec *svec, guint i, guint p)
{
  g_assert (svec->save_x);
  g_assert (i < svec->nitens);
  return ncm_vector_get (g_ptr_array_index (svec->saved_x, i), p);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_STATS_VEC_INLINE_H_ */

