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

#include <math.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_VEC             (ncm_stats_vec_get_type ())
#define NCM_STATS_VEC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_VEC, NcmStatsVec))
#define NCM_STATS_VEC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_VEC, NcmStatsVecClass))
#define NCM_IS_STATS_VEC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_VEC))
#define NCM_IS_STATS_VEC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_VEC))
#define NCM_STATS_VEC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_VEC, NcmStatsVecClass))

typedef struct _NcmStatsVecClass NcmStatsVecClass;
typedef struct _NcmStatsVec NcmStatsVec;

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
  NCM_STATS_VEC_COV,       /*< private >*/
  NCM_STATS_VEC_TYPES_LEN, /*< skip >*/
} NcmStatsVecType;

struct _NcmStatsVec
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsVecType t;
  guint len;
  gboolean save_x;
  gdouble weight;
  gdouble weight2;
  gdouble bias_wt;
  gdouble mean_inc;
  gdouble var_inc;
  gdouble cov_inc;
  guint nitens;
  NcmVector *x;
  NcmVector *mean;
  NcmVector *var;
  NcmMatrix *cov;
  GPtrArray *saved_x;
};

struct _NcmStatsVecClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_stats_vec_get_type (void) G_GNUC_CONST;

NcmStatsVec *ncm_stats_vec_new (guint len, NcmStatsVecType t, gboolean save_x);
NcmStatsVec *ncm_stats_vec_ref (NcmStatsVec *svec);
void ncm_stats_vec_free (NcmStatsVec *svec);
void ncm_stats_vec_clear (NcmStatsVec **svec);

void ncm_stats_vec_reset (NcmStatsVec *svec);
void ncm_stats_vec_update_weight (NcmStatsVec *svec, gdouble w);

G_INLINE_FUNC NcmVector *ncm_stats_vec_peek_x (NcmStatsVec *svec);
G_INLINE_FUNC void ncm_stats_vec_set (NcmStatsVec *svec, guint i, gdouble x_i);
G_INLINE_FUNC gdouble ncm_stats_vec_get (NcmStatsVec *svec, guint i);
G_INLINE_FUNC void ncm_stats_vec_update (NcmStatsVec *svec);
G_INLINE_FUNC gdouble ncm_stats_vec_get_mean (NcmStatsVec *svec, guint i);
G_INLINE_FUNC gdouble ncm_stats_vec_get_var (NcmStatsVec *svec, guint i);
G_INLINE_FUNC gdouble ncm_stats_vec_get_sd (NcmStatsVec *svec, guint i);
G_INLINE_FUNC gdouble ncm_stats_vec_get_cov (NcmStatsVec *svec, guint i, guint j);
G_INLINE_FUNC gdouble ncm_stats_vec_get_cor (NcmStatsVec *svec, guint i, guint j);
G_INLINE_FUNC void ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *mean);
G_INLINE_FUNC void ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m);

G_END_DECLS

#endif /* _NCM_STATS_VEC_H_ */

#ifndef _NCM_STATS_VEC_INLINE_H_
#define _NCM_STATS_VEC_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC NcmVector *
ncm_stats_vec_peek_x (NcmStatsVec *svec)
{
  return svec->x;
}

G_INLINE_FUNC void 
ncm_stats_vec_set (NcmStatsVec *svec, guint i, gdouble x_i)
{
  ncm_vector_fast_set (svec->x, i, x_i);
}

G_INLINE_FUNC gdouble 
ncm_stats_vec_get (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->x, i);
}

G_INLINE_FUNC void 
ncm_stats_vec_update (NcmStatsVec *svec)
{
  ncm_stats_vec_update_weight (svec, 1.0);
}

G_INLINE_FUNC gdouble 
ncm_stats_vec_get_mean (NcmStatsVec *svec, guint i)
{
  return ncm_vector_fast_get (svec->mean, i);
}

G_INLINE_FUNC gdouble 
ncm_stats_vec_get_var (NcmStatsVec *svec, guint i)
{
  g_assert (svec->t == NCM_STATS_VEC_VAR || svec->t == NCM_STATS_VEC_COV);
  return ncm_vector_fast_get (svec->var, i) * svec->bias_wt;
}

G_INLINE_FUNC gdouble 
ncm_stats_vec_get_sd (NcmStatsVec *svec, guint i)
{
  return sqrt (ncm_stats_vec_get_var (svec, i));
}

G_INLINE_FUNC gdouble 
ncm_stats_vec_get_cov (NcmStatsVec *svec, guint i, guint j)
{
  g_assert (svec->t == NCM_STATS_VEC_COV);
  if (i == j)
    return ncm_stats_vec_get_var (svec, i);
  else
    return ncm_matrix_get (svec->cov, i, j) * svec->bias_wt;
}

G_INLINE_FUNC gdouble 
ncm_stats_vec_get_cor (NcmStatsVec *svec, guint i, guint j)
{
  if (i == j)
    return 1.0;
  else
    return ncm_stats_vec_get_cov (svec, i, j) / (ncm_stats_vec_get_sd (svec, i) * ncm_stats_vec_get_sd (svec, j));
}

G_INLINE_FUNC void 
ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *x)
{
  g_assert (x != NULL);
  ncm_vector_memcpy (x, svec->mean);
}

G_INLINE_FUNC void 
ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m)
{
  guint i;
  g_assert (m != NULL);
  
  ncm_matrix_memcpy (m, svec->cov);
  for (i = 0; i < svec->len; i++)
    ncm_matrix_set (m, i, i, ncm_vector_fast_get (svec->var, i));

  ncm_matrix_scale (m, svec->bias_wt);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_STATS_VEC_INLINE_H_ */

