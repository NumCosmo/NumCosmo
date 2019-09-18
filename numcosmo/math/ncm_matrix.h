/***************************************************************************
 *            ncm_matrix.h
 *
 *  Thu January 05 20:18:45 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_MATRIX_H_
#define _NCM_MATRIX_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_cfg.h>
#include <numcosmo/math/ncm_util.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_rng.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_matrix.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_MATRIX             (ncm_matrix_get_type ())
#define NCM_MATRIX(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MATRIX, NcmMatrix))
#define NCM_MATRIX_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MATRIX, NcmMatrixClass))
#define NCM_IS_MATRIX(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MATRIX))
#define NCM_IS_MATRIX_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MATRIX))
#define NCM_MATRIX_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MATRIX, NcmMatrixClass))

typedef struct _NcmMatrixClass NcmMatrixClass;
typedef struct _NcmMatrix NcmMatrix;

/**
 * NcmMatrixInternal:
 * @NCM_MATRIX_SLICE: FIXME
 * @NCM_MATRIX_GSL_MATRIX: FIXME
 * @NCM_MATRIX_MALLOC: FIXME
 * @NCM_MATRIX_GARRAY: FIXME
 * @NCM_MATRIX_DERIVED: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcmMatrixInternal
{
  NCM_MATRIX_SLICE = 0,
  NCM_MATRIX_GSL_MATRIX,
  NCM_MATRIX_MALLOC,
  NCM_MATRIX_GARRAY,
  NCM_MATRIX_DERIVED,
} NcmMatrixInternal;

struct _NcmMatrixClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmMatrix:
 *
 * FIXME
 */
struct _NcmMatrix
{
  /*< private >*/
  GObject parent_instance;
  gsl_matrix_view mv;
  gpointer pdata;
  GDestroyNotify pfree;
  NcmMatrixInternal type;
};

GType ncm_matrix_get_type (void) G_GNUC_CONST;

NcmMatrix *ncm_matrix_new (const guint nrows, const guint ncols);
NcmMatrix *ncm_matrix_new0 (const guint nrows, const guint ncols);
NcmMatrix *ncm_matrix_new_full (gdouble *d, guint nrows, guint ncols, guint tda, gpointer pdata, GDestroyNotify pfree);
NcmMatrix *ncm_matrix_new_gsl (gsl_matrix *gm);
NcmMatrix *ncm_matrix_new_gsl_static (gsl_matrix *gm);
NcmMatrix *ncm_matrix_new_array (GArray *a, const guint ncols);
NcmMatrix *ncm_matrix_new_data_slice (gdouble *d, const guint nrows, const guint ncols);
NcmMatrix *ncm_matrix_new_data_malloc (gdouble *d, const guint nrows, const guint ncols);
NcmMatrix *ncm_matrix_new_data_static (gdouble *d, const guint nrows, const guint ncols);
NcmMatrix *ncm_matrix_new_data_static_tda (gdouble *d, const guint nrows, const guint ncols, const guint tda);
NcmMatrix *ncm_matrix_new_variant (GVariant *var);

NcmMatrix *ncm_matrix_ref (NcmMatrix *cm);

const NcmMatrix *ncm_matrix_const_new_data (const gdouble *d, guint nrows, guint ncols);
const NcmMatrix *ncm_matrix_const_new_variant (GVariant *var);

NcmMatrix *ncm_matrix_get_submatrix (NcmMatrix *cm, const guint k1, const guint k2, const guint nrows, const guint ncols);
NcmVector *ncm_matrix_get_col (NcmMatrix *cm, const guint col);
NcmVector *ncm_matrix_get_row (NcmMatrix *cm, const guint row);

void ncm_matrix_set_from_variant (NcmMatrix *cm, GVariant *var);
GVariant *ncm_matrix_get_variant (NcmMatrix *cm);
GVariant *ncm_matrix_peek_variant (NcmMatrix *cm);

void ncm_matrix_set_from_data (NcmMatrix *cm, gdouble *data);
void ncm_matrix_set_from_array (NcmMatrix *cm, GArray *a);

NCM_INLINE const NcmMatrix *ncm_matrix_new_gsl_const (gsl_matrix *m);
NCM_INLINE gdouble ncm_matrix_get (const NcmMatrix *cm, const guint i, const guint j);
NCM_INLINE gdouble *ncm_matrix_ptr (NcmMatrix *cm, const guint i, const guint j);
NCM_INLINE const gdouble *ncm_matrix_const_ptr (const NcmMatrix *cm, const guint i, const guint j);
NCM_INLINE GArray *ncm_matrix_get_array (NcmMatrix *cm);
NCM_INLINE void ncm_matrix_set (NcmMatrix *cm, const guint i, const guint j, const gdouble val);
NCM_INLINE void ncm_matrix_set_colmajor (NcmMatrix *cm, const guint i, const guint j, gdouble val);
NCM_INLINE void ncm_matrix_addto (NcmMatrix *cm, const guint i, const guint j, const gdouble val);
NCM_INLINE void ncm_matrix_transpose (NcmMatrix *cm);
NCM_INLINE void ncm_matrix_set_identity (NcmMatrix *cm);
NCM_INLINE void ncm_matrix_set_zero (NcmMatrix *cm);
NCM_INLINE void ncm_matrix_set_all (NcmMatrix *cm, const gdouble val);

NCM_INLINE void ncm_matrix_add (NcmMatrix *cm1, const NcmMatrix *cm2);
NCM_INLINE void ncm_matrix_sub (NcmMatrix *cm1, const NcmMatrix *cm2);
NCM_INLINE void ncm_matrix_mul_elements (NcmMatrix *cm1, const NcmMatrix *cm2);
NCM_INLINE void ncm_matrix_div_elements (NcmMatrix *cm1, const NcmMatrix *cm2);
NCM_INLINE void ncm_matrix_scale (NcmMatrix *cm, const gdouble val);
NCM_INLINE void ncm_matrix_add_constant (NcmMatrix *cm, const gdouble val);

NCM_INLINE void ncm_matrix_mul_row (NcmMatrix *cm, const guint row_i, const gdouble val);
NCM_INLINE void ncm_matrix_mul_col (NcmMatrix *cm, const guint col_i, const gdouble val);
NCM_INLINE void ncm_matrix_get_diag (NcmMatrix *cm, NcmVector *diag);
NCM_INLINE void ncm_matrix_set_diag (NcmMatrix *cm, NcmVector *diag);

NCM_INLINE void ncm_matrix_memcpy (NcmMatrix *cm1, const NcmMatrix *cm2);
NCM_INLINE void ncm_matrix_set_col (NcmMatrix *cm, const guint n, const NcmVector *cv);
NCM_INLINE void ncm_matrix_set_row (NcmMatrix *cm, const guint n, const NcmVector *cv);
NCM_INLINE gdouble ncm_matrix_fast_get (NcmMatrix *cm, const guint ij);
NCM_INLINE void ncm_matrix_fast_set (NcmMatrix *cm, const guint ij, const gdouble val);

NCM_INLINE gsl_matrix *ncm_matrix_gsl (NcmMatrix *cm);
NCM_INLINE const gsl_matrix *ncm_matrix_const_gsl (const NcmMatrix *cm);
NCM_INLINE guint ncm_matrix_col_len (const NcmMatrix *cm);
NCM_INLINE guint ncm_matrix_row_len (const NcmMatrix *cm);
NCM_INLINE guint ncm_matrix_nrows (const NcmMatrix *cm);
NCM_INLINE guint ncm_matrix_ncols (const NcmMatrix *cm);
NCM_INLINE guint ncm_matrix_size (const NcmMatrix *cm);
NCM_INLINE guint ncm_matrix_tda (const NcmMatrix *cm);
NCM_INLINE gdouble *ncm_matrix_data (NcmMatrix *cm);
NCM_INLINE const gdouble *ncm_matrix_const_data (const NcmMatrix *cm);

NcmMatrix *ncm_matrix_dup (const NcmMatrix *cm);
void ncm_matrix_substitute (NcmMatrix **cm, NcmMatrix *nm, gboolean check_size);
void ncm_matrix_add_mul (NcmMatrix *cm, const gdouble alpha, NcmMatrix *b);

gdouble ncm_matrix_cmp (const NcmMatrix *cm1, const NcmMatrix *cm2, const gdouble scale);
gdouble ncm_matrix_cmp_diag (const NcmMatrix *cm1, const NcmMatrix *cm2, const gdouble scale);

NcmMatrix *ncm_matrix_norma_diag (const NcmMatrix *cm1, NcmMatrix *cm2);

void ncm_matrix_free (NcmMatrix *cm);
void ncm_matrix_clear (NcmMatrix **cm);
void ncm_matrix_const_free (const NcmMatrix *cm);

void ncm_matrix_copy_triangle (NcmMatrix *cm, gchar UL);
void ncm_matrix_dsymm (NcmMatrix *cm, gchar UL, const gdouble alpha, NcmMatrix *A, NcmMatrix *B, const gdouble beta);
void ncm_matrix_dgemm (NcmMatrix *cm, gchar TransA, gchar TransB, const gdouble alpha, NcmMatrix *A, NcmMatrix *B, const gdouble beta);

gint ncm_matrix_cholesky_decomp (NcmMatrix *cm, gchar UL);
gint ncm_matrix_cholesky_inverse (NcmMatrix *cm, gchar UL);
gdouble ncm_matrix_cholesky_lndet (NcmMatrix *cm);
gint ncm_matrix_cholesky_solve (NcmMatrix *cm, NcmVector *b, gchar UL);
gint ncm_matrix_cholesky_solve2 (NcmMatrix *cm, NcmVector *b, gchar UL);
gint ncm_matrix_nearPD (NcmMatrix *cm, gchar UL, gboolean cholesky_decomp, const guint maxiter);
void ncm_matrix_sym_exp_cholesky (NcmMatrix *cm, gchar UL, NcmMatrix *exp_cm_dec);
void ncm_matrix_sym_posdef_log (NcmMatrix *cm, gchar UL, NcmMatrix *ln_cm);
void ncm_matrix_triang_to_sym (NcmMatrix *cm, gchar UL, gboolean zero, NcmMatrix *sym);
void ncm_matrix_log_vals (NcmMatrix *cm, gchar *prefix, gchar *format);

void ncm_matrix_fill_rand_cor (NcmMatrix *cm, const gdouble cor_level, NcmRNG *rng);
void ncm_matrix_fill_rand_cov (NcmMatrix *cm, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, NcmRNG *rng);
void ncm_matrix_fill_rand_cov2 (NcmMatrix *cm, NcmVector *mu, const gdouble reltol_min, const gdouble reltol_max, const gdouble cor_level, NcmRNG *rng);

G_END_DECLS

#endif /* _NCM_MATRIX_H_ */

#ifndef _NCM_MATRIX_INLINE_H_
#define _NCM_MATRIX_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE const NcmMatrix *
ncm_matrix_new_gsl_const (gsl_matrix *m)
{
  return ncm_matrix_new_data_static_tda ((m)->data, (m)->size1, (m)->size2, (m)->tda);
}

NCM_INLINE gdouble
ncm_matrix_get (const NcmMatrix *cm, const guint i, const guint j)
{
  return gsl_matrix_get (ncm_matrix_const_gsl (cm), i, j);
}

NCM_INLINE gdouble *
ncm_matrix_ptr (NcmMatrix *cm, guint i, guint j)
{
  return gsl_matrix_ptr (ncm_matrix_gsl (cm), i, j);
}

NCM_INLINE const gdouble *
ncm_matrix_const_ptr (const NcmMatrix *cm, guint i, guint j)
{
  return gsl_matrix_const_ptr (ncm_matrix_const_gsl (cm), i, j);
}

NCM_INLINE void
ncm_matrix_set (NcmMatrix *cm, const guint i, const guint j, gdouble val)
{
  gsl_matrix_set (ncm_matrix_gsl (cm), i, j, val);
}

NCM_INLINE void
ncm_matrix_set_colmajor (NcmMatrix *cm, const guint i, const guint j, gdouble val)
{
  ncm_matrix_data (cm) [i + ncm_matrix_nrows (cm) * j] = val;
}

NCM_INLINE void
ncm_matrix_addto (NcmMatrix *cm, guint i, guint j, gdouble val)
{
  gdouble *m = gsl_matrix_ptr (ncm_matrix_gsl (cm), i, j);
  *m += val;
}

NCM_INLINE void
ncm_matrix_transpose (NcmMatrix *cm)
{
  const gint ret = gsl_matrix_transpose (ncm_matrix_gsl (cm));
  NCM_TEST_GSL_RESULT ("gsl_matrix_transpose", ret);
}

NCM_INLINE void
ncm_matrix_set_identity (NcmMatrix *cm)
{
  gsl_matrix_set_identity (ncm_matrix_gsl (cm));
}

NCM_INLINE void
ncm_matrix_set_zero (NcmMatrix *cm)
{
  gsl_matrix_set_zero (ncm_matrix_gsl (cm));
}

NCM_INLINE void
ncm_matrix_set_all (NcmMatrix *cm, const gdouble val)
{
  gsl_matrix_set_all (ncm_matrix_gsl (cm), val);
}

NCM_INLINE void
ncm_matrix_add (NcmMatrix *cm1, const NcmMatrix *cm2)
{
  gsl_matrix_add (ncm_matrix_gsl (cm1), ncm_matrix_const_gsl (cm2));
}

NCM_INLINE void
ncm_matrix_sub (NcmMatrix *cm1, const NcmMatrix *cm2)
{
  gsl_matrix_sub (ncm_matrix_gsl (cm1), ncm_matrix_const_gsl (cm2));
}

NCM_INLINE void
ncm_matrix_mul_elements (NcmMatrix *cm1, const NcmMatrix *cm2)
{
  gsl_matrix_mul_elements (ncm_matrix_gsl (cm1), ncm_matrix_const_gsl (cm2));
}

NCM_INLINE void
ncm_matrix_div_elements (NcmMatrix *cm1, const NcmMatrix *cm2)
{
  gsl_matrix_div_elements (ncm_matrix_gsl (cm1), ncm_matrix_const_gsl (cm2));
}

NCM_INLINE void
ncm_matrix_scale (NcmMatrix *cm, const gdouble val)
{
  gsl_matrix_scale (ncm_matrix_gsl (cm), val);
}

NCM_INLINE void
ncm_matrix_add_constant (NcmMatrix *cm, const gdouble val)
{
  gsl_matrix_add_constant (ncm_matrix_gsl (cm), val);
}

NCM_INLINE void 
ncm_matrix_mul_row (NcmMatrix *cm, const guint row_i, const gdouble val)
{
  const guint ncols = ncm_matrix_ncols (cm);
  guint i;

  for (i = 0; i < ncols; i++)
  {
    ncm_matrix_ptr (cm, row_i, i)[0] *= val;
  }
}

NCM_INLINE void 
ncm_matrix_mul_col (NcmMatrix *cm, const guint col_i, const gdouble val)
{
  const guint nrows = ncm_matrix_nrows (cm);
  guint i;

  for (i = 0; i < nrows; i++)
  {
    ncm_matrix_ptr (cm, i, col_i)[0] *= val;
  }
}

NCM_INLINE void 
ncm_matrix_get_diag (NcmMatrix *cm, NcmVector *diag)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  const guint n     = MIN (nrows, ncols);
  guint i;

  g_assert_cmpuint (ncm_vector_len (diag), >=, n);
  
  for (i = 0; i < n; i++)
  {
    ncm_vector_set (diag, i, ncm_matrix_get (cm, i, i));
  }
}

NCM_INLINE void 
ncm_matrix_set_diag (NcmMatrix *cm, NcmVector *diag)
{
  const guint nrows = ncm_matrix_nrows (cm);
  const guint ncols = ncm_matrix_ncols (cm);
  const guint n     = MIN (nrows, ncols);
  guint i;

  g_assert_cmpuint (ncm_vector_len (diag), >=, n);
  
  for (i = 0; i < n; i++)
  {
    ncm_matrix_set (cm, i, i, ncm_vector_get (diag, i));
  }
}

NCM_INLINE void
ncm_matrix_memcpy (NcmMatrix *cm1, const NcmMatrix *cm2)
{
  const guint nrows = ncm_matrix_nrows (cm1);
  const guint ncols = ncm_matrix_ncols (cm1);
  const guint total = nrows * ncols;

  g_assert_cmpuint (nrows, ==, ncm_matrix_nrows (cm2));
  g_assert_cmpuint (ncols, ==, ncm_matrix_ncols (cm2));
  
  if ((ncm_matrix_tda (cm1) != ncols) || (ncm_matrix_tda (cm2) != ncm_matrix_ncols (cm2)))
  {
    register guint i;
    const guint ncols_bytes = sizeof (gdouble) * ncols;
    for (i = 0; i < nrows; i++)
      memcpy (ncm_matrix_ptr (cm1, i, 0), ncm_matrix_const_ptr (cm2, i, 0), ncols_bytes);  
  }
  else
    memcpy (ncm_matrix_data (cm1), ncm_matrix_const_data (cm2), sizeof (gdouble) * total);
}

NCM_INLINE void
ncm_matrix_set_col (NcmMatrix *cm, const guint n, const NcmVector *cv)
{
  gint ret = gsl_matrix_set_col (ncm_matrix_gsl (cm), n, ncm_vector_const_gsl (cv));
  g_assert (ret == GSL_SUCCESS);
}

NCM_INLINE void
ncm_matrix_set_row (NcmMatrix *cm, const guint n, const NcmVector *cv)
{
  gint ret = gsl_matrix_set_row (ncm_matrix_gsl (cm), n, ncm_vector_const_gsl (cv));
  g_assert (ret == GSL_SUCCESS);
}

NCM_INLINE GArray *
ncm_matrix_get_array (NcmMatrix *cm)
{
  g_assert (cm->type == NCM_MATRIX_GARRAY);
  return g_array_ref (cm->pdata);
}

NCM_INLINE gdouble 
ncm_matrix_fast_get (NcmMatrix *cm, const guint ij)
{
  return ncm_matrix_data (cm)[ij];
}

NCM_INLINE void 
ncm_matrix_fast_set (NcmMatrix *cm, const guint ij, const gdouble val)
{
  ncm_matrix_data (cm)[ij] = val;
}

NCM_INLINE gsl_matrix *
ncm_matrix_gsl (NcmMatrix *cm)
{ 
  return &cm->mv.matrix; 
}

NCM_INLINE const gsl_matrix *
ncm_matrix_const_gsl (const NcmMatrix *cm)
{
  return &(cm->mv.matrix);
}

NCM_INLINE guint 
ncm_matrix_col_len (const NcmMatrix *cm)
{
  return cm->mv.matrix.size1;
}

NCM_INLINE guint 
ncm_matrix_row_len (const NcmMatrix *cm)
{
  return cm->mv.matrix.size2;
}

NCM_INLINE guint 
ncm_matrix_nrows (const NcmMatrix *cm)
{
  return cm->mv.matrix.size1;
}

NCM_INLINE guint 
ncm_matrix_ncols (const NcmMatrix *cm)
{
  return cm->mv.matrix.size2;
}

NCM_INLINE guint 
ncm_matrix_size (const NcmMatrix *cm)
{
  return cm->mv.matrix.size1 * cm->mv.matrix.size2;
}

NCM_INLINE guint 
ncm_matrix_tda (const NcmMatrix *cm)
{
  return cm->mv.matrix.tda;
}

NCM_INLINE gdouble *
ncm_matrix_data (NcmMatrix *cm)
{
  return cm->mv.matrix.data;
}

NCM_INLINE const gdouble *
ncm_matrix_const_data (const NcmMatrix *cm)
{
  return cm->mv.matrix.data;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_MATRIX_INLINE_H_ */
