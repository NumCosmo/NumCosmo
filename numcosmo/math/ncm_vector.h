/***************************************************************************
 *            ncm_vector.h
 *
 *  Thu Nov  6 12:19:03 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NCM_VECTOR_H_
#define _NCM_VECTOR_H_

#include <string.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_vector.h>
#include <sundials/sundials_nvector.h>

G_BEGIN_DECLS

#define NCM_TYPE_VECTOR             (ncm_vector_get_type ())
#define NCM_VECTOR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_VECTOR, NcmVector))
#define NCM_VECTOR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_VECTOR, NcmVectorClass))
#define NCM_IS_VECTOR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_VECTOR))
#define NCM_IS_VECTOR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_VECTOR))
#define NCM_VECTOR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_VECTOR, NcmVectorClass))

typedef struct _NcmVectorClass NcmVectorClass;
typedef struct _NcmVector NcmVector;

/**
 * NcmVectorInternal:
 * @NCM_VECTOR_SLICE: FIXME
 * @NCM_VECTOR_GSL_VECTOR: FIXME
 * @NCM_VECTOR_MALLOC: FIXME
 * @NCM_VECTOR_ARRAY: FIXME
 * @NCM_VECTOR_DERIVED: FIXME
 *
 * FIXME
 *
 */
typedef enum _NcmVectorInternal
{
  NCM_VECTOR_SLICE = 0,
  NCM_VECTOR_GSL_VECTOR,
  NCM_VECTOR_MALLOC,
  NCM_VECTOR_ARRAY,
  NCM_VECTOR_DERIVED,
} NcmVectorInternal;

struct _NcmVectorClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmVector
{
  /*< private >*/
  GObject parent_instance;
  GArray *a;
  gsl_vector *gv;
  gsl_vector_view vv;
  GObject *pobj;
  NcmVectorInternal type;
};

GType ncm_vector_get_type (void) G_GNUC_CONST;

#define NCM_N2VECTOR(v) ((NcmVector *)((v)->content))
#define NCM_VECTOR_DATA(cv) ((cv)->vv.vector.data)

NcmVector *ncm_vector_new (gsize n);
NcmVector *ncm_vector_new_gsl (gsl_vector *gv);
NcmVector *ncm_vector_new_gsl_static (gsl_vector *gv);
NcmVector *ncm_vector_new_array (GArray *a);
NcmVector *ncm_vector_new_data_slice (gdouble *d, const gsize size, const gsize stride);
NcmVector *ncm_vector_new_data_malloc (gdouble *d, const gsize size, const gsize stride);
NcmVector *ncm_vector_new_data_static (gdouble *d, const gsize size, const gsize stride);
NcmVector *ncm_vector_ref (NcmVector *cv);
const NcmVector *ncm_vector_new_data_const (const gdouble *d, const gsize size, const gsize stride);

NcmVector *ncm_vector_get_subvector (NcmVector *cv, const gsize k, const gsize size);

G_INLINE_FUNC const NcmVector *ncm_vector_new_gsl_const (const gsl_vector *v);
G_INLINE_FUNC gdouble ncm_vector_get (const NcmVector *cv, const guint i);
G_INLINE_FUNC gdouble ncm_vector_fast_get (const NcmVector *cv, const guint i);
G_INLINE_FUNC gdouble *ncm_vector_ptr (NcmVector *cv, const guint i);
G_INLINE_FUNC void ncm_vector_set (NcmVector *cv, const guint i, const gdouble val);
G_INLINE_FUNC void ncm_vector_fast_set (NcmVector *cv, const guint i, const gdouble val);
G_INLINE_FUNC void ncm_vector_addto (NcmVector *cv, const guint i, const gdouble val);
G_INLINE_FUNC void ncm_vector_subfrom (NcmVector *cv, const guint i, const gdouble val);
G_INLINE_FUNC void ncm_vector_fast_subfrom (NcmVector *cv, const guint i, const gdouble val);
G_INLINE_FUNC void ncm_vector_set_all (NcmVector *cv, const gdouble val);
G_INLINE_FUNC void ncm_vector_scale (NcmVector *cv, const gdouble val);
G_INLINE_FUNC void ncm_vector_div (NcmVector *cv1, const NcmVector *cv2);
G_INLINE_FUNC void ncm_vector_add (NcmVector *cv1, const NcmVector *cv2);
G_INLINE_FUNC void ncm_vector_sub (NcmVector *cv1, const NcmVector *cv2);
G_INLINE_FUNC void ncm_vector_set_zero (NcmVector *cv);
G_INLINE_FUNC void ncm_vector_memcpy (NcmVector *cv1, const NcmVector *cv2);
G_INLINE_FUNC void ncm_vector_memcpy2 (NcmVector *cv1, const NcmVector *cv2, const guint cv1_start, const guint cv2_start, const guint size);
G_INLINE_FUNC GArray *ncm_vector_get_array (NcmVector *cv);
G_INLINE_FUNC GArray *ncm_vector_dup_array (NcmVector *cv);

G_INLINE_FUNC gsl_vector *ncm_vector_gsl (NcmVector *cv);
G_INLINE_FUNC const gsl_vector *ncm_vector_const_gsl (const NcmVector *cv);
G_INLINE_FUNC guint ncm_vector_len (const NcmVector *cv);
G_INLINE_FUNC guint ncm_vector_stride (const NcmVector *cv);

NcmVector *ncm_vector_dup (const NcmVector *cv);
void ncm_vector_free (NcmVector *cv);
void ncm_vector_clear (NcmVector **cv);
void ncm_vector_const_free (const NcmVector *cv);

N_Vector ncm_vector_nvector (NcmVector *cv);

G_END_DECLS

#endif /* _NCM_VECTOR_H_ */

#ifndef _NCM_VECTOR_INLINE_H_
#define _NCM_VECTOR_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC const NcmVector *
ncm_vector_new_gsl_const (const gsl_vector *v)
{
  return ncm_vector_new_data_static ((v)->data, (v)->size, (v)->stride);
}

G_INLINE_FUNC gdouble
ncm_vector_get (const NcmVector *cv, const guint i)
{
  return gsl_vector_get (ncm_vector_const_gsl (cv), i);
}

G_INLINE_FUNC gdouble
ncm_vector_fast_get (const NcmVector *cv, const guint i)
{
  return ncm_vector_const_gsl (cv)->data[i];
}

G_INLINE_FUNC gdouble *
ncm_vector_ptr (NcmVector *cv, const guint i)
{
  return gsl_vector_ptr (ncm_vector_gsl (cv), i);
}

G_INLINE_FUNC void
ncm_vector_set (NcmVector *cv, const guint i, const gdouble val)
{
  gsl_vector_set (ncm_vector_gsl (cv), i, val);
}

G_INLINE_FUNC void
ncm_vector_fast_set (NcmVector *cv, const guint i, const gdouble val)
{
  ncm_vector_gsl (cv)->data[i] = val;
}

G_INLINE_FUNC void
ncm_vector_addto (NcmVector *cv, const guint i, const gdouble val)
{
  gsl_vector_set (ncm_vector_gsl (cv), i, gsl_vector_get (ncm_vector_gsl (cv), i) + val);
}

G_INLINE_FUNC void
ncm_vector_subfrom (NcmVector *cv, const guint i, const gdouble val)
{
  gsl_vector_set (ncm_vector_gsl (cv), i, gsl_vector_get (ncm_vector_gsl (cv), i) - val);
}

G_INLINE_FUNC void
ncm_vector_fast_subfrom (NcmVector *cv, const guint i, const gdouble val)
{
  ncm_vector_gsl (cv)->data[i] = ncm_vector_gsl (cv)->data[i] - val;
}

G_INLINE_FUNC void
ncm_vector_set_all (NcmVector *cv, const gdouble val)
{
  gsl_vector_set_all (ncm_vector_gsl (cv), val);
}

G_INLINE_FUNC void
ncm_vector_scale (NcmVector *cv, const gdouble val)
{
  gsl_vector_scale (ncm_vector_gsl (cv), val);
}

G_INLINE_FUNC void
ncm_vector_div (NcmVector *cv1, const NcmVector *cv2)
{
  gsl_vector_div (ncm_vector_gsl(cv1), ncm_vector_const_gsl(cv2));
}

G_INLINE_FUNC void
ncm_vector_add (NcmVector *cv1, const NcmVector *cv2)
{
  gsl_vector_add (ncm_vector_gsl (cv1), ncm_vector_const_gsl (cv2));
}

G_INLINE_FUNC void
ncm_vector_sub (NcmVector *cv1, const NcmVector *cv2)
{
  gsl_vector_sub (ncm_vector_gsl (cv1), ncm_vector_const_gsl (cv2));
}

G_INLINE_FUNC void
ncm_vector_set_zero (NcmVector *cv)
{
  gsl_vector_set_zero (ncm_vector_gsl (cv));
}

G_INLINE_FUNC void
ncm_vector_memcpy (NcmVector *cv1, const NcmVector *cv2)
{
  gsl_vector_memcpy (ncm_vector_gsl (cv1), ncm_vector_const_gsl (cv2));
}

G_INLINE_FUNC void
ncm_vector_memcpy2 (NcmVector *cv1, const NcmVector *cv2, const guint cv1_start, const guint cv2_start, const guint size)
{
  memcpy (&NCM_VECTOR_DATA (cv1)[cv1_start], &NCM_VECTOR_DATA (cv2)[cv2_start], sizeof (gdouble) * size);
}

G_INLINE_FUNC GArray *
ncm_vector_get_array (NcmVector *cv)
{
  g_assert (cv->a != NULL);
  return g_array_ref (cv->a); 
}

G_INLINE_FUNC GArray *
ncm_vector_dup_array (NcmVector *cv)
{
  const guint len = ncm_vector_len (cv);
  GArray *a = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);
  g_array_append_vals (a, &ncm_vector_gsl (cv)->data[0], len);
  return a;
}

G_INLINE_FUNC gsl_vector *
ncm_vector_gsl (NcmVector *cv)
{
  return &(cv->vv.vector);
}

G_INLINE_FUNC const gsl_vector *
ncm_vector_const_gsl (const NcmVector *cv)
{
  return &(cv->vv.vector);
}

G_INLINE_FUNC guint
ncm_vector_len (const NcmVector *cv)
{
  return cv->vv.vector.size;
}

G_INLINE_FUNC guint
ncm_vector_stride (const NcmVector *cv)
{
  return cv->vv.vector.stride;
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_VECTOR_INLINE_H_ */
