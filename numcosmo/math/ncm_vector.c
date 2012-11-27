/***************************************************************************
 *            ncm_vector.c
 *
 *  Tue Jul  8 15:05:41 2008
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

/**
 * SECTION:ncm_vector
 * @title: Vector
 * @short_description: Allocation, access and operations.
 *
 * This object defines the functions for allocating and accessing vectors.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_vector.h"

G_DEFINE_TYPE (NcmVector, ncm_vector, G_TYPE_INITIALLY_UNOWNED);

/**
 * ncm_vector_new:
 * @n: defines the size of the vector.
 *
 * This function allocates memory for a new #NcmVector of double
 * with @n components.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new (gsize n)
{
  gdouble *d = g_slice_alloc (sizeof (gdouble) * n);
  return ncm_vector_new_data_slice (d, n, 1);
}

/**
 * ncm_vector_new_sunk:
 * @n: defines the size of the vector.
 *
 * This function allocates memory for a new #NcmVector of double
 * with @n components. Returns a sunk reference.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_sunk (gsize n)
{
  return ncm_vector_ref (ncm_vector_new (n));
}

/**
 * ncm_vector_new_gsl: (skip)
 * @gv: vector from GNU Scientific Library (GSL) to be converted into a #NcmVector.
 * 
 * This function saves @gv internally and frees it when it is no longer necessary.
 * The @gv vector must not be freed.
 * 
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_gsl (gsl_vector *gv)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  cv->gv = gv;
  cv->vv = gsl_vector_subvector (gv, 0, gv->size);
  cv->pobj = NULL;
  cv->a = NULL;
  cv->type = NCM_VECTOR_GSL_VECTOR;

  return cv;
}

/**
 * ncm_vector_new_gsl_static: (skip)
 * @gv: vector from GNU Scientific Library (GSL) to be converted into a #NcmVector.
 * 
 * This function saves @gv internally and does not frees.
 * The @gv vector must be valid during the life of the created #NcmVector.
 * 
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_gsl_static (gsl_vector *gv)
{
  NcmVector *cv = ncm_vector_new_gsl (gv);
  cv->type = NCM_VECTOR_DERIVED;
  return cv;
}

/**
 * ncm_vector_new_array:
 * @a: (array) (element-type double) (transfer full): array of doubles to be converted into a #NcmVector.
 *
 * This function saves @a internally and frees it when it is no longer necessary.
 * The @a array must not be freed.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_array (GArray *a)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  cv->vv = gsl_vector_view_array (&g_array_index (a, gdouble, 0), a->len);
  cv->a = g_array_ref (a);
  cv->pobj = NULL;
  cv->type = NCM_VECTOR_ARRAY;
  
  return cv;
}

/**
 * ncm_vector_new_data_slice:
 * @d: pointer to the first double allocated.
 * @size: number of doubles allocated.
 * @stride: the step-size from one element to the next in physical memory, measured in units of double.
 *
 * This function returns a #NcmVector of the array @d allocated using g_slice function.
 * This function saves @a internally and frees it when it is no longer necessary.
 * The @a vector must not be freed.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_slice (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  if (stride != 1)
	cv->vv = gsl_vector_view_array_with_stride (d, stride, size);
  else
	cv->vv = gsl_vector_view_array (d, size);
  cv->a = NULL;
  cv->pobj = NULL;
  cv->type = NCM_VECTOR_SLICE;
  return cv;
}

/**
 * ncm_vector_new_data_malloc:
 * @d: pointer to the first double allocated.
 * @size: number of doubles allocated.
 * @stride: the step-size from one element to the next in physical memory, measured in units of double.
 *
 * This function returns a #NcmVector of the array @d allocated using malloc.
 * It saves @d internally and frees it when it is no longer necessary.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_malloc (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_data_slice (d, size, stride);
  cv->type = NCM_VECTOR_MALLOC;
  return cv;
}

/**
 * ncm_vector_new_data_static:
 * @d: pointer to the first double allocated.
 * @size: number of doubles allocated.
 * @stride: the step-size from one element to the next in physical memory, measured in units of double.
 *
 * This function returns a #NcmVector of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new #NcmVector.
 */
NcmVector *
ncm_vector_new_data_static (gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = ncm_vector_new_data_slice (d, size, stride);
  cv->type = NCM_VECTOR_DERIVED;
  return cv;
}

/**
 * ncm_vector_new_data_const:
 * @d: pointer to the first double allocated.
 * @size: number of doubles allocated.
 * @stride: the step-size from one element to the next in physical memory, measured in units of double.
 *
 * This function returns a constant #NcmVector of the array @d.
 * The memory allocated is kept during all time life of the object and
 * must not be freed during this period.
 *
 * Returns: A new constant #NcmVector.
 */
const NcmVector *
ncm_vector_new_data_const (const gdouble *d, gsize size, gsize stride)
{
  NcmVector *cv = g_object_new (NCM_TYPE_VECTOR, NULL);
  if (stride != 1)
	cv->vv = gsl_vector_view_array_with_stride ((gdouble *)d, stride, size);
  else
	cv->vv = gsl_vector_view_array ((gdouble *)d, size);
  cv->a = NULL;
  cv->pobj = NULL;

  cv->type = NCM_VECTOR_DERIVED;
  return cv;
}

/**
 * ncm_vector_ref:
 * @cv: a NcmVector.
 *
 * This function increses the reference count of @cv or sink
 * the object.
 *
 * Returns: (transfer full): @cv
 */
NcmVector *
ncm_vector_ref (NcmVector *cv)
{
  return g_object_ref_sink (cv);
}

/**
 * ncm_vector_dup:
 * @cv: a constant #NcmVector.
 *
 * This function copies the elements of the constant vector @cv into a new #NcmVector.
 *
 * Returns: (transfer full): A #NcmVector.
   */
NcmVector *
ncm_vector_dup (const NcmVector *cv)
{
  NcmVector *cv_cp = ncm_vector_new (ncm_vector_len(cv));
  gsl_vector_memcpy (ncm_vector_gsl(cv_cp), ncm_vector_const_gsl(cv));
  return cv_cp;
}

/**
 * ncm_vector_get_subvector:
 * @cv: a #NcmVector.
 * @k: component index of the original vector.
 * @size: number of components of the subvector.
 *
 * This function returns a #NcmVector which is a subvector of the vector @cv.
 * The start of the new vector is the @k-th component from the original vector @cv.
 * The new vector has @size elements.
 *
 * Returns: (transfer full): A #NcmVector.
 */
NcmVector *
ncm_vector_get_subvector (NcmVector *cv, gsize k, gsize size)
{
  NcmVector *scv = g_object_new (NCM_TYPE_VECTOR, NULL);

  scv->vv = gsl_vector_subvector (ncm_vector_gsl (cv), k, size);
  scv->type = NCM_VECTOR_DERIVED;
  scv->pobj = G_OBJECT (ncm_vector_ref (cv));

  return scv;
}

/**
 * ncm_vector_new_gsl_const: (skip)
 * @v: vector from GNU Scientific Library (GSL).
 *
 * This function converts @v into a constant #NcmVector.
 *
 * Returns: A new constant #NcmVector.
 */
/**
 * ncm_vector_get:
 * @cv: a constant #NcmVector.
 * @i: component index.
 *
 * Returns: The @i-th component of the vector @cv.
 */
/**
 * ncm_vector_fast_get:
 * @cv: a constant #NcmVector.
 * @i: component index.
 *
 * Returns: The @i-th component of the vector @cv assuming stride == 1.
 */
/**
 * ncm_vector_ptr:
 * @cv: a #NcmVector.
 * @i: component index.
 *
 * Returns: A pointer to the @i-th component of the vector @cv.
 */
/**
 * ncm_vector_set:
 * @cv: a #NcmVector.
 * @i: component index.
 * @val: a constant double.
 *
 * This function sets the value of the @i-th component of the vector @cv to @val.
 */
/**
 * ncm_vector_fast_set:
 * @cv: a #NcmVector.
 * @i: component index.
 * @val: a constant double.
 *
 * This function sets the value of the @i-th component of the vector @cv to @val assuming stride == 1.
 */
/**
 * ncm_vector_addto:
 * @cv: a #NcmVector.
 * @i: component index.
 * @val: a constant double.
 *
 * This function adds @val to the value of the @i-th component of @cv.
 *
 */
/**
 * ncm_vector_subfrom:
 * @cv: a #NcmVector.
 * @i: component index.
 * @val: a cosntant double.
 *
 * This function subtracts @val from the value of the @i-th component of @cv.
 *
 */
/**
 * ncm_vector_fast_subfrom:
 * @cv: a #NcmVector.
 * @i: component index.
 * @val: a cosntant double.
 * 
 * This function subtracts @val from the value of the @i-th component of @cv assuming stride == 1.
 * 
 */
/**
 * ncm_vector_set_all:
 * @cv: a #NcmVector.
 * @val: a cosntant double.
 *
 * This function sets all the components of the vector @cv to the value @val.
 *
 */
/**
 * ncm_vector_scale:
 * @cv: a #NcmVector.
 * @val: a cosntant double.
 *
 * This function multiplies the components of the vector @cv by the constant factor @val.
 *
 */
/**
 * ncm_vector_div:
 * @cv1: a #NcmVector, numerator.
 * @cv2: a #NcmVector, denominator.
 *
 * This function divides the components of the vector @cv1 by the components of the vector @cv2.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_add:
 * @cv1: a #NcmVector.
 * @cv2: a #NcmVector.
 *
 * This function adds the components of the vector @cv2 to the components of the vector @cv1.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_sub:
 * @cv1: a #NcmVector.
 * @cv2: a #NcmVector.
 *
 * This function subtracts the components of the vector @cv2 to the components of the vector @cv1.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_set_zero:
 * @cv: a #NcmVector.
 *
 * This function sets all the components of the vector @cv to zero.
 *
 */
/**
 * ncm_vector_memcpy:
 * @cv1: a #NcmVector.
 * @cv2: a #NcmVector.
 *
 * This function copies the components of the vector @cv2 into the vector @cv1.
 * The two vectors must have the same length.
 *
 */
/**
 * ncm_vector_memcpy2:
 * @cv1: a #NcmVector.
 * @cv2: a #NcmVector.
 * @cv1_start: component of @cv1.
 * @cv2_start: component of @cv2.
 * @size: number of components.
 *
 * This function copies @size components of @cv2, counting from @cv2_start,
 * to the vector @cv1, starting from the @cv1_start component.
 * It is useful for vectors with different sizes.
 *
 */
/**
 * ncm_vector_get_array:
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: (transfer container) (element-type double): FIXME
 */
/**
 * ncm_vector_dup_array:
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: (transfer full) (element-type double): FIXME
 */
/**
 * ncm_vector_gsl: (skip)
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_const_gsl: (skip)
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_len:
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_vector_stride:
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: FIXME
 */

static void
_ncm_vector_dispose (GObject *object)
{
  NcmVector *cv = NCM_VECTOR (object);

  if (cv->a != NULL)
  {
    g_array_unref (cv->a);
    cv->a = NULL;
  }
  
  g_clear_object (&cv->pobj);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_vector_parent_class)->dispose (object);
}

static void
_ncm_vector_finalize (GObject *object)
{
  NcmVector *cv = NCM_VECTOR (object);
  switch (cv->type)
  {
	case NCM_VECTOR_SLICE:
	  g_slice_free1 (sizeof(gdouble) * ncm_vector_len (cv) * ncm_vector_stride (cv), NCM_VECTOR_DATA (cv));
	  NCM_VECTOR_DATA (cv) = NULL;
	  break;
	case NCM_VECTOR_ARRAY:
	  break;
	case NCM_VECTOR_MALLOC:
	  g_free (ncm_vector_gsl (cv)->data);
	  NCM_VECTOR_DATA (cv) = NULL;
	  break;
	case NCM_VECTOR_GSL_VECTOR:
	  gsl_vector_free (cv->gv);
	  cv->gv = NULL;
	  break;
	case NCM_VECTOR_DERIVED:
	  NCM_VECTOR_DATA (cv) = NULL;
	  break;
  }
  G_OBJECT_CLASS (ncm_vector_parent_class)->finalize (object);
}

/**
 * ncm_vector_free:
 * @cv: a #NcmVector.
 *
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released.
 *
 */
void
ncm_vector_free (NcmVector *cv)
{
  if (g_object_is_floating (cv))
    g_object_ref_sink (cv);
  g_object_unref (cv);
}

/**
 * ncm_vector_clear:
 * @cv: a #NcmVector.
 *
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released. The pointer is set to NULL.
 *
 */
void 
ncm_vector_clear (NcmVector **cv)
{
  if (*cv != NULL && g_object_is_floating (*cv))
    g_object_ref_sink (*cv);
  g_clear_object (cv);  
}

/**
 * ncm_vector_const_free:
 * @cv: a constant #NcmVector.
 *
 * Atomically decrements the reference count of @cv by one. If the reference count drops to 0,
 * all memory allocated by @cv is released.
 *
 */
void
ncm_vector_const_free (const NcmVector *cv)
{
  ncm_vector_free (NCM_VECTOR (cv));
}

static N_Vector
_ncm_nvclone(N_Vector nv)
{
  NcmVector *v = ncm_vector_new (ncm_vector_len(NCM_N2VECTOR(nv)));
  return ncm_vector_nvector (v);
}

static N_Vector
_ncm_nvcloneempty(N_Vector nv)
{
  return ncm_vector_nvector (NULL);
}

static void
_ncm_nvspace(N_Vector nv, glong *lrw, glong *liw)
{
  *lrw = ncm_vector_len(NCM_N2VECTOR(nv));
  *liw = (sizeof(NcmVector) % 4 == 0) ? sizeof(NcmVector)/4 : sizeof(NcmVector)/4 + 1;
}

static realtype *
_ncm_nvgetarraypointer(N_Vector nv)
{
  return NCM_VECTOR_DATA (NCM_N2VECTOR(nv));
}

static void
_ncm_nvsetarraypointer(realtype *data, N_Vector nv)
{
  NCM_VECTOR_DATA (NCM_N2VECTOR(nv)) = data;
}

static void
_ncm_nvlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{

}

static void
_ncm_nvconst(realtype a, N_Vector nv)
{
  gsl_vector_set_all (ncm_vector_gsl(NCM_N2VECTOR(nv)), a);
}

/**
 * ncm_vector_nvector: (skip)
 * @cv: a #NcmVector.
 *
 * FIXME
 *
 * Returns: FIXME
 */
N_Vector
ncm_vector_nvector (NcmVector *cv)
{
  struct _generic_N_Vector *nv = g_slice_new (struct _generic_N_Vector);
  g_object_ref_sink (cv);
  nv->content = cv;
  nv->ops = NULL;
  return nv;
}

static void
_ncm_vector_nvector_free (N_Vector nv)
{
  ncm_vector_free (NCM_N2VECTOR(nv));
  g_slice_free (struct _generic_N_Vector, nv);
}

static const struct _generic_N_Vector_Ops _ncm_ops =
{
  &_ncm_nvclone,
  &_ncm_nvcloneempty,
  &_ncm_vector_nvector_free,
  &_ncm_nvspace,
  &_ncm_nvgetarraypointer,
  &_ncm_nvsetarraypointer,
  &_ncm_nvlinearsum,
  &_ncm_nvconst
};

static void
ncm_vector_init (NcmVector *v)
{
  v->gv = NULL;
  v->a = NULL;
  v->pobj = NULL;
  v->type = 0;
  memset (&v->vv, 0, sizeof (gsl_vector_view));
}

static void
ncm_vector_class_init (NcmVectorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose = &_ncm_vector_dispose;
  object_class->finalize = &_ncm_vector_finalize;
}
