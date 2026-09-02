/***************************************************************************
 *            ncm_function_sample_set.c
 *
 *  Mon Mar 17 2026
 *  Copyright  2026 Sandro Dias Pinto Vitenti
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
 * NcmFunctionSampleSet:
 *
 * Ordered sample set for vector-valued functions $\vec{F}: \mathbb{R} \to \mathbb{R}^n$.
 *
 * Stores samples $(x_i, \vec{y}_i)$ in ascending $x_i$ order. Each interval
 * has an `interval_ok` flag for refinement status.
 *
 * The intended use is iterative refinement of splines built from
 * vector-valued functions:
 *
 * 1. Add initial samples with ncm_function_sample_set_add() or
 *    ncm_function_sample_set_add_func().
 * 2. Convert to #NcmSplineVec and test the interpolation error with
 *    ncm_function_sample_set_refine().
 * 3. Insert new samples where the error exceeds the tolerance, using the
 *    iterators.
 * 4. Mark an interval as `interval_ok` once it passes the error test.
 * 5. Repeat until ncm_function_sample_set_all_intervals_ok() returns TRUE.
 *
 * Iterators provide traversal, interval access, and insertion operations, with
 * O(1) access to sample data once positioned.
 *
 * |[<!-- language="C" -->
 * // Create iterator and traverse all samples (stack-allocated - no free needed)
 * NcmFunctionSampleSetIter iter_s;
 * NcmFunctionSampleSetIter *iter = &iter_s;
 * ncm_function_sample_set_iter_begin (fss, &iter);
 * while (ncm_function_sample_set_iter_is_valid (iter))
 * {
 *   gdouble x = ncm_function_sample_set_iter_get_x (iter);
 *   NcmVector *y = ncm_function_sample_set_iter_get_y (iter);
 *   // Process sample...
 *   ncm_function_sample_set_iter_next (iter);
 * }
 *
 * // Iterate over intervals using iter_next_pair (stack-allocated iterators)
 * NcmFunctionSampleSetIter it_s, next_it_s;
 * NcmFunctionSampleSetIter *it = &it_s;
 * NcmFunctionSampleSetIter *next_it = &next_it_s;
 * ncm_function_sample_set_iter_begin (fss, &it);
 * ncm_function_sample_set_iter_copy (it, &next_it);
 * while (ncm_function_sample_set_iter_next_pair (it, next_it))
 * {
 *   if (ncm_function_sample_set_iter_get_interval_ok (it) < threshold)
 *   {
 *     gdouble x_mid = 0.5 * (ncm_function_sample_set_iter_get_x (it) +
 *                            ncm_function_sample_set_iter_get_x (next_it));
 *     NcmFunctionSampleSetIter new_it_s;
 *     NcmFunctionSampleSetIter *new_it = &new_it_s;
 *     ncm_function_sample_set_iter_insert_after_func (fss, it, x_mid, func, NULL, &new_it);
 *   }
 *   ncm_function_sample_set_iter_next (it);
 * }
 * ]|
 *
 * Conversion to #NcmSplineVec reuses internal arrays and invalidates the
 * previously returned spline. Duplicate the spline to retain it.
 *
 * The vector dimension $n$ is fixed at creation and checked on insertion.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_function_sample_set.h"
#include "ncm/core/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_LEN,
  PROP_TRACK_RESIDUAL,
};

typedef struct _NcmFunctionSamplePoint
{
  gdouble x;
  NcmVector *y;
  NcmVector *residual;
  gint interval_ok;
  gboolean new_point;
} NcmFunctionSamplePoint;

struct _NcmFunctionSampleSet
{
  /*< private >*/
  GObject parent_instance;
  guint len;
  GList *samples;
  GArray *x_array_cache;
  GPtrArray *y_arrays_cache;
  GArray *absmaxF;   /* Maximum absolute value for each component */
  GArray *absmaxF_x; /* The x position of the maximum absolute value for each component */
  gdouble x_min;
  gdouble x_max;
  gboolean track_residual;
};

G_DEFINE_TYPE (NcmFunctionSampleSet, ncm_function_sample_set, G_TYPE_OBJECT)

static void
_ncm_function_sample_point_free (gpointer data)
{
  NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) data;

  ncm_vector_clear (&sp->y);
  ncm_vector_clear (&sp->residual);
  g_slice_free (NcmFunctionSamplePoint, sp);
}

static NcmFunctionSamplePoint *
_ncm_function_sample_point_new (const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp = g_slice_new (NcmFunctionSamplePoint);

  sp->x           = x;
  sp->y           = ncm_vector_ref (y);
  sp->residual    = NULL;
  sp->interval_ok = 0;
  sp->new_point   = TRUE;

  return sp;
}

static NcmFunctionSamplePoint *
_ncm_function_sample_point_new_old (const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp = g_slice_new (NcmFunctionSamplePoint);

  sp->x           = x;
  sp->y           = ncm_vector_ref (y);
  sp->residual    = NULL;
  sp->interval_ok = 0;
  sp->new_point   = FALSE;

  return sp;
}

static gint
_ncm_function_sample_point_cmp (gconstpointer a, gconstpointer b)
{
  const NcmFunctionSamplePoint *spa = (const NcmFunctionSamplePoint *) a;
  const NcmFunctionSamplePoint *spb = (const NcmFunctionSamplePoint *) b;

  if (spa->x < spb->x)
    return -1;
  else if (spa->x > spb->x)
    return 1;
  else
    return 0;
}

static void
ncm_function_sample_set_init (NcmFunctionSampleSet *fss)
{
  fss->len            = 0;
  fss->samples        = NULL;
  fss->x_array_cache  = NULL;
  fss->y_arrays_cache = NULL;
  fss->absmaxF        = NULL;
  fss->absmaxF_x      = NULL;
  fss->x_min          = GSL_POSINF;
  fss->x_max          = GSL_NEGINF;
  fss->track_residual = FALSE;
}

static void
_ncm_function_sample_set_dispose (GObject *object)
{
  NcmFunctionSampleSet *fss = NCM_FUNCTION_SAMPLE_SET (object);

  g_list_free_full (fss->samples, _ncm_function_sample_point_free);
  fss->samples = NULL;

  g_clear_pointer (&fss->x_array_cache, g_array_unref);
  g_clear_pointer (&fss->y_arrays_cache, g_ptr_array_unref);
  g_clear_pointer (&fss->absmaxF, g_array_unref);
  g_clear_pointer (&fss->absmaxF_x, g_array_unref);
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_function_sample_set_parent_class)->dispose (object);
}

static void
_ncm_function_sample_set_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_function_sample_set_parent_class)->finalize (object);
}

static void
_ncm_function_sample_set_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFunctionSampleSet *fss = NCM_FUNCTION_SAMPLE_SET (object);

  g_return_if_fail (NCM_IS_FUNCTION_SAMPLE_SET (object));

  switch (prop_id)
  {
    case PROP_LEN:
      fss->len = g_value_get_uint (value);
      break;
    case PROP_TRACK_RESIDUAL:
      ncm_function_sample_set_set_track_residual (fss, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_function_sample_set_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFunctionSampleSet *fss = NCM_FUNCTION_SAMPLE_SET (object);

  g_return_if_fail (NCM_IS_FUNCTION_SAMPLE_SET (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, fss->len);
      break;
    case PROP_TRACK_RESIDUAL:
      g_value_set_boolean (value, ncm_function_sample_set_get_track_residual (fss));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_function_sample_set_constructed (GObject *object)
{
  NcmFunctionSampleSet *fss = NCM_FUNCTION_SAMPLE_SET (object);
  guint i;

  /* Chain up : start */
  G_OBJECT_CLASS (ncm_function_sample_set_parent_class)->constructed (object);

  /* Initialize cached arrays for efficient spline building */
  fss->x_array_cache  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 100);
  fss->y_arrays_cache = g_ptr_array_new_full (fss->len, (GDestroyNotify) g_array_unref);
  fss->absmaxF        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), fss->len);
  fss->absmaxF_x      = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), fss->len);

  g_array_set_size (fss->absmaxF, fss->len);
  g_array_set_size (fss->absmaxF_x, fss->len);

  /* Initialize absmaxF to 0.0 for all components */
  for (i = 0; i < fss->len; i++)
  {
    GArray *y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 100);

    g_ptr_array_add (fss->y_arrays_cache, y_array);
    g_array_index (fss->absmaxF, gdouble, i)   = 0.0;
    g_array_index (fss->absmaxF_x, gdouble, i) = GSL_NAN;
  }
}

static void
ncm_function_sample_set_class_init (NcmFunctionSampleSetClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_function_sample_set_constructed;
  object_class->dispose      = &_ncm_function_sample_set_dispose;
  object_class->finalize     = &_ncm_function_sample_set_finalize;
  object_class->set_property = &_ncm_function_sample_set_set_property;
  object_class->get_property = &_ncm_function_sample_set_get_property;

  /**
   * NcmFunctionSampleSet:len:
   *
   * The dimension of the vector-valued function output.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      NULL,
                                                      "Vector dimension",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFunctionSampleSet:track-residual:
   *
   * Whether ncm_function_sample_set_refine() records the residual it actually
   * achieved on each interval, instead of only the pass/fail bit counted by
   * interval_ok. Off by default: the record costs one extra vector of length
   * #NcmFunctionSampleSet:len per sample, and only a caller that reports an
   * error estimate needs it. See ncm_function_sample_set_get_residuals().
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TRACK_RESIDUAL,
                                   g_param_spec_boolean ("track-residual",
                                                         NULL,
                                                         "Whether to record the achieved refinement residual",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_function_sample_set_new:
 * @len: dimension of the vector-valued function output
 *
 * Creates a new #NcmFunctionSampleSet for storing samples of a vector-valued
 * function with output dimension @len.
 *
 * Returns: (transfer full): a new #NcmFunctionSampleSet
 */
NcmFunctionSampleSet *
ncm_function_sample_set_new (const guint len)
{
  NcmFunctionSampleSet *fss = g_object_new (NCM_TYPE_FUNCTION_SAMPLE_SET,
                                            "len", len,
                                            NULL);

  return fss;
}

/**
 * ncm_function_sample_set_ref:
 * @fss: a #NcmFunctionSampleSet
 *
 * Increases the reference count of @fss.
 *
 * Returns: (transfer full): @fss
 */
NcmFunctionSampleSet *
ncm_function_sample_set_ref (NcmFunctionSampleSet *fss)
{
  return g_object_ref (fss);
}

/**
 * ncm_function_sample_set_free:
 * @fss: a #NcmFunctionSampleSet
 *
 * Decreases the reference count of @fss. If the reference count reaches zero,
 * @fss is freed.
 *
 */
void
ncm_function_sample_set_free (NcmFunctionSampleSet *fss)
{
  g_object_unref (fss);
}

/**
 * ncm_function_sample_set_clear:
 * @fss: a #NcmFunctionSampleSet
 *
 * Decreases the reference count of *@fss and sets the pointer *@fss to NULL.
 *
 */
void
ncm_function_sample_set_clear (NcmFunctionSampleSet **fss)
{
  g_clear_object (fss);
}

/* ============================================================================
 * Iterator implementation
 * ============================================================================ */

static NcmFunctionSampleSetIter *
_ncm_function_sample_set_iter_copy (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSampleSetIter *copy = g_slice_new (NcmFunctionSampleSetIter);

  *copy = *iter;

  return copy;
}

static void
_ncm_function_sample_set_iter_free (NcmFunctionSampleSetIter *iter)
{
  g_slice_free (NcmFunctionSampleSetIter, iter);
}

G_DEFINE_BOXED_TYPE (NcmFunctionSampleSetIter, ncm_function_sample_set_iter,
                     _ncm_function_sample_set_iter_copy, _ncm_function_sample_set_iter_free)

/**
 * ncm_function_sample_set_iter_begin:
 * @fss: a #NcmFunctionSampleSet
 * @iter_out: (out callee-allocates): iterator set to the first sample; if *@iter_out
 *   is %NULL a new iterator is heap-allocated, otherwise the existing memory is reused
 *
 * Positions an iterator at the first sample in @fss.
 * If @fss is empty, the iterator will be invalid.
 * When @iter_out points to a %NULL pointer the callee allocates a new iterator
 * that must be freed with ncm_function_sample_set_iter_free().
 * When @iter_out points to an already-allocated iterator (e.g. a stack variable)
 * no allocation occurs and no free is required.
 *
 */
void
ncm_function_sample_set_iter_begin (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter **iter_out)
{
  if (*iter_out == NULL)
    *iter_out = g_slice_new (NcmFunctionSampleSetIter);

  (*iter_out)->node  = fss->samples;
  (*iter_out)->owner = fss;
}

/**
 * ncm_function_sample_set_iter_end:
 * @fss: a #NcmFunctionSampleSet
 * @iter_out: (out callee-allocates): iterator set to the last sample; if *@iter_out
 *   is %NULL a new iterator is heap-allocated, otherwise the existing memory is reused
 *
 * Positions an iterator at the last sample in @fss.
 * If @fss is empty, the iterator will be invalid.
 * When @iter_out points to a %NULL pointer the callee allocates a new iterator
 * that must be freed with ncm_function_sample_set_iter_free().
 * When @iter_out points to an already-allocated iterator (e.g. a stack variable)
 * no allocation occurs and no free is required.
 *
 */
void
ncm_function_sample_set_iter_end (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter **iter_out)
{
  if (*iter_out == NULL)
    *iter_out = g_slice_new (NcmFunctionSampleSetIter);

  (*iter_out)->node  = g_list_last (fss->samples);
  (*iter_out)->owner = fss;
}

/**
 * ncm_function_sample_set_iter_copy:
 * @iter: a #NcmFunctionSampleSetIter
 * @iter_out: (out callee-allocates): copy of @iter pointing to the same position;
 *   if *@iter_out is %NULL a new iterator is heap-allocated, otherwise the
 *   existing memory is reused
 *
 * Copies the position of @iter into @iter_out.
 * When @iter_out points to a %NULL pointer the callee allocates a new iterator
 * that must be freed with ncm_function_sample_set_iter_free().
 * When @iter_out points to an already-allocated iterator (e.g. a stack variable)
 * no allocation occurs and no free is required.
 *
 */
void
ncm_function_sample_set_iter_copy (NcmFunctionSampleSetIter *iter, NcmFunctionSampleSetIter **iter_out)
{
  if (*iter_out == NULL)
    *iter_out = g_slice_new (NcmFunctionSampleSetIter);

  **iter_out = *iter;
}

/**
 * ncm_function_sample_set_iter_free:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Frees a heap-allocated iterator. Only call this for iterators allocated by the
 * callee (i.e. when a %NULL pointer was passed as the @iter_out argument to
 * ncm_function_sample_set_iter_begin(), ncm_function_sample_set_iter_end(),
 * ncm_function_sample_set_iter_copy(), or the insert functions). Stack-allocated
 * iterators must not be freed with this function.
 *
 */
void
ncm_function_sample_set_iter_free (NcmFunctionSampleSetIter *iter)
{
  _ncm_function_sample_set_iter_free (iter);
}

/**
 * ncm_function_sample_set_iter_is_valid:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Checks if @iter points to a valid sample.
 *
 * Returns: TRUE if @iter is valid, FALSE otherwise
 */
gboolean
ncm_function_sample_set_iter_is_valid (NcmFunctionSampleSetIter *iter)
{
  return iter->node != NULL;
}

/**
 * ncm_function_sample_set_iter_has_next:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Checks if there is a next sample after @iter.
 *
 * Returns: TRUE if there is a next sample, FALSE otherwise
 */
gboolean
ncm_function_sample_set_iter_has_next (NcmFunctionSampleSetIter *iter)
{
  return iter->node != NULL && iter->node->next != NULL;
}

/**
 * ncm_function_sample_set_iter_has_prev:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Checks if there is a previous sample before @iter.
 *
 * Returns: TRUE if there is a previous sample, FALSE otherwise
 */
gboolean
ncm_function_sample_set_iter_has_prev (NcmFunctionSampleSetIter *iter)
{
  return iter->node != NULL && iter->node->prev != NULL;
}

/**
 * ncm_function_sample_set_iter_next:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Moves @iter to the next sample. If @iter is at the last sample or invalid,
 * @iter becomes invalid.
 *
 */
void
ncm_function_sample_set_iter_next (NcmFunctionSampleSetIter *iter)
{
  g_assert (iter->node != NULL);
  iter->node = iter->node->next;
}

/**
 * ncm_function_sample_set_iter_prev:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Moves @iter to the previous sample. If @iter is at the first sample or invalid,
 * behavior is undefined.
 *
 */
void
ncm_function_sample_set_iter_prev (NcmFunctionSampleSetIter *iter)
{
  g_assert (iter->node != NULL);
  iter->node = iter->node->prev;
}

/**
 * ncm_function_sample_set_iter_get_x:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Gets the x value of the sample pointed to by @iter.
 *
 * Returns: the x value
 */
gdouble
ncm_function_sample_set_iter_get_x (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  return sp->x;
}

/**
 * ncm_function_sample_set_iter_get_y:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Gets the y vector of the sample pointed to by @iter.
 *
 * Returns: (transfer none): the y vector
 */
NcmVector *
ncm_function_sample_set_iter_get_y (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  return sp->y;
}

/**
 * ncm_function_sample_set_iter_get_interval_ok:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Gets the interval_ok flag of the sample pointed to by @iter.
 * This indicates whether the interval from this sample to the next has
 * passed refinement tests.
 *
 * Returns: the interval_ok value
 */
gint
ncm_function_sample_set_iter_get_interval_ok (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  return sp->interval_ok;
}

/**
 * ncm_function_sample_set_iter_get_residual:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Gets the residual recorded for the interval running from the sample pointed
 * to by @iter to the next one, or %NULL when that interval carries no record.
 * See ncm_function_sample_set_get_residuals() for what the value means.
 *
 * Returns: (transfer none) (nullable): the residual vector
 */
NcmVector *
ncm_function_sample_set_iter_get_residual (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  return sp->residual;
}

/**
 * ncm_function_sample_set_iter_get_new_point:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Gets the new_point flag of the sample pointed to by @iter.
 *
 * Returns: TRUE if the sample is marked as new, FALSE otherwise
 */
gboolean
ncm_function_sample_set_iter_get_new_point (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  return sp->new_point;
}

/**
 * ncm_function_sample_set_iter_set_interval_ok:
 * @iter: a #NcmFunctionSampleSetIter
 * @interval_ok: new interval_ok value
 *
 * Sets the interval_ok flag of the sample pointed to by @iter.
 *
 */
void
ncm_function_sample_set_iter_set_interval_ok (NcmFunctionSampleSetIter *iter, const gint interval_ok)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  sp->interval_ok = interval_ok;
}

/**
 * ncm_function_sample_set_iter_inc_interval_ok:
 * @iter: a #NcmFunctionSampleSetIter
 *
 * Increments the interval_ok flag of the sample pointed to by @iter.
 *
 */
void
ncm_function_sample_set_iter_inc_interval_ok (NcmFunctionSampleSetIter *iter)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  sp->interval_ok++;
}

/**
 * ncm_function_sample_set_iter_set_new_point:
 * @iter: a #NcmFunctionSampleSetIter
 * @new_point: new new_point flag value
 *
 * Sets the new_point flag of the sample pointed to by @iter.
 *
 */
void
ncm_function_sample_set_iter_set_new_point (NcmFunctionSampleSetIter *iter, const gboolean new_point)
{
  NcmFunctionSamplePoint *sp;

  g_assert (iter->node != NULL);
  sp = (NcmFunctionSamplePoint *) iter->node->data;

  sp->new_point = new_point;
}

/**
 * ncm_function_sample_set_iter_insert_after:
 * @fss: a #NcmFunctionSampleSet
 * @iter: a #NcmFunctionSampleSetIter
 * @x: knot position
 * @y: vector value at @x
 * @iter_out: (out callee-allocates): iterator pointing to the newly inserted sample;
 *   if *@iter_out is %NULL a new iterator is heap-allocated, otherwise the
 *   existing memory is reused
 *
 * Inserts a new sample after the position of @iter. The vector @y must have dimension
 * matching @fss:len property. The new sample's interval_ok is initialized to 0 and
 * new_point is set to TRUE. When @iter_out points to a %NULL pointer the callee
 * allocates a new iterator that must be freed with
 * ncm_function_sample_set_iter_free(). When @iter_out points to an already-allocated
 * iterator (e.g. a stack variable) no allocation occurs and no free is required.
 *
 * Warning: This does not check if the insertion maintains x-order. Use with care.
 *
 */
void
ncm_function_sample_set_iter_insert_after (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmVector *y, NcmFunctionSampleSetIter **iter_out)
{
  NcmFunctionSamplePoint *sp;
  guint i;

  g_assert (iter->node != NULL);
  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  sp = _ncm_function_sample_point_new (x, y);

  fss->samples = g_list_insert_before (fss->samples, iter->node->next, sp);

  if (*iter_out == NULL)
    *iter_out = g_slice_new (NcmFunctionSampleSetIter);

  (*iter_out)->node  = iter->node->next;
  (*iter_out)->owner = fss;

  /* Update x_min and x_max tracking */
  fss->x_min = MIN (fss->x_min, x);
  fss->x_max = MAX (fss->x_max, x);

  /* Update absmaxF tracking for each component */
  for (i = 0; i < fss->len; i++)
  {
    gdouble y_i        = ncm_vector_get (y, i);
    gdouble abs_y_i    = fabs (y_i);
    gdouble *absmaxF_i = &g_array_index (fss->absmaxF, gdouble, i);

    if (abs_y_i > *absmaxF_i)
    {
      *absmaxF_i                                 = abs_y_i;
      g_array_index (fss->absmaxF_x, gdouble, i) = x;
    }
  }
}

/**
 * ncm_function_sample_set_iter_insert_after_func:
 * @fss: a #NcmFunctionSampleSet
 * @iter: a #NcmFunctionSampleSetIter
 * @x: knot position
 * @f: (scope call): function to evaluate at @x
 * @user_data: user data to pass to @f
 * @iter_out: (out callee-allocates): iterator pointing to the newly inserted sample;
 *   if *@iter_out is %NULL a new iterator is heap-allocated, otherwise the
 *   existing memory is reused
 *
 * Evaluates @f at @x and inserts the result as a new sample after @iter. The new
 * sample's interval_ok is initialized to 0 and new_point is set to TRUE. When
 * @iter_out points to a %NULL pointer the callee allocates a new iterator that must be
 * freed with ncm_function_sample_set_iter_free(). When @iter_out points to an
 * already-allocated iterator (e.g. a stack variable) no allocation occurs and no free
 * is required.
 *
 * Warning: This does not check if the insertion maintains x-order. Use with care.
 *
 */
void
ncm_function_sample_set_iter_insert_after_func (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data, NcmFunctionSampleSetIter **iter_out)
{
  NcmVector *y = ncm_vector_new (fss->len);

  f (x, y, user_data);
  ncm_function_sample_set_iter_insert_after (fss, iter, x, y, iter_out);
  ncm_vector_free (y);
}

/**
 * ncm_function_sample_set_iter_insert_before:
 * @fss: a #NcmFunctionSampleSet
 * @iter: a #NcmFunctionSampleSetIter
 * @x: knot position
 * @y: vector value at @x
 * @iter_out: (out callee-allocates): iterator pointing to the newly inserted sample;
 *   if *@iter_out is %NULL a new iterator is heap-allocated, otherwise the
 *   existing memory is reused
 *
 * Inserts a new sample before the position of @iter. The vector @y must have dimension
 * matching @fss:len property. The new sample's interval_ok is initialized to 0 and
 * new_point is set to TRUE. When @iter_out points to a %NULL pointer the callee
 * allocates a new iterator that must be freed with
 * ncm_function_sample_set_iter_free(). When @iter_out points to an already-allocated
 * iterator (e.g. a stack variable) no allocation occurs and no free is required.
 *
 * Warning: This does not check if the insertion maintains x-order. Use with care.
 *
 */
void
ncm_function_sample_set_iter_insert_before (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmVector *y, NcmFunctionSampleSetIter **iter_out)
{
  NcmFunctionSamplePoint *sp;
  guint i;

  g_assert (iter->node != NULL);
  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  sp = _ncm_function_sample_point_new (x, y);

  fss->samples = g_list_insert_before (fss->samples, iter->node, sp);

  if (*iter_out == NULL)
    *iter_out = g_slice_new (NcmFunctionSampleSetIter);

  (*iter_out)->node  = iter->node->prev;
  (*iter_out)->owner = fss;

  /* Update x_min and x_max tracking */
  fss->x_min = MIN (fss->x_min, x);
  fss->x_max = MAX (fss->x_max, x);

  /* Update absmaxF tracking for each component */
  for (i = 0; i < fss->len; i++)
  {
    gdouble y_i        = ncm_vector_get (y, i);
    gdouble abs_y_i    = fabs (y_i);
    gdouble *absmaxF_i = &g_array_index (fss->absmaxF, gdouble, i);

    if (abs_y_i > *absmaxF_i)
    {
      *absmaxF_i                                 = abs_y_i;
      g_array_index (fss->absmaxF_x, gdouble, i) = x;
    }
  }
}

/**
 * ncm_function_sample_set_iter_insert_before_func:
 * @fss: a #NcmFunctionSampleSet
 * @iter: a #NcmFunctionSampleSetIter
 * @x: knot position
 * @f: (scope call): function to evaluate at @x
 * @user_data: user data to pass to @f
 * @iter_out: (out callee-allocates): iterator pointing to the newly inserted sample;
 *   if *@iter_out is %NULL a new iterator is heap-allocated, otherwise the
 *   existing memory is reused
 *
 * Evaluates @f at @x and inserts the result as a new sample before @iter. The new
 * sample's interval_ok is initialized to 0 and new_point is set to TRUE. When
 * @iter_out points to a %NULL pointer the callee allocates a new iterator that must be
 * freed with ncm_function_sample_set_iter_free(). When @iter_out points to an
 * already-allocated iterator (e.g. a stack variable) no allocation occurs and no free
 * is required.
 *
 * Warning: This does not check if the insertion maintains x-order. Use with care.
 *
 */
void
ncm_function_sample_set_iter_insert_before_func (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data, NcmFunctionSampleSetIter **iter_out)
{
  NcmVector *y = ncm_vector_new (fss->len);

  f (x, y, user_data);
  ncm_function_sample_set_iter_insert_before (fss, iter, x, y, iter_out);
  ncm_vector_free (y);
}

/**
 * ncm_function_sample_set_iter_next_pair:
 * @iter: a #NcmFunctionSampleSetIter
 * @next_iter: (out): a #NcmFunctionSampleSetIter to be set to the next position
 *
 * Helper for interval operations. Sets @next_iter to point to the sample after @iter.
 * This is useful for iterating over intervals where you need both endpoints.
 *
 * Returns: TRUE if both iterators are valid (there is an interval), FALSE otherwise
 */
gboolean
ncm_function_sample_set_iter_next_pair (NcmFunctionSampleSetIter *iter, NcmFunctionSampleSetIter *next_iter)
{
  if ((iter->node == NULL) || (iter->node->next == NULL))
    return FALSE;

  next_iter->node  = iter->node->next;
  next_iter->owner = iter->owner;

  return TRUE;
}

/**
 * ncm_function_sample_set_add:
 * @fss: a #NcmFunctionSampleSet
 * @x: knot position
 * @y: vector value at @x
 *
 * Adds a new sample point to @fss with position @x and vector value @y. The sample is
 * inserted in the correct position to maintain ascending x-order. The vector @y is
 * copied and must have dimension matching the @fss:len property. The interval_ok flag
 * for the new sample is initialized to 0. The sample is marked as a new point.
 *
 */
void
ncm_function_sample_set_add (NcmFunctionSampleSet *fss, const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp;
  guint i;

  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  sp           = _ncm_function_sample_point_new (x, y);
  fss->samples = g_list_insert_sorted (fss->samples, sp, _ncm_function_sample_point_cmp);

  /* Update x_min and x_max */
  fss->x_min = MIN (fss->x_min, x);
  fss->x_max = MAX (fss->x_max, x);

  /* Update absmaxF for each component */
  for (i = 0; i < fss->len; i++)
  {
    const gdouble abs_yi = fabs (ncm_vector_get (y, i));
    gdouble current_max  = g_array_index (fss->absmaxF, gdouble, i);

    if (abs_yi > current_max)
    {
      g_array_index (fss->absmaxF, gdouble, i)   = abs_yi;
      g_array_index (fss->absmaxF_x, gdouble, i) = x;
    }
  }
}

/**
 * ncm_function_sample_set_add_func:
 * @fss: a #NcmFunctionSampleSet
 * @x: knot position
 * @f: (scope call): function to evaluate at @x
 * @user_data: user data to pass to @f
 *
 * Evaluates the vector-valued function @f at @x and adds the result as a new sample
 * point. The sample is inserted in the correct position to maintain ascending x-order.
 * The interval_ok flag for the new sample is initialized to 0. The sample is marked as
 * a new point.
 *
 */
void
ncm_function_sample_set_add_func (NcmFunctionSampleSet *fss, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data)
{
  NcmVector *y = ncm_vector_new (fss->len);

  f (x, y, user_data);
  ncm_function_sample_set_add (fss, x, y);
  ncm_vector_free (y);
}

/**
 * ncm_function_sample_set_add_old:
 * @fss: a #NcmFunctionSampleSet
 * @x: knot position
 * @y: vector value at @x
 *
 * Adds a new sample point to @fss with position @x and vector value @y, marked as OLD.
 * This function is useful for boundary extensions where the added point should be
 * considered part of the base spline rather than a refinement target. The sample is
 * inserted in the correct position to maintain ascending x-order. The vector @y is
 * copied and must have dimension matching the @fss:len property. The interval_ok flag
 * for the new sample is initialized to 0. The sample is marked as an old point
 * (new_point = FALSE).
 *
 */
void
ncm_function_sample_set_add_old (NcmFunctionSampleSet *fss, const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp;
  guint i;

  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  sp           = _ncm_function_sample_point_new_old (x, y);
  fss->samples = g_list_insert_sorted (fss->samples, sp, _ncm_function_sample_point_cmp);

  /* Update x_min and x_max */
  fss->x_min = MIN (fss->x_min, x);
  fss->x_max = MAX (fss->x_max, x);

  /* Update absmaxF for each component */
  for (i = 0; i < fss->len; i++)
  {
    const gdouble abs_yi = fabs (ncm_vector_get (y, i));
    gdouble current_max  = g_array_index (fss->absmaxF, gdouble, i);

    if (abs_yi > current_max)
    {
      g_array_index (fss->absmaxF, gdouble, i)   = abs_yi;
      g_array_index (fss->absmaxF_x, gdouble, i) = x;
    }
  }
}

/**
 * ncm_function_sample_set_add_old_func:
 * @fss: a #NcmFunctionSampleSet
 * @x: knot position
 * @f: (scope call): function to evaluate at @x
 * @user_data: user data to pass to @f
 *
 * Evaluates the vector-valued function @f at @x and adds the result as an old sample
 * point. This function is useful for boundary extensions where the added point should
 * be considered part of the base spline rather than a refinement target. The sample is
 * inserted in the correct position to maintain ascending x-order. The interval_ok flag
 * for the new sample is initialized to 0. The sample is marked as an old point
 * (new_point = FALSE).
 *
 */
void
ncm_function_sample_set_add_old_func (NcmFunctionSampleSet *fss, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data)
{
  NcmVector *y = ncm_vector_new (fss->len);

  f (x, y, user_data);
  ncm_function_sample_set_add_old (fss, x, y);
  ncm_vector_free (y);
}

/**
 * ncm_function_sample_set_get_len:
 * @fss: a #NcmFunctionSampleSet
 *
 * Gets the dimension of the vector-valued function output.
 *
 * Returns: the dimension $n$
 */
guint
ncm_function_sample_set_get_len (NcmFunctionSampleSet *fss)
{
  return fss->len;
}

/**
 * ncm_function_sample_set_set_track_residual:
 * @fss: a #NcmFunctionSampleSet
 * @track_residual: whether to record the achieved residual
 *
 * Sets #NcmFunctionSampleSet:track-residual. Set it before refining: only the
 * passes that run while it is on leave a record behind.
 *
 */
void
ncm_function_sample_set_set_track_residual (NcmFunctionSampleSet *fss, const gboolean track_residual)
{
  fss->track_residual = track_residual;
}

/**
 * ncm_function_sample_set_get_track_residual:
 * @fss: a #NcmFunctionSampleSet
 *
 * Gets #NcmFunctionSampleSet:track-residual.
 *
 * Returns: whether the achieved residual is being recorded
 */
gboolean
ncm_function_sample_set_get_track_residual (NcmFunctionSampleSet *fss)
{
  return fss->track_residual;
}

/**
 * ncm_function_sample_set_get_nsamples:
 * @fss: a #NcmFunctionSampleSet
 *
 * Gets the number of samples currently stored in @fss.
 *
 * Returns: the number of samples
 */
guint
ncm_function_sample_set_get_nsamples (NcmFunctionSampleSet *fss)
{
  return g_list_length (fss->samples);
}

/**
 * ncm_function_sample_set_get_x_min:
 * @fss: a #NcmFunctionSampleSet
 *
 * Gets the minimum x value in the sample set.
 *
 * Returns: the minimum x value, or GSL_POSINF if no samples exist
 */
gdouble
ncm_function_sample_set_get_x_min (NcmFunctionSampleSet *fss)
{
  return fss->x_min;
}

/**
 * ncm_function_sample_set_get_x_max:
 * @fss: a #NcmFunctionSampleSet
 *
 * Gets the maximum x value in the sample set.
 *
 * Returns: the maximum x value, or GSL_NEGINF if no samples exist
 */
gdouble
ncm_function_sample_set_get_x_max (NcmFunctionSampleSet *fss)
{
  return fss->x_max;
}

/**
 * ncm_function_sample_set_get_absmaxF:
 * @fss: a #NcmFunctionSampleSet
 * @i: component index
 * @x: (out): x value where the maximum absolute value for component @i occurs
 *
 * Gets the maximum absolute value observed for component @i across all samples. This
 * is useful for determining appropriate tolerances in refinement algorithms.
 *
 * Returns: the maximum absolute value for component @i
 */
gdouble
ncm_function_sample_set_get_absmaxF (NcmFunctionSampleSet *fss, const guint i, gdouble *x)
{
  g_assert_cmpuint (i, <, fss->len);

  if (x != NULL)
    *x = g_array_index (fss->absmaxF_x, gdouble, i);

  return g_array_index (fss->absmaxF, gdouble, i);
}

/**
 * ncm_function_sample_set_get_absmaxF_l2_norm:
 * @fss: a #NcmFunctionSampleSet
 *
 * Computes the $L_2$ norm (Euclidean norm) of the maximum absolute values
 * across all components:
 * $$\|\vec{F}\|_2 = \sqrt{\sum_{i=0}^{n-1} (\max_x |F_i(x)|)^2}$$
 *
 * This is useful for determining a reference scale when all components
 * should contribute equally to error estimation.
 *
 * Returns: the $L_2$ norm of the component-wise maximum absolute values
 */
gdouble
ncm_function_sample_set_get_absmaxF_l2_norm (NcmFunctionSampleSet *fss)
{
  gdouble sum_sq = 0.0;
  guint i;

  for (i = 0; i < fss->len; i++)
  {
    const gdouble absmaxF_i = g_array_index (fss->absmaxF, gdouble, i);

    sum_sq += gsl_pow_2 (absmaxF_i);
  }

  return sqrt (sum_sq);
}

/**
 * ncm_function_sample_set_get_absmaxF_linf_norm:
 * @fss: a #NcmFunctionSampleSet
 *
 * Computes the $L_\infty$ norm (maximum norm) of the maximum absolute values
 * across all components:
 * $$\|\vec{F}\|_\infty = \max_{i=0}^{n-1} (\max_x |F_i(x)|)$$
 *
 * This represents the largest absolute value observed across all components
 * and samples, useful for conservative error bounds.
 *
 * Returns: the $L_\infty$ norm of the component-wise maximum absolute values
 */
gdouble
ncm_function_sample_set_get_absmaxF_linf_norm (NcmFunctionSampleSet *fss)
{
  gdouble max_val = 0.0;
  guint i;

  for (i = 0; i < fss->len; i++)
  {
    const gdouble absmaxF_i = g_array_index (fss->absmaxF, gdouble, i);

    max_val = GSL_MAX (max_val, absmaxF_i);
  }

  return max_val;
}

/**
 * ncm_function_sample_set_get_absmaxF_min:
 * @fss: a #NcmFunctionSampleSet
 *
 * Computes the minimum of the maximum absolute values across all components
 * that are not identically zero:
 * $$\min_{i : \max_x |F_i(x)| > 0} (\max_x |F_i(x)|)$$
 *
 * This returns the smallest (nonzero) peak value among all components, which
 * is useful for setting conservative tolerances that ensure even the weakest
 * component is adequately resolved in adaptive refinement algorithms.
 *
 * A component whose peak is exactly zero (e.g. a spin-2 field's ℓ=0,1
 * multipoles, which vanish identically) is excluded from the minimum: it is
 * already exactly represented by any spline through zero-valued samples and
 * needs no tolerance budget, so letting it collapse the tolerance to zero
 * for every *other*, genuinely nonzero component would be wrong -- that
 * tolerance would then be unsatisfiable by any finite refinement (see
 * ncm_function_sample_set_refine()).
 *
 * Returns: the minimum of the nonzero component-wise maximum absolute
 *          values, or 0.0 if every component is identically zero (or there
 *          are no components)
 */
gdouble
ncm_function_sample_set_get_absmaxF_min (NcmFunctionSampleSet *fss)
{
  gdouble min_val = GSL_POSINF;
  guint i;

  for (i = 0; i < fss->len; i++)
  {
    const gdouble absmaxF_i = g_array_index (fss->absmaxF, gdouble, i);

    if (absmaxF_i > 0.0)
      min_val = GSL_MIN (min_val, absmaxF_i);
  }

  return (min_val == GSL_POSINF) ? 0.0 : min_val;
}

/**
 * ncm_function_sample_set_expand_domain:
 * @fss: a #NcmFunctionSampleSet
 * @f: (scope call): callback function to evaluate at new points
 * @x_min_hard: minimum allowed x value (hard limit)
 * @x_max_hard: maximum allowed x value (hard limit)
 * @expansion_factor: multiplicative factor for domain expansion (e.g., 0.2 for 20%)
 * @epsilon: convergence threshold relative to reference scale
 * @max_iterations: maximum number of expansion iterations
 * @consecutive_tries: number of consecutive converged points required on each side
 * @user_data: user data passed to callback function
 *
 * Expands the domain of the function sample set by alternating between left and
 * right boundary expansion until convergence or limits are reached. At each step:
 *
 * - Left expansion proposes: $x_{\text{new}} = x_{\text{min}} \times (1 - \alpha)$
 * - Right expansion proposes: $x_{\text{new}} = x_{\text{max}} \times (1 + \alpha)$
 * - where $\alpha$ is the @expansion_factor
 *
 * The function recomputes the reference scale $F_{\text{ref}} = \|\text{absmaxF}\|_2$
 * after each point insertion, ensuring proper convergence detection even when the
 * function is still growing. Expansion on each side stops when:
 *
 * - $\|F(x)\| < \epsilon \times F_{\text{ref}}$ for @consecutive_tries consecutive points, OR
 * - A hard limit is reached
 *
 * The interleaved expansion pattern (left, right, left, right, ...) ensures the
 * reference scale tracks the global function behavior correctly.
 */
void
ncm_function_sample_set_expand_domain (NcmFunctionSampleSet     *fss,
                                       NcmFunctionSampleSetFunc f,
                                       const gdouble            x_min_hard,
                                       const gdouble            x_max_hard,
                                       const gdouble            expansion_factor,
                                       const gdouble            epsilon,
                                       const guint              max_iterations,
                                       const guint              consecutive_tries,
                                       gpointer                 user_data)
{
  gboolean expanding_left  = TRUE;
  gboolean expanding_right = TRUE;
  guint left_converged     = 0;
  guint right_converged    = 0;
  guint i;

  for (i = 0; i < max_iterations; i++)
  {
    /* Try left expansion */
    if (expanding_left)
    {
      const gdouble x_min         = ncm_function_sample_set_get_x_min (fss);
      const gdouble new_x         = x_min * (1.0 - expansion_factor);
      const gdouble new_x_clamped = GSL_MAX (new_x, x_min_hard);
      NcmVector *y                = ncm_vector_new (fss->len);
      gdouble norm_y;

      /* Evaluate function at new point (clamped to hard limit) */
      f (new_x_clamped, y, user_data);

      /* Compute norm of the function value */
      norm_y = ncm_vector_dnrm2 (y);

      /* Add point to sample set */
      ncm_function_sample_set_add_old (fss, new_x_clamped, y);

      /* Check if we hit the hard limit */
      if (new_x < x_min_hard)
      {
        expanding_left = FALSE;
      }
      else
      {
        /* Get updated reference scale */
        const gdouble F_ref = ncm_function_sample_set_get_absmaxF_l2_norm (fss);

        /* Check convergence */
        if (norm_y < epsilon * F_ref)
        {
          left_converged++;

          if (left_converged >= consecutive_tries)
            expanding_left = FALSE;
        }
        else
        {
          left_converged = 0;
        }
      }

      ncm_vector_free (y);
    }

    /* Try right expansion */
    if (expanding_right)
    {
      const gdouble x_max         = ncm_function_sample_set_get_x_max (fss);
      const gdouble new_x         = x_max * (1.0 + expansion_factor);
      const gdouble new_x_clamped = GSL_MIN (new_x, x_max_hard);
      NcmVector *y                = ncm_vector_new (fss->len);
      gdouble norm_y;

      /* Evaluate function at new point (clamped to hard limit) */
      f (new_x_clamped, y, user_data);

      /* Compute norm of the function value */
      norm_y = ncm_vector_dnrm2 (y);

      /* Add point to sample set */
      ncm_function_sample_set_add_old (fss, new_x_clamped, y);

      /* Check if we hit the hard limit */
      if (new_x > x_max_hard)
      {
        expanding_right = FALSE;
      }
      else
      {
        /* Get updated reference scale */
        const gdouble F_ref = ncm_function_sample_set_get_absmaxF_l2_norm (fss);

        /* Check convergence */
        if (norm_y < epsilon * F_ref)
        {
          right_converged++;

          if (right_converged >= consecutive_tries)
            expanding_right = FALSE;
        }
        else
        {
          right_converged = 0;
        }
      }

      ncm_vector_free (y);
    }

    /* Stop if both sides converged or hit limits */
    if (!expanding_left && !expanding_right)
      break;
  }
}

/**
 * ncm_function_sample_set_reset_interval_ok:
 * @fss: a #NcmFunctionSampleSet
 *
 * Resets all interval_ok flags to 0. This is useful when starting a new refinement pass.
 *
 */
void
ncm_function_sample_set_reset_interval_ok (NcmFunctionSampleSet *fss)
{
  GList *node;

  for (node = fss->samples; node != NULL; node = node->next)
  {
    ((NcmFunctionSamplePoint *) node->data)->interval_ok = 0;
  }
}

/**
 * ncm_function_sample_set_all_intervals_ok:
 * @fss: a #NcmFunctionSampleSet
 * @threshold: minimum interval_ok value required
 *
 * Checks if all intervals have interval_ok >= @threshold. This is useful for
 * determining convergence in refinement algorithms - when all intervals have passed
 * the refinement test enough times.
 *
 * Note: The last sample point is excluded since it doesn't define an interval.
 *
 * Returns: TRUE if all intervals have interval_ok >= @threshold, FALSE otherwise
 */
gboolean
ncm_function_sample_set_all_intervals_ok (NcmFunctionSampleSet *fss, const gint threshold)
{
  const guint nsamples = ncm_function_sample_set_get_nsamples (fss);
  GList *node;
  guint i;

  /* Need at least 2 points to have intervals */
  if (nsamples < 2)
    return FALSE;

  /* Check all intervals (all points except the last) */
  for (node = fss->samples, i = 0; node != NULL && i < nsamples - 1; node = node->next, i++)
  {
    NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) node->data;

    if (sp->interval_ok < threshold)
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_function_sample_set_mark_all_old:
 * @fss: a #NcmFunctionSampleSet
 *
 * Marks all points in @fss as old. This is typically called after a refinement pass to
 * indicate that all current points should be used in the next spline construction.
 *
 */
void
ncm_function_sample_set_mark_all_old (NcmFunctionSampleSet *fss)
{
  GList *node;

  for (node = fss->samples; node != NULL; node = node->next)
  {
    ((NcmFunctionSamplePoint *) node->data)->new_point = FALSE;
  }
}

/**
 * ncm_function_sample_set_to_spline_vec:
 * @fss: a #NcmFunctionSampleSet
 * @base_spline: a #NcmSpline to use as the base spline type
 *
 * Converts the sample set to a #NcmSplineVec. This reuses cached internal arrays for
 * efficiency, which means that:
 *
 * - **The returned #NcmSplineVec is invalidated by subsequent calls** to this function
 *   or ncm_function_sample_set_to_spline_vec_old() on the same @fss object.
 * - If you need to keep multiple #NcmSplineVec objects from the same sample set, you
 *   must call ncm_spline_vec_dup() on the returned object before calling this function
 *   again.
 *
 * The sample set itself is not modified and can continue to be used for further
 * refinement.
 *
 * Returns: (transfer full): a new #NcmSplineVec containing the interpolated function
 */
NcmSplineVec *
ncm_function_sample_set_to_spline_vec (NcmFunctionSampleSet *fss, NcmSpline *base_spline)
{
  const guint nsamples = ncm_function_sample_set_get_nsamples (fss);
  NcmVector *xv;
  GPtrArray *yv = g_ptr_array_new_full (fss->len, (GDestroyNotify) ncm_vector_free);
  GList *node;
  guint i, j;

  /* Resize cached arrays */
  g_array_set_size (fss->x_array_cache, nsamples);

  for (i = 0; i < fss->len; i++)
  {
    GArray *y_array = g_ptr_array_index (fss->y_arrays_cache, i);

    g_array_set_size (y_array, nsamples);
  }

  /* Fill cached arrays from samples */
  for (node = fss->samples, i = 0; node != NULL; node = node->next, i++)
  {
    NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) node->data;

    g_array_index (fss->x_array_cache, gdouble, i) = sp->x;

    for (j = 0; j < fss->len; j++)
    {
      GArray *y_array = g_ptr_array_index (fss->y_arrays_cache, j);

      g_array_index (y_array, gdouble, i) = ncm_vector_get (sp->y, j);
    }
  }

  /* Create NcmVectors from cached arrays */
  xv = ncm_vector_new_array (fss->x_array_cache);

  for (i = 0; i < fss->len; i++)
  {
    GArray *y_array   = g_ptr_array_index (fss->y_arrays_cache, i);
    NcmVector *y_comp = ncm_vector_new_array (y_array);

    g_ptr_array_add (yv, y_comp);
  }

  {
    NcmSplineVec *sv = ncm_spline_vec_new_gpa (base_spline, xv, yv, TRUE);

    ncm_vector_free (xv);
    g_ptr_array_unref (yv);

    return sv;
  }
}

/**
 * ncm_function_sample_set_to_spline_vec_old:
 * @fss: a #NcmFunctionSampleSet
 * @base_spline: a #NcmSpline to use as the base spline type
 *
 * Converts only the OLD sample points to a #NcmSplineVec. This reuses cached internal
 * arrays for efficiency, which means that:
 *
 * - **The returned #NcmSplineVec is invalidated by subsequent calls** to this function
 *   or ncm_function_sample_set_to_spline_vec() on the same @fss object.
 * - If you need to keep multiple #NcmSplineVec objects from the same sample set, you
 *   must call ncm_spline_vec_dup() on the returned object before calling this function
 *   again.
 *
 * This creates arrays from the samples where new_point is FALSE. The sample set is not
 * modified. This is useful for building a spline to test against NEW points during
 * refinement.
 *
 * Returns: (transfer full): a new #NcmSplineVec containing the interpolated function
 */
NcmSplineVec *
ncm_function_sample_set_to_spline_vec_old (NcmFunctionSampleSet *fss, NcmSpline *base_spline)
{
  GList *node;
  guint nold = 0;
  NcmVector *xv;
  GPtrArray *yv;
  guint i, j;

  /* Count old points */
  for (node = fss->samples; node != NULL; node = node->next)
  {
    NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) node->data;

    if (!sp->new_point)
      nold++;
  }

  g_assert (nold >= 2); /* Need at least 2 points for a spline */

  /* Resize cached arrays */
  g_array_set_size (fss->x_array_cache, nold);

  for (i = 0; i < fss->len; i++)
  {
    GArray *y_array = g_ptr_array_index (fss->y_arrays_cache, i);

    g_array_set_size (y_array, nold);
  }

  /* Fill cached arrays from old samples only */
  for (node = fss->samples, i = 0; node != NULL; node = node->next)
  {
    NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) node->data;

    if (!sp->new_point)
    {
      g_array_index (fss->x_array_cache, gdouble, i) = sp->x;

      for (j = 0; j < fss->len; j++)
      {
        GArray *y_array = g_ptr_array_index (fss->y_arrays_cache, j);

        g_array_index (y_array, gdouble, i) = ncm_vector_get (sp->y, j);
      }

      i++;
    }
  }

  /* Create NcmVectors from cached arrays */
  xv = ncm_vector_new_array (fss->x_array_cache);
  yv = g_ptr_array_new_full (fss->len, (GDestroyNotify) ncm_vector_free);

  for (i = 0; i < fss->len; i++)
  {
    GArray *y_array   = g_ptr_array_index (fss->y_arrays_cache, i);
    NcmVector *y_comp = ncm_vector_new_array (y_array);

    g_ptr_array_add (yv, y_comp);
  }

  {
    NcmSplineVec *sv = ncm_spline_vec_new_gpa (base_spline, xv, yv, TRUE);

    ncm_vector_free (xv);
    g_ptr_array_unref (yv);

    return sv;
  }
}

/*
 * Records @diff as the residual of the interval owned by @sp, keeping the
 * larger of the two whenever an interval is measured more than once. What is
 * stored is the componentwise |f - spline| of the pass that accepted the
 * interval, so it describes the interval as it stood *before* its own
 * subdivision -- see the note in ncm_function_sample_set_get_residuals().
 */
static void
_ncm_function_sample_point_update_residual (NcmFunctionSamplePoint *sp, NcmVector *diff, const guint len)
{
  guint j;

  if (sp->residual == NULL)
  {
    sp->residual = ncm_vector_new (len);
    ncm_vector_set_all (sp->residual, 0.0);
  }

  for (j = 0; j < len; j++)
  {
    const gdouble d = fabs (ncm_vector_get (diff, j));

    if (d > ncm_vector_get (sp->residual, j))
      ncm_vector_set (sp->residual, j, d);
  }
}

/**
 * ncm_function_sample_set_get_residuals:
 * @fss: a #NcmFunctionSampleSet
 *
 * Builds the per-interval residuals recorded while
 * #NcmFunctionSampleSet:track-residual was on, as a matrix with one row per
 * sample, in ascending $x$ order, and one column per component. Row $i$ holds
 * the residual of the interval $[x_i, x_{i+1}]$, matching the convention
 * ncm_function_sample_set_iter_get_interval_ok() uses; the last row is the
 * trailing edge and is always NaN.
 *
 * An entry is the componentwise $|f - \tilde f|$ measured by
 * ncm_function_sample_set_refine() on the pass that accepted the interval, not
 * the tolerance that pass was asked for -- refinement usually beats its own
 * threshold by orders, and by an amount that varies with the function. An
 * interval that was never accepted (refinement stopped at
 * @max_iter, or tracking was off for the pass that accepted it) reads NaN, and
 * the caller has to decide what to do with it.
 *
 * Two limits worth stating. The measurement is a *midpoint sample*, not a
 * supremum over the interval. And it is taken before the interval is split, so
 * for a cubic base spline it overstates each half by roughly $2^4$.
 *
 * Returns: (transfer full) (nullable): a #NcmMatrix of residuals, or %NULL if
 * nothing was recorded
 */
NcmMatrix *
ncm_function_sample_set_get_residuals (NcmFunctionSampleSet *fss)
{
  const guint nsamples = g_list_length (fss->samples);
  NcmMatrix *residuals;
  GList *node;
  guint i;

  if ((nsamples == 0) || !fss->track_residual)
    return NULL;

  residuals = ncm_matrix_new (nsamples, fss->len);

  for (node = fss->samples, i = 0; node != NULL; node = node->next, i++)
  {
    const NcmFunctionSamplePoint *sp = (const NcmFunctionSamplePoint *) node->data;
    guint j;

    for (j = 0; j < fss->len; j++)
      ncm_matrix_set (residuals, i, j, (sp->residual != NULL) ? ncm_vector_get (sp->residual, j) : GSL_NAN);
  }

  /* The last sample owns no interval to its right. */
  for (i = 0; i < fss->len; i++)
    ncm_matrix_set (residuals, nsamples - 1, i, GSL_NAN);

  return residuals;
}

/**
 * ncm_function_sample_set_estimate_residuals:
 * @fss: a #NcmFunctionSampleSet
 * @base_spline: the #NcmSpline whose error is wanted
 * @ref_spline: a higher-order #NcmSpline to measure it against
 *
 * Estimates @base_spline's interpolation error from the samples already held,
 * **without evaluating the function anywhere**: both splines are fitted to the
 * same data and differenced at each interval's midpoint, which is where the
 * samples say least. An embedded pair, in the sense a Runge-Kutta pair is one.
 *
 * The layout matches ncm_function_sample_set_get_residuals(): one row per
 * sample in ascending $x$ order, one column per component, row $i$ owning the
 * interval $[x_i, x_{i+1}]$, and a trailing row of NaN.
 *
 * The two differ in what they measure, and it is worth being clear which is
 * wanted. get_residuals() reports what refinement *observed*, against the true
 * function, but before the interval was split and at a point that then became a
 * knot -- so it describes the parent interval and overstates the final grid.
 * This reports the error of the finished fit on the final intervals, and costs
 * one extra banded solve rather than one function evaluation per interval.
 *
 * It inherits the usual embedded-pair assumption: the estimate is only as good
 * as @ref_spline being the better fit. Where that fails the difference measures
 * @ref_spline's error instead, which overstates rather than hides.
 *
 * Returns: (transfer full) (nullable): a #NcmMatrix of estimates, or %NULL if
 * there are too few samples
 */
NcmMatrix *
ncm_function_sample_set_estimate_residuals (NcmFunctionSampleSet *fss, NcmSpline *base_spline, NcmSpline *ref_spline)
{
  const guint nsamples = g_list_length (fss->samples);
  NcmSplineVec *sv_base, *sv_ref;
  NcmVector *y_base, *y_ref;
  NcmMatrix *residuals;
  GList *node;
  guint i;

  if (nsamples < 2)
    return NULL;

  sv_base = ncm_function_sample_set_to_spline_vec (fss, base_spline);
  sv_ref  = ncm_function_sample_set_to_spline_vec (fss, ref_spline);

  residuals = ncm_matrix_new (nsamples, fss->len);
  y_base    = ncm_vector_new (fss->len);
  y_ref     = ncm_vector_new (fss->len);

  for (node = fss->samples, i = 0; node->next != NULL; node = node->next, i++)
  {
    const NcmFunctionSamplePoint *sp   = (const NcmFunctionSamplePoint *) node->data;
    const NcmFunctionSamplePoint *next = (const NcmFunctionSamplePoint *) node->next->data;
    const gdouble xm                   = 0.5 * (sp->x + next->x);
    guint j;

    ncm_spline_vec_eval (sv_base, xm, y_base);
    ncm_spline_vec_eval (sv_ref, xm, y_ref);

    for (j = 0; j < fss->len; j++)
      ncm_matrix_set (residuals, i, j, fabs (ncm_vector_get (y_base, j) - ncm_vector_get (y_ref, j)));
  }

  /* The last sample owns no interval to its right. */
  for (i = 0; i < fss->len; i++)
    ncm_matrix_set (residuals, nsamples - 1, i, GSL_NAN);

  ncm_vector_free (y_base);
  ncm_vector_free (y_ref);
  ncm_spline_vec_free (sv_base);
  ncm_spline_vec_free (sv_ref);

  return residuals;
}

/**
 * ncm_function_sample_set_refine:
 * @fss: a #NcmFunctionSampleSet
 * @reltol: relative tolerance for refinement test
 * @abstol: absolute tolerance for refinement test
 * @base_spline: a #NcmSpline to use as the base spline type
 *
 * Performs a refinement pass on all NEW points. For each NEW point, this function:
 * 1. Creates a spline using OLD points only
 * 2. Evaluates the spline at the NEW point position
 * 3. Computes the error: ||f(x) - spline_f(x)||_2 <= reltol * ||f(x)||_2 + abstol
 * 4. If the test passes, increments interval_ok for both the NEW point and its left neighbor
 * 5. Marks all NEW points as OLD
 *
 * When #NcmFunctionSampleSet:track-residual is set, step 4 also records the
 * residual it just measured on both intervals, retrievable afterwards with
 * ncm_function_sample_set_get_residuals().
 *
 * The interval_ok counter at node i indicates how many times the interval [i, i+1] has
 * passed the refinement test. After refinement, all points are marked as OLD for the
 * next iteration.
 *
 */
void
ncm_function_sample_set_refine (NcmFunctionSampleSet *fss, const gdouble reltol, const gdouble abstol, NcmSpline *base_spline)
{
  NcmSplineVec *sv_old;
  NcmVector *y_spline;
  NcmFunctionSampleSetIter iter_s;
  NcmFunctionSampleSetIter *iter = &iter_s;

  /* Create spline from OLD points only */
  sv_old = ncm_function_sample_set_to_spline_vec_old (fss, base_spline);

  y_spline = ncm_vector_new (fss->len);

  /* Test each NEW point using iterator */
  ncm_function_sample_set_iter_begin (fss, &iter);

  while (ncm_function_sample_set_iter_is_valid (iter))
  {
    if (ncm_function_sample_set_iter_get_new_point (iter))
    {
      gdouble x, norm_f, norm_diff, threshold;
      NcmVector *y;

      x = ncm_function_sample_set_iter_get_x (iter);
      y = ncm_function_sample_set_iter_get_y (iter);

      /* Evaluate spline at new point */
      ncm_spline_vec_eval (sv_old, x, y_spline);

      /* Compute norms */
      norm_f = ncm_vector_dnrm2 (y); /* ||f(x)||_2 */

      /* Compute difference vector in place */
      {
        guint j;

        for (j = 0; j < fss->len; j++)
        {
          const gdouble diff = ncm_vector_get (y, j) - ncm_vector_get (y_spline, j);

          ncm_vector_set (y_spline, j, diff);
        }
      }

      norm_diff = ncm_vector_dnrm2 (y_spline); /* ||f(x) - spline_f(x)||_2 */
      threshold = reltol * norm_f + abstol;

      /* If test passes, increment interval_ok for this point and its left neighbor.
       * Uses <= (not <) so an exactly converged point (norm_diff == threshold == 0,
       * e.g. every sampled component is identically zero) still passes instead of
       * failing forever on a strict comparison it can never satisfy. */
      if (norm_diff <= threshold)
      {
        ncm_function_sample_set_iter_inc_interval_ok (iter);

        if (fss->track_residual)
          _ncm_function_sample_point_update_residual ((NcmFunctionSamplePoint *) iter->node->data, y_spline, fss->len);

        /* Also increment left neighbor if it exists */
        if (ncm_function_sample_set_iter_has_prev (iter))
        {
          NcmFunctionSampleSetIter prev = *iter;

          ncm_function_sample_set_iter_prev (&prev);
          ncm_function_sample_set_iter_inc_interval_ok (&prev);

          /* Both halves of the interval this point just split are covered by
           * the one measurement taken at the point itself. */
          if (fss->track_residual)
            _ncm_function_sample_point_update_residual ((NcmFunctionSamplePoint *) prev.node->data, y_spline, fss->len);
        }
      }
    }

    ncm_function_sample_set_iter_next (iter);
  }

  /* iter is stack-allocated; no free needed */

  /* Mark all points as OLD for next iteration */
  ncm_function_sample_set_mark_all_old (fss);

  ncm_vector_free (y_spline);
  ncm_spline_vec_free (sv_old);
}

/**
 * ncm_function_sample_set_adaptive_midpoint:
 * @fss: a #NcmFunctionSampleSet
 * @f: (scope call): function used to evaluate new midpoints
 * @user_data: user data passed to @f
 * @reltol: relative tolerance for refinement test
 * @abstol: absolute tolerance for refinement test
 * @max_iter: maximum number of refinement iterations
 * @min_pass_threshold: minimum interval_ok threshold to consider interval passed
 * @base_spline: base spline used for refinement tests
 *
 * Performs an iterative midpoint-based adaptive refinement similar to the
 * Python test harness. On each iteration intervals with interval_ok <
 * @min_pass_threshold receive their midpoint inserted (evaluated via @f). After
 * insertion, a refinement pass is performed via ncm_function_sample_set_refine(). The
 * process stops when all intervals reach @min_pass_threshold or when @max_iter is
 * reached.
 */
void
ncm_function_sample_set_adaptive_midpoint (NcmFunctionSampleSet     *fss,
                                           NcmFunctionSampleSetFunc f,
                                           const gdouble            reltol,
                                           const gdouble            abstol,
                                           const guint              max_iter,
                                           const gint               min_pass_threshold,
                                           NcmSpline                *base_spline,
                                           gpointer                 user_data)
{
  guint iteration;

  /* The Cubic not-a-knot spline requires at least 6 samples.
   * If we have fewer but at least 2, add midpoints until we reach 6.
   */
  {
    const guint nsamples = ncm_function_sample_set_get_nsamples (fss);

    if (nsamples < 2)
      g_error ("Adaptive ncm_function_sample_set_adaptive_midpoint: "
               "need at least 2 samples, got %u", nsamples);

    if (nsamples < 6)
    {
      /* Add midpoints to reach at least 6 samples */
      while (ncm_function_sample_set_get_nsamples (fss) < 6)
      {
        NcmFunctionSampleSetIter it_s;
        NcmFunctionSampleSetIter *it = &it_s;

        ncm_function_sample_set_iter_begin (fss, &it);

        while (ncm_function_sample_set_iter_has_next (it))
        {
          NcmFunctionSampleSetIter iter_new_s;
          NcmFunctionSampleSetIter iter_right = *it;
          NcmFunctionSampleSetIter *iter_new  = &iter_new_s;
          gdouble x_left, x_right, x_mid;

          ncm_function_sample_set_iter_next (&iter_right);

          x_left  = ncm_function_sample_set_iter_get_x (it);
          x_right = ncm_function_sample_set_iter_get_x (&iter_right);
          x_mid   = 0.5 * (x_left + x_right);

          /* Add midpoint as an old point (part of base spline) */
          ncm_function_sample_set_iter_insert_after_func (fss, it, x_mid, f, user_data, &iter_new);

          /* Mark as old immediately */
          ncm_function_sample_set_iter_set_new_point (iter_new, FALSE);

          /* Skip the newly inserted point to move to the next interval */
          *it = *iter_new;
          ncm_function_sample_set_iter_next (it);
        }
      }
    }
  }

  for (iteration = 0; iteration < max_iter; iteration++)
  {
    if (ncm_function_sample_set_all_intervals_ok (fss, min_pass_threshold))
      break;

    {
      NcmFunctionSampleSetIter it_s;
      NcmFunctionSampleSetIter *it = &it_s;

      ncm_function_sample_set_iter_begin (fss, &it);

      while (ncm_function_sample_set_iter_has_next (it))
      {
        if (ncm_function_sample_set_iter_get_interval_ok (it) < min_pass_threshold)
        {
          NcmFunctionSampleSetIter iter_new_s;
          NcmFunctionSampleSetIter iter_right = *it;
          NcmFunctionSampleSetIter *iter_new  = &iter_new_s;

          ncm_function_sample_set_iter_next (&iter_right);

          {
            gdouble x_left  = ncm_function_sample_set_iter_get_x (it);
            gdouble x_right = ncm_function_sample_set_iter_get_x (&iter_right);
            gdouble x_mid   = 0.5 * (x_left + x_right);

            if (fabs (x_left - x_right) < GSL_MACH_EPS * (fabs (x_left) + fabs (x_right)))
            {
              ncm_function_sample_set_iter_set_interval_ok (it, min_pass_threshold); /* Mark as passed to avoid infinite loop */
              continue;
            }

            ncm_function_sample_set_iter_insert_after_func (fss, it, x_mid, f, user_data, &iter_new);
            *it = *iter_new; /* Move iterator to new point for next loop */
          }
        }

        ncm_function_sample_set_iter_next (it);
      }

      ncm_function_sample_set_refine (fss, reltol, abstol, base_spline);
    }
  }

  if (iteration == max_iter)
    g_message ("# ncm_function_sample_set_adaptive_midpoint: Max iterations (%u) reached with %u knots\n",
               max_iter, ncm_function_sample_set_get_nsamples (fss));
}

/**
 * ncm_function_sample_set_log_vals:
 * @fss: a #NcmFunctionSampleSet
 *
 * Logs all sample values in @fss for debugging purposes. This prints the x position,
 * vector components, interval_ok flag, and new_point flag for each sample.
 *
 */
void
ncm_function_sample_set_log_vals (NcmFunctionSampleSet *fss)
{
  GList *node;
  guint i;

  g_message ("# NcmFunctionSampleSet: %u samples, dimension %u\n",
             ncm_function_sample_set_get_nsamples (fss),
             fss->len);

  for (node = fss->samples, i = 0; node != NULL; node = node->next, i++)
  {
    NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) node->data;
    GString *str               = g_string_new ("");
    guint j;

    g_string_append_printf (str, "# [%4u] x = % 22.15e  y = [", i, sp->x);

    for (j = 0; j < fss->len; j++)
    {
      g_string_append_printf (str, "% 22.15e", ncm_vector_get (sp->y, j));

      if (j < fss->len - 1)
        g_string_append (str, ", ");
    }

    g_string_append_printf (str, "]  interval_ok = %d  new = %d", sp->interval_ok, sp->new_point);

    g_message ("%s\n", str->str);
    g_string_free (str, TRUE);
  }
}

