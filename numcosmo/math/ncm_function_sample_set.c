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
 * This object stores an ordered set of samples $(x_i, \vec{y}_i)$ where each sample
 * consists of a knot position $x_i$ and a vector value $\vec{y}_i \in \mathbb{R}^n$.
 * Each sample also has an associated "OK" flag that can be used to mark whether
 * the sample point is acceptable for interpolation purposes.
 *
 * The primary use case is for iterative refinement algorithms that build splines
 * from vector-valued functions. The typical workflow is:
 * 1. Add initial samples to the set
 * 2. Convert to NcmSplineVec and test interpolation error
 * 3. Insert new samples where error exceeds tolerance
 * 4. Mark samples as OK when bins pass error tests
 * 5. Repeat until convergence
 *
 * Samples are maintained in ascending x-order using a GList internally, which
 * provides efficient insertion operations during the building phase. When
 * converting to NcmSplineVec, the samples are extracted to arrays.
 *
 * The dimension $n$ of the vector values is fixed at creation time and validated
 * on every insertion.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_function_sample_set.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_LEN,
};

typedef struct _NcmFunctionSamplePoint
{
  gdouble x;
  NcmVector *y;
  gint ok;
} NcmFunctionSamplePoint;

struct _NcmFunctionSampleSet
{
  /*< private >*/
  GObject parent_instance;
  guint len;
  GList *samples;
};

G_DEFINE_TYPE (NcmFunctionSampleSet, ncm_function_sample_set, G_TYPE_OBJECT)

static void
_ncm_function_sample_point_free (gpointer data)
{
  NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) data;

  ncm_vector_clear (&sp->y);
  g_slice_free (NcmFunctionSamplePoint, sp);
}

static NcmFunctionSamplePoint *
_ncm_function_sample_point_new (const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp = g_slice_new (NcmFunctionSamplePoint);

  sp->x  = x;
  sp->y  = ncm_vector_ref (y);
  sp->ok = 0;

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
  fss->len     = 0;
  fss->samples = NULL;
}

static void
_ncm_function_sample_set_dispose (GObject *object)
{
  NcmFunctionSampleSet *fss = NCM_FUNCTION_SAMPLE_SET (object);

  g_list_free_full (fss->samples, _ncm_function_sample_point_free);
  fss->samples = NULL;

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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_function_sample_set_class_init (NcmFunctionSampleSetClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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

/**
 * ncm_function_sample_set_add:
 * @fss: a #NcmFunctionSampleSet
 * @x: knot position
 * @y: vector value at @x
 *
 * Adds a new sample point to @fss with position @x and vector value @y.
 * The sample is inserted in the correct position to maintain ascending x-order.
 * The vector @y is copied and must have dimension matching the @fss:len property.
 * The OK flag for the new sample is initialized to 0.
 *
 */
void
ncm_function_sample_set_add (NcmFunctionSampleSet *fss, const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp;

  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  sp           = _ncm_function_sample_point_new (x, y);
  fss->samples = g_list_insert_sorted (fss->samples, sp, _ncm_function_sample_point_cmp);
}

/**
 * ncm_function_sample_set_insert_before:
 * @fss: a #NcmFunctionSampleSet
 * @index: index of the sample before which to insert
 * @x: knot position
 * @y: vector value at @x
 *
 * Inserts a new sample point before the sample at position @index.
 * The vector @y is copied and must have dimension matching the @fss:len property.
 * The OK flag for the new sample is initialized to 0.
 *
 * Warning: This does not check if the insertion maintains x-order. Use with care.
 *
 */
void
ncm_function_sample_set_insert_before (NcmFunctionSampleSet *fss, const guint index, const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp;
  GList *node;

  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  node = g_list_nth (fss->samples, index);
  g_assert (node != NULL);

  sp           = _ncm_function_sample_point_new (x, y);
  fss->samples = g_list_insert_before (fss->samples, node, sp);
}

/**
 * ncm_function_sample_set_insert_after:
 * @fss: a #NcmFunctionSampleSet
 * @index: index of the sample after which to insert
 * @x: knot position
 * @y: vector value at @x
 *
 * Inserts a new sample point after the sample at position @index.
 * The vector @y is copied and must have dimension matching the @fss:len property.
 * The OK flag for the new sample is initialized to 0.
 *
 * Warning: This does not check if the insertion maintains x-order. Use with care.
 *
 */
void
ncm_function_sample_set_insert_after (NcmFunctionSampleSet *fss, const guint index, const gdouble x, NcmVector *y)
{
  NcmFunctionSamplePoint *sp;
  GList *node;

  g_assert_cmpuint (ncm_vector_len (y), ==, fss->len);

  node = g_list_nth (fss->samples, index);
  g_assert (node != NULL);

  sp           = _ncm_function_sample_point_new (x, y);
  fss->samples = g_list_insert_before (fss->samples, node->next, sp);
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
 * ncm_function_sample_set_peek_x:
 * @fss: a #NcmFunctionSampleSet
 * @index: sample index
 *
 * Gets the knot position $x$ of the sample at @index.
 *
 * Returns: the x value
 */
gdouble
ncm_function_sample_set_peek_x (NcmFunctionSampleSet *fss, const guint index)
{
  GList *node = g_list_nth (fss->samples, index);

  g_assert (node != NULL);

  return ((NcmFunctionSamplePoint *) node->data)->x;
}

/**
 * ncm_function_sample_set_peek_y:
 * @fss: a #NcmFunctionSampleSet
 * @index: sample index
 *
 * Gets a borrowed reference to the vector value $\vec{y}$ of the sample at @index.
 * The returned vector should not be modified or freed.
 *
 * Returns: (transfer none): the vector value
 */
NcmVector *
ncm_function_sample_set_peek_y (NcmFunctionSampleSet *fss, const guint index)
{
  GList *node = g_list_nth (fss->samples, index);

  g_assert (node != NULL);

  return ((NcmFunctionSamplePoint *) node->data)->y;
}

/**
 * ncm_function_sample_set_get_sample:
 * @fss: a #NcmFunctionSampleSet
 * @index: sample index
 * @x: (out) (nullable): location to store the x value
 * @y: (out) (transfer none) (nullable): location to store the y vector
 * @ok: (out) (nullable): location to store the OK flag
 *
 * Gets the complete sample data at @index. Any of the output parameters
 * can be NULL if that value is not needed.
 *
 */
void
ncm_function_sample_set_get_sample (NcmFunctionSampleSet *fss, const guint index, gdouble *x, NcmVector **y, gint *ok)
{
  GList *node = g_list_nth (fss->samples, index);
  NcmFunctionSamplePoint *sp;

  g_assert (node != NULL);

  sp = (NcmFunctionSamplePoint *) node->data;

  if (x != NULL)
    *x = sp->x;

  if (y != NULL)
    *y = sp->y;

  if (ok != NULL)
    *ok = sp->ok;
}

/**
 * ncm_function_sample_set_get_ok:
 * @fss: a #NcmFunctionSampleSet
 * @index: sample index
 *
 * Gets the OK flag value of the sample at @index.
 *
 * Returns: the OK flag value
 */
gint
ncm_function_sample_set_get_ok (NcmFunctionSampleSet *fss, const guint index)
{
  GList *node = g_list_nth (fss->samples, index);

  g_assert (node != NULL);

  return ((NcmFunctionSamplePoint *) node->data)->ok;
}

/**
 * ncm_function_sample_set_set_ok:
 * @fss: a #NcmFunctionSampleSet
 * @index: sample index
 * @ok: new OK flag value
 *
 * Sets the OK flag of the sample at @index to @ok.
 *
 */
void
ncm_function_sample_set_set_ok (NcmFunctionSampleSet *fss, const guint index, const gint ok)
{
  GList *node = g_list_nth (fss->samples, index);

  g_assert (node != NULL);

  ((NcmFunctionSamplePoint *) node->data)->ok = ok;
}

/**
 * ncm_function_sample_set_inc_ok:
 * @fss: a #NcmFunctionSampleSet
 * @index: sample index
 *
 * Increments the OK flag of the sample at @index by 1.
 *
 */
void
ncm_function_sample_set_inc_ok (NcmFunctionSampleSet *fss, const guint index)
{
  GList *node = g_list_nth (fss->samples, index);

  g_assert (node != NULL);

  ((NcmFunctionSamplePoint *) node->data)->ok++;
}

/**
 * ncm_function_sample_set_reset_ok:
 * @fss: a #NcmFunctionSampleSet
 *
 * Resets all OK flags to 0. This is useful when starting a new refinement pass.
 *
 */
void
ncm_function_sample_set_reset_ok (NcmFunctionSampleSet *fss)
{
  GList *node;

  for (node = fss->samples; node != NULL; node = node->next)
  {
    ((NcmFunctionSamplePoint *) node->data)->ok = 0;
  }
}

/**
 * ncm_function_sample_set_to_spline_vec:
 * @fss: a #NcmFunctionSampleSet
 * @base_spline: a #NcmSpline to use as the base spline type
 *
 * Converts the sample set to a #NcmSplineVec. This creates arrays from the
 * current samples and constructs a new NcmSplineVec object. The sample set
 * is not modified and can continue to be used for further refinement.
 *
 * Returns: (transfer full): a new #NcmSplineVec containing the interpolated function
 */
NcmSplineVec *
ncm_function_sample_set_to_spline_vec (NcmFunctionSampleSet *fss, NcmSpline *base_spline)
{
  const guint nsamples = ncm_function_sample_set_get_nsamples (fss);
  NcmVector *xv        = ncm_vector_new (nsamples);
  GPtrArray *yv        = g_ptr_array_new_full (fss->len, (GDestroyNotify) ncm_vector_free);
  GList *node;
  guint i, j;

  /* Allocate y vectors */
  for (i = 0; i < fss->len; i++)
  {
    NcmVector *y_comp = ncm_vector_new (nsamples);

    g_ptr_array_add (yv, y_comp);
  }

  /* Fill arrays from samples */
  for (node = fss->samples, i = 0; node != NULL; node = node->next, i++)
  {
    NcmFunctionSamplePoint *sp = (NcmFunctionSamplePoint *) node->data;

    ncm_vector_set (xv, i, sp->x);

    for (j = 0; j < fss->len; j++)
    {
      NcmVector *y_comp = g_ptr_array_index (yv, j);

      ncm_vector_set (y_comp, i, ncm_vector_get (sp->y, j));
    }
  }

  {
    NcmSplineVec *sv = ncm_spline_vec_new_gpa (base_spline, xv, yv, TRUE);

    ncm_vector_free (xv);
    g_ptr_array_unref (yv);

    return sv;
  }
}

/**
 * ncm_function_sample_set_log_vals:
 * @fss: a #NcmFunctionSampleSet
 *
 * Logs all sample values in @fss for debugging purposes. This prints the
 * x position, vector components, and OK flag for each sample.
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

    g_string_append_printf (str, "]  ok = %d", sp->ok);

    g_message ("%s\n", str->str);
    g_string_free (str, TRUE);
  }
}

