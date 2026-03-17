/***************************************************************************
 *            ncm_function_sample_set.h
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

#ifndef _NCM_FUNCTION_SAMPLE_SET_H_
#define _NCM_FUNCTION_SAMPLE_SET_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_vec.h>

G_BEGIN_DECLS

#define NCM_TYPE_FUNCTION_SAMPLE_SET (ncm_function_sample_set_get_type ())

G_DECLARE_FINAL_TYPE (NcmFunctionSampleSet, ncm_function_sample_set, NCM, FUNCTION_SAMPLE_SET, GObject)

/**
 * NcmFunctionSampleSetFunc:
 * @x: input value
 * @y: output vector to be filled
 * @user_data: user data
 *
 * Vector-valued function type $\vec{F}: \mathbb{R} \to \mathbb{R}^n$.
 * The function should evaluate at @x and write the result into @y.
 *
 */
typedef void (*NcmFunctionSampleSetFunc) (const gdouble x, NcmVector *y, gpointer user_data);

/**
 * NcmFunctionSampleSetIter:
 * @node: (skip): internal GList node pointer
 * @owner: the #NcmFunctionSampleSet being iterated
 *
 * Iterator for traversing samples in a #NcmFunctionSampleSet.
 * Provides O(1) access to sample data once positioned, making it efficient
 * for sequential traversal and interval operations.
 *
 */
typedef struct _NcmFunctionSampleSetIter
{
  /*< private >*/
  GList *node;
  /*< public >*/
  NcmFunctionSampleSet *owner;
} NcmFunctionSampleSetIter;

#define NCM_TYPE_FUNCTION_SAMPLE_SET_ITER (ncm_function_sample_set_iter_get_type ())
GType ncm_function_sample_set_iter_get_type (void) G_GNUC_CONST;

/* Container creation and management */
NcmFunctionSampleSet *ncm_function_sample_set_new (const guint len);
NcmFunctionSampleSet *ncm_function_sample_set_ref (NcmFunctionSampleSet *fss);

void ncm_function_sample_set_free (NcmFunctionSampleSet *fss);
void ncm_function_sample_set_clear (NcmFunctionSampleSet **fss);

/* Append operations (kept for initial population) */
void ncm_function_sample_set_add (NcmFunctionSampleSet *fss, const gdouble x, NcmVector *y);
void ncm_function_sample_set_add_func (NcmFunctionSampleSet *fss, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data);
void ncm_function_sample_set_add_old (NcmFunctionSampleSet *fss, const gdouble x, NcmVector *y);
void ncm_function_sample_set_add_old_func (NcmFunctionSampleSet *fss, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data);

/* Container-level properties */
guint ncm_function_sample_set_get_len (NcmFunctionSampleSet *fss);
guint ncm_function_sample_set_get_nsamples (NcmFunctionSampleSet *fss);
gdouble ncm_function_sample_set_get_x_min (NcmFunctionSampleSet *fss);
gdouble ncm_function_sample_set_get_x_max (NcmFunctionSampleSet *fss);
gdouble ncm_function_sample_set_get_absmaxF (NcmFunctionSampleSet *fss, const guint i);

/* Container-level operations */
void ncm_function_sample_set_reset_interval_ok (NcmFunctionSampleSet *fss);
void ncm_function_sample_set_mark_all_old (NcmFunctionSampleSet *fss);
gboolean ncm_function_sample_set_all_intervals_ok (NcmFunctionSampleSet *fss, const gint threshold);

/* Spline conversion */
NcmSplineVec *ncm_function_sample_set_to_spline_vec (NcmFunctionSampleSet *fss, NcmSpline *base_spline);
NcmSplineVec *ncm_function_sample_set_to_spline_vec_old (NcmFunctionSampleSet *fss, NcmSpline *base_spline);

/* Refinement */
void ncm_function_sample_set_refine (NcmFunctionSampleSet *fss, const gdouble reltol, const gdouble abstol, NcmSpline *base_spline);

/* Debug */
void ncm_function_sample_set_log_vals (NcmFunctionSampleSet *fss);

/* ============================================================================
 * Iterator API - O(1) access for efficient traversal
 * ============================================================================ */

/* Iterator creation */
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_begin (NcmFunctionSampleSet *fss);
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_end (NcmFunctionSampleSet *fss);
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_copy (NcmFunctionSampleSetIter *iter);
void ncm_function_sample_set_iter_free (NcmFunctionSampleSetIter *iter);

/* Iterator state checks */
gboolean ncm_function_sample_set_iter_is_valid (NcmFunctionSampleSetIter *iter);
gboolean ncm_function_sample_set_iter_has_next (NcmFunctionSampleSetIter *iter);
gboolean ncm_function_sample_set_iter_has_prev (NcmFunctionSampleSetIter *iter);

/* Iterator movement */
void ncm_function_sample_set_iter_next (NcmFunctionSampleSetIter *iter);
void ncm_function_sample_set_iter_prev (NcmFunctionSampleSetIter *iter);

/* Iterator accessors */
gdouble ncm_function_sample_set_iter_get_x (NcmFunctionSampleSetIter *iter);
NcmVector *ncm_function_sample_set_iter_get_y (NcmFunctionSampleSetIter *iter);
gint ncm_function_sample_set_iter_get_interval_ok (NcmFunctionSampleSetIter *iter);
gboolean ncm_function_sample_set_iter_get_new_point (NcmFunctionSampleSetIter *iter);

/* Iterator mutators */
void ncm_function_sample_set_iter_set_interval_ok (NcmFunctionSampleSetIter *iter, const gint interval_ok);
void ncm_function_sample_set_iter_inc_interval_ok (NcmFunctionSampleSetIter *iter);
void ncm_function_sample_set_iter_set_new_point (NcmFunctionSampleSetIter *iter, const gboolean new_point);

/* Iterator-based insertion */
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_insert_after (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmVector *y);
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_insert_after_func (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data);
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_insert_before (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmVector *y);
NcmFunctionSampleSetIter *ncm_function_sample_set_iter_insert_before_func (NcmFunctionSampleSet *fss, NcmFunctionSampleSetIter *iter, const gdouble x, NcmFunctionSampleSetFunc f, gpointer user_data);

/* Iterator helpers for interval operations */
gboolean ncm_function_sample_set_iter_next_pair (NcmFunctionSampleSetIter *iter, NcmFunctionSampleSetIter *next_iter);

G_END_DECLS

#endif /* _NCM_FUNCTION_SAMPLE_SET_H_ */

