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

NcmFunctionSampleSet *ncm_function_sample_set_new (const guint len);
NcmFunctionSampleSet *ncm_function_sample_set_ref (NcmFunctionSampleSet *fss);

void ncm_function_sample_set_free (NcmFunctionSampleSet *fss);
void ncm_function_sample_set_clear (NcmFunctionSampleSet **fss);

void ncm_function_sample_set_add (NcmFunctionSampleSet *fss, const gdouble x, NcmVector *y);
void ncm_function_sample_set_insert_before (NcmFunctionSampleSet *fss, const guint index, const gdouble x, NcmVector *y);
void ncm_function_sample_set_insert_after (NcmFunctionSampleSet *fss, const guint index, const gdouble x, NcmVector *y);

guint ncm_function_sample_set_get_len (NcmFunctionSampleSet *fss);
guint ncm_function_sample_set_get_nsamples (NcmFunctionSampleSet *fss);

gdouble ncm_function_sample_set_peek_x (NcmFunctionSampleSet *fss, const guint index);
NcmVector *ncm_function_sample_set_peek_y (NcmFunctionSampleSet *fss, const guint index);
void ncm_function_sample_set_get_sample (NcmFunctionSampleSet *fss, const guint index, gdouble *x, NcmVector **y, gint *ok);

gint ncm_function_sample_set_get_ok (NcmFunctionSampleSet *fss, const guint index);
void ncm_function_sample_set_set_ok (NcmFunctionSampleSet *fss, const guint index, const gint ok);
void ncm_function_sample_set_inc_ok (NcmFunctionSampleSet *fss, const guint index);
void ncm_function_sample_set_reset_ok (NcmFunctionSampleSet *fss);

NcmSplineVec *ncm_function_sample_set_to_spline_vec (NcmFunctionSampleSet *fss, NcmSpline *base_spline);

void ncm_function_sample_set_log_vals (NcmFunctionSampleSet *fss);

G_END_DECLS

#endif /* _NCM_FUNCTION_SAMPLE_SET_H_ */

