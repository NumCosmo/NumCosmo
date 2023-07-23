/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_diff.h
 *
 *  Fri July 21 12:59:15 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_diff.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_DIFF_H_
#define _NCM_DIFF_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_DIFF (ncm_diff_get_type ())

G_DECLARE_FINAL_TYPE (NcmDiff, ncm_diff, NCM, DIFF, GObject);

/**
 * NcmDiffFuncNtoM:
 * @x: function argument
 * @y: function value
 * @user_data: (nullable): user data
 *
 * Functon $f(x)$ call back.
 *
 */
typedef void (*NcmDiffFuncNtoM) (NcmVector *x, NcmVector *y, gpointer user_data);

/**
 * NcmDiffFunc1toM:
 * @x: function argument
 * @y: function value
 * @user_data: (nullable): user data
 *
 * Functon $f(x)$ call back.
 *
 */
typedef void (*NcmDiffFunc1toM) (const gdouble x, NcmVector *y, gpointer user_data);

/**
 * NcmDiffFuncNto1:
 * @x: function argument
 * @user_data: (nullable): user data
 *
 * Functon $f(x)$ call back.
 *
 */
typedef gdouble (*NcmDiffFuncNto1) (NcmVector *x, gpointer user_data);

/**
 * NcmDiffFunc1to1:
 * @x: function argument
 * @user_data: (nullable): user data
 *
 * Functon $f(x)$ call back.
 *
 */
typedef gdouble (*NcmDiffFunc1to1) (const gdouble x, gpointer user_data);

NcmDiff *ncm_diff_new (void);
NcmDiff *ncm_diff_ref (NcmDiff *diff);

void ncm_diff_free (NcmDiff *diff);
void ncm_diff_clear (NcmDiff **diff);

guint ncm_diff_get_max_order (NcmDiff *diff);
gdouble ncm_diff_get_richardson_step (NcmDiff *diff);
gdouble ncm_diff_get_round_off_pad (NcmDiff *diff);
gdouble ncm_diff_get_trunc_error_pad (NcmDiff *diff);
gdouble ncm_diff_get_ini_h (NcmDiff *diff);

void ncm_diff_set_max_order (NcmDiff *diff, const guint maxorder);
void ncm_diff_set_richardson_step (NcmDiff *diff, const gdouble rs);
void ncm_diff_set_round_off_pad (NcmDiff *diff, const gdouble roff_pad);
void ncm_diff_set_trunc_error_pad (NcmDiff *diff, const gdouble terr_pad);
void ncm_diff_set_ini_h (NcmDiff *diff, const gdouble ini_h);

void ncm_diff_log_central_tables (NcmDiff *diff);
void ncm_diff_log_forward_tables (NcmDiff *diff);
void ncm_diff_log_backward_tables (NcmDiff *diff);

GArray *ncm_diff_rf_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr);
GArray *ncm_diff_rc_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr);
GArray *ncm_diff_rc_d2_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr);

GArray *ncm_diff_rf_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr);
GArray *ncm_diff_rc_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr);
GArray *ncm_diff_rc_d2_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr);

GArray *ncm_diff_rf_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr);
GArray *ncm_diff_rc_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr);
GArray *ncm_diff_rc_d2_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr);

GArray *ncm_diff_rf_Hessian_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr);

gdouble ncm_diff_rf_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err);
gdouble ncm_diff_rc_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err);
gdouble ncm_diff_rc_d2_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err);

G_END_DECLS

#endif /* _NCM_DIFF_H_ */

