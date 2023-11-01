/***************************************************************************
 *            ncm_fit_state.h
 *
 *  Thu November 29 15:27:17 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_FIT_STATE_H_
#define _NCM_FIT_STATE_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_FIT_STATE (ncm_fit_state_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitState, ncm_fit_state, NCM, FIT_STATE, GObject)

NcmFitState *ncm_fit_state_new (guint data_len, guint fparam_len, gint dof, gboolean is_least_squares);
NcmFitState *ncm_fit_state_ref (NcmFitState *fstate);
void ncm_fit_state_free (NcmFitState *fstate);
void ncm_fit_state_clear (NcmFitState **fstate);

void ncm_fit_state_set_all (NcmFitState *fstate, guint data_len, guint fparam_len, gint dof, gboolean is_least_squares);
void ncm_fit_state_reset (NcmFitState *fstate);

void ncm_fit_state_set_ls (NcmFitState *fstate, NcmVector *f, NcmMatrix *J);

void ncm_fit_state_set_fparam_len (NcmFitState *fstate, guint fparam_len);
guint ncm_fit_state_get_fparam_len (NcmFitState *fstate);

void ncm_fit_state_set_data_len (NcmFitState *fstate, guint data_len);
guint ncm_fit_state_get_data_len (NcmFitState *fstate);

void ncm_fit_state_set_dof (NcmFitState *fstate, gint dof);
gint ncm_fit_state_get_dof (NcmFitState *fstate);

void ncm_fit_state_add_iter (NcmFitState *fstate, guint niter);
void ncm_fit_state_set_niter (NcmFitState *fstate, guint niter);
guint ncm_fit_state_get_niter (NcmFitState *fstate);

void ncm_fit_state_add_func_eval (NcmFitState *fstate, guint func_eval);
void ncm_fit_state_set_func_eval (NcmFitState *fstate, guint func_eval);
guint ncm_fit_state_get_func_eval (NcmFitState *fstate);

void ncm_fit_state_add_grad_eval (NcmFitState *fstate, guint grad_eval);
void ncm_fit_state_set_grad_eval (NcmFitState *fstate, guint grad_eval);
guint ncm_fit_state_get_grad_eval (NcmFitState *fstate);

void ncm_fit_state_set_m2lnL_prec (NcmFitState *fstate, gdouble prec);
gdouble ncm_fit_state_get_m2lnL_prec (NcmFitState *fstate);

void ncm_fit_state_set_m2lnL_curval (NcmFitState *fstate, gdouble m2lnL_curval);
gdouble ncm_fit_state_get_m2lnL_curval (NcmFitState *fstate);

void ncm_fit_state_set_params_prec (NcmFitState *fstate, gdouble prec);
gdouble ncm_fit_state_get_params_prec (NcmFitState *fstate);

void ncm_fit_state_set_elapsed_time (NcmFitState *fstate, gdouble elapsed_time);
gdouble ncm_fit_state_get_elapsed_time (NcmFitState *fstate);

void ncm_fit_state_set_has_covar (NcmFitState *fstate, gboolean has_covar);
gboolean ncm_fit_state_has_covar (NcmFitState *fstate);

void ncm_fit_state_set_is_best_fit (NcmFitState *fstate, gboolean is_best_fit);
gboolean ncm_fit_state_is_best_fit (NcmFitState *fstate);

void ncm_fit_state_set_is_least_squares (NcmFitState *fstate, gboolean is_least_squares);
gboolean ncm_fit_state_is_least_squares (NcmFitState *fstate);

NcmVector *ncm_fit_state_peek_fparams (NcmFitState *fstate);
NcmMatrix *ncm_fit_state_peek_hessian (NcmFitState *fstate);
NcmMatrix *ncm_fit_state_peek_covar (NcmFitState *fstate);
NcmVector *ncm_fit_state_peek_f (NcmFitState *fstate);
NcmMatrix *ncm_fit_state_peek_J (NcmFitState *fstate);

G_END_DECLS

#endif /* _NCM_FIT_STATE_H_ */

