/***************************************************************************
 *            ncm_fit.h
 *
 *  Fri Aug 15 15:26:31 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_FIT_H_
#define _NCM_FIT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_diff.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/math/ncm_fit_state.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT (ncm_fit_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmFit, ncm_fit, NCM, FIT, GObject)

/**
 * NcmFitType:
 * @NCM_FIT_TYPE_GSL_LS: GSL Least Squares
 * @NCM_FIT_TYPE_GSL_MM: GSL Multidimensional Minimization
 * @NCM_FIT_TYPE_GSL_MMS: GSL Multidimensional Minimization (simplex)
 * @NCM_FIT_TYPE_LEVMAR: Levmar Least Squares Library
 * @NCM_FIT_TYPE_NLOPT: Non Linear Optimization (NLOpt)
 *
 * Defines the subclasse of NcmFit to be used.
 *
 */
typedef enum _NcmFitType
{
  NCM_FIT_TYPE_GSL_LS = 0,
  NCM_FIT_TYPE_GSL_MM,
  NCM_FIT_TYPE_GSL_MMS,
  NCM_FIT_TYPE_LEVMAR,
  NCM_FIT_TYPE_NLOPT,
} NcmFitType;

/**
 * NcmFitGradType:
 * @NCM_FIT_GRAD_NUMDIFF_FORWARD: Numerical gradient (forward)
 * @NCM_FIT_GRAD_NUMDIFF_CENTRAL: Numerical gradient (central)
 * @NCM_FIT_GRAD_NUMDIFF_ACCURATE: Numerical gradient (accurate)
 *
 * Defines the type of gradient calculation.
 *
 */
typedef enum _NcmFitGradType /*< enum,prefix=NCM_FIT_GRAD >*/
{
  NCM_FIT_GRAD_NUMDIFF_FORWARD = 0,
  NCM_FIT_GRAD_NUMDIFF_CENTRAL,
  NCM_FIT_GRAD_NUMDIFF_ACCURATE,
} NcmFitGradType;

typedef void (*_NcmFitLSJ) (NcmFit *fit, NcmMatrix *J);
typedef void (*_NcmFitLSFJ) (NcmFit *fit, NcmVector *f, NcmMatrix *J);
typedef void (*_NcmFitM2lnLGrad) (NcmFit *fit, NcmVector *grad);
typedef void (*_NcmFitM2lnLValGrad) (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);

/**
 * NcmFitGrad:
 *
 * Container for gradient functions.
 *
 */
typedef struct _NcmFitGrad
{
  /*< private >*/
  NcmFitGradType gtype;
  const gchar *diff_name;
  _NcmFitLSJ ls_J;
  _NcmFitLSFJ ls_f_J;
  _NcmFitM2lnLGrad m2lnL_grad;
  _NcmFitM2lnLValGrad m2lnL_val_grad;
} NcmFitGrad;

/**
 * NcmFitRunMsgs:
 * @NCM_FIT_RUN_MSGS_NONE: Messages disabled
 * @NCM_FIT_RUN_MSGS_SIMPLE: Messages enabled
 * @NCM_FIT_RUN_MSGS_FULL: Messages enabled (full)
 *
 * Defines the type of messages to be printed during the fit.
 *
 */
typedef enum _NcmFitRunMsgs
{
  NCM_FIT_RUN_MSGS_NONE = 0,
  NCM_FIT_RUN_MSGS_SIMPLE,
  NCM_FIT_RUN_MSGS_FULL,
} NcmFitRunMsgs;

struct _NcmFitClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmFit *(*copy_new) (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
  void (*reset) (NcmFit *fit);
  gboolean (*run) (NcmFit *fit, NcmFitRunMsgs mtype);
  const gchar *(*get_desc) (NcmFit *fit);
  gboolean is_least_squares;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[13];
};

NcmFit *ncm_fit_new (NcmFitType ftype, gchar *algo_name, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_ref (NcmFit *fit);
NcmFit *ncm_fit_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_dup (NcmFit *fit, NcmSerialize *ser);
void ncm_fit_free (NcmFit *fit);
void ncm_fit_clear (NcmFit **fit);

void ncm_fit_set_sub_fit (NcmFit *fit, NcmFit *sub_fit);
NcmFit *ncm_fit_get_sub_fit (NcmFit *fit);

void ncm_fit_set_grad_type (NcmFit *fit, NcmFitGradType gtype);
void ncm_fit_set_maxiter (NcmFit *fit, guint maxiter);
void ncm_fit_set_m2lnL_reltol (NcmFit *fit, gdouble tol);
void ncm_fit_set_m2lnL_abstol (NcmFit *fit, gdouble tol);
void ncm_fit_set_params_reltol (NcmFit *fit, gdouble tol);
void ncm_fit_set_messages (NcmFit *fit, NcmFitRunMsgs mtype);

NcmFitGradType ncm_fit_get_grad_type (NcmFit *fit);
guint ncm_fit_get_maxiter (NcmFit *fit);
gdouble ncm_fit_get_m2lnL_reltol (NcmFit *fit);
gdouble ncm_fit_get_m2lnL_abstol (NcmFit *fit);
gdouble ncm_fit_get_params_reltol (NcmFit *fit);
NcmFitRunMsgs ncm_fit_get_messages (NcmFit *fit);
gboolean ncm_fit_is_least_squares (NcmFit *fit);

NcmMSet *ncm_fit_peek_mset (NcmFit *fit);
NcmFitState *ncm_fit_peek_state (NcmFit *fit);
NcmLikelihood *ncm_fit_peek_likelihood (NcmFit *fit);
NcmDiff *ncm_fit_peek_diff (NcmFit *fit);

void ncm_fit_params_set (NcmFit *fit, guint i, const gdouble x);
void ncm_fit_params_set_vector (NcmFit *fit, NcmVector *x);
void ncm_fit_params_set_vector_offset (NcmFit *fit, NcmVector *x, guint offset);
void ncm_fit_params_set_array (NcmFit *fit, const gdouble *x);
void ncm_fit_params_set_gsl_vector (NcmFit *fit, const gsl_vector *x);
void ncm_fit_params_update (NcmFit *fit);

void ncm_fit_add_equality_constraint (NcmFit *fit, NcmMSetFunc *func, const gdouble tot);
void ncm_fit_add_inequality_constraint (NcmFit *fit, NcmMSetFunc *func, const gdouble tot);
void ncm_fit_remove_equality_constraints (NcmFit *fit);
void ncm_fit_remove_inequality_constraints (NcmFit *fit);
guint ncm_fit_equality_constraints_len (NcmFit *fit);
guint ncm_fit_inequality_constraints_len (NcmFit *fit);
void ncm_fit_get_equality_constraint (NcmFit *fit, guint i, NcmMSetFunc **func, gdouble *tot);
void ncm_fit_get_inequality_constraint (NcmFit *fit, guint i, NcmMSetFunc **func, gdouble *tot);

const gchar *ncm_fit_get_desc (NcmFit *fit);
void ncm_fit_log_info (NcmFit *fit);
void ncm_fit_log_covar (NcmFit *fit);
void ncm_fit_log_start (NcmFit *fit);
void ncm_fit_log_state (NcmFit *fit);
void ncm_fit_log_step (NcmFit *fit);
void ncm_fit_log_step_error (NcmFit *fit, const gchar *strerror, ...);
void ncm_fit_log_end (NcmFit *fit);

void ncm_fit_data_m2lnL_val (NcmFit *fit, gdouble *data_m2lnL);
void ncm_fit_priors_m2lnL_val (NcmFit *fit, gdouble *priors_m2lnL);

void ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL);
void ncm_fit_ls_f (NcmFit *fit, NcmVector *f);

void ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df);
void ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df);

void ncm_fit_ls_J (NcmFit *fit, NcmMatrix *J);
void ncm_fit_ls_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J);

void ncm_fit_obs_fisher (NcmFit *fit);
void ncm_fit_ls_fisher (NcmFit *fit);
void ncm_fit_fisher (NcmFit *fit);
void ncm_fit_numdiff_m2lnL_covar (NcmFit *fit);

gdouble ncm_fit_numdiff_m2lnL_lndet_covar (NcmFit *fit);
NcmMatrix *ncm_fit_get_covar (NcmFit *fit);

gdouble ncm_fit_covar_var (NcmFit *fit, NcmModelID mid, guint pid);
gdouble ncm_fit_covar_sd (NcmFit *fit, NcmModelID mid, guint pid);
gdouble ncm_fit_covar_cov (NcmFit *fit, NcmModelID mid1, guint pid1, NcmModelID mid2, guint pid2);
gdouble ncm_fit_covar_cor (NcmFit *fit, NcmModelID mid1, guint pid1, NcmModelID mid2, guint pid2);

gdouble ncm_fit_covar_fparam_var (NcmFit *fit, guint fpi);
gdouble ncm_fit_covar_fparam_sd (NcmFit *fit, guint fpi);
gdouble ncm_fit_covar_fparam_cov (NcmFit *fit, guint fpi1, guint fpi2);
gdouble ncm_fit_covar_fparam_cor (NcmFit *fit, guint fpi1, guint fpi2);

void ncm_fit_lr_test_range (NcmFit *fit, NcmModelID mid, guint pid, gdouble start, gdouble stop, gdouble step);
void ncm_fit_dprob (NcmFit *fit, NcmModelID mid, guint pid, gdouble a, gdouble b, gdouble step, gdouble norm);
gdouble ncm_fit_lr_test (NcmFit *fit, NcmModelID mid, guint pid, gdouble val, gint dof);
gdouble ncm_fit_prob (NcmFit *fit, NcmModelID mid, guint pid, gdouble a, gdouble b);
gdouble ncm_fit_chisq_test (NcmFit *fit, size_t bins);

void ncm_fit_reset (NcmFit *fit);
gboolean ncm_fit_run (NcmFit *fit, NcmFitRunMsgs mtype);
void ncm_fit_run_restart (NcmFit *fit, NcmFitRunMsgs mtype, const gdouble abstol, const gdouble reltol, NcmMSet *save_mset, const gchar *mset_file);

gdouble ncm_fit_type_constrain_error (NcmFit *fit, gdouble p, gint nu, gdouble dir, NcmMSetFunc *func, gdouble z, gboolean walk);
void ncm_fit_function_error (NcmFit *fit, NcmMSetFunc *func, gdouble *x, gboolean pretty_print, gdouble *f, gdouble *sigma_f);
gdouble ncm_fit_function_cov (NcmFit *fit, NcmMSetFunc *func1, gdouble z1, NcmMSetFunc *func2, gdouble z2, gboolean pretty_print);

#define NCM_FIT_DEFAULT_M2LNL_RELTOL (1e-8)
#define NCM_FIT_DEFAULT_M2LNL_ABSTOL (0.0)
#define NCM_FIT_DEFAULT_PARAMS_RELTOL (1e-5)
#define NCM_FIT_DEFAULT_MAXITER 100000

G_END_DECLS

#endif /* _NCM_FIT_H_ */

