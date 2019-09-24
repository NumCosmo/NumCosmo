/***************************************************************************
 *            ncm_fit.h
 *
 *  Fri Aug 15 15:26:31 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_NLOPT_2_2
#include <nlopt.h>
#endif /* HAVE_NLOPT_2_2 */

#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_FIT             (ncm_fit_get_type ())
#define NCM_FIT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT, NcmFit))
#define NCM_FIT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT, NcmFitClass))
#define NCM_IS_FIT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT))
#define NCM_IS_FIT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT))
#define NCM_FIT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT, NcmFitClass))

typedef struct _NcmFitClass NcmFitClass;
typedef struct _NcmFit NcmFit;

/**
 * NcmFitType:
 * @NCM_FIT_TYPE_GSL_LS: GSL Least Squares
 * @NCM_FIT_TYPE_GSL_MM: GSL Multidimensional Minimization
 * @NCM_FIT_TYPE_GSL_MMS: GSL Multidimensional Minimization (simplex)
 * @NCM_FIT_TYPE_LEVMAR: Levmar Least Squares Library
 * @NCM_FIT_TYPE_NLOPT: Non Linear Optimization (NLOpt)
 *
 * FIXME
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
 * @NCM_FIT_GRAD_ANALYTICAL: FIXME
 * @NCM_FIT_GRAD_NUMDIFF_FORWARD: FIXME
 * @NCM_FIT_GRAD_NUMDIFF_CENTRAL: FIXME
 * @NCM_FIT_GRAD_NUMDIFF_ACCURATE: FIXME
 *
 * FIXME
 */
typedef enum _NcmFitGradType
{
  NCM_FIT_GRAD_ANALYTICAL = 0,
  NCM_FIT_GRAD_NUMDIFF_FORWARD,
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
 * FIXME
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
 * @NCM_FIT_RUN_MSGS_NONE: FIXME
 * @NCM_FIT_RUN_MSGS_SIMPLE: FIXME
 * @NCM_FIT_RUN_MSGS_FULL: FIXME
 *
 * FIXME
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
};

struct _NcmFit
{
  /*< private >*/
  GObject parent_instance;
  NcmLikelihood *lh;
  NcmMSet *mset;
  NcmFitState *fstate;
  NcmFitRunMsgs mtype;
  NcmFitGrad grad;
  guint maxiter;
  gdouble m2lnL_reltol;
  gdouble m2lnL_abstol;
  gdouble params_reltol;
  GTimer *timer;
  NcmObjArray *equality_constraints;
	GArray *equality_constraints_tot;
  NcmObjArray *inequality_constraints;
	GArray *inequality_constraints_tot;
  NcmFit *sub_fit;
  NcmDiff *diff;
};

GType ncm_fit_get_type (void) G_GNUC_CONST;

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
guint ncm_fit_get_maxiter (NcmFit *fit);
gdouble ncm_fit_get_m2lnL_reltol (NcmFit *fit);
gdouble ncm_fit_get_m2lnL_abstol (NcmFit *fit);
gdouble ncm_fit_get_params_reltol (NcmFit *fit);

NcmMSet *ncm_fit_peek_mset (NcmFit *fit);

NCM_INLINE void ncm_fit_params_set (NcmFit *fit, guint i, const gdouble x);
NCM_INLINE void ncm_fit_params_set_vector (NcmFit *fit, NcmVector *x);
NCM_INLINE void ncm_fit_params_set_vector_offset (NcmFit *fit, NcmVector *x, guint offset);
NCM_INLINE void ncm_fit_params_set_array (NcmFit *fit, const gdouble *x);
NCM_INLINE void ncm_fit_params_set_gsl_vector (NcmFit *fit, const gsl_vector *x);
NCM_INLINE void ncm_fit_params_update (NcmFit *fit);

void ncm_fit_add_equality_constraint (NcmFit *fit, NcmMSetFunc *func, const gdouble tot);
void ncm_fit_add_inequality_constraint (NcmFit *fit, NcmMSetFunc *func, const gdouble tot);
void ncm_fit_remove_equality_constraints (NcmFit *fit);
void ncm_fit_remove_inequality_constraints (NcmFit *fit);
guint ncm_fit_has_equality_constraints (NcmFit *fit);
guint ncm_fit_has_inequality_constraints (NcmFit *fit);

gboolean ncm_fit_is_least_squares (NcmFit *fit);

const gchar *ncm_fit_get_desc (NcmFit *fit);
void ncm_fit_log_info (NcmFit *fit);
void ncm_fit_log_covar (NcmFit *fit);
void ncm_fit_log_start (NcmFit *fit);
void ncm_fit_log_state (NcmFit *fit);
void ncm_fit_log_step (NcmFit *fit);
void ncm_fit_log_step_error (NcmFit *fit, const gchar *strerror, ...);
void ncm_fit_log_end (NcmFit *fit);

void ncm_fit_fishermatrix_print (NcmFit *fit, FILE *out, gchar *header);

NCM_INLINE void ncm_fit_data_m2lnL_val (NcmFit *fit, gdouble *data_m2lnL);
NCM_INLINE void ncm_fit_priors_m2lnL_val (NcmFit *fit, gdouble *priors_m2lnL);

NCM_INLINE void ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL);
NCM_INLINE void ncm_fit_ls_f (NcmFit *fit, NcmVector *f);

NCM_INLINE void ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df);
void ncm_fit_m2lnL_grad_an (NcmFit *fit, NcmVector *df);
void ncm_fit_m2lnL_grad_nd_fo (NcmFit *fit, NcmVector *grad);
void ncm_fit_m2lnL_grad_nd_ce (NcmFit *fit, NcmVector *grad);
void ncm_fit_m2lnL_grad_nd_ac (NcmFit *fit, NcmVector *grad);
void ncm_fit_m2lnL_hessian_nd_ce (NcmFit *fit, NcmMatrix *hessian);

NCM_INLINE void ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df);
void ncm_fit_m2lnL_val_grad_an (NcmFit *fit, gdouble *result, NcmVector *df);
void ncm_fit_m2lnL_val_grad_nd_fo (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
void ncm_fit_m2lnL_val_grad_nd_ce (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
void ncm_fit_m2lnL_val_grad_nd_ac (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);

NCM_INLINE void ncm_fit_ls_J (NcmFit *fit, NcmMatrix *J);
void ncm_fit_ls_J_an (NcmFit *fit, NcmMatrix *J);
void ncm_fit_ls_J_nd_fo (NcmFit *fit, NcmMatrix *J);
void ncm_fit_ls_J_nd_ce (NcmFit *fit, NcmMatrix *J);
void ncm_fit_ls_J_nd_ac (NcmFit *fit, NcmMatrix *J);

NCM_INLINE void ncm_fit_ls_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J);
void ncm_fit_ls_f_J_an (NcmFit *fit, NcmVector *f, NcmMatrix *J);
void ncm_fit_ls_f_J_nd_fo (NcmFit *fit, NcmVector *f, NcmMatrix *J);
void ncm_fit_ls_f_J_nd_ce (NcmFit *fit, NcmVector *f, NcmMatrix *J);
void ncm_fit_ls_f_J_nd_ac (NcmFit *fit, NcmVector *f, NcmMatrix *J);

void ncm_fit_fisher_to_covar (NcmFit *fit, NcmMatrix *fisher);
void ncm_fit_obs_fisher (NcmFit *fit);
void ncm_fit_fisher (NcmFit *fit);

void ncm_fit_numdiff_m2lnL_hessian (NcmFit *fit, NcmMatrix *H, gdouble reltol);
void ncm_fit_numdiff_m2lnL_covar (NcmFit *fit);
void ncm_fit_ls_covar (NcmFit *fit);
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

gdouble ncm_fit_residual_ks_test (NcmFit *fit, gdouble *o_mean, gdouble *o_sd, gdouble *o_skew, gdouble *o_kurtosis, gdouble *o_max);
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

#ifndef _NCM_FIT_INLINE_H_
#define _NCM_FIT_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void 
ncm_fit_params_set (NcmFit *fit, guint i, const gdouble x)
{
  ncm_mset_fparam_set (fit->mset, i, x);
  ncm_fit_params_update (fit);
}

NCM_INLINE void 
ncm_fit_params_set_vector (NcmFit *fit, NcmVector *x)
{
  ncm_mset_fparams_set_vector (fit->mset, x);
  ncm_fit_params_update (fit);
}

NCM_INLINE void 
ncm_fit_params_set_vector_offset (NcmFit *fit, NcmVector *x, guint offset)
{
  ncm_mset_fparams_set_vector_offset (fit->mset, x, offset);
  ncm_fit_params_update (fit);
}

NCM_INLINE void 
ncm_fit_params_set_array (NcmFit *fit, const gdouble *x)
{
  ncm_mset_fparams_set_array (fit->mset, x);
  ncm_fit_params_update (fit);
}

NCM_INLINE void 
ncm_fit_params_set_gsl_vector (NcmFit *fit, const gsl_vector *x)
{
  ncm_mset_fparams_set_gsl_vector (fit->mset, x);
  ncm_fit_params_update (fit);
}

NCM_INLINE void 
ncm_fit_params_update (NcmFit *fit)
{
  if (fit->sub_fit != NULL)
  {
    NcmFitRunMsgs mlevel = fit->mtype == NCM_FIT_RUN_MSGS_FULL ? NCM_FIT_RUN_MSGS_SIMPLE : NCM_FIT_RUN_MSGS_NONE;    
    ncm_fit_run (fit->sub_fit, mlevel);
  }
}

NCM_INLINE void
ncm_fit_data_m2lnL_val (NcmFit *fit, gdouble *data_m2lnL)
{
  ncm_dataset_m2lnL_val (fit->lh->dset, fit->mset, data_m2lnL);
}

NCM_INLINE void
ncm_fit_priors_m2lnL_val (NcmFit *fit, gdouble *priors_m2lnL)
{
  ncm_likelihood_priors_m2lnL_val (fit->lh, fit->mset, priors_m2lnL);
}

NCM_INLINE void
ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL)
{
  ncm_likelihood_m2lnL_val (fit->lh, fit->mset, m2lnL);
  fit->fstate->func_eval++;
}

NCM_INLINE void
ncm_fit_ls_f (NcmFit *fit, NcmVector *f)
{
  ncm_likelihood_leastsquares_f (fit->lh, fit->mset, f);
  fit->fstate->func_eval++;
}

NCM_INLINE void
ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df)
{
  fit->grad.m2lnL_grad (fit, df);
}

NCM_INLINE void
ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df)
{
  fit->grad.m2lnL_val_grad (fit, result, df);
}

NCM_INLINE void
ncm_fit_ls_J (NcmFit *fit, NcmMatrix *J)
{
  fit->grad.ls_J (fit, J);
}

NCM_INLINE void
ncm_fit_ls_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J)
{
  fit->grad.ls_f_J (fit, f, J);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_FIT_INLINE_H_ */
