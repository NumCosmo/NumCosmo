/***************************************************************************
 *            fit.h
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
#ifdef HAVE_NLOPT_2_2
#include <nlopt.h>
#endif /* HAVE_NLOPT_2_2 */

#include <glib-object.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT             (ncm_fit_get_type ())
#define NCM_FIT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT, NcmFit))
#define NCM_FIT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT, NcmFitClass))
#define NCM_IS_FIT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT))
#define NCM_IS_FIT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT))
#define NCM_FIT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT, NcmFitClass))

typedef struct _NcmFitClass NcmFitClass;
typedef struct _NcmFit NcmFit;

typedef void (*NcmFitLsJ) (NcmFit *fit, NcmMatrix *J);
typedef void (*NcmFitLsFJ) (NcmFit *fit, NcmVector *f, NcmMatrix *J);
typedef void (*NcmFitM2lnLGrad) (NcmFit *fit, NcmVector *grad);
typedef void (*NcmFitM2lnLValGrad) (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);

/**
 * NcmFitType:
 * @NCM_FIT_TYPE_LEAST_SQUARES: FIXME
 * @NCM_FIT_TYPE_MULTIMIN: FIXME
 * @NCM_FIT_TYPE_SIMPLEX: FIXME
 * @NCM_FIT_TYPE_LEVMAR_DER: FIXME
 * @NCM_FIT_TYPE_LEVMAR_DIF: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_COBYLA: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_BOBYQA: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_NEWUOA: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_NEWUOA_BOUND: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_PRAXIS: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_NELDERMEAD: FIXME
 * @NCM_FIT_TYPE_NLOPT_LN_SBPLX: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_MMA: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_LBFGS: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND_RESTART: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_TNEWTON_RESTART: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_TNEWTON: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_VAR1: FIXME
 * @NCM_FIT_TYPE_NLOPT_LD_VAR2: FIXME
 *
 * FIXME
 */
typedef enum _NcmFitType
{
  NCM_FIT_TYPE_LEAST_SQUARES = 0,
  NCM_FIT_TYPE_MULTIMIN,
  NCM_FIT_TYPE_SIMPLEX,
  NCM_FIT_TYPE_LEVMAR_DER,
  NCM_FIT_TYPE_LEVMAR_DIF,
  NCM_FIT_TYPE_NLOPT_LN_COBYLA,
  NCM_FIT_TYPE_NLOPT_LN_BOBYQA,
  NCM_FIT_TYPE_NLOPT_LN_NEWUOA,
  NCM_FIT_TYPE_NLOPT_LN_NEWUOA_BOUND,
  NCM_FIT_TYPE_NLOPT_LN_PRAXIS,
  NCM_FIT_TYPE_NLOPT_LN_NELDERMEAD,
  NCM_FIT_TYPE_NLOPT_LN_SBPLX,
  NCM_FIT_TYPE_NLOPT_LD_MMA,
  NCM_FIT_TYPE_NLOPT_LD_LBFGS,
  NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND_RESTART,
  NCM_FIT_TYPE_NLOPT_LD_TNEWTON_PRECOND,
  NCM_FIT_TYPE_NLOPT_LD_TNEWTON_RESTART,
  NCM_FIT_TYPE_NLOPT_LD_TNEWTON,
  NCM_FIT_TYPE_NLOPT_LD_VAR1,
  NCM_FIT_TYPE_NLOPT_LD_VAR2,
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
};

struct _NcmFit
{
  /*< private >*/
  GObject parent_instance;
  NcLikelihood *lh;
  NcmMSet *mset;
  NcmFitType type;
  NcmFitGradType gtype;
  NcmFitRunMsgs mtype;
  const gchar *solver_name;
  const gchar *diff_name;
  guint data_len;
  guint fparam_len;
  guint dof;
  guint niter;
  guint n_func_eval;
  guint n_grad_eval;
  gdouble m2lnL_prec;
  gdouble m2lnL_prec_target;
  gdouble params_prec;
  gdouble params_prec_target;
  GTimer *timer;
  gsl_multifit_fdfsolver *least_squares;
  gsl_multimin_fdfminimizer *multimin;
  gsl_multimin_fminimizer *simplex;
#ifdef HAVE_NLOPT_2_2
  nlopt_opt nlopt;
  nlopt_algorithm nlopt_algo;
#endif /* HAVE_NLOPT_2_2 */
  gpointer levmar;
  gboolean bootstrap;
  gboolean jackknife;
  NcmVector *x;
  NcmVector *f;
  NcmVector *df;
  NcmMatrix *J;
  NcmVector *full_df;
  NcmMatrix *full_J;
  NcmMatrix *covar;
  NcmMatrix *hessian;
  NcmVector *jackpv;
  gsl_vector_uint *bs;
  gdouble sqrt_m2lnL;
  gdouble m2lnL;
  gdouble m2lnL_dof;
  NcmFitLsJ ls_J;
  NcmFitLsFJ ls_f_J;
  NcmFitM2lnLGrad m2lnL_grad;
  NcmFitM2lnLValGrad m2lnL_val_grad;
};

GType ncm_fit_get_type (void) G_GNUC_CONST;

NcmFit *ncm_fit_new (NcLikelihood *lh, NcmMSet *mset, NcmFitType type, NcmFitGradType gtype);
NcmFit *ncm_fit_ref (NcmFit *fit);
NcmFit *ncm_fit_copy (NcmFit *fit);
void ncm_fit_free (NcmFit *fit);

void ncm_fit_save (NcmFit *fit, NcmVector *x, NcmVector *f, NcmMatrix *J);

void ncm_fit_log_info (NcmFit *fit);
void ncm_fit_log_covar (NcmFit *fit);
void ncm_fit_log_start (NcmFit *fit);
void ncm_fit_log_state (NcmFit *fit, gdouble m2lnL);
void ncm_fit_log_step (NcmFit *fit, gdouble m2lnL);
void ncm_fit_log_step_error (NcmFit *fit, const gchar *strerror, ...);
void ncm_fit_log_end (NcmFit *fit);

void ncm_fit_fishermatrix_print (NcmFit *fit, FILE *out, gchar *header);

void ncm_fit_gen_bootstrap (NcmFit *fit);

void ncm_fit_leastsquares_f (NcmFit *fit, NcmVector *f);
void ncm_fit_leastsquares_J (NcmFit *fit, NcmMatrix *J);
void ncm_fit_leastsquares_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J);

void ncm_fit_numdiff_forward_leastsquares_J (NcmFit *fit, NcmMatrix *J);
void ncm_fit_numdiff_forward_leastsquares_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J);
void ncm_fit_numdiff_central_leastsquares_J (NcmFit *fit, NcmMatrix *J);
void ncm_fit_numdiff_central_leastsquares_f_J (NcmFit *fit, NcmVector *f, NcmMatrix *J);

void ncm_fit_m2lnL_val (NcmFit *fit, gdouble *m2lnL);
void ncm_fit_data_m2lnL_val (NcmFit *fit, gdouble *data_m2lnL);
void ncm_fit_priors_m2lnL_val (NcmFit *fit, gdouble *priors_m2lnL);
void ncm_fit_m2lnL_grad (NcmFit *fit, NcmVector *df);
void ncm_fit_m2lnL_val_grad (NcmFit *fit, gdouble *result, NcmVector *df);

void ncm_fit_numdiff_accurate_m2lnL_grad (NcmFit *fit, NcmVector *grad);
void ncm_fit_numdiff_accurate_m2lnL_val_grad (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
void ncm_fit_numdiff_forward_m2lnL_grad (NcmFit *fit, NcmVector *grad);
void ncm_fit_numdiff_forward_m2lnL_val_grad (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);
void ncm_fit_numdiff_central_m2lnL_grad (NcmFit *fit, NcmVector *grad);
void ncm_fit_numdiff_central_m2lnL_val_grad (NcmFit *fit, gdouble *m2lnL, NcmVector *grad);

void ncm_fit_numdiff_m2lnL_hessian (NcmFit *fit, NcmMatrix *H);
void ncm_fit_numdiff_m2lnL_covar (NcmFit *fit);
void ncm_fit_covar_leastsquares_calc (NcmFit *fit);

gdouble ncm_fit_covar_var (NcmFit *fit, NcmModelID gmid, guint pid);
gdouble ncm_fit_covar_sd (NcmFit *fit, NcmModelID gmid, guint pid);
gdouble ncm_fit_covar_cov (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2);
gdouble ncm_fit_covar_cor (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2);

gdouble ncm_fit_covar_fparam_var (NcmFit *fit, guint fpi);
gdouble ncm_fit_covar_fparam_sd (NcmFit *fit, guint fpi);
gdouble ncm_fit_covar_fparam_cov (NcmFit *fit, guint fpi1, guint fpi2);
gdouble ncm_fit_covar_fparam_cor (NcmFit *fit, guint fpi1, guint fpi2);

void ncm_fit_lr_test_range (NcmFit *fit, NcmModelID gmid, guint pid, gdouble start, gdouble stop, gdouble step);
void ncm_fit_dprob (NcmFit *fit, NcmModelID gmid, guint pid, gdouble a, gdouble b, gdouble step, gdouble norm);
gdouble ncm_fit_lr_test (NcmFit *fit, NcmModelID gmid, guint pid, gdouble val, gint dof);
gdouble ncm_fit_prob (NcmFit *fit, NcmModelID gmid, guint pid, gdouble a, gdouble b);
gdouble ncm_fit_chisq_test (NcmFit *fit, size_t bins);

gboolean ncm_fit_run (NcmFit *fit, gint niters, NcmFitRunMsgs mtype);
gboolean ncm_fit_cr (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2, gdouble p);
GList *ncm_fit_cr2 (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2, gdouble p);
GList *ncm_fit_cr2_fisher (NcmFit *fit, NcmModelID gmid1, guint pid1, NcmModelID gmid2, guint pid2, gdouble p);
gboolean ncm_fit_cr_1dim (NcmFit *fit, NcmModelID gmid, guint pid, gdouble p, gint nu, gdouble *err_inf, gdouble *err_sup);
gboolean ncm_fit_cr_points_print (GList *points, FILE *out);
gboolean ncm_fit_cr_points_free (GList *points);

gdouble ncm_fit_type_constrain_error (NcmFit *fit, gdouble p, gint nu, gdouble dir, NcmMSetFunc *func, gdouble z, gboolean walk);
gdouble ncm_fit_function_error (NcmFit *fit, NcmMSetFunc *func, gdouble z, gboolean pretty_print);
gdouble ncm_fit_function_cov (NcmFit *fit, NcmMSetFunc *func1, gdouble z1, NcmMSetFunc *func2, gdouble z2, gboolean pretty_print);

NcmMatrix *ncm_fit_montecarlo_matrix (NcmFit *fit, NcmMSet *mset, guint maxiter, guint ni, guint nf, NcmFitRunMsgs mtype);
void ncm_fit_montecarlo_matrix_print (NcmFit *fit, NcmMatrix *param_matrix);
void ncm_fit_montecarlo_matrix_mean_covar (NcmFit *fit, NcmMatrix *param_matrix);

#define NCM_FIT_NUMDIFF_SCALE (1.0e-4)
#define NCM_FIT_NPARAM(fit) ((fit)->pt->nfree)
#define NCM_FIT_DEFAULT_M2LNL_ABSTOL (0.0)
#define NCM_FIT_DEFAULT_M2LNL_RELTOL (1e-13)

G_END_DECLS

#endif /* _NCM_FIT_H_ */
