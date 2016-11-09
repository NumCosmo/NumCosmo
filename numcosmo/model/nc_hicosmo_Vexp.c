/***************************************************************************
 *            nc_hicosmo_Vexp.c
 *
 *  Fri October 28 13:27:53 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <sandro@isoftware.com.br>
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
 * SECTION:nc_hicosmo_Vexp
 * @title: NcHICosmoVexp
 * @short_description: $\Lambda$CDM model.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_Vexp.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

G_DEFINE_TYPE (NcHICosmoVexp, nc_hicosmo_Vexp, NC_TYPE_HICOSMO);

enum {
  PROP_0,
  PROP_SIZE,
};

typedef struct _NcHICosmoVexpState
{
  gdouble alpha;
  gdouble E;
} NcHICosmoVexpState;

static void
nc_hicosmo_Vexp_init (NcHICosmoVexp *Vexp)
{
  Vexp->cvode_qt  = CVodeCreate (CV_BDF, CV_NEWTON);
  Vexp->cvode_clp = CVodeCreate (CV_BDF, CV_NEWTON);
  Vexp->cvode_clm = CVodeCreate (CV_BDF, CV_NEWTON);
  Vexp->qt_init   = FALSE;
  Vexp->clp_init  = FALSE;
  Vexp->clm_init  = FALSE;
  Vexp->glue_de   = TRUE;
  Vexp->y_qt      = N_VNew_Serial (2);
  Vexp->ydot_qt   = N_VNew_Serial (2);
  Vexp->y_cl      = N_VNew_Serial (2);
  Vexp->RH_lp     = 0.0;
  Vexp->alpha_b   = 0.0;
  Vexp->a_0de     = 0.0;
  Vexp->a_0c      = 0.0;
  Vexp->a_0e      = 0.0;
  Vexp->qc        = 0.0;
  Vexp->qe        = 0.0;
  Vexp->Ec        = 0.0;
  Vexp->Ee        = 0.0;
  Vexp->alpha_qc  = 0.0;
  Vexp->alpha_qe  = 0.0;
  Vexp->c1c       = 0.0;
  Vexp->c1e       = 0.0;
  Vexp->c2c       = 0.0;
  Vexp->c2e       = 0.0;
  Vexp->evol_c    = g_array_sized_new (TRUE, TRUE, sizeof (NcHICosmoVexpState), 1000);
  Vexp->evol_e    = g_array_sized_new (TRUE, TRUE, sizeof (NcHICosmoVexpState), 1000);
  Vexp->E_alpha_s = ncm_spline_cubic_notaknot_new ();
}

static void
_nc_hicosmo_Vexp_dispose (GObject *object)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (object);

  ncm_spline_clear (&Vexp->E_alpha_s);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_Vexp_parent_class)->dispose (object);
}


static void
_nc_hicosmo_Vexp_finalize (GObject *object)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (object);
  
  if (Vexp->cvode_qt != NULL)
  {
    CVodeFree (&Vexp->cvode_qt);
    Vexp->cvode_qt = NULL;
  }

  if (Vexp->cvode_clp != NULL)
  {
    CVodeFree (&Vexp->cvode_clp);
    Vexp->cvode_clp = NULL;
  }

  if (Vexp->cvode_clm != NULL)
  {
    CVodeFree (&Vexp->cvode_clm);
    Vexp->cvode_clm = NULL;
  }

  if (Vexp->y_qt != NULL)
  {
    N_VDestroy (Vexp->y_qt);
    Vexp->y_qt = NULL;
  }

  if (Vexp->ydot_qt != NULL)
  {
    N_VDestroy (Vexp->ydot_qt);
    Vexp->ydot_qt = NULL;
  }

  if (Vexp->y_cl != NULL)
  {
    N_VDestroy (Vexp->y_cl);
    Vexp->y_cl = NULL;
  }

  g_clear_pointer (&Vexp->evol_c, g_array_unref);
  g_clear_pointer (&Vexp->evol_e, g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_Vexp_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_Vexp_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Vexp_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Vexp_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Vexp_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_Vexp_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_Vexp_bgp_cs2 (NcHICosmo *cosmo, gdouble z);

static void
nc_hicosmo_Vexp_class_init (NcHICosmoVexpClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->dispose      = &_nc_hicosmo_Vexp_dispose;
  object_class->finalize     = &_nc_hicosmo_Vexp_finalize;

  ncm_model_class_set_name_nick (model_class, "V_\\exp", "Vexp");
  ncm_model_class_add_params (model_class, NC_HICOSMO_VEXP_SPARAM_LEN, 0, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_c0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_OMEGA_C, "\\Omega_{c0}", "Omegac",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_OMEGA_C,
                               NCM_PARAM_TYPE_FREE);
  /* Set Omega_L0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_OMEGA_L, "\\Omega_{\\Lambda0}", "OmegaL",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_OMEGA_L,
                               NCM_PARAM_TYPE_FREE);
  /* Set sigmaphi param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_SIGMA_PHI, "\\sigma_{\\phi}", "sigmaphi",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_SIGMA_PHI,
                               NCM_PARAM_TYPE_FREE);
  /* Set cphi param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_D_PHI, "d_\\phi", "dphi",
                               -10.0, 10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_D_PHI,
                               NCM_PARAM_TYPE_FIXED);
  /* Set alpha_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_ALPHA_B, "\\alpha_b", "alphab",
                              0.0,  G_MAXDOUBLE, 1.0e-5,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_ALPHA_B,
                              NCM_PARAM_TYPE_FIXED);
  /* Set xb param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_X_B, "x_b", "xb",
                               1.0e10,  1.0e60, 1.0e20,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_X_B,
                               NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_Vexp_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_Vexp_E2);
  nc_hicosmo_set_Omega_t0_impl  (parent_class, &_nc_hicosmo_Vexp_Omega_t0);

  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_Vexp_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_Vexp_d2E2_dz2);

  nc_hicosmo_set_bgp_cs2_impl   (parent_class, &_nc_hicosmo_Vexp_bgp_cs2);
}

#define VECTOR    (NCM_MODEL (cosmo)->params)
#define MACRO_H0  (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_H0))
#define OMEGA_C   (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_OMEGA_C))
#define OMEGA_L   (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_OMEGA_L))
#define SIGMA_PHI (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_SIGMA_PHI))
#define D_PHI     (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_D_PHI))
#define ALPHA_B   (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_ALPHA_B))
#define X_B       (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_X_B))

#define LAMBDAp (1.0 + 1.0 / M_SQRT2)
#define LAMBDAm (1.0 - 1.0 / M_SQRT2)

static gdouble 
_nc_hicosmo_Vexp_q (const gdouble epsilon, const gint cl_b)
{
  if (cl_b > 0)
    return fabs (epsilon) / (LAMBDAm - fabs (epsilon));
  else
    return fabs (epsilon) / (LAMBDAp - fabs (epsilon));
}

static gdouble 
_nc_hicosmo_Vexp_x_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return 1.0 - LAMBDAm * q / (1.0 + q);
  else
    return -1.0 + LAMBDAp * q / (1.0 + q);  
}

static gdouble 
_nc_hicosmo_Vexp_1mx_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return LAMBDAm * q / (1.0 + q);
  else
    return 2.0 - LAMBDAp * q / (1.0 + q);
}

static gdouble 
_nc_hicosmo_Vexp_1px_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return 2.0 - LAMBDAm * q / (1.0 + q);
  else
    return LAMBDAp * q / (1.0 + q);
}

static gdouble 
_nc_hicosmo_Vexp_1sqrt2_x_q (const gdouble q, const gint cl_b)
{
  if (cl_b > 0)
    return - LAMBDAm * (1.0 / (1.0 + q));
  else
    return + LAMBDAp * (1.0 / (1.0 + q));

}

static gdouble
_nc_hicosmo_Vexp_ac_a0c_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  const gdouble onemx      = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx      = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);
  
  return cbrt (Vexp->c2c * pow (onemx, LAMBDAp) * pow (onepx, LAMBDAm) / gsl_pow_2 (onesqrt2mx));
}

static gdouble
_nc_hicosmo_Vexp_ae_a0e_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  const gdouble onemx      = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx      = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);
  
  return cbrt (Vexp->c2e * pow (onemx, LAMBDAp) * pow (onepx, LAMBDAm) / gsl_pow_2 (onesqrt2mx));
}

static gdouble
_nc_hicosmo_Vexp_Hc_H0_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  const gdouble onemx      = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx      = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);
  
  return - sqrt (Vexp->c1c) * pow (onemx, - LAMBDAp) * pow (onepx, - LAMBDAm) * fabs (onesqrt2mx) / Vexp->c2c;
}

static gdouble
_nc_hicosmo_Vexp_He_H0_q (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  const gdouble onemx      = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx      = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);
  
  return + sqrt (Vexp->c1e) * pow (onemx, - LAMBDAp) * pow (onepx, - LAMBDAm) * fabs (onesqrt2mx) / Vexp->c2e;
}

static void 
_nc_hicosmo_Vexp_dalpha_dtau_dphi_dtau (NcHICosmoVexp *Vexp, const gdouble t, const gdouble alpha, const gdouble phi, gdouble *dalpha, gdouble *dphi)
{  
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
  
	const gdouble sigma        = SIGMA_PHI;
	const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num     = - alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
	const gdouble dalpha_num   = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
	const gdouble den          = 2.0 * (cos_cosh + 1.0);
  
  dphi[0]   = dphi_num / den;
	dalpha[0] = dalpha_num / den;
  
	return;
}

static void 
_nc_hicosmo_Vexp_H_x (NcHICosmoVexp *Vexp, const gdouble t, const gdouble alpha, const gdouble phi, gdouble *H, gdouble *x)
{  
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
  
	const gdouble sigma        = SIGMA_PHI;
	const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num     = - alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
	const gdouble dalpha_num   = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
	const gdouble den          = 2.0 * (cos_cosh + 1.0);
  
  H[0] = dalpha_num / (den * exp (3.0 * alpha));
	x[0] = dphi_num / dalpha_num;

	return;
}

static gdouble 
_nc_hicosmo_Vexp_dphi_dalpha_m1_tgh (NcHICosmoVexp *Vexp, const gdouble t, const gdouble alpha, const gdouble phi)
{  
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
  
	const gdouble sigma        = SIGMA_PHI;
	const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble num          = (- alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh) * tanh_arg_hyp - phi * sigma2 * sin_cosh;
	const gdouble dem          = (phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp) * tanh_arg_hyp;
  
	return num / dem;
}

static gdouble 
_nc_hicosmo_Vexp_a_0_by_ea (NcHICosmoVexp *Vexp, const gdouble t, const gdouble epsilon, const gdouble a, const gint cl_b)
{  
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
  const gdouble d     = D_PHI;
  const gdouble Oc    = OMEGA_C;
  const gdouble LpLm  = pow (LAMBDAp, LAMBDAm);
  const gdouble LmLp  = pow (LAMBDAm, LAMBDAp);

  if (cl_b > 0)
  {
    const gdouble twoLm = pow (2.0, LAMBDAm);
    const gdouble a_0   = (pow (fabs (epsilon), LAMBDAp / 3.0) / a) * cbrt (gsl_pow_2 (d * Vexp->RH_lp) * twoLm / (LpLm * LmLp * Oc));
    return a_0;
  }
  else
  {
    if (Vexp->glue_de)
    {
      const gdouble twoLp = pow (2.0, LAMBDAp);
      const gdouble a_0   = (a / pow (fabs (epsilon), LAMBDAm / 3.0)) * cbrt (2.0 * gsl_pow_2 (LAMBDAp) / twoLp);
      return a_0;
    }
    else
    {
      const gdouble twoLp = pow (2.0, LAMBDAp);
      const gdouble a_0   = (pow (fabs (epsilon), LAMBDAm / 3.0) / a) * cbrt (gsl_pow_2 (d * Vexp->RH_lp) * twoLp / (LpLm * LmLp * Oc));
      return a_0;
    }
  }
}

static gint 
_nc_hicosmo_Vexp_qt_f (realtype tau, N_Vector y_qt, N_Vector ydot_qt, void *f_data)
{
	NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (f_data);
	const gdouble alpha = NV_Ith_S (y_qt, 0) + Vexp->alpha_b;
	const gdouble phi	  = NV_Ith_S (y_qt, 1);
  gdouble dalpha, dphi;

  _nc_hicosmo_Vexp_dalpha_dtau_dphi_dtau (Vexp, tau, alpha, phi, &dalpha, &dphi);
  
	NV_Ith_S (ydot_qt, 0) = dalpha;
	NV_Ith_S (ydot_qt, 1) = dphi;
	
	return 0;
}

static gint
_nc_hicosmo_Vexp_qt_J (_NCM_SUNDIALS_INT_TYPE N, realtype tau, N_Vector y_qt, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (jac_data);
  NcHICosmo *cosmo    = NC_HICOSMO (Vexp);

	const gdouble alpha = NV_Ith_S (y_qt, 0) + Vexp->alpha_b;
	const gdouble phi	  = NV_Ith_S (y_qt, 1);

	const gdouble sigma        = SIGMA_PHI;
	const gdouble sigma2       = sigma * sigma;
  const gdouble d            = D_PHI;
  const gdouble arg_tri      = 2.0 * alpha * d;
  const gdouble arg_hyp      = sigma2 * alpha * phi;
  const gdouble cosh_arg_hyp = cosh (arg_hyp);
  const gdouble tanh_arg_hyp = tanh (arg_hyp);
  const gdouble sin_arg_tri  = sin (arg_tri);
  const gdouble cos_arg_tri  = cos (arg_tri);
  const gdouble sin_cosh     = sin_arg_tri / cosh_arg_hyp;
  const gdouble cos_cosh     = cos_arg_tri / cosh_arg_hyp;

  const gdouble dphi_num     = - alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
	const gdouble dalpha_num   = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
	const gdouble den          = 2.0 * (cos_cosh + 1.0);

  const gdouble dalpha       = dalpha_num / den;
  const gdouble dphi         = dphi_num / den;
  
  DENSE_ELEM (J, 0, 0) = d * sigma2 * phi - dalpha * (- 2.0 * d * sin_cosh + sigma2 * phi * tanh_arg_hyp) / (cos_cosh + 1.0);
  DENSE_ELEM (J, 0, 1) = ((sigma2 * sin_cosh + 2.0 * d * alpha * sigma2) * cos_cosh + sigma2 * sin_cosh - phi * sigma2 * sin_cosh * sigma2 * alpha * tanh_arg_hyp + 2.0 * d * alpha * sigma2 / gsl_pow_2 (cosh_arg_hyp)) / (den * (cos_cosh + 1.0));
    /*(sigma2 * sin_cosh + 2.0 * d * alpha * sigma2) / den - dalpha * sigma2 * alpha * tanh_arg_hyp / (cos_cosh + 1.0);*/

  DENSE_ELEM (J, 1, 0) = - dphi * (-2.0 * d * sin_cosh + sigma2 * phi * tanh_arg_hyp) / (cos_cosh + 1.0) + (-sigma2 * sin_cosh - 4.0 * d * d * sin_cosh - 2.0 * d * alpha * sigma2 * cos_cosh + 2.0 * d * sigma2 * phi * tanh_arg_hyp) / den;
  DENSE_ELEM (J, 1, 1) = sigma2 * sigma2 * alpha * alpha * sin_cosh * tanh_arg_hyp / (2.0 * gsl_pow_2 (cos_cosh + 1.0));

/*  printf ("JAC % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", alpha, phi, dalpha, dphi, DENSE_ELEM (J, 0, 0), DENSE_ELEM (J, 0, 1), DENSE_ELEM (J, 1, 0), DENSE_ELEM (J, 1, 1));*/
  
  return 0;
}

static gint 
_nc_hicosmo_Vexp_qt_root (realtype tau, N_Vector y_qt, realtype *gout, void *user_data)
{
	NcHICosmoVexp *Vexp   = NC_HICOSMO_VEXP (user_data);
  NcHICosmo *cosmo      = NC_HICOSMO (Vexp);
	const gdouble alpha   = NV_Ith_S (y_qt, 0) + Vexp->alpha_b;
	const gdouble phi	    = NV_Ith_S (y_qt, 1);
  const gdouble epsilon = _nc_hicosmo_Vexp_dphi_dalpha_m1_tgh (Vexp, tau, alpha, phi);
  const gdouble x_b     = X_B;
  const gdouble a_b     = exp (Vexp->alpha_b);
  const gdouble a_0p    = x_b * a_b;
  const gdouble a_0m    = Vexp->glue_de ? Vexp->a_0de : a_0p;

  gdouble H_lp, x;

  _nc_hicosmo_Vexp_H_x (Vexp, tau, alpha, phi, &H_lp, &x);
  
	gout[0] = log (fabs (epsilon * x) * 1.0e5);
  gout[1] = log (_nc_hicosmo_Vexp_a_0_by_ea (Vexp, tau, epsilon, exp (alpha), +1) / a_0p);
  gout[2] = log (_nc_hicosmo_Vexp_a_0_by_ea (Vexp, tau, epsilon, exp (alpha), -1) / a_0m);

/*printf ("% 21.15g % 21.15g\n", gout[0], gout[1]);*/
	return 0;
} 

static gint 
_nc_hicosmo_Vexp_clp_f (realtype tau, N_Vector y_clp, N_Vector ydot_clp, void *f_data)
{
	const gdouble q = NV_Ith_S (y_clp, 0);
  const gdouble E = NV_Ith_S (y_clp, 1);
  const gdouble x = _nc_hicosmo_Vexp_x_q (q, +1);
  
	NV_Ith_S (ydot_clp, 0) = + 6.0 * q * (LAMBDAm + 0.25 * q) / (1.0 + q);
  NV_Ith_S (ydot_clp, 1) = - 3.0 * x * x * E;
	
	return 0;
}

static gint 
_nc_hicosmo_Vexp_clm_f (realtype tau, N_Vector y_clm, N_Vector ydot_clm, void *f_data)
{
	const gdouble q = NV_Ith_S (y_clm, 0);
  const gdouble E = NV_Ith_S (y_clm, 1);
  const gdouble x = _nc_hicosmo_Vexp_x_q (q, -1);
  
	NV_Ith_S (ydot_clm, 0) = + 6.0 * q * (LAMBDAp + 0.25 * q) / (1.0 + q);
  NV_Ith_S (ydot_clm, 1) = - 3.0 * x * x * E;
	
	return 0;
}

#define ABSTOL (1.0e-14)

static void
_nc_hicosmo_Vexp_init_qt (NcHICosmoVexp *Vexp)
{
  gint flag;

  NV_Ith_S (Vexp->y_qt, 0) = 0.0;
  NV_Ith_S (Vexp->y_qt, 1) = 0.0;
  
  if (!Vexp->qt_init)
  {
    flag = CVodeInit (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_f, 0.0, Vexp->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

		flag = CVodeSetMaxStep (Vexp->cvode_qt, G_MAXUINT32);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (Vexp->cvode_qt, ABSTOL, 1.0e-150);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (Vexp->cvode_qt, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVDense (Vexp->cvode_qt, 2);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVodeRootInit (Vexp->cvode_qt, 3, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

    flag = CVDlsSetDenseJacFn (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );    
  }
  else
  {
    flag = CVodeReInit (Vexp->cvode_qt, 0.0, Vexp->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeRootInit (Vexp->cvode_qt, 3, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );    

    flag = CVDlsSetDenseJacFn (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );    
  }
}

static void
_nc_hicosmo_Vexp_init_clp (NcHICosmoVexp *Vexp, gdouble alpha_q, gdouble q0, gdouble E0)
{
  gint flag;

  NV_Ith_S (Vexp->y_cl, 0) = q0;
  NV_Ith_S (Vexp->y_cl, 1) = E0;
  
  if (!Vexp->clp_init)
  {
    flag = CVodeInit (Vexp->cvode_clp, &_nc_hicosmo_Vexp_clp_f, alpha_q, Vexp->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

		flag = CVodeSetMaxStep (Vexp->cvode_clp, G_MAXUINT32);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (Vexp->cvode_clp, ABSTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (Vexp->cvode_clp, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVDense (Vexp->cvode_clp, 1);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );
  }
  else
  {
    flag = CVodeReInit (Vexp->cvode_clp, alpha_q, Vexp->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }
}

static void
_nc_hicosmo_Vexp_init_clm (NcHICosmoVexp *Vexp, gdouble alpha_q, gdouble q0, gdouble E0)
{
  gint flag;

  NV_Ith_S (Vexp->y_cl, 0) = q0;
  NV_Ith_S (Vexp->y_cl, 1) = E0;
  
  if (!Vexp->clm_init)
  {
    flag = CVodeInit (Vexp->cvode_clm, &_nc_hicosmo_Vexp_clm_f, alpha_q, Vexp->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

		flag = CVodeSetMaxStep (Vexp->cvode_clm, G_MAXUINT32);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (Vexp->cvode_clm, ABSTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (Vexp->cvode_clm, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVDense (Vexp->cvode_clm, 1);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );
  }
  else
  {
    flag = CVodeReInit (Vexp->cvode_clm, alpha_q, Vexp->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }
}

static void
_nc_hicosmo_Vexp_evolve_qt (NcHICosmoVexp *Vexp, gdouble tau_f)
{
  NcHICosmo *cosmo      = NC_HICOSMO (Vexp);
  GArray *earray        = tau_f > 0.0 ? Vexp->evol_e : Vexp->evol_c;
  gint flag;

  _nc_hicosmo_Vexp_init_qt (Vexp);

  flag = CVodeSetStopTime (Vexp->cvode_qt, tau_f);
  NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

  while (TRUE)
  {
    gdouble tau;
    gboolean root_found;

    flag = CVode (Vexp->cvode_qt, tau_f, Vexp->y_qt, &tau, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    root_found = (flag == CV_ROOT_RETURN);

    {
      const gdouble ralpha = NV_Ith_S (Vexp->y_qt, 0);
      const gdouble alpha  = ralpha + Vexp->alpha_b;
      const gdouble phi    = NV_Ith_S (Vexp->y_qt, 1);
      gdouble H_lp, x;

      _nc_hicosmo_Vexp_H_x (Vexp, tau, alpha, phi, &H_lp, &x);

      {
        const gdouble E            = H_lp * Vexp->RH_lp;
        const NcHICosmoVexpState s = {ralpha, E};

        g_array_append_val (earray, s);
        
        if (TRUE)
        {
          const gdouble a3 = exp (3.0 * alpha);
          const gdouble ep = _nc_hicosmo_Vexp_dphi_dalpha_m1_tgh (Vexp, tau, alpha, phi);
          
          printf ("% 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %s\n", 
                  tau, 
                  alpha, 
                  phi,
                  E,
                  x,
                  a3 * H_lp,
                  ep,
                  cbrt ((a3 / pow (fabs (ep), LAMBDAm))),
                  cbrt ((pow (fabs (ep), LAMBDAp) / a3) * gsl_pow_2 (Vexp->RH_lp)),
                  x * SIGMA_PHI * sqrt (2.0 * Vexp->alpha_b * fabs (ralpha)),
                  root_found ? "ROOT" : "NOT-FOUND"
                  );
        }
      }
    }

    if (root_found)
    {
      const gdouble alpha = NV_Ith_S (Vexp->y_qt, 0) + Vexp->alpha_b;
      const gdouble phi   = NV_Ith_S (Vexp->y_qt, 1);
      gint rootsfound[3]  = {0, 0, 0};
      gdouble H_lp, x;

      _nc_hicosmo_Vexp_H_x (Vexp, tau, alpha, phi, &H_lp, &x);

      flag = CVodeGetRootInfo (Vexp->cvode_qt, rootsfound);
      NCM_CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, );

      if (((x > 0.0) && (rootsfound[1] != 0)) || ((x < 0.0) && (rootsfound[2] != 0)))
      {
        const gint cl_b       = GSL_SIGN (x);
        const gdouble epsilon = _nc_hicosmo_Vexp_dphi_dalpha_m1_tgh (Vexp, tau, alpha, phi);
        const gdouble a_q     = exp (alpha);
        const gdouble LpLm    = pow (LAMBDAp, LAMBDAm);
        const gdouble LmLp    = pow (LAMBDAm, LAMBDAp);
        const gdouble d       = D_PHI;

        const gdouble a_0     = _nc_hicosmo_Vexp_a_0_by_ea (Vexp, tau, epsilon, a_q, cl_b);
        const gdouble alpha_0 = log (a_0);
        const gdouble alpha_q = alpha;
        const gdouble q       = _nc_hicosmo_Vexp_q (epsilon, cl_b);
        const gdouble E       = H_lp * Vexp->RH_lp;
        const gdouble c1      = (cl_b > 0) ? gsl_pow_2 (LAMBDAm * d * Vexp->RH_lp) / gsl_pow_6 (a_0) : gsl_pow_2 (LAMBDAp * d * Vexp->RH_lp) / gsl_pow_6 (a_0);
        const gdouble c2      = c1 / (OMEGA_C * LpLm * LmLp);

        if (tau_f < 0.0)
        {
          Vexp->cl_bc    = cl_b;
          Vexp->a_0c     = a_0;
          Vexp->alpha_0c = alpha_0;
          Vexp->alpha_qc = alpha_q;
          Vexp->qc       = q;
          Vexp->Ec       = E;
          Vexp->c1c      = c1;
          Vexp->c2c      = c2;

          if ((cl_b < 0) && Vexp->glue_de)
          {
            Vexp->c1c      = 0.5 * OMEGA_L;
            Vexp->c2c      = 0.5;
/*            
            Vexp->alpha_0c = log (2.0 * gsl_pow_2 (LAMBDAp * d * Vexp->RH_lp) / OMEGA_L) / 6.0;
            Vexp->a_0c     = exp (Vexp->alpha_0c);
*/
          }
        }
        else
        {
          Vexp->cl_be    = cl_b;
          Vexp->a_0e     = a_0;
          Vexp->alpha_0e = alpha_0;
          Vexp->alpha_qe = alpha_q;
          Vexp->qe       = q;
          Vexp->Ee       = E;
          Vexp->c1e      = c1;
          Vexp->c2e      = c2;

          if ((cl_b < 0) && Vexp->glue_de)
          {
            Vexp->c1e      = 0.5 * OMEGA_L;
            Vexp->c2e      = 0.5;
/*
            Vexp->alpha_0e = log (2.0 * gsl_pow_2 (LAMBDAp * d * Vexp->RH_lp) / OMEGA_L) / 6.0;
            Vexp->a_0e     = exp (Vexp->alpha_0c);
*/          
          }
        }
        break;
      }
    }

    if (tau == tau_f)
      break;
  }
}

static void
_nc_hicosmo_Vexp_evolve_cl_c (NcHICosmoVexp *Vexp)
{
  const gdouble alpha_f = 100.0 + Vexp->alpha_0c;
  gpointer cvode;
  gint flag;

  if (Vexp->cl_bc > 0)
  {
    _nc_hicosmo_Vexp_init_clp (Vexp, Vexp->alpha_qc, Vexp->qc, Vexp->Ec);
    cvode = Vexp->cvode_clp;
  }
  else
  {
    _nc_hicosmo_Vexp_init_clm (Vexp, Vexp->alpha_qc, Vexp->qc, Vexp->Ec);
    cvode = Vexp->cvode_clm;
  }

  flag = CVodeSetStopTime (cvode, alpha_f);
  NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

  while (TRUE)
  {
    gdouble alpha;

    flag = CVode (cvode, alpha_f, Vexp->y_cl, &alpha, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    {
      const gdouble q    = NV_Ith_S (Vexp->y_cl, 0);
      const gdouble En   = NV_Ith_S (Vexp->y_cl, 1);
      const NcHICosmoVexpState s = {alpha - Vexp->alpha_b, En};

      g_array_append_val (Vexp->evol_c, s);
      if (FALSE)
      {
        const gdouble a_a0 = _nc_hicosmo_Vexp_ac_a0c_q (Vexp, q, Vexp->cl_bc);
        const gdouble E    = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, Vexp->cl_bc);
        printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                alpha, q, exp (alpha - Vexp->alpha_0c), a_a0, exp (alpha - Vexp->alpha_0c) / a_a0 - 1.0, E, En, E / En - 1.0);
      }
    }

    if (alpha == alpha_f)
      break;      
  }
}

static void
_nc_hicosmo_Vexp_evolve_cl_e (NcHICosmoVexp *Vexp)
{
  const gdouble alpha_f = 100.0 + Vexp->alpha_0e;
  gpointer cvode;
  gint flag;

  if (Vexp->cl_be > 0)
  {
    _nc_hicosmo_Vexp_init_clp (Vexp, Vexp->alpha_qe, Vexp->qe, Vexp->Ee);
    cvode = Vexp->cvode_clp;
  }
  else
  {
    _nc_hicosmo_Vexp_init_clm (Vexp, Vexp->alpha_qe, Vexp->qe, Vexp->Ee);
    cvode = Vexp->cvode_clm;
  }

  flag = CVodeSetStopTime (cvode, alpha_f);
  NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

  while (TRUE)
  {
    gdouble alpha;

    flag = CVode (cvode, alpha_f, Vexp->y_cl, &alpha, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    {
      const gdouble q    = NV_Ith_S (Vexp->y_cl, 0);
      const gdouble En   = NV_Ith_S (Vexp->y_cl, 1);
      const NcHICosmoVexpState s = {alpha - Vexp->alpha_b, En};

      g_array_append_val (Vexp->evol_e, s);
      if (FALSE)
      {
        const gdouble a_a0 = _nc_hicosmo_Vexp_ae_a0e_q (Vexp, q, Vexp->cl_be);
        const gdouble E    = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, Vexp->cl_be);
        printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                alpha, q, exp (alpha - Vexp->alpha_0e), a_a0, exp (alpha - Vexp->alpha_0e) / a_a0 - 1.0, E, En, E / En - 1.0);
      }
    }

    if (alpha == alpha_f)
      break;      
  }
}

static void
_nc_hicosmo_Vexp_prepare (NcHICosmoVexp *Vexp)
{
  NcHICosmo *cosmo      = NC_HICOSMO (Vexp);  
  const gdouble tau_f   = 1.0e10;
  
  if (!ncm_model_state_is_update (NCM_MODEL (Vexp)))
  {
    const gdouble d = D_PHI;
    Vexp->RH_lp     = nc_hicosmo_RH_planck (cosmo);
    Vexp->alpha_b   = ALPHA_B;
    Vexp->a_0de     = pow (2.0 * gsl_pow_2 (d * LAMBDAp * Vexp->RH_lp) / OMEGA_L, 1.0 / 6.0);

    g_array_set_size (Vexp->evol_c, 0);
    g_array_set_size (Vexp->evol_e, 0);
    
    _nc_hicosmo_Vexp_evolve_qt (Vexp, -tau_f);
    _nc_hicosmo_Vexp_evolve_qt (Vexp, +tau_f);
    _nc_hicosmo_Vexp_evolve_cl_c (Vexp);
    _nc_hicosmo_Vexp_evolve_cl_e (Vexp);

    printf ("# Contracting branch %u; Expanding branch %u;\n", Vexp->evol_c->len, Vexp->evol_e->len);
    {
      const guint len    = Vexp->evol_c->len + Vexp->evol_e->len;
      NcmVector *alpha_v = ncm_vector_new (len);
      NcmVector *E_v     = ncm_vector_new (len);
      glong i, j = 0;

      for (i = Vexp->evol_c->len - 1; i >= 0; i--)
      {
        const NcHICosmoVexpState s = g_array_index (Vexp->evol_c, NcHICosmoVexpState, i);
        ncm_vector_fast_set (alpha_v, j, -s.alpha);
        ncm_vector_fast_set (E_v, j,      s.E);
        printf ("%6ld %6ld % 21.15e % 21.15e % 21.15e\n", i, j, -s.alpha, s.E, s.alpha - Vexp->alpha_0c);
        j++;
      }

      for (i = 0; i < Vexp->evol_e->len; i++)
      {
        const NcHICosmoVexpState s = g_array_index (Vexp->evol_e, NcHICosmoVexpState, i);
        ncm_vector_fast_set (alpha_v, j, s.alpha);
        ncm_vector_fast_set (E_v,     j, s.E);
        printf ("%6ld %6ld % 21.15e % 21.15e % 21.15e\n", i, j, s.alpha, s.E, s.alpha - Vexp->alpha_0e);
        j++;
      }

      ncm_spline_set (Vexp->E_alpha_s, alpha_v, E_v, TRUE);
      
    }
    
    ncm_model_state_set_update (NCM_MODEL (Vexp));
  }
  else
    return;
}

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_Vexp_E2 (NcHICosmo *cosmo, gdouble z)
{
  _nc_hicosmo_Vexp_prepare (NC_HICOSMO_VEXP (cosmo)); 
  return 1.0;
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_Vexp_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  _nc_hicosmo_Vexp_prepare (NC_HICOSMO_VEXP (cosmo)); 
  return 0.0;
}

static gdouble
_nc_hicosmo_Vexp_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  _nc_hicosmo_Vexp_prepare (NC_HICOSMO_VEXP (cosmo)); 
  return 0.0;
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_Vexp_H0 (NcHICosmo *cosmo) { return MACRO_H0; }
static gdouble _nc_hicosmo_Vexp_Omega_t0 (NcHICosmo *cosmo) { return OMEGA_C; }
static gdouble
_nc_hicosmo_Vexp_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
{
  _nc_hicosmo_Vexp_prepare (NC_HICOSMO_VEXP (cosmo)); 
  return 0.0;
}

/**
 * nc_hicosmo_Vexp_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoVexp *
nc_hicosmo_Vexp_new (void)
{
  NcHICosmoVexp *Vexp = g_object_new (NC_TYPE_HICOSMO_VEXP, NULL);
  return Vexp;
}
