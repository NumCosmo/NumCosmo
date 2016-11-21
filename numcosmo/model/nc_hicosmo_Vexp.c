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
#include "perturbations/nc_hipert_adiab.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

static void nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface);

G_DEFINE_TYPE_WITH_CODE (NcHICosmoVexp, nc_hicosmo_Vexp, NC_TYPE_HICOSMO,
                         G_IMPLEMENT_INTERFACE (NC_TYPE_HIPERT_IADIAB,
                                                nc_hipert_iadiab_interface_init)
                         );

enum {
  PROP_0,
  PROP_SIZE,
};

typedef struct _NcHICosmoVexpState
{
  gdouble tau;
  gdouble N;
  gdouble xtau;
  gdouble dlnx_dalpha;
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
  Vexp->alpha_0c  = 0.0;
  Vexp->alpha_0e  = 0.0;
  Vexp->tau_x0    = 0.0;
  Vexp->c1c       = 0.0;
  Vexp->c1e       = 0.0;
  Vexp->c2c       = 0.0;
  Vexp->c2e       = 0.0;
  Vexp->evol_c    = g_array_sized_new (TRUE, TRUE, sizeof (NcHICosmoVexpState), 1000);
  Vexp->evol_e    = g_array_sized_new (TRUE, TRUE, sizeof (NcHICosmoVexpState), 1000);

  Vexp->xtau_s           = ncm_spline_cubic_notaknot_new ();
  Vexp->lnN_s            = ncm_spline_cubic_notaknot_new ();
  Vexp->tau_dlnx_dtau_s  = ncm_spline_cubic_notaknot_new ();
  Vexp->E2_s             = ncm_spline_cubic_notaknot_new ();
}

static void
_nc_hicosmo_Vexp_dispose (GObject *object)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (object);

  ncm_spline_clear (&Vexp->xtau_s);
  ncm_spline_clear (&Vexp->lnN_s);
  ncm_spline_clear (&Vexp->tau_dlnx_dtau_s);
  ncm_spline_clear (&Vexp->E2_s);
  
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
}

void _nc_hicosmo_Vexp_init_x_alpha_recur (NcmVector *alpha_x_series, NcmVector *x_alpha_series, const gdouble pf, guint nmax);

static gdouble 
_nc_hicosmo_Vexp_init_x_alpha_series (const gdouble dalpha)
{
  const gdouble an[25] = {
    +2.1213203435596425732,
    -3.1819805153394638598,
    +0.0,
    +14.318912319027587369,
    -38.661063261374485897,
    +30.069715869957933475,
    +133.47271840236429655,
    -562.31135761409854934,
    +852.95970820407459509,
    +879.60822680923629430,
    -7833.1017094058379035,
    +17969.469996896215007,
    -5132.5309707529267568,
    -97212.641401182025488,
    +325940.66821809623956,
    -370774.67123243092981,
    -949452.26167748074539,
    +5.2556176004654939455e6,
    -9.9203647608405914131e6,
    -3.5322017625266955972e6,
    +7.4749084074478288480e7,
    -2.0578628366558384105e8,
    +1.4591930291931796152e8,
    +8.8839574285216097773e8,
    -3.6743882793601701195e9
  };
  register gint i;
  register gdouble res = an[24];
  
  for (i = 23; i >= 0; i--)
  {
#ifdef HAVE_FMA
    res = fma (dalpha, res, an[i]);
#else
    res = dalpha * res + an[i];
#endif
  }

  return dalpha * res;
}

static gdouble 
_nc_hicosmo_Vexp_init_dx_alpha_series (const gdouble dalpha)
{
  const gdouble an[25] = {
    +2.1213203435596425732   *  1.0,
    -3.1819805153394638598   *  2.0,
    +0.0                     *  3.0,
    +14.318912319027587369   *  4.0,
    -38.661063261374485897   *  5.0,
    +30.069715869957933475   *  6.0,
    +133.47271840236429655   *  7.0,
    -562.31135761409854934   *  8.0,
    +852.95970820407459509   *  9.0,
    +879.60822680923629430   * 10.0,
    -7833.1017094058379035   * 11.0,
    +17969.469996896215007   * 12.0,
    -5132.5309707529267568   * 13.0,
    -97212.641401182025488   * 14.0,
    +325940.66821809623956   * 15.0,
    -370774.67123243092981   * 16.0,
    -949452.26167748074539   * 17.0,
    +5.2556176004654939455e6 * 18.0,
    -9.9203647608405914131e6 * 19.0,
    -3.5322017625266955972e6 * 20.0,
    +7.4749084074478288480e7 * 21.0,
    -2.0578628366558384105e8 * 22.0,
    +1.4591930291931796152e8 * 23.0,
    +8.8839574285216097773e8 * 24.0,
    -3.6743882793601701195e9 * 25.0
  };
  register gint i;
  register gdouble res = an[24];
  
  for (i = 23; i >= 0; i--)
  {
#ifdef HAVE_FMA
    res = fma (dalpha, res, an[i]);
#else
    res = dalpha * res + an[i];
#endif
  }

  return res;
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

static gdouble _nc_hicosmo_Vexp_zeta_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_zeta_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static gdouble _nc_hicosmo_Vexp_zeta_eval_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k);
static void _nc_hicosmo_Vexp_zeta_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu);

static guint _nc_hicosmo_Vexp_zeta_nsing (NcHIPertIAdiab *iad, const gdouble k);
static void _nc_hicosmo_Vexp_zeta_get_sing_info (NcHIPertIAdiab *iad, const gdouble k, const guint sing, gdouble *ts, NcmHOAASingType *st);
static gdouble _nc_hicosmo_Vexp_zeta_eval_sing_mnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, guint sing);
static gdouble _nc_hicosmo_Vexp_zeta_eval_sing_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, guint sing);
static void _nc_hicosmo_Vexp_zeta_eval_sing_system (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, guint sing, gdouble *nu, gdouble *dlnmnu);

static void
nc_hipert_iadiab_interface_init (NcHIPertIAdiabInterface *iface)
{
  iface->eval_mnu         = &_nc_hicosmo_Vexp_zeta_eval_mnu;
  iface->eval_nu          = &_nc_hicosmo_Vexp_zeta_eval_nu;
  iface->eval_dlnmnu      = &_nc_hicosmo_Vexp_zeta_eval_dlnmnu;
  iface->eval_system      = &_nc_hicosmo_Vexp_zeta_eval_system;

  iface->nsing            = &_nc_hicosmo_Vexp_zeta_nsing;
  iface->get_sing_info    = &_nc_hicosmo_Vexp_zeta_get_sing_info;
  iface->eval_sing_mnu    = &_nc_hicosmo_Vexp_zeta_eval_sing_mnu;
  iface->eval_sing_dlnmnu = &_nc_hicosmo_Vexp_zeta_eval_sing_dlnmnu;
  iface->eval_sing_system = &_nc_hicosmo_Vexp_zeta_eval_sing_system;
}

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

static gdouble 
_nc_hicosmo_Vexp_cl_dlnx_dalpha (NcHICosmoVexp *Vexp, const gdouble q, const gint cl_b)
{
  const gdouble x          = _nc_hicosmo_Vexp_x_q (q, cl_b);
  const gdouble onemx      = _nc_hicosmo_Vexp_1mx_q (q, cl_b);
  const gdouble onepx      = _nc_hicosmo_Vexp_1px_q (q, cl_b);
  const gdouble onesqrt2mx = _nc_hicosmo_Vexp_1sqrt2_x_q (q, cl_b);

  return 3.0 * onepx * onemx * onesqrt2mx / x;
}

static void 
_nc_hicosmo_Vexp_dalpha_dtau_dphi_dtQ (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi, gdouble *dalpha, gdouble *dphi)
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
_nc_hicosmo_Vexp_H_x (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi, gdouble *H, gdouble *x)
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
_nc_hicosmo_Vexp_epsilon (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi)
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
  
	const gdouble x            = dphi_num / dalpha_num;
  const gdouble x2           = x * x;
  gdouble x2m1;
  
  if (x2 > 1.1)
  {
    x2m1 = x2 - 1.0;
  }
  else
  {
    const gdouble epsilon_num = (- alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh) * tanh_arg_hyp - phi * sigma2 * sin_cosh;
    const gdouble epsilon_dem = (phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp) * tanh_arg_hyp;
    const gdouble epsilon     = epsilon_num / epsilon_dem;

    x2m1 = (epsilon * epsilon + 2.0 * epsilon / tanh_arg_hyp + 1.0 / gsl_pow_2 (tanh_arg_hyp * cosh_arg_hyp));
  }

  return - ncm_util_sqrt1px_m1 (x2m1);
}

static gdouble 
_nc_hicosmo_Vexp_qt_dlnx_dalpha (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble alpha, const gdouble phi)
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
  
	const gdouble x            = dphi_num / dalpha_num;
  const gdouble x2           = x * x;
  gdouble onemx2;
  
  if (x2 > 1.1)
  {
    onemx2 = 1.0 - x2;
  }
  else
  {
    const gdouble e_num = (- alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh) * tanh_arg_hyp - phi * sigma2 * sin_cosh;
    const gdouble e_dem = (phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp) * tanh_arg_hyp;
    const gdouble e     = e_num / e_dem;

    onemx2 = - (e * e + 2.0 * e / tanh_arg_hyp + 1.0 / gsl_pow_2 (tanh_arg_hyp * cosh_arg_hyp));
  }

  {
    const gdouble alpha2 = alpha * alpha;
    const gdouble f1     = + onemx2 * sigma2 * (alpha + phi / x);
    const gdouble f2     = - 2.0 * d * sigma2 * alpha * cos_cosh * onemx2 / dphi_num;
    const gdouble f3     = - sin_cosh * (4.0 * d * d + sigma2 * (1.0 + x2 * (1.0 + alpha2 * sigma2)) + sigma2 * sigma2 * phi * (2.0 * x * alpha + phi)) / dphi_num;
    
    return (f1 + f2 + f3);
  }
}

static gdouble 
_nc_hicosmo_Vexp_a_0_by_ea (NcHICosmoVexp *Vexp, const gdouble tQ, const gdouble epsilon, const gdouble a, const gint cl_b)
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
_nc_hicosmo_Vexp_qt_f (realtype tQ, N_Vector y_qt, N_Vector ydot_qt, void *f_data)
{
	NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (f_data);
	const gdouble alpha = NV_Ith_S (y_qt, 0) + Vexp->alpha_b;
	const gdouble phi	  = NV_Ith_S (y_qt, 1);
  gdouble dalpha, dphi;

  _nc_hicosmo_Vexp_dalpha_dtau_dphi_dtQ (Vexp, tQ, alpha, phi, &dalpha, &dphi);
  
	NV_Ith_S (ydot_qt, 0) = dalpha;
	NV_Ith_S (ydot_qt, 1) = dphi;

  /*printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", tQ, alpha, dalpha, phi, dphi);*/
  
	return 0;
}

static gint
_nc_hicosmo_Vexp_qt_J (_NCM_SUNDIALS_INT_TYPE N, realtype tQ, N_Vector y_qt, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
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
_nc_hicosmo_Vexp_qt_root (realtype tQ, N_Vector y_qt, realtype *gout, gpointer user_data)
{
	NcHICosmoVexp *Vexp   = NC_HICOSMO_VEXP (user_data);
  NcHICosmo *cosmo      = NC_HICOSMO (Vexp);
	const gdouble alpha   = NV_Ith_S (y_qt, 0) + Vexp->alpha_b;
	const gdouble phi	    = NV_Ith_S (y_qt, 1);
  const gdouble epsilon = _nc_hicosmo_Vexp_epsilon (Vexp, tQ, alpha, phi);
  const gdouble x_b     = X_B;
  const gdouble a_b     = exp (Vexp->alpha_b);
  const gdouble a_0p    = x_b * a_b;
  const gdouble a_0m    = Vexp->glue_de ? Vexp->a_0de : a_0p;

  gdouble H_lp, x;

  _nc_hicosmo_Vexp_H_x (Vexp, tQ, alpha, phi, &H_lp, &x);
  
	gout[0] = log (fabs (epsilon * x) * 1.0e4);
  gout[1] = log (_nc_hicosmo_Vexp_a_0_by_ea (Vexp, tQ, epsilon, exp (alpha), +1) / a_0p);
  gout[2] = log (_nc_hicosmo_Vexp_a_0_by_ea (Vexp, tQ, epsilon, exp (alpha), -1) / a_0m);
  gout[3] = log (fabs (H_lp * Vexp->RH_lp * x));

  /*printf ("# gout % 21.15g % 21.15g % 21.15g\n", gout[0], gout[1], gout[1]);*/
	return 0;
} 

static gint 
_nc_hicosmo_Vexp_clp_f (realtype tQ, N_Vector y_clp, N_Vector ydot_clp, void *f_data)
{
	const gdouble q = NV_Ith_S (y_clp, 0);
  const gdouble E = NV_Ith_S (y_clp, 1);
  const gdouble x = _nc_hicosmo_Vexp_x_q (q, +1);
  
	NV_Ith_S (ydot_clp, 0) = + 6.0 * q * (LAMBDAm + 0.25 * q) / (1.0 + q);
  NV_Ith_S (ydot_clp, 1) = - 3.0 * x * x * E;
	
	return 0;
}

static gint 
_nc_hicosmo_Vexp_clm_f (realtype tQ, N_Vector y_clm, N_Vector ydot_clm, void *f_data)
{
	const gdouble q = NV_Ith_S (y_clm, 0);
  const gdouble E = NV_Ith_S (y_clm, 1);
  const gdouble x = _nc_hicosmo_Vexp_x_q (q, -1);
  
	NV_Ith_S (ydot_clm, 0) = + 6.0 * q * (LAMBDAp + 0.25 * q) / (1.0 + q);
  NV_Ith_S (ydot_clm, 1) = - 3.0 * x * x * E;
	
	return 0;
}

static gint 
_nc_hicosmo_Vexp_clm_root (realtype tQ, N_Vector y_cl, realtype *gout, gpointer user_data)
{
	const gdouble q = NV_Ith_S (y_cl, 0);

	gout[0] = 1.0e4 * _nc_hicosmo_Vexp_x_q (q, -1);

  return 0;
} 

#define RELTOL (1.0e-14)

static void
_nc_hicosmo_Vexp_init_qt (NcHICosmoVexp *Vexp, const gdouble direction)
{
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
	const gdouble sigma    = SIGMA_PHI;
	const gdouble sigma2   = sigma * sigma;
	const gdouble d        = D_PHI;
  const gdouble tQ_i     = 1.0e-60 * direction;
  const gdouble tQ_i2    = tQ_i * tQ_i;
  const gdouble arg_trig = d * Vexp->alpha_b;
  gint flag;

  NV_Ith_S (Vexp->y_qt, 0) = 0.125 * tQ_i2 * sigma2 * (d * Vexp->alpha_b / gsl_pow_2 (cos (arg_trig)) + tan (arg_trig)) * (2.0 * d - Vexp->alpha_b * sigma2 *  tan (arg_trig));
  NV_Ith_S (Vexp->y_qt, 1) = 0.5 * tQ_i * (2.0 * d - Vexp->alpha_b * sigma2 *  tan (arg_trig));
  
  if (!Vexp->qt_init)
  {
    flag = CVodeInit (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_f, tQ_i, Vexp->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

		flag = CVodeSetMaxStep (Vexp->cvode_qt, G_MAXUINT32);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (Vexp->cvode_qt, RELTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (Vexp->cvode_qt, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVDense (Vexp->cvode_qt, 2);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVodeRootInit (Vexp->cvode_qt, 4, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

    flag = CVDlsSetDenseJacFn (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );    

    flag = CVodeSetInitStep (Vexp->cvode_qt, tQ_i * 1.0e-4);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
  }
  else
  {
    flag = CVodeReInit (Vexp->cvode_qt, tQ_i, Vexp->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeRootInit (Vexp->cvode_qt, 4, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );    

    flag = CVDlsSetDenseJacFn (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );    

    flag = CVodeSetInitStep (Vexp->cvode_qt, tQ_i * 1.0e-4);
    NCM_CVODE_CHECK (&flag, "CVodeSetInitStep", 1, );
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

    flag = CVodeSStolerances (Vexp->cvode_clp, RELTOL, 0.0);
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

    flag = CVodeSStolerances (Vexp->cvode_clm, RELTOL, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (Vexp->cvode_clm, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVDense (Vexp->cvode_clm, 1);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVodeRootInit (Vexp->cvode_clm, 1, &_nc_hicosmo_Vexp_clm_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );
  }
  else
  {
    flag = CVodeReInit (Vexp->cvode_clm, alpha_q, Vexp->y_cl);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  }
}

static void
_nc_hicosmo_Vexp_evolve_qt (NcHICosmoVexp *Vexp, gdouble tQ_f)
{
  NcHICosmo *cosmo     = NC_HICOSMO (Vexp);
  GArray *earray       = tQ_f > 0.0 ? Vexp->evol_e : Vexp->evol_c;
  gboolean root1_found = FALSE;
  gint flag;

  _nc_hicosmo_Vexp_init_qt (Vexp, GSL_SIGN (tQ_f));

  flag = CVodeSetStopTime (Vexp->cvode_qt, tQ_f);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  while (TRUE)
  {
    gdouble tQ;
    gboolean root_found;

    flag = CVode (Vexp->cvode_qt, tQ_f, Vexp->y_qt, &tQ, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    root_found = (flag == CV_ROOT_RETURN);

    {
      const gdouble ralpha = NV_Ith_S (Vexp->y_qt, 0);
      const gdouble alpha  = ralpha + Vexp->alpha_b;
      const gdouble phi    = NV_Ith_S (Vexp->y_qt, 1);
      gdouble H_lp, x;

      _nc_hicosmo_Vexp_H_x (Vexp, tQ, alpha, phi, &H_lp, &x);

      {
        const gdouble E            = H_lp * Vexp->RH_lp;
        const gdouble tau          = GSL_SIGN (E) * sqrt (2.0 * ralpha);
        const gdouble dlnx_dalpha  = _nc_hicosmo_Vexp_qt_dlnx_dalpha (Vexp, tQ, alpha, phi);
        const NcHICosmoVexpState s = {tau, tau / E, x * tau, dlnx_dalpha};

        g_array_append_val (earray, s);
        
        if (FALSE)
        {
          const gdouble a3 = exp (3.0 * alpha);
          const gdouble ep = _nc_hicosmo_Vexp_epsilon (Vexp, tQ, alpha, phi);
          
          printf ("% 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e % 21.15e %s\n", 
                  tQ, 
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
      gint rootsfound[4]  = {0, 0, 0, 0};
      gdouble H_lp, x;

      _nc_hicosmo_Vexp_H_x (Vexp, tQ, alpha, phi, &H_lp, &x);

      flag = CVodeGetRootInfo (Vexp->cvode_qt, rootsfound);
      NCM_CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, );

      if (rootsfound[3] != 0)
        g_error ("_nc_hicosmo_Vexp_evolve_qt: quantum phase is too long, unable to match with the classical phase.");

      if (rootsfound[0] != 0)
        root1_found = TRUE;

      if ((((x > 0.0) && (rootsfound[1] != 0)) || ((x < 0.0) && (rootsfound[2] != 0))))
      {
        const gint cl_b       = GSL_SIGN (x);
        const gdouble epsilon = _nc_hicosmo_Vexp_epsilon (Vexp, tQ, alpha, phi);
        const gdouble a_q     = exp (alpha);
        const gdouble LpLm    = pow (LAMBDAp, LAMBDAm);
        const gdouble LmLp    = pow (LAMBDAm, LAMBDAp);
        const gdouble d       = D_PHI;

        const gdouble a_0     = _nc_hicosmo_Vexp_a_0_by_ea (Vexp, tQ, epsilon, a_q, cl_b);
        const gdouble alpha_0 = log (a_0);
        const gdouble alpha_q = alpha;
        const gdouble q       = _nc_hicosmo_Vexp_q (epsilon, cl_b);
        const gdouble E       = H_lp * Vexp->RH_lp;
        const gdouble c1      = (cl_b > 0) ? gsl_pow_2 (LAMBDAm * d * Vexp->RH_lp) / gsl_pow_6 (a_0) : gsl_pow_2 (LAMBDAp * d * Vexp->RH_lp) / gsl_pow_6 (a_0);
        const gdouble c2      = c1 / (OMEGA_C * LpLm * LmLp);

        if (!root1_found)
          g_warning ("_nc_hicosmo_Vexp_evolve_qt: imperfect match of the classical solution.");

        if (tQ_f < 0.0)
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

    if (tQ == tQ_f)
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
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  while (TRUE)
  {
    gdouble alpha;

    flag = CVode (cvode, alpha_f, Vexp->y_cl, &alpha, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );
    
    {
      const gdouble q            = NV_Ith_S (Vexp->y_cl, 0);
      const gdouble En           = NV_Ith_S (Vexp->y_cl, 1);
      const gdouble tau          = GSL_SIGN (En) * sqrt (2.0 * (alpha - Vexp->alpha_b));      
      const gdouble x            = _nc_hicosmo_Vexp_x_q (q, Vexp->cl_bc);
      const gdouble dlnx_dalpha  = _nc_hicosmo_Vexp_cl_dlnx_dalpha (Vexp, q, Vexp->cl_bc);
      const NcHICosmoVexpState s = {tau, tau / En, x * tau, dlnx_dalpha};

      if (flag == CV_ROOT_RETURN)
        Vexp->tau_x0 = tau;

      g_array_append_val (Vexp->evol_c, s);
      if (TRUE)
      {
        const gdouble a_a0 = _nc_hicosmo_Vexp_ac_a0c_q (Vexp, q, Vexp->cl_bc);
        const gdouble E    = _nc_hicosmo_Vexp_Hc_H0_q (Vexp, q, Vexp->cl_bc);
        printf ("# C % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                alpha, q, x, exp (alpha - Vexp->alpha_0c), a_a0, exp (alpha - Vexp->alpha_0c) / a_a0 - 1.0, E, En, E / En - 1.0);
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
      const gdouble q            = NV_Ith_S (Vexp->y_cl, 0);
      const gdouble En           = NV_Ith_S (Vexp->y_cl, 1);
      const gdouble tau          = GSL_SIGN (En) * sqrt (2.0 * (alpha - Vexp->alpha_b));      
      const gdouble x            = _nc_hicosmo_Vexp_x_q (q, Vexp->cl_be);
      const gdouble dlnx_dalpha  = _nc_hicosmo_Vexp_cl_dlnx_dalpha (Vexp, q, Vexp->cl_be);
      const NcHICosmoVexpState s = {tau, tau / En, x * tau, dlnx_dalpha};

      if (flag == CV_ROOT_RETURN)
        Vexp->tau_x0 = tau;

      g_array_append_val (Vexp->evol_e, s);
      if (TRUE)
      {
        const gdouble a_a0 = _nc_hicosmo_Vexp_ae_a0e_q (Vexp, q, Vexp->cl_be);
        const gdouble E    = _nc_hicosmo_Vexp_He_H0_q (Vexp, q, Vexp->cl_be);
        printf ("# E % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                alpha, q, x, exp (alpha - Vexp->alpha_0e), a_a0, exp (alpha - Vexp->alpha_0e) / a_a0 - 1.0, E, En, E / En - 1.0);
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
  const gdouble tQ_f   = 1.0e20;
  
  if (!ncm_model_state_is_update (NCM_MODEL (Vexp)))
  {
    const gdouble d = D_PHI;
    Vexp->RH_lp     = nc_hicosmo_RH_planck (cosmo);
    Vexp->alpha_b   = ALPHA_B;
    Vexp->a_0de     = pow (2.0 * gsl_pow_2 (d * LAMBDAp * Vexp->RH_lp) / OMEGA_L, 1.0 / 6.0);

    g_array_set_size (Vexp->evol_c, 0);
    g_array_set_size (Vexp->evol_e, 0);
    
    _nc_hicosmo_Vexp_evolve_qt (Vexp, -tQ_f);
    _nc_hicosmo_Vexp_evolve_qt (Vexp, +tQ_f);
    _nc_hicosmo_Vexp_evolve_cl_c (Vexp);
    _nc_hicosmo_Vexp_evolve_cl_e (Vexp);

    /*printf ("# Contracting branch %u; Expanding branch %u;\n", Vexp->evol_c->len, Vexp->evol_e->len);*/

    {
      const guint len            = Vexp->evol_c->len + Vexp->evol_e->len;
      NcmVector *tau_v           = ncm_vector_new (len);
      NcmVector *xtau_v          = ncm_vector_new (len);
      NcmVector *lnN_v           = ncm_vector_new (len);
      NcmVector *tau_dlnx_dtau_v = ncm_vector_new (len);
      NcmVector *z_v             = ncm_vector_new (Vexp->evol_e->len);
      NcmVector *E2_v            = ncm_vector_new (Vexp->evol_e->len);
      gdouble last_z             = GSL_POSINF;
      glong zi                   = Vexp->evol_e->len;
      glong i, j = 0;

      for (i = Vexp->evol_c->len - 1; i >= 0; i--)
      {
        const NcHICosmoVexpState s = g_array_index (Vexp->evol_c, NcHICosmoVexpState, i);

        if (FALSE)
        {
          const gdouble alpha_alpha0 = 0.5 * s.tau * s.tau + Vexp->alpha_b - Vexp->alpha_0c;
          printf ("%6ld %6u % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                  j + 1, len, s.tau, s.N * exp (-alpha_alpha0), s.xtau, s.dlnx_dalpha * s.tau * s.tau, alpha_alpha0);
        }

        ncm_vector_fast_set (tau_v,           j, s.tau);
        ncm_vector_fast_set (xtau_v,          j, s.xtau);
        ncm_vector_fast_set (lnN_v,           j, log (s.N));
        ncm_vector_fast_set (tau_dlnx_dtau_v, j, s.dlnx_dalpha * s.tau * s.tau);
        
        j++;
      }

      for (i = 0; i < Vexp->evol_e->len; i++)
      {
        const NcHICosmoVexpState s = g_array_index (Vexp->evol_e, NcHICosmoVexpState, i);
        const gdouble alpha_alpha0 = 0.5 * s.tau * s.tau + Vexp->alpha_b - Vexp->alpha_0e;
        const gdouble z            = expm1 (-alpha_alpha0);

        if (FALSE)
        {
          printf ("%6ld %6u % 21.15g % 21.15e % 21.15e % 21.15e % 21.15e\n", 
                  j + 1, len, s.tau, s.N * exp (-alpha_alpha0), s.xtau, s.dlnx_dalpha * s.tau * s.tau, alpha_alpha0);
        }
        
        if ((z * 1.000001 < last_z) && z >= -0.1)
        {
          zi--;
          ncm_vector_fast_set (z_v,  zi, z);
          ncm_vector_fast_set (E2_v, zi, s.tau / s.N);
          last_z = z;
          /*printf ("% 21.15g % 21.15g\n", z, s.E * s.E);*/
        }
        
        ncm_vector_fast_set (tau_v,           j, s.tau);
        ncm_vector_fast_set (xtau_v,          j, s.xtau);
        ncm_vector_fast_set (lnN_v,           j, log (s.N));
        ncm_vector_fast_set (tau_dlnx_dtau_v, j, s.dlnx_dalpha * s.tau * s.tau);

        j++;
      }

      if (zi != 0)
      {
        NcmVector *tmp_z_v  = z_v;
        NcmVector *tmp_E2_v = E2_v;

        z_v  = ncm_vector_get_subvector (z_v,  zi, Vexp->evol_e->len - zi);
        E2_v = ncm_vector_get_subvector (E2_v, zi, Vexp->evol_e->len - zi);

        ncm_vector_free (tmp_z_v);
        ncm_vector_free (tmp_E2_v);
      }
      ncm_spline_set (Vexp->E2_s,              z_v, E2_v, TRUE);
      
      ncm_spline_set (Vexp->xtau_s,          tau_v, xtau_v, TRUE);
      ncm_spline_set (Vexp->lnN_s,           tau_v, lnN_v, TRUE);
      ncm_spline_set (Vexp->tau_dlnx_dtau_s, tau_v, tau_dlnx_dtau_v, TRUE);
      
      ncm_vector_free (z_v);
      ncm_vector_free (E2_v);
      
      ncm_vector_free (tau_dlnx_dtau_v);
      ncm_vector_free (xtau_v);
      ncm_vector_free (lnN_v);
      ncm_vector_free (tau_v);
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
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (cosmo);
  
  _nc_hicosmo_Vexp_prepare (Vexp); 

  return ncm_spline_eval (Vexp->E2_s, z);
}

/****************************************************************************
 * Normalized Hubble function redshift derivative
 ****************************************************************************/

static gdouble
_nc_hicosmo_Vexp_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (cosmo);
  
  _nc_hicosmo_Vexp_prepare (Vexp); 

  return ncm_spline_eval_deriv (Vexp->E2_s, z);
}

static gdouble
_nc_hicosmo_Vexp_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (cosmo);
  
  _nc_hicosmo_Vexp_prepare (Vexp); 

  return ncm_spline_eval_deriv2 (Vexp->E2_s, z);
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_Vexp_H0 (NcHICosmo *cosmo) { return MACRO_H0; }
static gdouble _nc_hicosmo_Vexp_Omega_t0 (NcHICosmo *cosmo) { return OMEGA_C; }

static gdouble 
_nc_hicosmo_Vexp_zeta_eval_mnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);
  
  _nc_hicosmo_Vexp_prepare (Vexp); 

  {
    const gdouble tau2   = tau * tau;
    const gdouble xtau   = ncm_spline_eval (Vexp->xtau_s, tau);
    const gdouble xtau2  = xtau * xtau;
    const gdouble alpha  = 0.5 * tau2 + Vexp->alpha_b;
    const gdouble lna_a0 = alpha - ((tau <= 0.0) ? Vexp->alpha_0c : Vexp->alpha_0e);
    
    if (fabs (xtau / tau) < 1.0e-1)
    {
      const gdouble dalpha = 0.5 * (tau - Vexp->tau_x0) * (tau + Vexp->tau_x0);
      const gdouble x      = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);
/*
      printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", 
              tau, Vexp->tau_x0, tau2, alpha, dalpha, xtau / tau, x, exp (2.0 * lna_a0) * x * x * k);
*/
      return exp (2.0 * lna_a0) * x * x * k;
    }
    
    return exp (2.0 * lna_a0) * xtau2 * k / tau2;
  }
}

static gdouble 
_nc_hicosmo_Vexp_zeta_eval_nu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);
  
  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau2   = tau * tau;
    const gdouble lnN    = ncm_spline_eval (Vexp->lnN_s, tau);
    const gdouble lna_a0 = 0.5 * tau2 + Vexp->alpha_b - ((tau <= 0.0) ? Vexp->alpha_0c : Vexp->alpha_0e);
    return exp (lnN - lna_a0) * k;
  }
}

static gdouble 
_nc_hicosmo_Vexp_zeta_eval_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k)
{
  NcHICosmoVexp *Vexp   = NC_HICOSMO_VEXP (iad);
  gdouble tau_dlnx_dtau; 
  
  _nc_hicosmo_Vexp_prepare (Vexp); 

  tau_dlnx_dtau = ncm_spline_eval (Vexp->tau_dlnx_dtau_s, tau);

  return 2.0 * tau_dlnx_dtau / tau + 2.0 * tau;
}

static void 
_nc_hicosmo_Vexp_zeta_eval_system (NcHIPertIAdiab *iad, const gdouble tau, const gdouble k, gdouble *nu, gdouble *dlnmnu)
{
  NcHICosmoVexp *Vexp   = NC_HICOSMO_VEXP (iad);
  gdouble tau_dlnx_dtau; 
  
  _nc_hicosmo_Vexp_prepare (Vexp); 

  tau_dlnx_dtau = ncm_spline_eval (Vexp->tau_dlnx_dtau_s, tau);
  
  {
    const gdouble tau2   = tau * tau;
    const gdouble lnN    = ncm_spline_eval (Vexp->lnN_s, tau);
    const gdouble lna_a0 = 0.5 * tau2 + Vexp->alpha_b - ((tau <= 0.0) ? Vexp->alpha_0c : Vexp->alpha_0e);

    nu[0]     = exp (lnN - lna_a0) * k;
    dlnmnu[0] = 2.0 * (tau_dlnx_dtau / tau + tau);
  }

  return;
}

static guint 
_nc_hicosmo_Vexp_zeta_nsing (NcHIPertIAdiab *iad, const gdouble k)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);

  _nc_hicosmo_Vexp_prepare (Vexp); 

  if (Vexp->tau_x0 != 0.0)
    return 1;
  else
    return 0;
}

static void 
_nc_hicosmo_Vexp_zeta_get_sing_info (NcHIPertIAdiab *iad, const gdouble k, const guint sing, gdouble *ts, NcmHOAASingType *st)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);

  _nc_hicosmo_Vexp_prepare (Vexp); 

  g_assert_cmpuint (sing, ==, 0);

  ts[0] = Vexp->tau_x0;
  st[0] = NCM_HOAA_SING_TYPE_ZERO;
}

static gdouble 
_nc_hicosmo_Vexp_zeta_eval_sing_mnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, guint sing)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau  = tau_m_taus + Vexp->tau_x0;
    const gdouble xtau = ncm_spline_eval (Vexp->xtau_s, tau);

    if (fabs (xtau / tau) > 0.1)
    {
      return _nc_hicosmo_Vexp_zeta_eval_mnu (iad, tau_m_taus + Vexp->tau_x0, k);
    }
    else
    {
      const gdouble tau2   = tau * tau;
      const gdouble alpha  = 0.5 * tau2 + Vexp->alpha_b;
      const gdouble lna_a0 = alpha - ((tau <= 0.0) ? Vexp->alpha_0c : Vexp->alpha_0e);

      const gdouble dalpha = 0.5 * tau_m_taus * (tau + Vexp->tau_x0);
      const gdouble x      = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);

      if (FALSE)
      {
        const gdouble dalpha0 = 0.5 * (tau - Vexp->tau_x0) * (tau + Vexp->tau_x0);
        const gdouble xtau    = ncm_spline_eval (Vexp->xtau_s, tau);
        const gdouble xx      = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha0);

        printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", 
                tau, Vexp->tau_x0, tau2, alpha, dalpha, x, xtau / tau, xx);
      }
      return exp (2.0 * lna_a0) * x * x * k;
    }
  }
}

static gdouble 
_nc_hicosmo_Vexp_zeta_eval_sing_dlnmnu (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, guint sing)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);

  _nc_hicosmo_Vexp_prepare (Vexp);

  {
    const gdouble tau  = tau_m_taus + Vexp->tau_x0;
    const gdouble xtau = ncm_spline_eval (Vexp->xtau_s, tau);

    if (fabs (xtau / tau) > 0.1)
    {
      return _nc_hicosmo_Vexp_zeta_eval_dlnmnu (iad, tau_m_taus + Vexp->tau_x0, k);
    }
    else
    {
      const gdouble dalpha    = 0.5 * tau_m_taus * (tau + Vexp->tau_x0);
      const gdouble x         = _nc_hicosmo_Vexp_init_x_alpha_series (dalpha);
      const gdouble dx        = _nc_hicosmo_Vexp_init_dx_alpha_series (dalpha);
      const gdouble dlnx      = dx / x;
      const gdouble dlnx_dtau = dlnx * tau;

      if (TRUE)
      {
        const gdouble tau_dlnx_dtau = ncm_spline_eval (Vexp->tau_dlnx_dtau_s, tau);
        printf ("% 21.15g % 21.15e % 21.15e\n", tau, 2.0 * (tau_dlnx_dtau / tau + tau) * tau_m_taus, 2.0 * (dlnx_dtau + tau) * tau_m_taus);        
      }

      return 2.0 * (dlnx_dtau + tau);
    }
  }
}

static void 
_nc_hicosmo_Vexp_zeta_eval_sing_system (NcHIPertIAdiab *iad, const gdouble tau_m_taus, const gdouble k, guint sing, gdouble *nu, gdouble *dlnmnu)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (iad);

  nu[0]     = _nc_hicosmo_Vexp_zeta_eval_nu (iad, tau_m_taus + Vexp->tau_x0, k);
  dlnmnu[0] = _nc_hicosmo_Vexp_zeta_eval_sing_dlnmnu (iad, tau_m_taus, k, sing);
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

/**
 * nc_hicosmo_Vexp_tau_min:
 * @Vexp: a #NcHICosmoVexp
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble 
nc_hicosmo_Vexp_tau_min (NcHICosmoVexp *Vexp)
{
  _nc_hicosmo_Vexp_prepare (Vexp);
  return ncm_vector_get (Vexp->xtau_s->xv, 0);
}

/**
 * nc_hicosmo_Vexp_tau_max:
 * @Vexp: a #NcHICosmoVexp
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble 
nc_hicosmo_Vexp_tau_max (NcHICosmoVexp *Vexp)
{
  _nc_hicosmo_Vexp_prepare (Vexp);
  
  return ncm_vector_get (Vexp->xtau_s->xv, ncm_vector_len (Vexp->xtau_s->xv) - 1);
}
