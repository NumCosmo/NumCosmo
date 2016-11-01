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
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

G_DEFINE_TYPE (NcHICosmoVexp, nc_hicosmo_Vexp, NC_TYPE_HICOSMO);

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_Vexp_init (NcHICosmoVexp *Vexp)
{
  Vexp->cvode_qt = CVodeCreate (CV_BDF, CV_NEWTON);
  Vexp->cvode_cl = CVodeCreate (CV_BDF, CV_NEWTON);
  Vexp->qt_init  = FALSE;
  Vexp->cl_init  = FALSE;
  Vexp->y_qt     = N_VNew_Serial (2);
  Vexp->ydot_qt  = N_VNew_Serial (2);
  Vexp->y_cl     = NULL;
}

static void
nc_hicosmo_Vexp_finalize (GObject *object)
{
  NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (object);
  
  if (Vexp->cvode_qt != NULL)
  {
    CVodeFree (&Vexp->cvode_qt);
    Vexp->cvode_qt = NULL;
  }

  if (Vexp->cvode_cl != NULL)
  {
    CVodeFree (&Vexp->cvode_cl);
    Vexp->cvode_cl = NULL;
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

  object_class->finalize     = &nc_hicosmo_Vexp_finalize;

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
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_VEXP_ALPHA_0, "\\alpha_0", "alpha0",
                              0.0,  G_MAXDOUBLE, 1.0e-5,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_VEXP_DEFAULT_ALPHA_0,
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
#define SIGMA_PHI (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_SIGMA_PHI))
#define D_PHI     (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_D_PHI))
#define ALPHA_0   (ncm_vector_get (VECTOR, NC_HICOSMO_VEXP_ALPHA_0))

static void 
_nc_hicosmo_Vexp_dalpha_dphi (NcHICosmoVexp *Vexp, const gdouble t, const gdouble alpha, const gdouble phi, gdouble *dalpha, gdouble *dphi)
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
  const gdouble exp3alpha    = exp (3.0 * alpha * 0.0);

  const gdouble dphi_num     = - alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
	const gdouble dalpha_num   = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
	const gdouble den          = 2.0 * exp3alpha * (cos_cosh + 1.0);
  
  dphi[0]   = dphi_num / den;
	dalpha[0] = dalpha_num / den ;

/*printf ("% 21.15g % 21.15g % 21.15g\n", dalpha_num, dphi_num, den);*/
  
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

static int 
_nc_hicosmo_Vexp_qt_f (realtype t, N_Vector y_qt, N_Vector ydot_qt, void *f_data)
{
	NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (f_data);
	const double alpha  = NV_Ith_S (y_qt, 0);
	const double phi	  = NV_Ith_S (y_qt, 1);
  gdouble dalpha, dphi;

  _nc_hicosmo_Vexp_dalpha_dphi (Vexp, t, alpha, phi, &dalpha, &dphi);
  
	NV_Ith_S (ydot_qt, 0) = dalpha;
	NV_Ith_S (ydot_qt, 1) = dphi;
	
	return 0;
}

static gint
_nc_hicosmo_Vexp_qt_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y_qt, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (jac_data);
  NcHICosmo *cosmo    = NC_HICOSMO (Vexp);

	const gdouble alpha = NV_Ith_S (y_qt, 0);
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
  const gdouble exp3alpha    = exp (3.0 * alpha);

  const gdouble dphi_num     = - alpha * sigma2 * sin_cosh + 2.0 * d * cos_cosh + 2.0 * d;
	const gdouble dalpha_num   = phi * sigma2 * sin_cosh + 2.0 * d * tanh_arg_hyp;
	const gdouble den          = 2.0 * exp3alpha * (cos_cosh + 1.0);

  const gdouble dalpha       = dalpha_num / den;
  const gdouble dphi         = dphi_num / den;
  
  DENSE_ELEM (J, 0, 0) = -3.0 * dalpha + d * sigma2 * phi / exp3alpha - dalpha * (- 2.0 * d * sin_cosh + sigma2 * phi * tanh_arg_hyp) / (cos_cosh + 1.0);
  DENSE_ELEM (J, 0, 1) = ((sigma2 * sin_cosh + 2.0 * d * alpha * sigma2) * cos_cosh + sigma2 * sin_cosh - phi * sigma2 * sin_cosh * sigma2 * alpha * tanh_arg_hyp + 2.0 * d * alpha * sigma2 / gsl_pow_2 (cosh_arg_hyp)) / (den * (cos_cosh + 1.0));
    /*(sigma2 * sin_cosh + 2.0 * d * alpha * sigma2) / den - dalpha * sigma2 * alpha * tanh_arg_hyp / (cos_cosh + 1.0);*/

  DENSE_ELEM (J, 1, 0) = -3.0 * dphi - dphi * (-2.0 * d * sin_cosh + sigma2 * phi * tanh_arg_hyp) / (cos_cosh + 1.0) + (-sigma2 * sin_cosh - 4.0 * d * d * sin_cosh - 2.0 * d * alpha * sigma2 * cos_cosh + 2.0 * d * sigma2 * phi * tanh_arg_hyp) / den;
  DENSE_ELEM (J, 1, 1) = sigma2 * sigma2 * alpha * alpha * sin_cosh * tanh_arg_hyp / (2.0 * exp3alpha * gsl_pow_2 (cos_cosh + 1.0));

  printf ("JAC % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g\n", alpha, phi, dalpha, dphi, DENSE_ELEM (J, 0, 0), DENSE_ELEM (J, 0, 1), DENSE_ELEM (J, 1, 0), DENSE_ELEM (J, 1, 1));
  
  return 0;
}

static gint 
_nc_hicosmo_Vexp_qt_root (realtype t, N_Vector y_qt, realtype *gout, void *user_data)
{
	NcHICosmoVexp *Vexp = NC_HICOSMO_VEXP (user_data);
	const gdouble alpha = NV_Ith_S (y_qt, 0);
	const gdouble phi	  = NV_Ith_S (y_qt, 1);

	gout[0] = log (fabs (_nc_hicosmo_Vexp_dphi_dalpha_m1_tgh (Vexp, t, alpha, phi)) * 1.0e150);
	
	return 0;
} 


static void
_nc_hicosmo_Vexp_init_qt (NcHICosmoVexp *Vexp)
{
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
  gint flag;

  NV_Ith_S (Vexp->y_qt, 0) = ALPHA_0;
  NV_Ith_S (Vexp->y_qt, 1) = 1.0e-100;
  
  if (!Vexp->qt_init)
  {
    flag = CVodeInit (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_f, 0.0, Vexp->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

		flag = CVodeSetMaxStep (Vexp->cvode_qt, G_MAXUINT32);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, );

    flag = CVodeSStolerances (Vexp->cvode_qt, 1.0e-14, 0.0);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetUserData (Vexp->cvode_qt, Vexp);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVDense (Vexp->cvode_qt, 2);
    NCM_CVODE_CHECK (&flag, "CVDense", 1, );

    flag = CVodeRootInit (Vexp->cvode_qt, 1, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );

//    flag = CVDlsSetDenseJacFn (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
//    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );    
  }
  else
  {
    flag = CVodeReInit (Vexp->cvode_qt, 0.0, Vexp->y_qt);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeRootInit (Vexp->cvode_qt, 1, &_nc_hicosmo_Vexp_qt_root);
    NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );    

//    flag = CVDlsSetDenseJacFn (Vexp->cvode_qt, &_nc_hicosmo_Vexp_qt_J);
//    NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );    
  }
}

static void
_nc_hicosmo_Vexp_prepare (NcHICosmoVexp *Vexp)
{
  NcHICosmo *cosmo = NC_HICOSMO (Vexp);
  
  if (!ncm_model_state_is_update (NCM_MODEL (Vexp)))
  {
    gint flag;
    printf ("Preparing!\n");

    _nc_hicosmo_Vexp_init_qt (Vexp);

    while (TRUE)
    {
      gdouble t;
      gboolean root_found;
      
      flag = CVode (Vexp->cvode_qt, 1.0e50, Vexp->y_qt, &t, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode", 1, );

      root_found = (flag == CV_ROOT_RETURN);
        
      flag = CVodeGetDky (Vexp->cvode_qt, t, 1, Vexp->ydot_qt);
      NCM_CVODE_CHECK (&flag, "CVodeGetDky", 1, );

      {
        const gdouble a3    = exp (3.0 * NV_Ith_S (Vexp->y_qt, 0));
        const gdouble ep    = _nc_hicosmo_Vexp_dphi_dalpha_m1_tgh (Vexp, t, NV_Ith_S (Vexp->y_qt, 0), NV_Ith_S (Vexp->y_qt, 1));
        const gdouble RH_lp = nc_hicosmo_RH_planck (cosmo);
        const gdouble H_lp  = NV_Ith_S (Vexp->ydot_qt, 0) / a3;
        const gdouble E     = H_lp * RH_lp;
        const gdouble x     = NV_Ith_S (Vexp->ydot_qt, 1) / NV_Ith_S (Vexp->ydot_qt, 0);
        
        printf ("% 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15g % 21.15e % 21.15e %s\n", 
                t, 
                NV_Ith_S (Vexp->y_qt, 0), 
                NV_Ith_S (Vexp->y_qt, 1),
                E,
                x,
                a3 * H_lp,
                ep,
                cbrt ((pow (fabs (ep), 2.0 - 2.0 / M_SQRT2) / (a3 * a3)) * RH_lp * RH_lp),
                cbrt ((pow (fabs (ep), 1.0 + 1.0 / M_SQRT2) / a3) * RH_lp * RH_lp),
                root_found ? "ROOT" : "NOT-FOUND"
                );
      }
      if (root_found)
        break;
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
