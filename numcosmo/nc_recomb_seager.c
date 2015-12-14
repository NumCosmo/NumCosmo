/***************************************************************************
 *            nc_recomb_seager.c
 *
 *  Mon November 05 18:28:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

/**
 * SECTION:nc_recomb_seager
 * @title: NcRecombSeager
 * @short_description: Cosmic recombination implementing Seager (1999).
 * @include: numcosmo/nc_recomb_seager.h
 *
 * Cosmic recobination as describe in [Seager (1999)][XSeager1999].
 * It uses nc_recomb_HeII_ion_saha_x_by_HeIII_He () to obtain the value of $\lambda$
 * where the numerical integration will start.
 * 
 * $
 *  \newcommand{\He}{\text{He}}
 *  \newcommand{\HeI}{\text{HeI}}
 *  \newcommand{\HeII}{\text{HeII}}
 *  \newcommand{\HeIII}{\text{HeIII}}
 *  \newcommand{\Hy}{\text{H}}
 *  \newcommand{\HyI}{\text{HI}}
 *  \newcommand{\HyII}{\text{HII}}
 *  \newcommand{\e}{{\text{e}^-}}
 * $
 * 
 * See #nc_recomb for symbol definitions.
 * 
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_recomb_seager.h"
#include "nc_recomb.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_util.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "nc_enum_types.h"

#include <nvector/nvector_serial.h>
#include <gsl/gsl_sf_exp.h>

G_DEFINE_TYPE (NcRecombSeager, nc_recomb_seager, NC_TYPE_RECOMB);

enum
{
  PROP_0,
  PROP_OPTS,
  PROP_SIZE,
};

static gint H_ion_full_f (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data);
static gint H_ion_full_J (_NCM_SUNDIALS_INT_TYPE N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void
nc_recomb_seager_init (NcRecombSeager *recomb_seager)
{
	NcRecomb *recomb = NC_RECOMB (recomb_seager);
  recomb_seager->cvode = CVodeCreate (CV_BDF, CV_NEWTON);
  NCM_CVODE_CHECK ((void *)recomb_seager->cvode, "CVodeCreate", 0, );
	recomb_seager->init    = FALSE;
  recomb_seager->opts    = 0;
  recomb_seager->H_fudge = 0.0;
  recomb_seager->AGauss1 = -0.140;
  recomb_seager->AGauss2 =  0.079;
  recomb_seager->zGauss1 =  7.280;
  recomb_seager->zGauss2 =  6.730;
  recomb_seager->wGauss1 =  0.180;
  recomb_seager->wGauss2 =  0.330;

	recomb_seager->ion   = &H_ion_full_f;
	recomb_seager->ion_J = &H_ion_full_J;

	recomb_seager->n = 3;

	recomb_seager->y0     = N_VNew_Serial(recomb_seager->n);
  recomb_seager->y      = N_VNew_Serial(recomb_seager->n);
  recomb_seager->abstol = N_VNew_Serial(recomb_seager->n);

	recomb->Xe_s = ncm_spline_cubic_notaknot_new ();
}

static void
nc_recomb_seager_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);
  g_return_if_fail (NC_IS_RECOMB_SEAGER (object));

  switch (prop_id)
  {
    case PROP_OPTS:
      nc_recomb_seager_set_options (recomb_seager, g_value_get_flags (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_recomb_seager_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);
  g_return_if_fail (NC_IS_RECOMB_SEAGER (object));

  switch (prop_id)
  {
    case PROP_OPTS:
      g_value_set_flags (value, nc_recomb_seager_get_options (recomb_seager));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_recomb_seager_constructed (GObject *object)
{
	/* Chain up : start */
  G_OBJECT_CLASS (nc_recomb_seager_parent_class)->constructed (object);
	{
	}
}

static void
nc_recomb_seager_finalize (GObject *object)
{
	NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (object);

	N_VDestroy (recomb_seager->y);
	N_VDestroy (recomb_seager->y0);
	N_VDestroy (recomb_seager->abstol);

	CVodeFree (&recomb_seager->cvode);

	/* Chain up : end */
  G_OBJECT_CLASS (nc_recomb_seager_parent_class)->finalize (object);
}

static void nc_recomb_seager_prepare (NcRecomb *recomb, NcHICosmo *cosmo);

static void
nc_recomb_seager_class_init (NcRecombSeagerClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcRecombClass *recomb_class = NC_RECOMB_CLASS (klass);

  object_class->set_property = &nc_recomb_seager_set_property;
  object_class->get_property = &nc_recomb_seager_get_property;
  object_class->constructed = nc_recomb_seager_constructed;
	object_class->finalize = nc_recomb_seager_finalize;

  /**
   * NcRecombSeager:options:
   *
   * Integration options.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OPTS,
                                   g_param_spec_flags ("options",
                                                       NULL,
                                                       "Integration options",
                                                       NC_TYPE_RECOMB_SEAGER_OPT, NC_RECOM_SEAGER_OPT_ALL,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	recomb_class->prepare = &nc_recomb_seager_prepare;
}

static gdouble nc_recomb_seager_HII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x);
static gdouble nc_recomb_seager_HII_ion_gcor_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x);
static gdouble nc_recomb_seager_HeII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x);
static gdouble nc_recomb_seager_Tm_dx (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x);

static void nc_recomb_seager_HII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad);
static void nc_recomb_seager_HII_ion_gcor_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad);
static void nc_recomb_seager_HeII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad);
static void nc_recomb_seager_Tm_dx_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad);

typedef struct _NcRecombSeagerParams
{
  NcRecombSeager *recomb_seager;
  NcHICosmo *cosmo;
} NcRecombSeagerParams;

static gint
H_ion_full_f (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcRecombSeagerParams *rsp = (NcRecombSeagerParams *) f_data;
	const gdouble x     = exp (-lambda);
  const gdouble XHII  = NV_Ith_S (y, 0);
  const gdouble Tm    = NV_Ith_S (y, 1);
  const gdouble XHeII = NV_Ith_S (y, 2);

  NV_Ith_S (ydot, 0) = -x * nc_recomb_seager_HII_ion_rate (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x);
  NV_Ith_S (ydot, 1) = -x * nc_recomb_seager_Tm_dx (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x);
  NV_Ith_S (ydot, 2) = -x * nc_recomb_seager_HeII_ion_rate (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x);

  return GSL_SUCCESS;
}

static gint
H_ion_full_J (_NCM_SUNDIALS_INT_TYPE N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcRecombSeagerParams *rsp = (NcRecombSeagerParams *) jac_data;
	const gdouble x     = exp (-lambda);
  gdouble XHII  = NV_Ith_S (y, 0);
  gdouble Tm    = NV_Ith_S (y, 1);
  gdouble XHeII = NV_Ith_S (y, 2);
  gdouble grad[3];

  NCM_UNUSED (N);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);

  nc_recomb_seager_HII_ion_rate_grad (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x, grad);
  DENSE_ELEM (J, 0, 0) = -x * grad[0];
  DENSE_ELEM (J, 0, 1) = -x * grad[1];
  DENSE_ELEM (J, 0, 2) = -x * grad[2];

  nc_recomb_seager_Tm_dx_grad (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x, grad);
  DENSE_ELEM (J, 1, 0) = -x * grad[0];
  DENSE_ELEM (J, 1, 1) = -x * grad[1];
  DENSE_ELEM (J, 1, 2) = -x * grad[2];

  nc_recomb_seager_HeII_ion_rate_grad (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x, grad);
  DENSE_ELEM (J, 2, 0) = -x * grad[0];
  DENSE_ELEM (J, 2, 1) = -x * grad[1];
  DENSE_ELEM (J, 2, 2) = -x * grad[2];

  return 0;
}

static gint
H_ion_full_gcor_f (realtype lambda, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcRecombSeagerParams *rsp = (NcRecombSeagerParams *) f_data;
	const gdouble x     = exp (-lambda);
  const gdouble XHII  = NV_Ith_S (y, 0);
  const gdouble Tm    = NV_Ith_S (y, 1);
  const gdouble XHeII = NV_Ith_S (y, 2);

  NV_Ith_S (ydot, 0) = -x * nc_recomb_seager_HII_ion_gcor_rate (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x);
  NV_Ith_S (ydot, 1) = -x * nc_recomb_seager_Tm_dx (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x);
  NV_Ith_S (ydot, 2) = -x * nc_recomb_seager_HeII_ion_rate (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x);

  return GSL_SUCCESS;
}

static gint
H_ion_full_gcor_J (_NCM_SUNDIALS_INT_TYPE N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcRecombSeagerParams *rsp = (NcRecombSeagerParams *) jac_data;
	const gdouble x     = exp (-lambda);
  gdouble XHII  = NV_Ith_S (y, 0);
  gdouble Tm    = NV_Ith_S (y, 1);
  gdouble XHeII = NV_Ith_S (y, 2);
  gdouble grad[3];

  NCM_UNUSED (N);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);

  nc_recomb_seager_HII_ion_gcor_rate_grad (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x, grad);
  DENSE_ELEM (J, 0, 0) = -x * grad[0];
  DENSE_ELEM (J, 0, 1) = -x * grad[1];
  DENSE_ELEM (J, 0, 2) = -x * grad[2];

  nc_recomb_seager_Tm_dx_grad (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x, grad);
  DENSE_ELEM (J, 1, 0) = -x * grad[0];
  DENSE_ELEM (J, 1, 1) = -x * grad[1];
  DENSE_ELEM (J, 1, 2) = -x * grad[2];

  nc_recomb_seager_HeII_ion_rate_grad (rsp->recomb_seager, rsp->cosmo, XHII, Tm, XHeII, x, grad);
  DENSE_ELEM (J, 2, 0) = -x * grad[0];
  DENSE_ELEM (J, 2, 1) = -x * grad[1];
  DENSE_ELEM (J, 2, 2) = -x * grad[2];

  return 0;
}

static gdouble
_nc_recomb_He_fully_ionized_Xe (gdouble lambda, gpointer p)
{
  NcHICosmo *cosmo = NC_HICOSMO (p);
	return nc_recomb_He_fully_ionized_Xe (cosmo, exp(-lambda));
}

static void
nc_recomb_seager_prepare (NcRecomb *recomb, NcHICosmo *cosmo)
{
	NcRecombSeager *recomb_seager = NC_RECOMB_SEAGER (recomb);
	const gdouble x_HeIII = nc_recomb_HeII_ion_saha_x_by_HeIII_He (cosmo, recomb->init_frac);
  const gdouble lambdai = recomb->lambdai;
	const gdouble lambda_HeIII = -log (x_HeIII);
	const gdouble lambdaf = -log (1.0);
  NcRecombSeagerParams pparams = {recomb_seager, cosmo};
	gsl_function F;

	F.function = &_nc_recomb_He_fully_ionized_Xe;
  F.params   = cosmo;

	ncm_spline_set_func (recomb->Xe_s, NCM_SPLINE_FUNCTION_SPLINE,
	                     &F, lambdai, lambda_HeIII, 0, recomb->prec);

  /*****************************************************************************
	 * Assuming hydrogen is completly ionized and no more double ionized helium
	 * i.e., $X_\HeIII = 0$.
	 ****************************************************************************/
	{
		const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
		const gdouble XHII = 1.0;
		gdouble XHeII, Tm, XeXHeII_XHeI;

		Tm = T0 * x_HeIII;
		XeXHeII_XHeI = nc_recomb_HeI_ion_saha (cosmo, x_HeIII);
    XHeII = (XeXHeII_XHeI + 1.0) * ncm_sqrt1px_m1 (4.0 * ncm_c_prim_XHe () * XeXHeII_XHeI / gsl_pow_2 (XeXHeII_XHeI + 1.0)) / 2.0;
    
		NV_Ith_S (recomb_seager->y0, 0) = XHII;
		NV_Ith_S (recomb_seager->y0, 1) = Tm;
		NV_Ith_S (recomb_seager->y0, 2) = XHeII;
	}

	if (!recomb_seager->init)
	{
		gint flag = CVodeInit (recomb_seager->cvode, recomb_seager->ion, lambda_HeIII, recomb_seager->y0);
		NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
		recomb_seager->init = TRUE;
	}
	else
	{
		gint flag = CVodeReInit (recomb_seager->cvode, lambda_HeIII, recomb_seager->y0);
		NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
	}

	{
		const gdouble reltol = GSL_MIN (recomb->prec, 1e-11);
		gint flag;

		NV_Ith_S (recomb_seager->abstol, 0) = GSL_MIN (reltol, recomb->prec * 1e-5);
		NV_Ith_S (recomb_seager->abstol, 2) = GSL_MIN (reltol, recomb->prec * 1e-5);
		NV_Ith_S (recomb_seager->abstol, 1) = 0.0;

		flag = CVodeSVtolerances (recomb_seager->cvode, reltol, recomb_seager->abstol);
		NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

		flag = CVodeSetUserData (recomb_seager->cvode, &pparams);
		NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

		flag = CVodeSetMaxNumSteps (recomb_seager->cvode, 100000);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

		flag = CVDense (recomb_seager->cvode, recomb_seager->n);
		NCM_CVODE_CHECK (&flag, "CVDense", 1, );

		flag = CVDlsSetDenseJacFn (recomb_seager->cvode, recomb_seager->ion_J);
		NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );

		flag = CVodeSetStopTime (recomb_seager->cvode, lambdaf);
		NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

		flag = CVodeSetMaxErrTestFails (recomb_seager->cvode, 14);
		NCM_CVODE_CHECK (&flag, "CVodeSetMaxErrTestFails", 1, );
	}

	{
		NcmVector *lambda_v = ncm_spline_get_xv (recomb->Xe_s);
		NcmVector *Xe_v     = ncm_spline_get_yv (recomb->Xe_s);
		GArray *lambda_a    = ncm_vector_get_array (lambda_v);
		GArray *Xe_a        = ncm_vector_get_array (Xe_v);
		gdouble lambda_last = 0.0;

		ncm_vector_free (lambda_v);
		ncm_vector_free (Xe_v);

		g_array_remove_range (lambda_a, lambda_a->len - 1, 1);
		g_array_remove_range (Xe_a,     Xe_a->len - 1,     1);

		while (TRUE)
		{
			gdouble lambda_i;
			gint flag = CVode (recomb_seager->cvode, lambdaf, recomb_seager->y, &lambda_i, CV_ONE_STEP);
			NCM_CVODE_CHECK (&flag, "CVode", 1, );
			{
				const gdouble XHII  = NV_Ith_S (recomb_seager->y, 0);
				const gdouble XHeII = NV_Ith_S (recomb_seager->y, 2);
				const gdouble Xe    = XHII + XHeII;
/*				printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", expm1 (-lambda_i), 
                NV_Ith_S (recomb_seager->y, 0), 
                NV_Ith_S (recomb_seager->y, 1),
                NV_Ith_S (recomb_seager->y, 2));
*/
        if (fabs ((lambda_last - lambda_i) / lambda_last) > 1e-7)
				{
					g_array_append_val (lambda_a, lambda_i);
					g_array_append_val (Xe_a,     Xe);
					lambda_last = lambda_i;
				}

				if (lambda_i == lambdaf)
					break;
			}
		}

    ncm_spline_set_array (recomb->Xe_s, lambda_a, Xe_a, TRUE);
		g_array_unref (lambda_a);
		g_array_unref (Xe_a);
	}

	if (FALSE)
	{
		guint i;
		printf ("# Xe spline len %u\n", ncm_vector_len (recomb->Xe_s->xv));

		for (i = 0; i < ncm_vector_len (recomb->Xe_s->xv); i++)
    {
			printf ("% 20.15g % 20.15g\n", ncm_vector_get (recomb->Xe_s->xv, i),
			        ncm_vector_get (recomb->Xe_s->yv, i));
		}

		for (i = 0; i < 30000; i++)
    {
			gdouble lambda = recomb->lambdai + (recomb->lambdaf - recomb->lambdai) / (30000.0 - 1.0) * i;
			printf ("% 20.15g % 20.15g\n", lambda,
			        ncm_spline_eval (recomb->Xe_s, lambda));
		}
	}

	recomb->tau_s          = ncm_spline_copy_empty (recomb->Xe_s);
	recomb->dtau_dlambda_s = ncm_spline_copy_empty (recomb->Xe_s);

	_nc_recomb_prepare_tau_splines (recomb, cosmo);
}

static gdouble
nc_recomb_seager_HII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x)
{
  const gdouble Xe                = XHII + XHeII;
  const gdouble kbTm              = ncm_c_kb () * Tm;
  const gdouble x3                = gsl_pow_3 (x);
  const gdouble h2                = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b           = nc_hicosmo_Omega_b (cosmo);
  const gdouble Tm3_2             = sqrt(gsl_pow_3 (Tm));
  const gdouble lambda_e3_T3_2    = gsl_pow_3 (ncm_c_thermal_wl_e ());
  const gdouble n_H               = ncm_c_prim_H_Yp () * Omega_b * x3 * (ncm_c_crit_density () * h2) / ncm_c_rest_energy_p ();
  const gdouble alpha_H           = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble beta_H            = gsl_sf_exp_mult (-ncm_c_H_bind_2s () / kbTm, alpha_H * Tm3_2 / lambda_e3_T3_2);
  const gdouble beta_H_exp_mE_kbT = gsl_sf_exp_mult (-ncm_c_H_Lyman_2s () / kbTm, beta_H);
  const gdouble f1                = (Xe * XHII * n_H * alpha_H - beta_H_exp_mE_kbT * (1.0 - XHII));

  const gdouble H                 = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble KH                = ncm_c_H_Lyman_2p_wl3_8pi () / H;
  const gdouble Lambda_H          = ncm_c_decay_H_rate_2s_1s ();
  const gdouble f2                = 1.0 + KH * Lambda_H * n_H * (1.0 - XHII);

  const gdouble f3                = H * x * (1.0 + KH * (Lambda_H + beta_H) * n_H * (1.0 - XHII));

  return f1 * f2 / f3;
}

static void
nc_recomb_seager_HII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad)
{
  const gdouble x2    = x * x;
  const gdouble x3    = x2 * x;
  const gdouble Xe    = XHII + XHeII;
  const gdouble Tm2   = Tm * Tm;
  const gdouble Tm3_2 = sqrt (Tm * Tm2);
  const gdouble alpha = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble beta  = ncm_c_boltzmann_factor_H_2s (Tm) * Tm3_2 * alpha;

  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
  const gdouble n_b0    = Omega_b * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0     = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H       = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();

  const gdouble f1 = alpha * n_0 * x2 / H;

  const gdouble f2na = ncm_c_decay_H_rate_2s_1s () * (1.0 - XHII);
  const gdouble f2nb = H / x3 / (n_0 * ncm_c_H_Lyman_2p_wl3_8pi ());
  const gdouble f2n = f2na + f2nb;

  const gdouble f2da = beta * (1.0 - XHII);
  const gdouble f2d = f2na + f2nb + f2da;

  const gdouble f2 = f2n / f2d;

  const gdouble S = ncm_c_boltzmann_factor_H_1s (Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = (XHII * Xe - (1.0 - XHII) * S);

  const gdouble ddX = f1 * (-ncm_c_decay_H_rate_2s_1s () / f2d) * f3 +
    f1 * ((ncm_c_decay_H_rate_2s_1s () + beta) * f2n / gsl_pow_2 (f2d)) * f3 +
    f1 * f2n / f2d * (Xe + XHII + S);

  const gdouble dalpha = nc_recomb_seager_pequignot_HI_case_B_dTm (recomb_seager, cosmo, Tm);
  const gdouble df1 = f1 * dalpha / alpha;
  const gdouble dbeta =
    (3.0 / 2.0 / Tm + ncm_c_H_bind_1s () / ncm_c_kb () / Tm2 + dalpha / alpha) * beta;
  const gdouble df2 = - f2n * (dbeta * (1.0 - XHII)) / gsl_pow_2 (f2d);
  const gdouble dS =
    (3.0 / 2.0 / Tm + ncm_c_H_bind_2s () / ncm_c_kb () / Tm2) * S;
  const gdouble df3 = -(1.0 - XHII) * dS;
  const gdouble ddTm = df1 * f2 * f3 + f1 * df2 * f3 + f1 * f2 * df3;

  const gdouble ddXHeII = f1 * f2 * XHII;

  grad[0] = ddX;
  grad[1] = ddTm;
  grad[2] = ddXHeII;

  return;
}

static gdouble
nc_recomb_seager_HII_ion_gcor_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x)
{
  const gdouble Xe                = XHII + XHeII;
  const gdouble kbTm              = ncm_c_kb () * Tm;
  const gdouble x3                = gsl_pow_3 (x);
  const gdouble h2                = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b           = nc_hicosmo_Omega_b (cosmo);
  const gdouble Tm3_2             = sqrt(gsl_pow_3 (Tm));
  const gdouble lambda_e3_T3_2    = gsl_pow_3 (ncm_c_thermal_wl_e ());
  const gdouble n_H               = ncm_c_prim_H_Yp () * Omega_b * x3 * (ncm_c_crit_density () * h2) / ncm_c_rest_energy_p ();
  const gdouble alpha_H           = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble beta_H            = gsl_sf_exp_mult (-ncm_c_H_bind_2s () / kbTm, alpha_H * Tm3_2 / lambda_e3_T3_2);
  const gdouble beta_H_exp_mE_kbT = gsl_sf_exp_mult (-ncm_c_H_Lyman_2s () / kbTm, beta_H);
  const gdouble f1                = (Xe * XHII * n_H * alpha_H - beta_H_exp_mE_kbT * (1.0 - XHII));

  const gdouble lnx               = log (x);
  const gdouble K_gcor            = 1.0 + 
    recomb_seager->AGauss1 * exp (-gsl_pow_2 ((lnx - recomb_seager->zGauss1) / recomb_seager->wGauss1)) +
    recomb_seager->AGauss2 * exp (-gsl_pow_2 ((lnx - recomb_seager->zGauss2) / recomb_seager->wGauss2));
  
  const gdouble H                 = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();  
  const gdouble KH                = ncm_c_H_Lyman_2p_wl3_8pi () * K_gcor / H;
  const gdouble Lambda_H          = ncm_c_decay_H_rate_2s_1s ();
  const gdouble f2                = 1.0 + KH * Lambda_H * n_H * (1.0 - XHII);

  const gdouble f3                = H * x * (1.0 + KH * (Lambda_H + beta_H) * n_H * (1.0 - XHII));

  return f1 * f2 / f3;
}

static void
nc_recomb_seager_HII_ion_gcor_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad)
{
  const gdouble x2    = x * x;
  const gdouble x3    = x2 * x;
  const gdouble Xe    = XHII + XHeII;
  const gdouble Tm2   = Tm * Tm;
  const gdouble Tm3_2 = sqrt (Tm * Tm2);
  const gdouble alpha = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble beta  = ncm_c_boltzmann_factor_H_2s (Tm) * Tm3_2 * alpha;

  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
  const gdouble n_b0    = Omega_b * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0     = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H       = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();

  const gdouble f1 = alpha * n_0 * x2 / H;

  const gdouble lnx               = log (x);
  const gdouble K_gcor            = 1.0 + 
    recomb_seager->AGauss1 * exp (-gsl_pow_2 ((lnx - recomb_seager->zGauss1) / recomb_seager->wGauss1)) +
    recomb_seager->AGauss2 * exp (-gsl_pow_2 ((lnx - recomb_seager->zGauss2) / recomb_seager->wGauss2));
  
  const gdouble f2na = ncm_c_decay_H_rate_2s_1s () * (1.0 - XHII);
  const gdouble f2nb = H / x3 / (n_0 * ncm_c_H_Lyman_2p_wl3_8pi () * K_gcor);
  const gdouble f2n = f2na + f2nb;

  const gdouble f2da = beta * (1.0 - XHII);
  const gdouble f2d = f2na + f2nb + f2da;

  const gdouble f2 = f2n / f2d;

  const gdouble S = ncm_c_boltzmann_factor_H_1s (Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = (XHII * Xe - (1.0 - XHII) * S);

  const gdouble ddX = f1 * (-ncm_c_decay_H_rate_2s_1s () / f2d) * f3 +
    f1 * ((ncm_c_decay_H_rate_2s_1s () + beta) * f2n / gsl_pow_2 (f2d)) * f3 +
    f1 * f2n / f2d * (Xe + XHII + S);

  const gdouble dalpha = nc_recomb_seager_pequignot_HI_case_B_dTm (recomb_seager, cosmo, Tm);
  const gdouble df1 = f1 * dalpha / alpha;
  const gdouble dbeta =
    (3.0 / 2.0 / Tm + ncm_c_H_bind_1s () / ncm_c_kb () / Tm2 + dalpha / alpha) * beta;
  const gdouble df2 = - f2n * (dbeta * (1.0 - XHII)) / gsl_pow_2 (f2d);
  const gdouble dS =
    (3.0 / 2.0 / Tm + ncm_c_H_bind_2s () / ncm_c_kb () / Tm2) * S;
  const gdouble df3 = -(1.0 - XHII) * dS;
  const gdouble ddTm = df1 * f2 * f3 + f1 * df2 * f3 + f1 * f2 * df3;

  const gdouble ddXHeII = f1 * f2 * XHII;

  grad[0] = ddX;
  grad[1] = ddTm;
  grad[2] = ddXHeII;

  return;
}

static gdouble
nc_recomb_seager_HeII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x)
{
  const gdouble Xe = XHII + XHeII;
	const gdouble x2 = x * x;
  const gdouble alpha = nc_recomb_seager_hummer_HeI_case_B (recomb_seager, cosmo, Tm);
  const gdouble h2 = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
  const gdouble n_b0 = Omega_b * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0 = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble f1 = alpha * n_0 * x2 / H;

  const gdouble x3 = x2 * x;
  const gdouble Tm3_2 = sqrt(gsl_pow_3 (Tm));
  const gdouble f2na = ncm_c_decay_He_rate_2s_1s () * (ncm_c_prim_XHe () - XHeII);
  const gdouble f2nb = gsl_sf_exp_mult (-ncm_c_HeI_2s_m_2p_kb () / Tm, H / x3 / (n_0 * ncm_c_HeI_Lyman_2p_wl3_8pi ()));
  const gdouble f2n = f2na + f2nb;
  const gdouble f2da = 4.0 * ncm_c_boltzmann_factor_HeI_2s (Tm) * Tm3_2 * alpha * (ncm_c_prim_XHe () - XHeII);
  const gdouble f2d = f2na + f2nb + f2da;
  const gdouble f2 = f2n / f2d;

  const gdouble S = 4.0 * ncm_c_boltzmann_factor_HeI_1s (Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = (XHeII * Xe - (ncm_c_prim_XHe () - XHeII) * S);

  return f1 * f2 * f3;
}

static void
nc_recomb_seager_HeII_ion_rate_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad)
{
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble Xe = XHII + XHeII;
  const gdouble Tm2 = Tm * Tm;
  const gdouble alpha = nc_recomb_seager_hummer_HeI_case_B (recomb_seager, cosmo, Tm);
  const gdouble h2 = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
  const gdouble n_b0 = Omega_b * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0 = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble f1 = alpha * n_0 * x2 / H;

  const gdouble Tm3_2 = sqrt (gsl_pow_3 (Tm));
  const gdouble f2na = ncm_c_decay_He_rate_2s_1s () * (ncm_c_prim_XHe () - XHeII);
  const gdouble f2nb = gsl_sf_exp_mult (-ncm_c_HeI_2s_m_2p_kb () / Tm, H / x3 / (n_0 * ncm_c_HeI_Lyman_2p_wl3_8pi ()));
  const gdouble f2n = f2na + f2nb;
  const gdouble beta = 4.0 * ncm_c_boltzmann_factor_HeI_2s (Tm) * Tm3_2 * alpha;
  const gdouble f2da = beta * (ncm_c_prim_XHe () - XHeII);
  const gdouble f2d = f2na + f2nb + f2da;
  const gdouble f2 = f2n / f2d;

  const gdouble S = 4.0 * ncm_c_boltzmann_factor_HeI_1s (Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = ( XHeII * Xe - (ncm_c_prim_XHe () - XHeII) * S);

  const gdouble ddXHeII = f1 * (-ncm_c_decay_He_rate_2s_1s () / f2d) * f3 +
    f1 * ((ncm_c_decay_He_rate_2s_1s () + beta) * f2n / gsl_pow_2 (f2d)) * f3 +
    f1 * f2n / f2d * (Xe + XHeII + S);

  const gdouble dalpha = nc_recomb_seager_hummer_HeI_case_B_dTm (recomb_seager, cosmo, Tm);
  const gdouble df1 = f1 * dalpha / alpha;
  const gdouble dbeta = (3.0 / 2.0 / Tm + ncm_c_HeI_bind_1s () / ncm_c_kb () / Tm2 + dalpha / alpha) * beta;
  const gdouble df2nb = ncm_c_HeI_2s_m_2p_kb () / Tm2 * f2nb;
  const gdouble df2 = df2nb / f2d - f2n * (dbeta * (ncm_c_prim_XHe () - XHeII) + df2nb) / gsl_pow_2 (f2d);
  const gdouble dS =
    (3.0 / 2.0 / Tm + ncm_c_HeI_bind_2s () / ncm_c_kb () / Tm2) * S;
  const gdouble df3 = -(ncm_c_prim_XHe () - XHeII) * dS;
  const gdouble ddTm = df1 * f2 * f3 + f1 * df2 * f3 + f1 * f2 * df3;

  const gdouble ddX = f1 * f2 * XHeII;

  grad[0] = ddX;
  grad[1] = ddTm;
  grad[2] = ddXHeII;
}

static gdouble
nc_recomb_seager_Tm_dx (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T = T0 * x;
	const gdouble x3 = gsl_pow_3 (x);
	const gdouble T04 = gsl_pow_4 (T0);
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();

  const gdouble f1 = (8.0 * ncm_c_thomson_cs () * ncm_c_AR () * T04 /
    (3.0 * ncm_c_c () * ncm_c_mass_e ())) * x3 / H;

  const gdouble Xe = XHII + XHeII;
  const gdouble f2 = Xe * (Tm - T) / (1.0 + ncm_c_prim_XHe () + Xe);

  const gdouble f3 = 2.0 * Tm / x;

  return f1 * f2 + f3;
}

static void
nc_recomb_seager_Tm_dx_grad (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, const gdouble XHII, const gdouble Tm, const gdouble XHeII, const gdouble x, gdouble *grad)
{
  const gdouble T0 = nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T = T0 * x;
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble f1 = (8.0  * ncm_c_thomson_cs () * ncm_c_AR () *
    T0 * T0 * T0 * T0 /
    (3.0 * ncm_c_c () * ncm_c_mass_e ())) *
    gsl_pow_3 (x) / H;

  const gdouble Xe = XHII + XHeII;
  const gdouble f2 = Xe * (Tm - T) / (1.0 + ncm_c_prim_XHe () + Xe);

  const gdouble df2_dX = (1.0 / Xe - 1.0 / (1.0 + ncm_c_prim_XHe () + Xe)) * f2;
  const gdouble ddXHII = f1 * df2_dX;

  const gdouble ddTm = f1 * Xe / (1.0 + ncm_c_prim_XHe () + Xe) + 2.0 / x;

  const gdouble ddXHeII = ddXHII;

  grad[0] = ddXHII;
  grad[1] = ddTm;
  grad[2] = ddXHeII;

  return;
}

/**
 * nc_recomb_seager_new:
 * 
 * Creates a new #NcRecombSeager using default properties.
 * 
 * Returns: (transfer full): a new #NcRecombSeager.
 */
NcRecombSeager *
nc_recomb_seager_new (void)
{
  return g_object_new (NC_TYPE_RECOMB_SEAGER,
                       NULL);
}

/**
 * nc_recomb_seager_new_full:
 * @init_frac: inital fraction of $X_{\HeIII}/X_{\He}$ where to start numerical integration
 * @zi: inital redshift
 * @prec: integration precision
 * 
 * Creates a new #NcRecombSeager using @init_frac, @zi and @prec.
 * 
 * Returns: (transfer full): a new #NcRecombSeager.
 */
NcRecombSeager *
nc_recomb_seager_new_full (gdouble init_frac, gdouble zi, gdouble prec)
{
  return g_object_new (NC_TYPE_RECOMB_SEAGER,
                       "init-frac", init_frac,
                       "zi", zi,
                       "prec", prec,
                       NULL);
}

/**
 * nc_recomb_seager_ref:
 * @recomb_seager: a #NcRecombSeager
 *
 * Increases the reference count of @recomb_seager.
 *
 * Returns: (transfer full): @recomb_seager.
 */
NcRecombSeager *
nc_recomb_seager_ref (NcRecombSeager *recomb_seager)
{
  return NC_RECOMB_SEAGER (g_object_ref (recomb_seager));
}

/**
 * nc_recomb_seager_free:
 * @recomb_seager: a #NcRecombSeager.
 *
 * Decreases the reference count of @recomb_seager.
 *
 */
void
nc_recomb_seager_free (NcRecombSeager *recomb_seager)
{
  g_object_unref (recomb_seager);
}

/**
 * nc_recomb_seager_clear:
 * @recomb_seager: a #NcRecombSeager.
 *
 * Decreases the reference count of *@recomb_seager if
 * *@recomb_seager is not NULL, then sets *@recomb_seager to NULL.
 *
 */
void
nc_recomb_seager_clear (NcRecombSeager **recomb_seager)
{
  g_clear_object (recomb_seager);
}

/**
 * nc_recomb_seager_set_options:
 * @recomb_seager: a #NcRecombSeager
 * @opts: a #NcRecombSeagerOpt
 * 
 * Sets integration options.
 *
 */
void 
nc_recomb_seager_set_options (NcRecombSeager *recomb_seager, NcRecombSeagerOpt opts)
{
  if (recomb_seager->opts != opts)
  {
    NcRecomb *recomb = NC_RECOMB (recomb_seager);
    
    if (opts & NC_RECOM_SEAGER_OPT_FUDGE_GAUSS_COR)
    {
      recomb_seager->H_fudge = 1.125;
      recomb_seager->ion     = &H_ion_full_gcor_f;
      recomb_seager->ion_J   = &H_ion_full_gcor_J;
    }
    else if (opts & NC_RECOM_SEAGER_OPT_FUDGE)
    {
      recomb_seager->H_fudge = 1.14;
      recomb_seager->ion     = &H_ion_full_f;
      recomb_seager->ion_J   = &H_ion_full_J;
    }
    else
    {
      recomb_seager->H_fudge = 1.0;
      recomb_seager->ion     = &H_ion_full_f;
      recomb_seager->ion_J   = &H_ion_full_J;
    }

    if (recomb_seager->init)
    {
      CVodeFree (&recomb_seager->cvode);
      recomb_seager->cvode = CVodeCreate (CV_BDF, CV_NEWTON);
      NCM_CVODE_CHECK ((void *)recomb_seager->cvode, "CVodeCreate", 0, );
      recomb_seager->init = FALSE;
    }

    ncm_model_ctrl_force_update (recomb->ctrl);
    recomb_seager->opts = opts;
  }
}

/**
 * nc_recomb_seager_get_options:
 * @recomb_seager: a #NcRecombSeager
 * 
 * Gets integration options.
 * 
 * Returns: currently used integration options.
 */
NcRecombSeagerOpt 
nc_recomb_seager_get_options (NcRecombSeager *recomb_seager)
{
  return recomb_seager->opts;
}


/**
 * nc_recomb_seager_pequignot_HI_case_B:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HyII$ recombination coefficient.
 *
 * The fitting formula of the case B recombination coefficient for $\HyII$ as
 * in [Pequignot (1991)][XPequignot1991].
 *
 * Returns: the value of the case B recombination coefficient for
 * $\HyII$, $\alpha_H$ .
 */
gdouble
nc_recomb_seager_pequignot_HI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble F =  recomb_seager->H_fudge;   /* fudge factor */
  const gdouble G =  1e-19;
  const gdouble a =  4.309;
  const gdouble b = -0.6166;
  const gdouble c =  0.6703;
  const gdouble d =  0.5300;
  const gdouble t = Tm * 1.0e-4;
  const gdouble res = F * G * a * pow (t, b) / (1.0 + c * pow (t, d));
  NCM_UNUSED (cosmo);
  return res;
}

/**
 * nc_recomb_seager_pequignot_HI_case_B_dTm:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HyII$ recombination coefficient derivative with respect to $T_m$.
 *
 * The derivative of the fitting formula of the case B recombination coefficient for $\HyII$
 * nc_recomb_seager_pequignot_HI_case_B ().
 *
 * Returns: the value of the case B recombination coefficient for $\HyII$, $d\alpha_H/dT_m$.
 */
gdouble
nc_recomb_seager_pequignot_HI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble F =  recomb_seager->H_fudge;   /* fudge factor */
  const gdouble G =  1e-19;
  const gdouble a =  4.309;
  const gdouble b = -0.6166;
  const gdouble c =  0.6703;
  const gdouble d =  0.5300;
  const gdouble t = Tm * 1e-4;
  const gdouble t_b = pow (t, b);
  const gdouble t_d = pow (t, d);
  const gdouble res = a * F * G * (b + c * (b - d) * t_d) * t_b / (Tm * gsl_pow_2 (1.0 + c * t_d));
  NCM_UNUSED (cosmo);
  return res;
}

/**
 * nc_recomb_seager_hummer_HeI_case_B:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HeII$ recombination coefficient.
 *
 * The fitting formula of the case B recombination coefficient for $\HeII$ as
 * in [Hummer (1998)][XHummer1998].
 *
 * Returns: the value of the case B recombination coefficient for $\HeII$, $\alpha_H$ .
 */
gdouble
nc_recomb_seager_hummer_HeI_case_B (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble sqrt_Tm = sqrt (Tm);
  const gdouble sqrt_T1 = pow (10.0, 5.114 / 2.0);
  const gdouble sqrt_T2 = sqrt (3.0);
  const gdouble p = 0.711;
  const gdouble q = pow (10.0, -16.744);
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble res = q / (sqrt_Tm_T2 * pow (1.0 + sqrt_Tm_T2, 1.0 - p) * pow (1.0 + sqrt_Tm_T1, 1.0 + p));
  NCM_UNUSED (cosmo);
  return res;
}

/**
 * nc_recomb_seager_hummer_HeI_case_B_dTm:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @Tm: the matter (baryons) temperature $T_m$
 *
 * The case B $\HeII$ recombination coefficient derivative with respect to Tm.
 *
 * The derivative of the fitting formula of the case B recombination coefficient for $\HeII$
 * nc_recomb_seager_hummer_HeI_case_B ().
 *
 * Returns: the value of the case B recombination coefficient for $\HeII$, $d\alpha_H/dT_m$.
 */
gdouble
nc_recomb_seager_hummer_HeI_case_B_dTm (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble sqrt_Tm = sqrt(Tm);
  const gdouble q = pow (10.0, -16.744);
  const gdouble T1 = pow (10.0, 5.114);
  const gdouble sqrt_T1 = sqrt (T1);
  const gdouble sqrt_T2 = sqrt (3.0);
  const gdouble p = 0.711;
  const gdouble T1_2 = T1 * T1;
  const gdouble sqrt_Tm_T2 = sqrt_Tm / sqrt_T2;
  const gdouble sqrt_Tm_T1 = sqrt_Tm / sqrt_T1;
  const gdouble Tm_T1_3_2 = gsl_pow_3 (sqrt_Tm_T1);
  const gdouble res = -q *
    (
      Tm * (2.0 + p + 3.0 * sqrt_Tm_T2) +
      T1 * sqrt_Tm_T1 * (1.0 + (2.0 - p) * sqrt_Tm_T2)
    ) /
    (
      2.0 * T1_2 * Tm_T1_3_2 * sqrt_Tm_T2 *
      pow (1.0 + sqrt_Tm_T2, 2.0 - p) *
      pow (1.0 + sqrt_Tm_T1, 2.0 + p)
      );
  NCM_UNUSED (cosmo);
  return res;
}

/**
 * nc_recomb_seager_weinberg_HII_ion_rate:
 * @recomb_seager: a #NcRecombSeager
 * @cosmo: a #NcHICosmo
 * @XHII: FIXME
 * @Tm: FIXME
 * @XHeII: FIXME
 * @x: normalized scale factor inverse $x = 1 + z = a_0/a$
 *
 * $dX_\e/dx$ implemented using Weinbergs book
 *
 * Returns: FIXME
 */
gdouble
nc_recomb_seager_weinberg_HII_ion_rate (NcRecombSeager *recomb_seager, NcHICosmo *cosmo, gdouble XHII, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble Xe = XHII + XHeII;
  const gdouble alpha = nc_recomb_seager_pequignot_HI_case_B (recomb_seager, cosmo, Tm);
  const gdouble h2 = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
  const gdouble n_b0 = Omega_b * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0 = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / ncm_c_kpc ();
  const gdouble f1 = alpha * n_0 * x * x / H;

  const gdouble x3 = gsl_pow_3 (x);
  const gdouble Tm3_2 = sqrt (gsl_pow_3 (Tm));
  const gdouble f2na = ncm_c_decay_H_rate_2s_1s () * (1.0 - XHII);
  const gdouble f2nb = H / x3 / (n_0 * ncm_c_H_Lyman_2p_wl3_8pi ());
  const gdouble f2n = f2na + f2nb;
  const gdouble f2da = ncm_c_boltzmann_factor_H_2s (Tm) * Tm3_2 * alpha * (1.0 - XHII);
  const gdouble f2d = f2na + f2nb + f2da;
  const gdouble f2 = f2n / f2d;

  const gdouble S = ncm_c_boltzmann_factor_H_1s  (Tm) * Tm3_2 / (n_0 * x3);
  const gdouble f3 = (XHII * Xe - (1.0 - XHII) * S);

  return f1 * f2 * f3;
}
