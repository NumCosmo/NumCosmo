/***************************************************************************
 *            nc_recomb.h
 *
 *  Sun Oct  5 20:40:46 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_RECOMB_H_
#define _NC_RECOMB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_util.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/function_cache.h>

#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#endif

G_BEGIN_DECLS

#define NC_TYPE_RECOMB             (nc_recomb_get_type ())
#define NC_RECOMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_RECOMB, NcRecomb))
#define NC_RECOMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_RECOMB, NcRecombClass))
#define NC_IS_RECOMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_RECOMB))
#define NC_IS_RECOMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_RECOMB))
#define NC_RECOMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_RECOMB, NcRecombClass))

typedef struct _NcRecombClass NcRecombClass;
typedef struct _NcRecomb NcRecomb;

struct _NcRecomb
{
  /*< private >*/
  GObject parent_instance;
  gdouble zi, lambdai, lambdaf, prec;
  gdouble init_frac;
  gsl_min_fminimizer *fmin;
  gsl_root_fsolver *fsol;
  NcmSpline *Xe_s;
  NcmSpline *dtau_dlambda_s;
  NcmSpline *tau_s;
  NcmModelCtrl *ctrl;
};

struct _NcRecombClass
{
  /*< private >*/
  GObjectClass parent_class;
  void (*prepare) (NcRecomb *recomb, NcHICosmo *cosmo);
};

GType nc_recomb_get_type (void) G_GNUC_CONST;

NcRecomb *nc_recomb_new_from_name (const gchar *recomb_name);
NcRecomb *nc_recomb_ref (NcRecomb *recomb);
void nc_recomb_free (NcRecomb *recomb);
void nc_recomb_clear (NcRecomb **recomb);
void nc_recomb_prepare (NcRecomb *recomb, NcHICosmo *cosmo);
G_INLINE_FUNC void nc_recomb_prepare_if_needed (NcRecomb *recomb, NcHICosmo *cosmo);

gdouble nc_recomb_HI_ion_saha (NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_HeI_ion_saha (NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_HeII_ion_saha (NcHICosmo *cosmo, const gdouble x);

gdouble nc_recomb_HeII_ion_saha_x (NcHICosmo *cosmo, const gdouble f);
gdouble nc_recomb_HeII_ion_saha_x_by_HeIII_He (NcHICosmo *cosmo, const gdouble f);
gdouble nc_recomb_He_fully_ionized_Xe (NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);

gdouble nc_recomb_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_dtau_dx (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_dtau_dlambda (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_d2tau_dlambda2 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_d3tau_dlambda3 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_tau_lambda0_lambda1 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda0, const gdouble lambda1);
gdouble nc_recomb_log_v_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_v_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_dv_tau_dlambda (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_d2v_tau_dlambda2 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);

gdouble nc_recomb_v_tau_lambda_mode (NcRecomb *recomb, NcHICosmo *cosmo);
void nc_recomb_v_tau_lambda_features (NcRecomb *recomb, NcHICosmo *cosmo, gdouble logref, gdouble *lambda_max, gdouble *lambda_l, gdouble *lambda_u);

gdouble nc_recomb_tau_zstar (NcRecomb *recomb, NcHICosmo *cosmo);
gdouble nc_recomb_tau_cutoff (NcRecomb *recomb, NcHICosmo *cosmo);
gdouble nc_recomb_tau_zdrag (NcRecomb *recomb, NcHICosmo *cosmo);

G_INLINE_FUNC gdouble nc_recomb_dtau_dlambda_Xe (NcHICosmo *cosmo, const gdouble lambda);
G_INLINE_FUNC gdouble nc_recomb_He_fully_ionized_dtau_dlambda (NcHICosmo *cosmo, const gdouble lambda);

G_INLINE_FUNC gdouble nc_recomb_pequignot_HI_case_B (NcHICosmo *cosmo, const gdouble Tm);
G_INLINE_FUNC gdouble nc_recomb_pequignot_HI_case_B_dTm (NcHICosmo *cosmo, const gdouble Tm);
G_INLINE_FUNC gdouble nc_recomb_hummer_HeI_case_B (NcHICosmo *cosmo, const gdouble Tm);
G_INLINE_FUNC gdouble nc_recomb_hummer_HeI_case_B_dTm (NcHICosmo *cosmo, const gdouble Tm);

/* Internal use */
void _nc_recomb_prepare_tau_splines (NcRecomb *recomb, NcHICosmo *cosmo);

G_END_DECLS

#endif /* _NC_RECOMB_H_ */

#ifndef _NC_RECOMB_INLINE_H_
#define _NC_RECOMB_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC void
nc_recomb_prepare_if_needed (NcRecomb *recomb, NcHICosmo *cosmo)
{
  if (ncm_model_ctrl_update (recomb->ctrl, NCM_MODEL (cosmo)))
	nc_recomb_prepare (recomb, cosmo);
}

G_INLINE_FUNC gdouble
nc_recomb_dtau_dlambda_Xe (NcHICosmo *cosmo, const gdouble lambda)
{
	const gdouble x = exp (-lambda);
  const gdouble x3 = gsl_pow_3 (x);
  const gdouble h2 = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
  const gdouble n_b0 = Omega_b * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0 = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / (ncm_c_kpc ());

  return -ncm_c_c () * ncm_c_thomson_cs () * n_0 * x3 / H;
}

G_INLINE_FUNC gdouble
nc_recomb_He_fully_ionized_dtau_dlambda (NcHICosmo *cosmo, const gdouble lambda)
{
	const gdouble x = exp (-lambda);
	const gdouble x3 = gsl_pow_3 (x);
	const gdouble h2 = nc_hicosmo_h2 (cosmo);
	const gdouble Omega_b = nc_hicosmo_Omega_b (cosmo);
	const gdouble n_b0 = Omega_b * ncm_c_crit_number_density_p () * h2;
	const gdouble n_0 = ncm_c_prim_H_Yp () * n_b0;
  const gdouble H = nc_hicosmo_H (cosmo, x - 1.0) / (ncm_c_kpc ());
  const gdouble Xe = nc_recomb_He_fully_ionized_Xe (cosmo, x);

	return -Xe * ncm_c_c () * ncm_c_thomson_cs () * n_0 * x3 / H;
}

G_INLINE_FUNC gdouble
nc_recomb_pequignot_HI_case_B (NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble F =  1.14;   /* fudge factor */
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

G_INLINE_FUNC gdouble
nc_recomb_pequignot_HI_case_B_dTm (NcHICosmo *cosmo, gdouble Tm)
{
  const gdouble F =  1.14;   /* fudge factor */
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

G_INLINE_FUNC gdouble
nc_recomb_hummer_HeI_case_B (NcHICosmo *cosmo, gdouble Tm)
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

G_INLINE_FUNC gdouble
nc_recomb_hummer_HeI_case_B_dTm (NcHICosmo *cosmo, gdouble Tm)
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

G_INLINE_FUNC gdouble
nc_recomb_weinberg_HII_ion_rate (NcHICosmo *cosmo, gdouble XHII, gdouble Tm, gdouble XHeII, gdouble x)
{
  const gdouble Xe = XHII + XHeII;
  const gdouble alpha = nc_recomb_pequignot_HI_case_B (cosmo, Tm);
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

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_RECOMB_INLINE_H_ */

