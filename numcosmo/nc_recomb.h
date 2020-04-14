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
#include <numcosmo/nc_hireion.h>
#include <numcosmo/math/ncm_util.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_ode_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_RECOMB             (nc_recomb_get_type ())
#define NC_RECOMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_RECOMB, NcRecomb))
#define NC_RECOMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_RECOMB, NcRecombClass))
#define NC_IS_RECOMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_RECOMB))
#define NC_IS_RECOMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_RECOMB))
#define NC_RECOMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_RECOMB, NcRecombClass))

typedef struct _NcRecombClass NcRecombClass;
typedef struct _NcRecomb NcRecomb;

struct _NcRecombClass
{
  /*< private >*/
  GObjectClass parent_class;
  void (*prepare) (NcRecomb *recomb, NcHICosmo *cosmo);
  gdouble (*Xe) (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
  gdouble (*XHII) (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
  gdouble (*XHeII) (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda); 
};

struct _NcRecomb
{
  /*< private >*/
  GObject parent_instance;
  gdouble zi, lambdai, lambdaf, prec;
  gdouble init_frac;
  gsl_min_fminimizer *fmin;
  gsl_root_fsolver *fsol;
  NcmSpline *dtau_dlambda_s;
  NcmSpline *tau_s;
  NcmOdeSpline *tau_ode_s;
  NcmOdeSpline *tau_drag_ode_s;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_reion;
  gdouble v_tau_max_z, v_tau_max_lambda;
  gdouble tau_z, tau_lambda;
  gdouble tau_drag_z, tau_drag_lambda;
  gdouble tau_cutoff_z, tau_cutoff_lambda;
};

GType nc_recomb_get_type (void) G_GNUC_CONST;

NcRecomb *nc_recomb_new_from_name (const gchar *recomb_name);
NcRecomb *nc_recomb_ref (NcRecomb *recomb);
void nc_recomb_free (NcRecomb *recomb);
void nc_recomb_clear (NcRecomb **recomb);

void nc_recomb_prepare (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE void nc_recomb_prepare_if_needed (NcRecomb *recomb, NcHICosmo *cosmo);

void nc_recomb_set_zi (NcRecomb *recomb, const gdouble zi);
void nc_recomb_require_zi (NcRecomb *recomb, const gdouble zi);
gdouble nc_recomb_get_zi (NcRecomb *recomb);

NCM_INLINE gdouble nc_recomb_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
NCM_INLINE gdouble nc_recomb_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
NCM_INLINE gdouble nc_recomb_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);

gdouble nc_recomb_HI_ion_saha (NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_HeI_ion_saha (NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_HeII_ion_saha (NcHICosmo *cosmo, const gdouble x);

gdouble nc_recomb_HeII_ion_saha_x (NcHICosmo *cosmo, const gdouble f);
gdouble nc_recomb_HeII_ion_saha_x_by_HeIII_He (NcHICosmo *cosmo, const gdouble f);
gdouble nc_recomb_He_fully_ionized_Xe (NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_XHI (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_XHeI (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);
gdouble nc_recomb_equilibrium_XHeIII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble x);

gdouble nc_recomb_dtau_dx (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_dtau_dlambda (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_d2tau_dlambda2 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_d3tau_dlambda3 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_tau_drag (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_tau_lambda0_lambda1 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda0, const gdouble lambda1);
gdouble nc_recomb_log_v_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_v_tau (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_dv_tau_dlambda (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);
gdouble nc_recomb_d2v_tau_dlambda2 (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda);

void nc_recomb_v_tau_lambda_features (NcRecomb *recomb, NcHICosmo *cosmo, gdouble logref, gdouble *lambda_max, gdouble *lambda_l, gdouble *lambda_u);

NCM_INLINE gdouble nc_recomb_get_v_tau_max_lambda (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE gdouble nc_recomb_get_tau_lambda (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE gdouble nc_recomb_get_tau_drag_lambda (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE gdouble nc_recomb_get_tau_cutoff_lambda (NcRecomb *recomb, NcHICosmo *cosmo);

NCM_INLINE gdouble nc_recomb_get_v_tau_max_z (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE gdouble nc_recomb_get_tau_z (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE gdouble nc_recomb_get_tau_drag_z (NcRecomb *recomb, NcHICosmo *cosmo);
NCM_INLINE gdouble nc_recomb_get_tau_cutoff_z (NcRecomb *recomb, NcHICosmo *cosmo);

NCM_INLINE gdouble nc_recomb_dtau_dlambda_Xe (NcHICosmo *cosmo, const gdouble lambda);
NCM_INLINE gdouble nc_recomb_He_fully_ionized_dtau_dlambda (NcHICosmo *cosmo, const gdouble lambda);

/* Internal use */
void _nc_recomb_prepare_tau_splines (NcRecomb *recomb, NcHICosmo *cosmo);
void _nc_recomb_prepare_redshifts (NcRecomb *recomb, NcHICosmo *cosmo);

#define NC_RECOMB_STARTING_X (1.0e12)

G_END_DECLS

#endif /* _NC_RECOMB_H_ */

#ifndef _NC_RECOMB_INLINE_H_
#define _NC_RECOMB_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
nc_recomb_prepare_if_needed (NcRecomb *recomb, NcHICosmo *cosmo)
{
  gboolean cosmo_up = ncm_model_ctrl_update (recomb->ctrl_cosmo, NCM_MODEL (cosmo));

  if (cosmo_up)
    nc_recomb_prepare (recomb, cosmo);
}

NCM_INLINE gdouble 
nc_recomb_Xe (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  return NC_RECOMB_GET_CLASS (recomb)->Xe (recomb, cosmo, lambda);
}

NCM_INLINE gdouble 
nc_recomb_XHII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  return NC_RECOMB_GET_CLASS (recomb)->XHII (recomb, cosmo, lambda);
}

NCM_INLINE gdouble 
nc_recomb_XHeII (NcRecomb *recomb, NcHICosmo *cosmo, const gdouble lambda)
{
  return NC_RECOMB_GET_CLASS (recomb)->XHeII (recomb, cosmo, lambda);
}

NCM_INLINE gdouble
nc_recomb_dtau_dlambda_Xe (NcHICosmo *cosmo, const gdouble lambda)
{
	const gdouble x        = exp (-lambda);
  const gdouble x3       = gsl_pow_3 (x);
  const gdouble h2       = nc_hicosmo_h2 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble n_b0     = Omega_b0 * ncm_c_crit_number_density_p () * h2;
  const gdouble n_0      = nc_hicosmo_Yp_1H (cosmo) * n_b0;
  const gdouble H        = nc_hicosmo_H (cosmo, x - 1.0) / (ncm_c_kpc ());

  return -ncm_c_c () * ncm_c_thomson_cs () * n_0 * x3 / H;
}

NCM_INLINE gdouble
nc_recomb_He_fully_ionized_dtau_dlambda (NcHICosmo *cosmo, const gdouble lambda)
{
	const gdouble x        = exp (-lambda);
	const gdouble x3       = gsl_pow_3 (x);
	const gdouble h2       = nc_hicosmo_h2 (cosmo);
	const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
	const gdouble n_b0     = Omega_b0 * ncm_c_crit_number_density_p () * h2;
	const gdouble n_0      = nc_hicosmo_Yp_1H (cosmo) * n_b0;
  const gdouble H        = nc_hicosmo_H (cosmo, x - 1.0) / (ncm_c_kpc ());
  const gdouble Xe       = nc_recomb_He_fully_ionized_Xe (cosmo, x);

	return -Xe * ncm_c_c () * ncm_c_thomson_cs () * n_0 * x3 / H;
}

NCM_INLINE gdouble 
nc_recomb_get_v_tau_max_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->v_tau_max_lambda;
}

NCM_INLINE gdouble 
nc_recomb_get_tau_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->tau_lambda;
}

NCM_INLINE gdouble 
nc_recomb_get_tau_drag_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->tau_drag_lambda;
}

NCM_INLINE gdouble 
nc_recomb_get_tau_cutoff_lambda (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->tau_cutoff_lambda;
}

NCM_INLINE gdouble 
nc_recomb_get_v_tau_max_z (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->v_tau_max_z;
}

NCM_INLINE gdouble 
nc_recomb_get_tau_z (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->tau_z;
}

NCM_INLINE gdouble 
nc_recomb_get_tau_drag_z (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->tau_drag_z;
}

NCM_INLINE gdouble 
nc_recomb_get_tau_cutoff_z (NcRecomb *recomb, NcHICosmo *cosmo)
{
  return recomb->tau_cutoff_z;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_RECOMB_INLINE_H_ */
