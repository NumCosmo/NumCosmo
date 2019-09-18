/***************************************************************************
 *            nc_hicosmo.h
 *
 *  Mon Jul 16 18:03:42 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_HICOSMO_H_
#define _NC_HICOSMO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/math/ncm_powspec_filter.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO             (nc_hicosmo_get_type ())
#define NC_HICOSMO(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO, NcHICosmo))
#define NC_HICOSMO_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO, NcHICosmoClass))
#define NC_IS_HICOSMO(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO))
#define NC_IS_HICOSMO_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO))
#define NC_HICOSMO_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO, NcHICosmoClass))

/**
 * NcHICosmoImpl:
 * @NC_HICOSMO_IMPL_H0: Hubble constant
 * @NC_HICOSMO_IMPL_Omega_b0: Baryonic density today $\Omega_{b0} = \rho_{b0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_c0: Cold dark matter density today $\Omega_{c0} = \rho_{c0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_g0: Photons density today $\Omega_{\gamma0} = \rho_{\gamma0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_nu0: Ultra-relativistic neutrinos density today $\Omega_{\nu0} = \rho_{\nu0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_mnu0: Massive neutrinos density today $\Omega_{m\nu0} = \rho_{m\nu0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Press_mnu0: Massive neutrinos dimensionless pressure today $P_{m\nu0} = p_{m\nu0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_mnu0_n: The n-th massive neutrinos density today $\Omega_{m\nu0,n} = \rho_{m\nu0,n} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Press_mnu0_n: The n-th massive neutrinos dimensionless pressure today $P_{m\nu0,n} = p_{m\nu0,n} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_m0: Total matter density today $\Omega_{m0} = \rho_{m0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_r0: Total radiation density today $\Omega_{r0} = \rho_{r0} / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_Omega_t0: Total density today $\Omega_{t0}$
 * @NC_HICOSMO_IMPL_T_gamma0: Radiation temperature today
 * @NC_HICOSMO_IMPL_Yp_4He: Primordial Helium mass fraction 
 * @NC_HICOSMO_IMPL_z_lss: Redshift of the last scatering surface
 * @NC_HICOSMO_IMPL_as_drag: Acoustic Scale at drag redshift
 * @NC_HICOSMO_IMPL_xb: Maximum redshift
 * @NC_HICOSMO_IMPL_E2Omega_b: Baryonic density $E^2\Omega_{b} = \rho_b(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_c: Cold dark matter density $E^2\Omega_{c} = \rho_c(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_g: Photons density $E^2\Omega_{\gamma} = \rho_\gamma(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_nu: Ultra-relativistic neutrinos density $E^2\Omega_{\nu} = \rho_\nu(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_mnu: Massive neutrinos density $E^2\Omega_{m\nu} = \rho_{m\nu}(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Press_mnu: Massive neutrinos pressure $E^2P_{m\nu} = p_{m\nu}(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_mnu_n: The n-th Massive neutrinos density $E^2\Omega_{m\nu,n} = \rho_{m\nu,n}(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Press_mnu_n: The n-th Massive neutrinos pressure $E^2P_{m\nu,n} = p_{m\nu,n}(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_m: Total matter density $E^2\Omega_{m} = \rho_m(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_r: Total radiation density $\Omega_{r} = \rho_r(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2Omega_t: Total density $E2\Omega_{t0} = \rho_t(z) / \rho_{\mathrm{crit}0}$
 * @NC_HICOSMO_IMPL_E2: Dimensionless Hubble function squared $H^2(z) / H_0^2$
 * @NC_HICOSMO_IMPL_dE2_dz: Derivative of the dimensionless Hubble function squared.
 * @NC_HICOSMO_IMPL_d2E2_dz2: Second derivative of the dimensionless Hubble function squared.
 * @NC_HICOSMO_IMPL_bgp_cs2: Baryon-photon plasma speed of sound squared $c_s^2$.
 * @NC_HICOSMO_IMPL_Dc: Comoving distance
 * @NC_HICOSMO_IMPL_NMassNu: Number of massive neutrinos
 * @NC_HICOSMO_IMPL_MassNuInfo: Massive neutrino info
 * @NC_HICOSMO_IMPL_get_bg_var: Background variables interface for perturbations
 *
 * Flags defining the implementation options of the NcHICosmo abstract object. 
 * 
 */
typedef enum /*< flags,underscore_name=NC_HICOSMO_IMPL >*/
{
  NC_HICOSMO_IMPL_H0 = 0, 
  NC_HICOSMO_IMPL_Omega_b0,
  NC_HICOSMO_IMPL_Omega_c0,
  NC_HICOSMO_IMPL_Omega_g0,
  NC_HICOSMO_IMPL_Omega_nu0,
  NC_HICOSMO_IMPL_Omega_mnu0,
  NC_HICOSMO_IMPL_Press_mnu0,
  NC_HICOSMO_IMPL_Omega_mnu0_n,
  NC_HICOSMO_IMPL_Press_mnu0_n,
  NC_HICOSMO_IMPL_Omega_m0,
  NC_HICOSMO_IMPL_Omega_r0,
  NC_HICOSMO_IMPL_Omega_t0,
  NC_HICOSMO_IMPL_T_gamma0,
  NC_HICOSMO_IMPL_Yp_4He,
  NC_HICOSMO_IMPL_z_lss,
  NC_HICOSMO_IMPL_as_drag,
  NC_HICOSMO_IMPL_xb,
  NC_HICOSMO_IMPL_E2Omega_b,
  NC_HICOSMO_IMPL_E2Omega_c,
  NC_HICOSMO_IMPL_E2Omega_g,
  NC_HICOSMO_IMPL_E2Omega_nu,
  NC_HICOSMO_IMPL_E2Omega_mnu,
  NC_HICOSMO_IMPL_E2Press_mnu,
  NC_HICOSMO_IMPL_E2Omega_mnu_n,
  NC_HICOSMO_IMPL_E2Press_mnu_n,
  NC_HICOSMO_IMPL_E2Omega_m,
  NC_HICOSMO_IMPL_E2Omega_r,
  NC_HICOSMO_IMPL_E2Omega_t,
  NC_HICOSMO_IMPL_E2,
  NC_HICOSMO_IMPL_dE2_dz,
  NC_HICOSMO_IMPL_d2E2_dz2,
  NC_HICOSMO_IMPL_bgp_cs2,
  NC_HICOSMO_IMPL_Dc,
  NC_HICOSMO_IMPL_NMassNu,
  NC_HICOSMO_IMPL_MassNuInfo,
  NC_HICOSMO_IMPL_get_bg_var,
  /* < private > */
  NC_HICOSMO_IMPL_LAST,       /*< skip >*/
} NcHICosmoImpl;

#define NC_HICOSMO_IMPL_FLAG_RH_Mpc       NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_H0)
#define NC_HICOSMO_IMPL_FLAG_RH_planck    NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_H0)
#define NC_HICOSMO_IMPL_FLAG_Omega_k0     NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_Omega_t0)
#define NC_HICOSMO_IMPL_FLAG_h            NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_H0)
#define NC_HICOSMO_IMPL_FLAG_h2           NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_H0)
#define NC_HICOSMO_IMPL_FLAG_Omega_b0h2   NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_b0,   NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_Omega_c0h2   NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_c0,   NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_Omega_g0h2   NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_g0,   NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_Omega_nu0h2  NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_nu0,  NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_Omega_mnu0h2 NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_mnu0, NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_Omega_m0h2   NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_m0,   NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_Omega_r0h2   NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_Omega_r0,   NC_HICOSMO_IMPL_h2)
#define NC_HICOSMO_IMPL_FLAG_H_Yp         NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_Yp_4He)
#define NC_HICOSMO_IMPL_FLAG_XHe          NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_Yp_4He)

#define NC_HICOSMO_IMPL_FLAG_H        NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_H0, NC_HICOSMO_IMPL_E2)
#define NC_HICOSMO_IMPL_FLAG_dH_dz    NCM_MODEL_3OPT2IMPL (NC_HICOSMO_IMPL_H0, NC_HICOSMO_IMPL_E2, NC_HICOSMO_IMPL_dE2_dz)
#define NC_HICOSMO_IMPL_FLAG_E        NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_E2)
#define NC_HICOSMO_IMPL_FLAG_Em2      NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_E2)
#define NC_HICOSMO_IMPL_FLAG_q        NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_E2, NC_HICOSMO_IMPL_dE2_dz)
#define NC_HICOSMO_IMPL_FLAG_j        NCM_MODEL_3OPT2IMPL (NC_HICOSMO_IMPL_E2, NC_HICOSMO_IMPL_dE2_dz, NC_HICOSMO_IMPL_d2E2_dz2)
#define NC_HICOSMO_IMPL_FLAG_Omega_k0 NCM_MODEL_OPT2IMPL  (NC_HICOSMO_IMPL_Omega_t0)
#define NC_HICOSMO_IMPL_FLAG_wec      NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_E2, NC_HICOSMO_IMPL_Omega_k0)
#define NC_HICOSMO_IMPL_FLAG_dec      NCM_MODEL_2OPT2IMPL (NC_HICOSMO_IMPL_E2, NC_HICOSMO_IMPL_Omega_k0)

typedef struct _NcHICosmoClass NcHICosmoClass;
typedef struct _NcHICosmo NcHICosmo;

typedef gdouble (*NcHICosmoFunc0) (NcHICosmo *cosmo);
typedef gdouble (*NcHICosmoFunc1Z) (NcHICosmo *cosmo, const gdouble z);
typedef gdouble (*NcHICosmoFunc1K) (NcHICosmo *cosmo, const gdouble k);

typedef gdouble (*NcHICosmoVFunc0) (NcHICosmo *cosmo, const guint n);
typedef gdouble (*NcHICosmoVFunc1Z) (NcHICosmo *cosmo, const guint n, const gdouble z);
typedef gdouble (*NcHICosmoVFunc1K) (NcHICosmo *cosmo, const guint n, const gdouble k);

typedef guint (*NcHICosmoFuncNMassNu) (NcHICosmo *cosmo);
typedef void (*NcHICosmoFuncMassNuInfo) (NcHICosmo *cosmo, const guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *xi, gdouble *g);

#ifndef __GTK_DOC_IGNORE__
typedef struct _NcHIPrim NcHIPrim;
typedef struct _NcHIReion NcHIReion;
typedef struct _NcHIPertBGVar NcHIPertBGVar;
#endif

typedef void (*NcHICosmoGetBGVar) (NcHICosmo *cosmo, const gdouble t, NcHIPertBGVar *bg_var);

typedef struct _NcHICosmoFuncZ
{
  const gchar *name;
  const gchar *desc;
  NcHICosmoFunc1Z f;
  NcHICosmoImpl impl;
} NcHICosmoFuncZ;

typedef struct _NcHICosmoFunc
{
  const gchar *name;
  const gchar *desc;
  NcHICosmoFunc0 f;
  NcHICosmoImpl impl;
} NcHICosmoFunc;

struct _NcHICosmoClass
{
  /*< private >*/
  NcmModelClass parent_class;
  NcHICosmoFunc0   H0;
  NcHICosmoFunc0   Omega_b0;
  NcHICosmoFunc0   Omega_c0;
  NcHICosmoFunc0   Omega_g0;
  NcHICosmoFunc0   Omega_nu0;
  NcHICosmoFunc0   Omega_mnu0;
  NcHICosmoFunc0   Press_mnu0;
  NcHICosmoFunc0   Omega_m0;
  NcHICosmoFunc0   Omega_r0;
  NcHICosmoFunc0   Omega_t0;
  NcHICosmoFunc0   T_gamma0;
  NcHICosmoFunc0   Yp_4He;
  NcHICosmoFunc0   z_lss;
  NcHICosmoFunc0   as_drag;
  NcHICosmoFunc0   xb;
  NcHICosmoVFunc0  Omega_mnu0_n;
  NcHICosmoVFunc0  Press_mnu0_n;
  NcHICosmoFunc1Z  E2Omega_b;
  NcHICosmoFunc1Z  E2Omega_c;
  NcHICosmoFunc1Z  E2Omega_g;
  NcHICosmoFunc1Z  E2Omega_nu;
  NcHICosmoFunc1Z  E2Omega_mnu;
  NcHICosmoFunc1Z  E2Press_mnu;
  NcHICosmoFunc1Z  E2Omega_m;
  NcHICosmoFunc1Z  E2Omega_r;
  NcHICosmoFunc1Z  E2Omega_t;
  NcHICosmoFunc1Z  E2;
  NcHICosmoFunc1Z  dE2_dz;
  NcHICosmoFunc1Z  d2E2_dz2;
  NcHICosmoFunc1Z  bgp_cs2;
  NcHICosmoFunc1Z  Dc;
  NcHICosmoVFunc1Z E2Omega_mnu_n;
  NcHICosmoVFunc1Z E2Press_mnu_n;
  NcHICosmoFuncNMassNu NMassNu;
  NcHICosmoFuncMassNuInfo MassNuInfo;
  NcHICosmoGetBGVar get_bg_var;
};

/**
 * NcHICosmo:
 *
 * FIXME
 *
 */
struct _NcHICosmo
{
  /*< private >*/
  NcmModel parent_instance;
  gboolean is_eternal;
  NcHIPrim *prim;
  NcHIReion *reion;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  const gsl_min_fminimizer_type *Tmin;
  gsl_min_fminimizer *smin;
};

GType nc_hicosmo_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_hicosmo);

void nc_hicosmo_set_H0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_b0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_c0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_g0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_nu0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_mnu0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Press_mnu0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_mnu0_n_impl (NcHICosmoClass *model_class, NcHICosmoVFunc0 f);
void nc_hicosmo_set_Press_mnu0_n_impl (NcHICosmoClass *model_class, NcHICosmoVFunc0 f);
void nc_hicosmo_set_Omega_m0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_r0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Omega_t0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_T_gamma0_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_Yp_4He_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_z_lss_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_as_drag_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);
void nc_hicosmo_set_xb_impl (NcHICosmoClass *model_class, NcHICosmoFunc0 f);

void nc_hicosmo_set_E2Omega_b_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_c_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_g_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_nu_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_mnu_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Press_mnu_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_mnu_n_impl (NcHICosmoClass *model_class, NcHICosmoVFunc1Z f);
void nc_hicosmo_set_E2Press_mnu_n_impl (NcHICosmoClass *model_class, NcHICosmoVFunc1Z f);
void nc_hicosmo_set_E2Omega_m_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_r_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_E2Omega_t_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);

void nc_hicosmo_set_E2_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_dE2_dz_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_d2E2_dz2_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_bgp_cs2_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_Dc_impl (NcHICosmoClass *model_class, NcHICosmoFunc1Z f);
void nc_hicosmo_set_NMassNu_impl (NcHICosmoClass *model_class, NcHICosmoFuncNMassNu f);
void nc_hicosmo_set_MassNuInfo_impl (NcHICosmoClass *model_class, NcHICosmoFuncMassNuInfo f);

void nc_hicosmo_set_get_bg_var_impl (NcHICosmoClass *model_class, NcHICosmoGetBGVar f);

NcHICosmo *nc_hicosmo_new_from_name (GType parent_type, gchar *cosmo_name);
NcHICosmo *nc_hicosmo_ref (NcHICosmo *cosmo);
void nc_hicosmo_free (NcHICosmo *cosmo);
void nc_hicosmo_clear (NcHICosmo **cosmo);

void nc_hicosmo_log_all_models (GType parent);
gdouble nc_hicosmo_zt (NcHICosmo *cosmo, const gdouble z_max);

void nc_hicosmo_mqE2_max (NcHICosmo *cosmo, const gdouble z_max, gdouble *zm, gdouble *mqE2m);
void nc_hicosmo_dec_min (NcHICosmo *cosmo, const gdouble z_max, gdouble *zm, gdouble *decm);
void nc_hicosmo_q_min (NcHICosmo *cosmo, const gdouble z_max, gdouble *zm, gdouble *qm);

/*
 * Cosmological model constant functions
 */
NCM_INLINE gdouble nc_hicosmo_H0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_RH_Mpc (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_RH_planck (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_h (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_b0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_c0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_g0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_nu0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_mnu0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Press_mnu0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_mnu0_n (NcHICosmo *cosmo, const guint n);
NCM_INLINE gdouble nc_hicosmo_Press_mnu0_n (NcHICosmo *cosmo, const guint n);
NCM_INLINE gdouble nc_hicosmo_Omega_m0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_r0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_t0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_k0 (NcHICosmo *cosmo);

NCM_INLINE gdouble nc_hicosmo_Omega_b0h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_c0h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_g0h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_nu0h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_mnu0h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_m0h2 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Omega_r0h2 (NcHICosmo *cosmo);

NCM_INLINE gdouble nc_hicosmo_T_gamma0 (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Yp_4He (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_Yp_1H (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_XHe (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_crit_density (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_baryon_density (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_He_number_density (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_H_number_density (NcHICosmo *cosmo);

NCM_INLINE gdouble nc_hicosmo_z_lss (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_as_drag (NcHICosmo *cosmo);
NCM_INLINE gdouble nc_hicosmo_xb (NcHICosmo *cosmo);

NCM_INLINE gdouble nc_hicosmo_E2Omega_b (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_c (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_g (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_nu (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_mnu (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Press_mnu (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Press_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_m (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_r (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2Omega_t (NcHICosmo *cosmo, const gdouble z);

NCM_INLINE gdouble nc_hicosmo_H (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_dH_dz (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_E2 (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_Em2 (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_dE2_dz (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_d2E2_dz2 (NcHICosmo *cosmo, const gdouble z);

NCM_INLINE gdouble nc_hicosmo_bgp_cs2 (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_Dc (NcHICosmo *cosmo, const gdouble z);

NCM_INLINE guint nc_hicosmo_NMassNu (NcHICosmo *cosmo);
NCM_INLINE void nc_hicosmo_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *xi, gdouble *g);
NCM_INLINE gdouble nc_hicosmo_Neff (NcHICosmo *cosmo);

NCM_INLINE void nc_hicosmo_get_bg_var (NcHICosmo *cosmo, const gdouble t, NcHIPertBGVar *bg_var);

NCM_INLINE gdouble nc_hicosmo_E2Omega_k (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_q (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_nec (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_dec (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_wec (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_qp (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_j (NcHICosmo *cosmo, const gdouble z);
NCM_INLINE gdouble nc_hicosmo_kinetic_w (NcHICosmo *cosmo, const gdouble z);

NCM_INLINE gdouble nc_hicosmo_mqE2 (NcHICosmo *cosmo, const gdouble z);

NCM_INLINE gdouble nc_hicosmo_abs_alpha (NcHICosmo *cosmo, gdouble x);
NCM_INLINE gdouble nc_hicosmo_x_alpha (NcHICosmo *cosmo, gdouble alpha);

NCM_INLINE NcHIPrim *nc_hicosmo_peek_prim (NcHICosmo *cosmo);
NCM_INLINE NcHIReion *nc_hicosmo_peek_reion (NcHICosmo *cosmo);

gdouble nc_hicosmo_sigma8 (NcHICosmo *cosmo, NcmPowspecFilter *psf);

#define NC_HICOSMO_DEFAULT_PARAMS_RELTOL (1e-7)
#define NC_HICOSMO_DEFAULT_PARAMS_ABSTOL (0.0)
#define NC_HICOSMO_OMEGA_K0_LIMIT (1.0e-13)

G_END_DECLS

#endif /* _NC_HICOSMO_H_ */

#ifndef _NC_HICOSMO_INLINE_H_
#define _NC_HICOSMO_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,H0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_b0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_c0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_g0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_nu0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_mnu0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Press_mnu0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_m0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_r0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_t0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,T_gamma0)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Yp_4He)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,z_lss)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,as_drag)
NCM_MODEL_FUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,xb)

NCM_MODEL_VFUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Omega_mnu0_n)
NCM_MODEL_VFUNC0_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Press_mnu0_n)

NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_b,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_c,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_g,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_nu,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_mnu,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Press_mnu,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_m,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_r,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_t,z)

NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,dE2_dz,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,d2E2_dz2,z)

NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,bgp_cs2,z)
NCM_MODEL_FUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,Dc,z)

NCM_MODEL_VFUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Omega_mnu_n,z)
NCM_MODEL_VFUNC1_IMPL (NC_HICOSMO,NcHICosmo,nc_hicosmo,E2Press_mnu_n,z)

NCM_INLINE guint 
nc_hicosmo_NMassNu (NcHICosmo *cosmo)
{
  return NC_HICOSMO_GET_CLASS (cosmo)->NMassNu (cosmo);
}

NCM_INLINE void
nc_hicosmo_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *xi, gdouble *g)
{
  NC_HICOSMO_GET_CLASS (cosmo)->MassNuInfo (cosmo, nu_i, mass_eV, T_0, xi, g);
}

NCM_INLINE void 
nc_hicosmo_get_bg_var (NcHICosmo *cosmo, const gdouble t, NcHIPertBGVar *bg_var)
{
  NC_HICOSMO_GET_CLASS (cosmo)->get_bg_var (cosmo, t, bg_var);
}

NCM_INLINE gdouble
nc_hicosmo_RH_Mpc (NcHICosmo *cosmo)
{
  return (ncm_c_c () / (1.0e3 * nc_hicosmo_H0 (cosmo)));
}

NCM_INLINE gdouble
nc_hicosmo_RH_planck (NcHICosmo *cosmo)
{
  return nc_hicosmo_RH_Mpc (cosmo) * ncm_c_Mpc () / ncm_c_planck_length ();
}

NCM_INLINE gdouble
nc_hicosmo_Omega_k0 (NcHICosmo *cosmo)
{
  const gdouble Omega_k0 = (1.0 - nc_hicosmo_Omega_t0 (cosmo));
  return fabs (Omega_k0) < NC_HICOSMO_OMEGA_K0_LIMIT ? 0.0 : Omega_k0;
}

NCM_INLINE gdouble
nc_hicosmo_H (NcHICosmo *cosmo, const gdouble z)
{
  return (nc_hicosmo_H0 (cosmo) * sqrt (nc_hicosmo_E2 (cosmo, z)));
}

NCM_INLINE gdouble
nc_hicosmo_h (NcHICosmo *cosmo)
{
  return nc_hicosmo_H0 (cosmo) / 100.0;
}

NCM_INLINE gdouble
nc_hicosmo_h2 (NcHICosmo *cosmo)
{
  return gsl_pow_2 (nc_hicosmo_H0 (cosmo) / 100.0);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_b0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_b0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_c0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_c0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_g0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_g0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_nu0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_nu0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_mnu0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_mnu0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_m0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_m0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Omega_r0h2 (NcHICosmo *cosmo)
{
  return nc_hicosmo_h2 (cosmo) * nc_hicosmo_Omega_r0 (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_Yp_1H (NcHICosmo *cosmo)
{
  return 1.0 - nc_hicosmo_Yp_4He (cosmo);
}

NCM_INLINE gdouble
nc_hicosmo_XHe (NcHICosmo *cosmo)
{
  return nc_hicosmo_Yp_4He (cosmo) / (ncm_c_mass_ratio_4He_1H () * nc_hicosmo_Yp_1H (cosmo));
}

NCM_INLINE gdouble 
nc_hicosmo_crit_density (NcHICosmo *cosmo)
{
  const gdouble h2 = nc_hicosmo_h2 (cosmo);
  return ncm_c_crit_density_h2 () * h2;
}

NCM_INLINE gdouble 
nc_hicosmo_baryon_density (NcHICosmo *cosmo)
{
  const gdouble rho_crit = nc_hicosmo_crit_density (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble rho_b0   = Omega_b0 * rho_crit;

  return rho_b0;
}

NCM_INLINE gdouble 
nc_hicosmo_He_number_density (NcHICosmo *cosmo)
{
  const gdouble rho_b0 = nc_hicosmo_baryon_density (cosmo);
  const gdouble Yp_4He = nc_hicosmo_Yp_4He (cosmo);
  const gdouble E_4He  = ncm_c_rest_energy_4He ();

  return Yp_4He * rho_b0 / E_4He;
}

NCM_INLINE gdouble 
nc_hicosmo_H_number_density (NcHICosmo *cosmo)
{
  const gdouble rho_b0 = nc_hicosmo_baryon_density (cosmo);
  const gdouble Yp_1H  = nc_hicosmo_Yp_1H (cosmo);
  const gdouble E_1H   = ncm_c_rest_energy_1H ();

  return Yp_1H * rho_b0 / E_1H;
}

NCM_INLINE gdouble
nc_hicosmo_E (NcHICosmo *cosmo, const gdouble z)
{
  return sqrt (nc_hicosmo_E2 (cosmo, z));
}

NCM_INLINE gdouble
nc_hicosmo_Em2 (NcHICosmo *cosmo, const gdouble z)
{
  return 1.0 / nc_hicosmo_E2 (cosmo, z);
}

NCM_INLINE gdouble
nc_hicosmo_dH_dz (NcHICosmo *cosmo, const gdouble z)
{
  return nc_hicosmo_H0 (cosmo) *
	nc_hicosmo_dE2_dz (cosmo, z) / (2.0 * sqrt (nc_hicosmo_E2 (cosmo, z)));
}

NCM_INLINE gdouble 
nc_hicosmo_Neff (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  const gdouble z    = (1.0e3 * ncm_c_eV () / ncm_c_kb ()) / nc_hicosmo_T_gamma0 (cosmo) - 1.0;
  const gdouble Neff = (nc_hicosmo_E2Omega_nu (cosmo, z) + 3.0 * nc_hicosmo_E2Press_mnu (cosmo, z)) / (nc_hicosmo_E2Omega_g (cosmo, z) * conv);

  return Neff;
}

NCM_INLINE gdouble 
nc_hicosmo_E2Omega_k (NcHICosmo *cosmo, const gdouble z)
{
  return nc_hicosmo_Omega_k0 (cosmo) * gsl_pow_2 (1.0 + z);
}

NCM_INLINE gdouble
nc_hicosmo_q (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble E2     = nc_hicosmo_E2 (cosmo, z);
  const gdouble dE2_dz = nc_hicosmo_dE2_dz (cosmo, z);
  
  return (dE2_dz * (1.0 + z) / (2.0 * E2) - 1.0);
}

NCM_INLINE gdouble
nc_hicosmo_nec (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x  = 1.0 + z;
  const gdouble q  = nc_hicosmo_q (cosmo, z);
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
  const gdouble Ok = nc_hicosmo_Omega_k0 (cosmo);

  return (1.0 + q) * E2 - Ok * x * x;
}

NCM_INLINE gdouble
nc_hicosmo_dec (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x  = 1.0 + z;
  const gdouble q  = nc_hicosmo_q (cosmo, z);
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
  const gdouble Ok = nc_hicosmo_Omega_k0 (cosmo);

  return ((2.0 - q) * E2 - 2.0 * Ok * x * x) / 3.0;
}

NCM_INLINE gdouble
nc_hicosmo_wec (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x  = 1.0 + z;
  const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
  const gdouble Ok = nc_hicosmo_Omega_k0 (cosmo);

  return (E2 - Ok * x * x);
}

NCM_INLINE gdouble
nc_hicosmo_qp (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble E2       = nc_hicosmo_E2 (cosmo, z);
  const gdouble dE2_dz   = nc_hicosmo_dE2_dz (cosmo, z);
  const gdouble d2E2_dz2 = nc_hicosmo_d2E2_dz2 (cosmo, z);

  return (1.0 + z) / (2.0 * E2) * (dE2_dz / (1.0 + z) + d2E2_dz2 - dE2_dz * dE2_dz / E2);
}

NCM_INLINE gdouble
nc_hicosmo_j (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble E2       = nc_hicosmo_E2 (cosmo, z);
  const gdouble dE2_dz   = nc_hicosmo_dE2_dz (cosmo, z);
  const gdouble d2E2_dz2 = nc_hicosmo_d2E2_dz2 (cosmo, z);

  return gsl_pow_2 (1.0 + z) * (d2E2_dz2 - 2.0 * dE2_dz / (1.0 + z)) / (2.0 * E2) + 1.0;
}

NCM_INLINE gdouble
nc_hicosmo_kinetic_w (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_k = nc_hicosmo_E2Omega_k (cosmo, z);
  const gdouble q       = nc_hicosmo_q (cosmo, z);
  const gdouble kw      = (2.0 * q / (1.0 + Omega_k / E2) - 1.0) / 3.0;

  return kw;
}

NCM_INLINE gdouble
nc_hicosmo_mqE2 (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble q        = nc_hicosmo_q (cosmo, z);
  const gdouble E2       = nc_hicosmo_E2 (cosmo, z);

  return -q * E2;
}

NCM_INLINE gdouble
nc_hicosmo_x_alpha (NcHICosmo *cosmo, gdouble alpha)
{
  const gdouble xb = nc_hicosmo_xb (cosmo);
  return xb * exp (- fabs (alpha));
}

NCM_INLINE gdouble
nc_hicosmo_abs_alpha (NcHICosmo *cosmo, gdouble x)
{
  const gdouble xb = nc_hicosmo_xb (cosmo);
  return log (xb / x);
}

NCM_INLINE NcHIPrim *
nc_hicosmo_peek_prim (NcHICosmo *cosmo)
{
  return cosmo->prim;
}

NCM_INLINE NcHIReion *
nc_hicosmo_peek_reion (NcHICosmo *cosmo)
{
  return cosmo->reion;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HICOSMO_INLINE_H_ */
