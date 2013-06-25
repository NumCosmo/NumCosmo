/***************************************************************************
 *            ncm_c.h
 *
 *  Wed Oct 15 17:31:25 2008
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

#ifndef _NCM_C_H_
#define _NCM_C_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <gsl/gsl_const_num.h>
#include <math.h>

G_BEGIN_DECLS

#define NCM_TYPE_C             (ncm_c_get_type ())
#define NCM_C(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_C, NcmC))
#define NCM_C_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_C, NcmCClass))
#define NCM_IS_C(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_C))
#define NCM_IS_C_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_C))
#define NCM_C_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_C, NcmCClass))

typedef struct _NcmCClass NcmCClass;
typedef struct _NcmC NcmC;

struct _NcmCClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmC
{
  /*< private >*/
  GObject parent_instance;
};

GType ncm_c_get_type (void) G_GNUC_CONST;

/*******************************************************************************
 * Mathematical constants
 *******************************************************************************/

G_INLINE_FUNC long double ncm_c_sqrt_1_4pi (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_sqrt_2pi (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_sqrt_3_4pi (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_lnpi_4 (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_ln2pi (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_pi (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_tan_1arcsec (void) G_GNUC_CONST;

G_INLINE_FUNC gdouble ncm_c_degree_to_radian (const gdouble d) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_radian_to_degree (const gdouble r) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_radian_0_2pi (const gdouble r) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_sign_sin (const gdouble r) G_GNUC_CONST;

/*******************************************************************************
 * START: 2006 CODATA recommended values (see end of file)
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_c (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_h (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_hbar (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_fine_struct (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_kb (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_G (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_planck_length (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thomson_cs (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_stefan_boltzmann (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_mass_e (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_mass_p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_mass_n (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_mass_ratio_alpha_p (void) G_GNUC_CONST;

/*******************************************************************************
 * END: 2006 CODATA recommended values
 *******************************************************************************/

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_hc (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_fine_struct_square (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_kpc (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_Mpc (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_AR (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_c2 (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_planck_length2 (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_rest_energy_e (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_rest_energy_p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_rest_energy_n (void) G_GNUC_CONST;

/*******************************************************************************
 * Constants from other places
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_decay_H_rate_2s_1s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_decay_He_rate_2s_1s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_bind_1s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeII_bind_1s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2s_wl (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2p_wl (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_bind_2s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_bind_2p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_2s_m_2p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_2s_m_2p_kb (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2s_wl3_8pi (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2p_wl3_8pi (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_reduced_mass (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_reduced_energy (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_bind (const gint n, const gint j) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_bind_1s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_bind_2s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_bind_2p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_series (const gint n, const gint j) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_2s (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_2p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_series_wl (const gint n, const gint j) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_2s_wl (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_2p_wl (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_2s_wl3_8pi (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_H_Lyman_2p_wl3_8pi (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thermal_wl_e (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thermal_wl_p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thermal_wl_n (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thermal_wn_e (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thermal_wn_p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_thermal_wn_n (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_H_1s (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_H_2s (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_H_2p (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_HeI_1s (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_HeI_2s (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_HeI_2p (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_AU (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_pc (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_mass_solar (void) G_GNUC_CONST;

/* Statistics */

G_INLINE_FUNC long double ncm_c_stats_1sigma (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_stats_2sigma (void) G_GNUC_CONST;
G_INLINE_FUNC long double ncm_c_stats_3sigma (void) G_GNUC_CONST;

/*******************************************************************************
 * Observational data
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_wmap3_cmb_z (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap3_cmb_R (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap3_cmb_sigma_R (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_cmb_z (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_cmb_R (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_cmb_sigma_R (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap7_cmb_z (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap7_cmb_R (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap7_cmb_sigma_R (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_K (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_Ka (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_Q (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_V (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_W (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_z (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_A (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_sigma_A (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_DV (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_sigma_DV (void) ;
G_INLINE_FUNC gdouble ncm_c_bao_percival2007_DV_DV (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_percival2007_sigma_DV_DV (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_percival2010_DV_DV (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_bao_percival2010_sigma_DV_DV (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_hubble_cte_wmap (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_hubble_cte_hst (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_hubble_cte_msa (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_neutrino_n_eff (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_prim_He_Yp (void) ;
G_INLINE_FUNC gdouble ncm_c_prim_H_Yp (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_prim_XHe (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_hubble_radius (void) ;
G_INLINE_FUNC gdouble ncm_c_hubble_radius_planck (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_crit_density (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_crit_mass_density (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_crit_mass_density_solar_Mpc (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_crit_number_density_p (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_crit_number_density_n (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_blackbody_energy_density (void) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_radiation_temp_to_h2omega_r (const gdouble T) G_GNUC_CONST;
G_INLINE_FUNC gdouble ncm_c_radiation_h2Omega_r_to_temp (const gdouble omr) G_GNUC_CONST;

G_END_DECLS

#endif /* _NCM_C_H_ */

#ifndef _NCM_C_INLINE_H_
#define _NCM_C_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

/*******************************************************************************
 * Mathematical constants
 *******************************************************************************/

G_INLINE_FUNC long double ncm_c_sqrt_1_4pi (void)
{ return 0.28209479177387814347403972578038630L; }

G_INLINE_FUNC long double ncm_c_sqrt_2pi (void)
{ return 2.5066282746310005024157652848110452L; }

G_INLINE_FUNC long double ncm_c_sqrt_3_4pi (void)
{ return 0.48860251190291992158638462283834700L; }

G_INLINE_FUNC long double ncm_c_lnpi_4 (void)
{ return 0.28618247146235004353585683783826468L; }

G_INLINE_FUNC long double ncm_c_ln2pi (void)
{ return 1.8378770664093454835606594728112353L; }

G_INLINE_FUNC long double ncm_c_pi (void)
{ return 3.1415926535897932384626433832795029L; }

G_INLINE_FUNC long double ncm_c_tan_1arcsec (void)
{ return 4.8481368111333441675396429478852853e-6L; }


G_INLINE_FUNC gdouble ncm_c_degree_to_radian (const gdouble d)
{ return d * M_PI / 180.0; }

G_INLINE_FUNC gdouble ncm_c_radian_to_degree (const gdouble r)
{ return r * 180.0 / M_PI; }

G_INLINE_FUNC gdouble ncm_c_radian_0_2pi (const gdouble r)
{ return r - 2.0 * M_PI * floor (r / (2.0 * M_PI)); }

G_INLINE_FUNC gdouble ncm_c_sign_sin (const gdouble r)
{ return ncm_c_radian_0_2pi (r) < M_PI ? 1.0 : -1.0; }

/*******************************************************************************
 * START: 2006 CODATA recommended values (see end of file)
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_c (void)
{ return 299792458.0; }

G_INLINE_FUNC gdouble ncm_c_h (void)
{ return 6.62606896e-34; }

G_INLINE_FUNC gdouble ncm_c_hbar (void)
{ return 1.054571628e-34; }

G_INLINE_FUNC gdouble ncm_c_fine_struct (void)
{ return 7.2973525376e-3; }

G_INLINE_FUNC gdouble ncm_c_kb (void)
{ return 1.3806504e-23; }

G_INLINE_FUNC gdouble ncm_c_G (void)
{ return 6.67428e-11; }

G_INLINE_FUNC gdouble ncm_c_planck_length (void)
{ return 1.616252e-35; }

G_INLINE_FUNC gdouble ncm_c_thomson_cs (void)
{ return 0.6652458558e-28; }

G_INLINE_FUNC gdouble ncm_c_stefan_boltzmann (void)
{ return 5.670400e-8; }

G_INLINE_FUNC gdouble ncm_c_mass_e (void)
{ return 9.10938215e-31; }

G_INLINE_FUNC gdouble ncm_c_mass_p (void)
{ return 1.672621637e-27; }

G_INLINE_FUNC gdouble ncm_c_mass_n (void)
{ return 1.674927211e-27; }

G_INLINE_FUNC gdouble ncm_c_mass_ratio_alpha_p (void)
{ return 3.97259968951; }

/*******************************************************************************
 * END: 2006 CODATA recommended values
 *******************************************************************************/

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_hc (void)
{ return ncm_c_h () * ncm_c_c (); }

G_INLINE_FUNC gdouble ncm_c_fine_struct_square (void)
{ return ncm_c_fine_struct () * ncm_c_fine_struct (); }

G_INLINE_FUNC gdouble ncm_c_kpc (void)
{ return GSL_CONST_NUM_KILO * ncm_c_pc (); }

G_INLINE_FUNC gdouble ncm_c_Mpc (void)
{ return GSL_CONST_NUM_MEGA * ncm_c_pc (); }

G_INLINE_FUNC gdouble ncm_c_AR (void)
{ return 4.0 * ncm_c_stefan_boltzmann () / ncm_c_c (); }

G_INLINE_FUNC gdouble ncm_c_c2 (void)
{ return ncm_c_c () * ncm_c_c (); }

G_INLINE_FUNC gdouble ncm_c_planck_length2 (void)
{ return ncm_c_planck_length () * ncm_c_planck_length (); }

G_INLINE_FUNC gdouble ncm_c_rest_energy_e (void)
{ return ncm_c_mass_e () * ncm_c_c2 (); }

G_INLINE_FUNC gdouble ncm_c_rest_energy_p (void)
{ return ncm_c_mass_p () * ncm_c_c2 (); }

G_INLINE_FUNC gdouble ncm_c_rest_energy_n (void)
{ return ncm_c_mass_n () * ncm_c_c2 (); }

/*******************************************************************************
 * Constants from other places
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_decay_H_rate_2s_1s (void)
{ return 8.22458; }

G_INLINE_FUNC gdouble ncm_c_decay_He_rate_2s_1s (void)
{ return 51.3; }

G_INLINE_FUNC gdouble ncm_c_HeI_bind_1s (void)
{ return 1.98310772e7 * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeII_bind_1s (void)
{ return 4.389088863e7 * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2s (void)
{ return 1.66277434e7 * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2p (void)
{ return 1.71134891e7 * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2s_wl (void)
{ return 1.0 / 1.66277434e7; }

G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2p_wl (void)
{ return 1.0 / 1.71134891e7; }

G_INLINE_FUNC gdouble ncm_c_HeI_bind_2s (void)
{ return 3.2033338e6 * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeI_bind_2p (void)
{ return 2.7175881e6 * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeI_2s_m_2p (void)
{ return (3.2033338e6 - 2.7175881e6) * ncm_c_hc (); }

G_INLINE_FUNC gdouble ncm_c_HeI_2s_m_2p_kb (void)
{ return (3.2033338e6 - 2.7175881e6) * ncm_c_hc () / ncm_c_kb (); }

G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2s_wl3_8pi (void)
{ return ncm_c_HeI_Lyman_2s_wl() * ncm_c_HeI_Lyman_2s_wl () *
	ncm_c_HeI_Lyman_2s_wl () / (8.0 * M_PI); }

G_INLINE_FUNC gdouble ncm_c_HeI_Lyman_2p_wl3_8pi (void)
{ return ncm_c_HeI_Lyman_2p_wl() * ncm_c_HeI_Lyman_2p_wl () *
	ncm_c_HeI_Lyman_2p_wl () / (8.0 * M_PI); }

G_INLINE_FUNC gdouble ncm_c_H_reduced_mass (void)
{ return ncm_c_mass_e () / (1.0 + ncm_c_mass_e () / ncm_c_mass_p ()); }

G_INLINE_FUNC gdouble ncm_c_H_reduced_energy (void)
{ return ncm_c_H_reduced_mass () * ncm_c_c2 (); }

G_INLINE_FUNC gdouble ncm_c_H_bind (const gint n, const gint j)
{ return ncm_c_H_reduced_energy () *
	(1.0 - 1.0 / sqrt (1.0 + ncm_c_fine_struct_square () /
	                   pow (n - j - 0.5 + sqrt(pow(j + 0.5, 2.0) -
	                                           ncm_c_fine_struct_square () ),
	                        2.0) )); }

G_INLINE_FUNC gdouble ncm_c_H_bind_1s (void)
{ return ncm_c_H_bind (1.0, 0.5); }

G_INLINE_FUNC gdouble ncm_c_H_bind_2s (void)
{ return ncm_c_H_bind (2.0, 0.5); }

G_INLINE_FUNC gdouble ncm_c_H_bind_2p (void)
{ return ncm_c_H_bind (2.0, 1.5); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_series (const gint n, const gint j)
{ return ncm_c_H_bind_1s () - ncm_c_H_bind (n,j); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_2s (void)
{ return ncm_c_H_Lyman_series (2, 0.5); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_2p (void)
{ return ncm_c_H_Lyman_series (2, 1.5); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_series_wl (const gint n, const gint j)
{ return ncm_c_hc () / ncm_c_H_Lyman_series (n,j); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_2s_wl (void)
{ return ncm_c_H_Lyman_series_wl (2, 0.5); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_2p_wl (void)
{ return ncm_c_H_Lyman_series_wl (2, 1.5); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_2s_wl3_8pi (void)
{ return ncm_c_H_Lyman_2s_wl () * ncm_c_H_Lyman_2s_wl () *
	ncm_c_H_Lyman_2s_wl () / (8.0 * M_PI); }

G_INLINE_FUNC gdouble ncm_c_H_Lyman_2p_wl3_8pi (void)
{ return ncm_c_H_Lyman_2p_wl () * ncm_c_H_Lyman_2p_wl () *
	ncm_c_H_Lyman_2p_wl () / (8.0 * M_PI); }

G_INLINE_FUNC gdouble ncm_c_thermal_wl_e (void)
{ return sqrt ((2.0 * M_PI * ncm_c_hbar () * ncm_c_hbar () ) /
               (ncm_c_mass_e () * ncm_c_kb ())); }

G_INLINE_FUNC gdouble ncm_c_thermal_wl_p (void)
{ return sqrt ((2.0 * M_PI * ncm_c_hbar () * ncm_c_hbar () ) /
               (ncm_c_mass_p () * ncm_c_kb ())); }

G_INLINE_FUNC gdouble ncm_c_thermal_wl_n (void)
{ return sqrt ((2.0 * M_PI * ncm_c_hbar () * ncm_c_hbar () ) /
               (ncm_c_mass_n () * ncm_c_kb ())); }

G_INLINE_FUNC gdouble ncm_c_thermal_wn_e (void)
{ return 1.0 / ncm_c_thermal_wl_e (); }

G_INLINE_FUNC gdouble ncm_c_thermal_wn_p (void)
{ return 1.0 / ncm_c_thermal_wl_p (); }

G_INLINE_FUNC gdouble ncm_c_thermal_wn_n (void)
{ return 1.0 / ncm_c_thermal_wl_n (); }

G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_H_1s (const gdouble T)
{ return pow(ncm_c_thermal_wn_e (), 3.0) *
 exp(- ncm_c_H_bind_1s () / (ncm_c_kb () * T)); }

G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_H_2s (const gdouble T)
{ return pow(ncm_c_thermal_wn_e (), 3.0) *
 exp(- ncm_c_H_bind_2s () / (ncm_c_kb () * T)); }

G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_H_2p (const gdouble T)
{ return pow(ncm_c_thermal_wn_e (), 3.0) *
 exp(- ncm_c_H_bind_2p () / (ncm_c_kb () * T)); }

G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_HeI_1s (const gdouble T)
{ return pow(ncm_c_thermal_wn_e (), 3.0) *
 exp(- ncm_c_HeI_bind_1s () / (ncm_c_kb () * T)); }

G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_HeI_2s (const gdouble T)
{ return pow(ncm_c_thermal_wn_e (), 3.0) *
 exp(- ncm_c_HeI_bind_2s () / (ncm_c_kb () * T)); }

G_INLINE_FUNC gdouble ncm_c_boltzmann_factor_HeI_2p (const gdouble T)
{ return pow(ncm_c_thermal_wn_e (), 3.0) *
 exp(- ncm_c_HeI_bind_2p () / (ncm_c_kb () * T)); }

G_INLINE_FUNC gdouble ncm_c_AU (void)
{ return 1.49597870691e11; }

G_INLINE_FUNC gdouble ncm_c_pc (void)
{ return 3.085678e16; }

G_INLINE_FUNC gdouble ncm_c_mass_solar (void)
{ return 1.98892e30; }

/* Statistics */

G_INLINE_FUNC long double ncm_c_stats_1sigma (void)
{ return 0.6826894921370858971704650912640758449558L; }

G_INLINE_FUNC long double ncm_c_stats_2sigma (void)
{ return 0.9544997361036415855994347256669331250564L; }

G_INLINE_FUNC long double ncm_c_stats_3sigma (void)
{ return 0.9973002039367398109466963704648100452443L; }

/*******************************************************************************
 * Observational data
 *******************************************************************************/

G_INLINE_FUNC gdouble ncm_c_wmap3_cmb_z (void)
{ return 1089.0; }

G_INLINE_FUNC gdouble ncm_c_wmap3_cmb_R (void)
{ return 1.70; }

G_INLINE_FUNC gdouble ncm_c_wmap3_cmb_sigma_R (void)
{ return 0.03; }

G_INLINE_FUNC gdouble ncm_c_wmap5_cmb_z (void)
{ return 1090.0; }

G_INLINE_FUNC gdouble ncm_c_wmap5_cmb_R (void)
{ return 1.71; }

G_INLINE_FUNC gdouble ncm_c_wmap5_cmb_sigma_R (void)
{ return 0.019; }

G_INLINE_FUNC gdouble ncm_c_wmap7_cmb_z (void)
{ return 1091.3; }

G_INLINE_FUNC gdouble ncm_c_wmap7_cmb_R (void)
{ return 1.725; }

G_INLINE_FUNC gdouble ncm_c_wmap7_cmb_sigma_R (void)
{ return 0.018; }

G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_K (void)
{ return 1.436; }

G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_Ka (void)
{ return 1.470; }

G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_Q (void)
{ return 2.197; }

G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_V (void)
{ return 3.133; }

G_INLINE_FUNC gdouble ncm_c_wmap5_coadded_I_W (void)
{ return 6.538; }

G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_z (void)
{ return 0.35; }

G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_A (void)
{ return 0.469; }

G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_sigma_A (void)
{ return 0.017; }

G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_DV (void)
{ return 1334.0; }

G_INLINE_FUNC gdouble ncm_c_bao_eisenstein_sigma_DV (void)
{ return 88.0; }

G_INLINE_FUNC gdouble ncm_c_bao_percival2007_DV_DV (void)
{ return 1.812; }

G_INLINE_FUNC gdouble ncm_c_bao_percival2007_sigma_DV_DV (void)
{ return 0.060; }

G_INLINE_FUNC gdouble ncm_c_bao_percival2010_DV_DV (void)
{ return 1.736; }

G_INLINE_FUNC gdouble ncm_c_bao_percival2010_sigma_DV_DV (void)
{ return 0.065; }

G_INLINE_FUNC gdouble ncm_c_hubble_cte_wmap (void)
{ return 73.0; }

G_INLINE_FUNC gdouble ncm_c_hubble_cte_hst (void)
{ return 72.0; }

G_INLINE_FUNC gdouble ncm_c_hubble_cte_msa (void)
{ return 68.0; }

G_INLINE_FUNC gdouble ncm_c_neutrino_n_eff (void)
{ return 3.04; }

G_INLINE_FUNC gdouble ncm_c_prim_He_Yp (void)
{ return 0.24; }

G_INLINE_FUNC gdouble ncm_c_prim_H_Yp (void)
{ return 1.0 - ncm_c_prim_He_Yp (); }

G_INLINE_FUNC gdouble ncm_c_prim_XHe (void)
{ return ncm_c_prim_He_Yp () / (ncm_c_mass_ratio_alpha_p () * ncm_c_prim_H_Yp ()); }

G_INLINE_FUNC gdouble ncm_c_hubble_radius (void)
{ return ncm_c_c () / (100.0e3); }

G_INLINE_FUNC gdouble ncm_c_hubble_radius_planck (void)
{ return ncm_c_hubble_radius () * ncm_c_Mpc () / ncm_c_planck_length (); }

G_INLINE_FUNC gdouble ncm_c_crit_density (void)
{ return 3.0 * pow (ncm_c_c () / (10.0 * ncm_c_pc ()), 2.0) / (8.0 * M_PI * ncm_c_G ()); }

G_INLINE_FUNC gdouble ncm_c_crit_mass_density (void)
{ return 3.0 * pow (1.0 / (10.0 * ncm_c_pc ()), 2.0) / (8.0 * M_PI * ncm_c_G ()); }

G_INLINE_FUNC gdouble ncm_c_crit_mass_density_solar_Mpc (void)
{ return ncm_c_crit_mass_density () / ncm_c_mass_solar () * pow (ncm_c_Mpc (), 3.0); }

G_INLINE_FUNC gdouble ncm_c_crit_number_density_p (void)
{ return ncm_c_crit_density () / ncm_c_rest_energy_p (); }

G_INLINE_FUNC gdouble ncm_c_crit_number_density_n (void)
{ return ncm_c_crit_density () / ncm_c_rest_energy_n (); }

G_INLINE_FUNC gdouble ncm_c_blackbody_energy_density (void)
{ return 4.0 * ncm_c_stefan_boltzmann () / ncm_c_c (); }

G_INLINE_FUNC gdouble ncm_c_radiation_temp_to_h2omega_r (const gdouble T)
{ return ncm_c_blackbody_energy_density () * ((T * T) * (T * T)) / ncm_c_crit_density (); }

G_INLINE_FUNC gdouble ncm_c_radiation_h2Omega_r_to_temp (const gdouble omr)
{ return pow(ncm_c_crit_density () * omr / ncm_c_blackbody_energy_density (), 0.25); }

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_C_INLINE_H_ */
