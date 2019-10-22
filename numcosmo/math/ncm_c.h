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
#include <numcosmo/math/ncm_util.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

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

NCM_INLINE long double ncm_c_sqrt_1_4pi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_sqrt_pi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_sqrt_2pi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_sqrt_pi_2 (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_sqrt_3_4pi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_ln2 (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_ln3 (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_lnpi_4 (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_ln2pi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_lnpi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_pi (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_2_pi_2 (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_tan_1arcsec (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_deg2_steradian (void) G_GNUC_CONST;

NCM_INLINE gdouble ncm_c_degree_to_radian (const gdouble d) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_radian_to_degree (const gdouble r) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_radian_0_2pi (const gdouble r) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_sign_sin (const gdouble r) G_GNUC_CONST;

/*******************************************************************************
 * START: 2018 CODATA recommended values (see end of file)
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_c (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_h (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_hbar (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_fine_struct (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_kb (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_G (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_planck_length (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thomson_cs (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_stefan_boltzmann (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_atomic (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_e (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_n (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_ratio_alpha_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_ratio_e_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_Rinf (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_Ry (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_eV (void) G_GNUC_CONST;

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_year (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_lightyear (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_lightyear_pc (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_Glightyear_Mpc (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_hc (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_fine_struct_square (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_AR (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_c2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_planck_length2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_atomic (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_e (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_n (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thermal_wl_e (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thermal_wl_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thermal_wl_n (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thermal_wn_e (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thermal_wn_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_thermal_wn_n (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_H_reduced_mass (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_H_reduced_energy (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_H_bind (const gdouble n, const gdouble j) G_GNUC_CONST;

/*******************************************************************************
 * END: 2018 CODATA recommended values
 *******************************************************************************/

/*******************************************************************************
 * START: IUPAC related constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_mass_1H_u (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_2H_u (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_3H_u (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_3He_u (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_4He_u (void) G_GNUC_CONST;

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_mass_1H (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_2H (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_3H (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_3He (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_4He (void) G_GNUC_CONST;

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_rest_energy_1H (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_2H (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_3H (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_3He (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_rest_energy_4He (void) G_GNUC_CONST;

NCM_INLINE gdouble ncm_c_mass_ratio_4He_1H (void) G_GNUC_CONST;

/*******************************************************************************
 * END: IUPAC related constants
 *******************************************************************************/

/*******************************************************************************
 * START: IAU related constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_au (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_pc (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_kpc (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_Mpc (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_G_mass_solar (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_mass_solar (void) G_GNUC_CONST;

/*******************************************************************************
 * END: IAU related constants
 *******************************************************************************/

/*******************************************************************************
 * START: NIST Atomic Spectra database
 *******************************************************************************/
/*******************************************************************************
 * -- START: Hydrogen I
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

NCM_INLINE gdouble ncm_c_HI_ion_wn_1s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_wn_2s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_wn_2p_2P0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_wn_2p_2P3_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_wn_2p_2Pmean (void) G_GNUC_CONST;

/* Ionization energy: E */

NCM_INLINE gdouble ncm_c_HI_ion_E_1s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_E_2s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_E_2p_2P0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_E_2p_2P3_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_ion_E_2p_2Pmean (void) G_GNUC_CONST;

/* Lyman series wavenumber: wn */

NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2p_2P0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2p_2P3_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2p_2Pmean (void) G_GNUC_CONST;

/* Lyman series wavelength: wl */

NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2p_2P0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2p_2P3_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2p_2Pmean (void) G_GNUC_CONST;

/* Lyman series factor: wl^3 / (8pi) */

NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2s_2S0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2p_2P0_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2p_2P3_5 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2p_2Pmean (void) G_GNUC_CONST;

/* Boltzmann factor */

NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_1s_2S0_5 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2s_2S0_5 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2p_2P0_5 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2p_2P3_5 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2p_2Pmean (const gdouble T) G_GNUC_CONST;

/*******************************************************************************
 * -- END: Hydrogen I
 *******************************************************************************/
/*******************************************************************************
 * -- START: Helium I
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeI_ion_wn_1s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2s_3S1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_1P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3P0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3P2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3Pmean (void) G_GNUC_CONST;

/* Ionization energy: E */

NCM_INLINE gdouble ncm_c_HeI_ion_E_1s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2s_3S1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_1P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3P0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3P2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3Pmean (void) G_GNUC_CONST;

/* Lyman series wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2s_3S1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_1P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3P0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3P2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3Pmean (void) G_GNUC_CONST;

/* Lyman series wavelength: wl */

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2s_3S1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_1P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3P0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3P2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3Pmean (void) G_GNUC_CONST;

/* Lyman series factor: wl^3 / (8pi) */

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2s_3S1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3P0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3P1 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3P2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean (void) G_GNUC_CONST;

/* Boltzmann factor */

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_1s_1S0 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2s_1S0 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2s_3S1 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_1P1 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3P0 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3P1 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3P2 (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3Pmean (const gdouble T) G_GNUC_CONST;

/* Balmer series wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeI_Balmer_wn_2p_1P1_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Balmer_wn_2p_3Pmean_2s_3S1 (void) G_GNUC_CONST;

/* Balmer series: E / k_B */

NCM_INLINE gdouble ncm_c_HeI_Balmer_E_kb_2p_1P1_2s_1S0 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_HeI_Balmer_E_kb_2p_3Pmean_2s_3S1 (void) G_GNUC_CONST;

/*******************************************************************************
 * -- END: Helium I
 *******************************************************************************/
/*******************************************************************************
 * -- START: Helium II
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeII_ion_wn_1s_2S0_5 (void) G_GNUC_CONST;

/* Ionization energy: E */

NCM_INLINE gdouble ncm_c_HeII_ion_E_1s_2S0_5 (void) G_GNUC_CONST;

/*******************************************************************************
 * -- END: Helium II
 *******************************************************************************/
/*******************************************************************************
 * END: NIST Atomic Spectra database
 *******************************************************************************/

/*******************************************************************************
 * Constants from other sources
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_decay_H_rate_2s_1s (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_decay_He_rate_2s_1s (void) G_GNUC_CONST;

/*******************************************************************************
 * START: Statistics
 *******************************************************************************/

NCM_INLINE long double ncm_c_stats_1sigma (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_stats_2sigma (void) G_GNUC_CONST;
NCM_INLINE long double ncm_c_stats_3sigma (void) G_GNUC_CONST;

/*******************************************************************************
 * END: Statistics
 *******************************************************************************/

/*******************************************************************************
 * START: Observational data
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_wmap5_coadded_I_K (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_wmap5_coadded_I_Ka (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_wmap5_coadded_I_Q (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_wmap5_coadded_I_V (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_wmap5_coadded_I_W (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_hubble_cte_wmap (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_hubble_cte_hst (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_hubble_cte_msa (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_hubble_radius_hm1_Mpc (void) ;
NCM_INLINE gdouble ncm_c_hubble_radius_hm1_planck (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_crit_density_h2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_crit_mass_density_h2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_crit_mass_density_h2_solar_mass_Mpc3 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_crit_number_density_p (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_crit_number_density_n (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_blackbody_energy_density (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_blackbody_per_crit_density_h2 (void) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_radiation_temp_to_h2omega_r (const gdouble T) G_GNUC_CONST;
NCM_INLINE gdouble ncm_c_radiation_h2Omega_r0_to_temp (const gdouble omr) G_GNUC_CONST;

/*******************************************************************************
 * END: Observational data
 *******************************************************************************/

G_END_DECLS

#endif /* _NCM_C_H_ */

#ifndef _NCM_C_INLINE_H_
#define _NCM_C_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

/*******************************************************************************
 * Mathematical constants
 *******************************************************************************/

NCM_INLINE long double ncm_c_sqrt_1_4pi (void)
{ return 0.28209479177387814347403972578038630L; }

NCM_INLINE long double ncm_c_sqrt_pi (void)
{ return 1.77245385090551602729816748334114518L; }

NCM_INLINE long double ncm_c_sqrt_2pi (void)
{ return 2.5066282746310005024157652848110452L; }

NCM_INLINE long double ncm_c_sqrt_pi_2 (void)
{ return 1.2533141373155002512078826424055226L; }

NCM_INLINE long double ncm_c_sqrt_3_4pi (void)
{ return 0.48860251190291992158638462283834700L; }

NCM_INLINE long double ncm_c_ln2 (void)
{ return 0.69314718055994530941723212145817657L; }

NCM_INLINE long double ncm_c_ln3 (void)
{ return 1.0986122886681096913952452369225257L; }

NCM_INLINE long double ncm_c_lnpi_4 (void)
{ return 0.28618247146235004353585683783826468L; }

NCM_INLINE long double ncm_c_ln2pi (void)
{ return 1.8378770664093454835606594728112353L; }

NCM_INLINE long double ncm_c_lnpi (void)
{ return 1.1447298858494001741434273513530587L; }

NCM_INLINE long double ncm_c_pi (void)
{ return 3.1415926535897932384626433832795029L; }

NCM_INLINE long double ncm_c_2_pi_2 (void)
{ return 19.739208802178717237668981999752302L; }

NCM_INLINE long double ncm_c_tan_1arcsec (void)
{ return 4.8481368111333441675396429478852853e-6L; }

NCM_INLINE long double ncm_c_deg2_steradian (void)
{ return 3.0461741978670859934674354937889355e-4L; }

NCM_INLINE gdouble ncm_c_degree_to_radian (const gdouble d)
{ return d * ncm_c_pi () / 180.0; }

NCM_INLINE gdouble ncm_c_radian_to_degree (const gdouble r)
{ return r * 180.0 / ncm_c_pi (); }

NCM_INLINE gdouble ncm_c_radian_0_2pi (const gdouble r)
{ return r - 2.0 * ncm_c_pi () * floor (r / (2.0 * ncm_c_pi ())); }

NCM_INLINE gdouble ncm_c_sign_sin (const gdouble r)
{ return ncm_c_radian_0_2pi (r) < ncm_c_pi () ? 1.0 : -1.0; }

/*******************************************************************************
 * START: 2018 CODATA recommended values (see constants.txt)
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_c (void)
{ return 299792458.0; }

NCM_INLINE gdouble ncm_c_h (void)
{ return 6.62607015e-34; }
          
NCM_INLINE gdouble ncm_c_hbar (void)
{ return 1.054571817e-34; }

NCM_INLINE gdouble ncm_c_fine_struct (void)
{ return 7.2973525693e-3; }

NCM_INLINE gdouble ncm_c_kb (void)
{ return 1.380649e-23; }

NCM_INLINE gdouble ncm_c_G (void)
{ return 6.67430e-11; }

NCM_INLINE gdouble ncm_c_planck_length (void)
{ return 1.616255e-35; }

NCM_INLINE gdouble ncm_c_thomson_cs (void)
{ return 0.66524587321e-28; }

NCM_INLINE gdouble ncm_c_stefan_boltzmann (void)
{ return 5.670374419e-8; }

NCM_INLINE gdouble ncm_c_mass_atomic (void)
{ return 1.66053906660e-27; }

NCM_INLINE gdouble ncm_c_mass_e (void)
{ return 9.1093837015e-31; }
         
NCM_INLINE gdouble ncm_c_mass_p (void)
{ return 1.67262192369e-27; }

NCM_INLINE gdouble ncm_c_mass_n (void)
{ return 1.67492749804e-27; }

NCM_INLINE gdouble ncm_c_mass_ratio_alpha_p (void)
{ return 3.97259969009; }

NCM_INLINE gdouble ncm_c_mass_ratio_e_p (void)
{ return 5.44617021487e-4; }

NCM_INLINE gdouble ncm_c_Rinf (void)
{ return 10973731.568160; }

NCM_INLINE gdouble ncm_c_Ry (void)
{ return 2.1798723611035e-18; }

NCM_INLINE gdouble ncm_c_eV (void)
{ return 1.602176634e-19; }

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_year (void)
{ return 365.25 * 24.0 * 60.0 * 60.0; }

NCM_INLINE gdouble ncm_c_lightyear (void)
{ return ncm_c_c () * ncm_c_year (); }

NCM_INLINE gdouble ncm_c_lightyear_pc (void)
{ return ncm_c_lightyear () / ncm_c_pc (); }

NCM_INLINE gdouble ncm_c_Glightyear_Mpc (void)
{ return 1.0e3 * ncm_c_lightyear () / ncm_c_pc (); }

NCM_INLINE gdouble ncm_c_hc (void)
{ return ncm_c_h () * ncm_c_c (); }

NCM_INLINE gdouble ncm_c_fine_struct_square (void)
{ return ncm_c_fine_struct () * ncm_c_fine_struct (); }

NCM_INLINE gdouble ncm_c_AR (void)
{ return 4.0 * ncm_c_stefan_boltzmann () / ncm_c_c (); }

NCM_INLINE gdouble ncm_c_c2 (void)
{ return ncm_c_c () * ncm_c_c (); }

NCM_INLINE gdouble ncm_c_planck_length2 (void)
{ return ncm_c_planck_length () * ncm_c_planck_length (); }

NCM_INLINE gdouble ncm_c_rest_energy_atomic (void)
{ return ncm_c_mass_atomic () * ncm_c_c2 (); }

NCM_INLINE gdouble ncm_c_rest_energy_e (void)
{ return ncm_c_mass_e () * ncm_c_c2 (); }

NCM_INLINE gdouble ncm_c_rest_energy_p (void)
{ return ncm_c_mass_p () * ncm_c_c2 (); }

NCM_INLINE gdouble ncm_c_rest_energy_n (void)
{ return ncm_c_mass_n () * ncm_c_c2 (); }

NCM_INLINE gdouble ncm_c_H_reduced_mass (void)
{ return ncm_c_mass_e () / (1.0 + ncm_c_mass_ratio_e_p ()); }

NCM_INLINE gdouble ncm_c_thermal_wl_e (void)
{ return sqrt ((2.0 * ncm_c_pi () * ncm_c_hbar () * ncm_c_hbar () ) /
               (ncm_c_mass_e () * ncm_c_kb ())); }

NCM_INLINE gdouble ncm_c_thermal_wl_p (void)
{ return sqrt ((2.0 * ncm_c_pi () * ncm_c_hbar () * ncm_c_hbar () ) /
               (ncm_c_mass_p () * ncm_c_kb ())); }

NCM_INLINE gdouble ncm_c_thermal_wl_n (void)
{ return sqrt ((2.0 * ncm_c_pi () * ncm_c_hbar () * ncm_c_hbar () ) /
               (ncm_c_mass_n () * ncm_c_kb ())); }

NCM_INLINE gdouble ncm_c_thermal_wn_e (void)
{ return 1.0 / ncm_c_thermal_wl_e (); }

NCM_INLINE gdouble ncm_c_thermal_wn_p (void)
{ return 1.0 / ncm_c_thermal_wl_p (); }

NCM_INLINE gdouble ncm_c_thermal_wn_n (void)
{ return 1.0 / ncm_c_thermal_wl_n (); }

NCM_INLINE gdouble ncm_c_H_reduced_energy (void)
{ return ncm_c_H_reduced_mass () * ncm_c_c2 (); }

NCM_INLINE gdouble ncm_c_H_bind (const gdouble n, const gdouble j)
{ 
  const gdouble j_1_2            = j + 0.5;
  const gdouble j_1_2_pow_2      = gsl_pow_2 (j_1_2);
  const gdouble alpha2           = ncm_c_fine_struct_square ();
  const gdouble d_j              = j_1_2 * (- ncm_util_sqrt1px_m1 (- alpha2 / j_1_2_pow_2)); /*(1.0 - sqrt (1.0 - alpha2 / j_1_2_pow_2))*/
  const gdouble f_arg            = ncm_c_fine_struct_square () / gsl_pow_2 (n - d_j);
  const gdouble f_nj_m_one       = -ncm_util_sqrt1px_m1 (f_arg) / sqrt (1.0 + f_arg); /*(1.0 - sqrt (1.0 + f_arg))*/
  const gdouble r                = ncm_c_mass_ratio_e_p ();
  const gdouble E_binding        = ncm_c_H_reduced_energy () * f_nj_m_one * (1.0 - 0.5 * f_nj_m_one * r / gsl_pow_2 (1.0 + r)); 
  return - E_binding;
}

/*******************************************************************************
 * END: 2018 CODATA recommended values
 *******************************************************************************/

/*******************************************************************************
 * START: IUPAC related constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_mass_1H_u (void)
{ return 1.00782503223; }

NCM_INLINE gdouble ncm_c_mass_2H_u (void)
{ return 2.01410177812; }

NCM_INLINE gdouble ncm_c_mass_3H_u (void)
{ return 3.0160492779; }

NCM_INLINE gdouble ncm_c_mass_3He_u (void)
{ return 3.0160293201; }

NCM_INLINE gdouble ncm_c_mass_4He_u (void)
{ return 4.00260325413; }

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_mass_1H (void)
{ return ncm_c_mass_1H_u () * ncm_c_mass_atomic (); }

NCM_INLINE gdouble ncm_c_mass_2H (void)
{ return ncm_c_mass_2H_u () * ncm_c_mass_atomic (); }

NCM_INLINE gdouble ncm_c_mass_3H (void)
{ return ncm_c_mass_3H_u () * ncm_c_mass_atomic (); }

NCM_INLINE gdouble ncm_c_mass_3He (void)
{ return ncm_c_mass_3He_u () * ncm_c_mass_atomic (); }

NCM_INLINE gdouble ncm_c_mass_4He (void)
{ return ncm_c_mass_4He_u () * ncm_c_mass_atomic ();}

NCM_INLINE gdouble ncm_c_rest_energy_1H (void)
{ return ncm_c_mass_1H_u () * ncm_c_rest_energy_atomic (); }

NCM_INLINE gdouble ncm_c_rest_energy_2H (void)
{ return ncm_c_mass_2H_u () * ncm_c_rest_energy_atomic (); }

NCM_INLINE gdouble ncm_c_rest_energy_3H (void)
{ return ncm_c_mass_3H_u () * ncm_c_rest_energy_atomic (); }

NCM_INLINE gdouble ncm_c_rest_energy_3He (void)
{ return ncm_c_mass_3He_u () * ncm_c_rest_energy_atomic (); }

NCM_INLINE gdouble ncm_c_rest_energy_4He (void)
{ return ncm_c_mass_4He_u () * ncm_c_rest_energy_atomic (); }

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_mass_ratio_4He_1H (void)
{ return ncm_c_mass_4He_u () / ncm_c_mass_1H_u (); }

/*******************************************************************************
 * END: IUPAC related constants
 *******************************************************************************/

/*******************************************************************************
 * START: IAU related constants
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_au (void)
{ return 1.49597870700e11; }

NCM_INLINE gdouble ncm_c_pc (void)
{ return 3.0856775814913672789139379577965e16L; }

NCM_INLINE gdouble ncm_c_kpc (void)
{ return 1.0e3 * ncm_c_pc (); }

NCM_INLINE gdouble ncm_c_Mpc (void)
{ return 1.0e6 * ncm_c_pc (); }

NCM_INLINE gdouble ncm_c_G_mass_solar (void)
{ return 1.3271244e20; }

NCM_INLINE gdouble ncm_c_mass_solar (void)
{ return ncm_c_G_mass_solar () / ncm_c_G (); }

/*******************************************************************************
 * END: IAU related constants
 *******************************************************************************/

/*******************************************************************************
 * START: NIST Atomic Spectra database
 *******************************************************************************/
/*******************************************************************************
 * -- START: Hydrogen I
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

NCM_INLINE gdouble ncm_c_HI_ion_wn_1s_2S0_5 (void)
{ return 1.0967877174307e7; }

NCM_INLINE gdouble ncm_c_HI_ion_wn_2s_2S0_5 (void)
{ return ncm_c_HI_ion_wn_1s_2S0_5 () - ncm_c_HI_Lyman_wn_2s_2S0_5 (); }

NCM_INLINE gdouble ncm_c_HI_ion_wn_2p_2P0_5 (void)
{ return ncm_c_HI_ion_wn_1s_2S0_5 () - ncm_c_HI_Lyman_wn_2p_2P0_5 (); }

NCM_INLINE gdouble ncm_c_HI_ion_wn_2p_2P3_5 (void)
{ return ncm_c_HI_ion_wn_1s_2S0_5 () - ncm_c_HI_Lyman_wn_2p_2P3_5 (); }

NCM_INLINE gdouble ncm_c_HI_ion_wn_2p_2Pmean (void)
{ return 0.5 * (ncm_c_HI_ion_wn_2p_2P0_5 () + ncm_c_HI_ion_wn_2p_2P3_5 ()); }

/* Ionization energy: E */

NCM_INLINE gdouble ncm_c_HI_ion_E_1s_2S0_5 (void)
{ return ncm_c_HI_ion_wn_1s_2S0_5 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HI_ion_E_2s_2S0_5 (void)
{ return ncm_c_HI_ion_wn_2s_2S0_5 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HI_ion_E_2p_2P0_5 (void)
{ return ncm_c_HI_ion_wn_2p_2P0_5 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HI_ion_E_2p_2P3_5 (void) 
{ return ncm_c_HI_ion_wn_2p_2P3_5 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HI_ion_E_2p_2Pmean (void) 
{ return ncm_c_HI_ion_wn_2p_2Pmean () * ncm_c_hc (); }

/* Lyman series wavenumber: wn */

NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2s_2S0_5 (void)
{ return 8.22589543992821e6; }

NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2p_2P0_5 (void)
{ return 8.22589191133e6; }

NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2p_2P3_5 (void)
{ return 8.22592850014e6; }

NCM_INLINE gdouble ncm_c_HI_Lyman_wn_2p_2Pmean (void)
{ return 0.5 * (ncm_c_HI_Lyman_wn_2p_2P0_5 () + ncm_c_HI_Lyman_wn_2p_2P3_5 ()); }

/* Lyman series wavelength: wl */

NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2s_2S0_5 (void)
{ return 1.0 / ncm_c_HI_Lyman_wn_2s_2S0_5 (); }

NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2p_2P0_5 (void)
{ return 1.0 / ncm_c_HI_Lyman_wn_2p_2P0_5 (); }

NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2p_2P3_5 (void)
{ return 1.0 / ncm_c_HI_Lyman_wn_2p_2P3_5 (); }

NCM_INLINE gdouble ncm_c_HI_Lyman_wl_2p_2Pmean (void)
{ return 1.0 / ncm_c_HI_Lyman_wn_2p_2Pmean (); }

/* Lyman series factor: wl^3 / (8pi) */

NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2s_2S0_5 (void)
{ return gsl_pow_3 (ncm_c_HI_Lyman_wl_2s_2S0_5 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2p_2P0_5 (void)
{ return gsl_pow_3 (ncm_c_HI_Lyman_wl_2p_2P0_5 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2p_2P3_5 (void)
{ return gsl_pow_3 (ncm_c_HI_Lyman_wl_2p_2P3_5 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HI_Lyman_wl3_8pi_2p_2Pmean (void)
{ return gsl_pow_3 (ncm_c_HI_Lyman_wl_2p_2Pmean ()) / (8.0 * ncm_c_pi ()); }

/* Boltzmann factor */

NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_1s_2S0_5 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp(- ncm_c_HI_ion_E_1s_2S0_5 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2s_2S0_5 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp(- ncm_c_HI_ion_E_2s_2S0_5 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2p_2P0_5 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp(- ncm_c_HI_ion_E_2p_2P0_5 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2p_2P3_5 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp(- ncm_c_HI_ion_E_2p_2P3_5 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HI_2p_2Pmean (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp(- ncm_c_HI_ion_E_2p_2Pmean () / (ncm_c_kb () * T)); }

/*******************************************************************************
 * -- END: Hydrogen I
 *******************************************************************************/
/*******************************************************************************
 * -- START: Helium I
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeI_ion_wn_1s_1S0 (void)
{ return 1.9831066637e7; }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2s_1S0 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2s_1S0 (); }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2s_3S1 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2s_3S1 (); }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_1P1 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2p_1P1 (); }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3P0 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2p_3P0 (); }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3P1 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2p_3P1 (); }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3P2 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2p_3P2 (); }

NCM_INLINE gdouble ncm_c_HeI_ion_wn_2p_3Pmean (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () - ncm_c_HeI_Lyman_wn_2p_3Pmean (); }

/* Ionization energy: E */

NCM_INLINE gdouble ncm_c_HeI_ion_E_1s_1S0 (void)
{ return ncm_c_HeI_ion_wn_1s_1S0 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2s_1S0 (void)
{ return  ncm_c_HeI_ion_wn_2s_1S0 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2s_3S1 (void)
{ return  ncm_c_HeI_ion_wn_2s_3S1 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_1P1 (void)
{ return  ncm_c_HeI_ion_wn_2p_1P1 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3P0 (void)
{ return  ncm_c_HeI_ion_wn_2p_3P0 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3P1 (void)
{ return  ncm_c_HeI_ion_wn_2p_3P1 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3P2 (void)
{ return  ncm_c_HeI_ion_wn_2p_3P2 () * ncm_c_hc (); }

NCM_INLINE gdouble ncm_c_HeI_ion_E_2p_3Pmean (void)
{ return  ncm_c_HeI_ion_wn_2p_3Pmean () * ncm_c_hc (); }

/* Lyman series wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2s_1S0 (void)
{ return 1.66277440141e7; }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2s_3S1 (void)
{ return 1.598559743297e7; }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_1P1 (void)
{ return 1.71134896946e7; }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3P0 (void)
{ return 1.690878308131e7; }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3P1 (void)
{ return 1.690868428979e7; }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3P2 (void)
{ return 1.690867664725e7; }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wn_2p_3Pmean (void)
{ return (ncm_c_HeI_Lyman_wn_2p_3P0 () + 
          ncm_c_HeI_Lyman_wn_2p_3P1 () + 
          ncm_c_HeI_Lyman_wn_2p_3P2 ()) / 3.0; }

/* Lyman series wavelength: wl */

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2s_1S0 (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2s_1S0 (); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2s_3S1 (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2s_3S1 (); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_1P1 (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2p_1P1 (); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3P0 (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2p_3P0 (); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3P1 (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2p_3P1 (); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3P2 (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2p_3P2 (); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl_2p_3Pmean (void)
{ return 1.0 / ncm_c_HeI_Lyman_wn_2p_3Pmean (); }

/* Lyman series factor: wl^3 / (8pi) */

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2s_1S0 (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2s_1S0 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2s_3S1 (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2s_3S1 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_1P1 (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2p_1P1 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3P0 (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2p_3P0 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3P1 (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2p_3P1 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3P2 (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2p_3P2 ()) / (8.0 * ncm_c_pi ()); }

NCM_INLINE gdouble ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean (void)
{ return gsl_pow_3 (ncm_c_HeI_Lyman_wl_2p_3Pmean ()) / (8.0 * ncm_c_pi ()); }

/* Boltzmann factor */

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_1s_1S0 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_1s_1S0 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2s_1S0 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2s_1S0 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2s_3S1 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2s_3S1 () / (ncm_c_kb () * T)); }

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_1P1 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2p_1P1 () / (ncm_c_kb () * T)); } 

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3P0 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2p_3P0 () / (ncm_c_kb () * T)); } 

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3P1 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2p_3P1 () / (ncm_c_kb () * T)); } 

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3P2 (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2p_3P2 () / (ncm_c_kb () * T)); } 

NCM_INLINE gdouble ncm_c_boltzmann_factor_HeI_2p_3Pmean (const gdouble T)
{ return gsl_pow_3 (ncm_c_thermal_wn_e ()) *
 exp (- ncm_c_HeI_ion_E_2p_3Pmean () / (ncm_c_kb () * T)); } 

/* Balmer series wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeI_Balmer_wn_2p_1P1_2s_1S0 (void)
{ return ncm_c_HeI_Lyman_wn_2p_1P1 () - ncm_c_HeI_Lyman_wn_2s_1S0 (); }

NCM_INLINE gdouble ncm_c_HeI_Balmer_wn_2p_3Pmean_2s_3S1 (void)
{ return ncm_c_HeI_Lyman_wn_2p_3Pmean () - ncm_c_HeI_Lyman_wn_2s_3S1 (); }

/* Balmer series: E / k_B */

NCM_INLINE gdouble ncm_c_HeI_Balmer_E_kb_2p_1P1_2s_1S0 (void)
{ return ncm_c_HeI_Balmer_wn_2p_1P1_2s_1S0 () * ncm_c_hc () / ncm_c_kb (); }

NCM_INLINE gdouble ncm_c_HeI_Balmer_E_kb_2p_3Pmean_2s_3S1 (void)
{ return ncm_c_HeI_Balmer_wn_2p_3Pmean_2s_3S1 () * ncm_c_hc () / ncm_c_kb (); }

/*******************************************************************************
 * -- END: Helium I
 *******************************************************************************/
/*******************************************************************************
 * -- START: Helium II
 *******************************************************************************/

/* Ionization energy wavenumber: wn */

NCM_INLINE gdouble ncm_c_HeII_ion_wn_1s_2S0_5 (void)
{ return 4.389088785e7; }

/* Ionization energy: E */

NCM_INLINE gdouble ncm_c_HeII_ion_E_1s_2S0_5 (void)
{ return ncm_c_HeII_ion_wn_1s_2S0_5 () * ncm_c_hc (); }

/*******************************************************************************
 * -- END: Helium II
 *******************************************************************************/
/*******************************************************************************
 * END: NIST Atomic Spectra database
 *******************************************************************************/

/*******************************************************************************
 * Constants from other sources
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_decay_H_rate_2s_1s (void)
{ return 8.2245809; }

NCM_INLINE gdouble ncm_c_decay_He_rate_2s_1s (void)
{ return 51.3; }

/*******************************************************************************
 * START: Statistics
 *******************************************************************************/

NCM_INLINE long double ncm_c_stats_1sigma (void)
{ return 0.6826894921370858971704650912640758449558L; }

NCM_INLINE long double ncm_c_stats_2sigma (void)
{ return 0.9544997361036415855994347256669331250564L; }

NCM_INLINE long double ncm_c_stats_3sigma (void)
{ return 0.9973002039367398109466963704648100452443L; }

/*******************************************************************************
 * END: Statistics
 *******************************************************************************/

/*******************************************************************************
 * START: Observational data
 *******************************************************************************/

NCM_INLINE gdouble ncm_c_wmap5_coadded_I_K (void)
{ return 1.436; }

NCM_INLINE gdouble ncm_c_wmap5_coadded_I_Ka (void)
{ return 1.470; }

NCM_INLINE gdouble ncm_c_wmap5_coadded_I_Q (void)
{ return 2.197; }

NCM_INLINE gdouble ncm_c_wmap5_coadded_I_V (void)
{ return 3.133; }

NCM_INLINE gdouble ncm_c_wmap5_coadded_I_W (void)
{ return 6.538; }

NCM_INLINE gdouble ncm_c_hubble_cte_wmap (void)
{ return 73.0; }

NCM_INLINE gdouble ncm_c_hubble_cte_hst (void)
{ return 72.0; }

NCM_INLINE gdouble ncm_c_hubble_cte_msa (void)
{ return 68.0; }

NCM_INLINE gdouble ncm_c_hubble_radius_hm1_Mpc (void)
{ return ncm_c_c () / (100.0e3); }

NCM_INLINE gdouble ncm_c_hubble_radius_hm1_planck (void)
{ return ncm_c_hubble_radius_hm1_Mpc () * ncm_c_Mpc () / ncm_c_planck_length (); }

NCM_INLINE gdouble ncm_c_crit_density_h2 (void)
{ return 3.0 * pow (ncm_c_c () / (10.0 * ncm_c_pc ()), 2.0) / (8.0 * ncm_c_pi () * ncm_c_G ()); }

NCM_INLINE gdouble ncm_c_crit_mass_density_h2 (void)
{ return 3.0 * pow (1.0 / (10.0 * ncm_c_pc ()), 2.0) / (8.0 * ncm_c_pi () * ncm_c_G ()); }

NCM_INLINE gdouble ncm_c_crit_mass_density_h2_solar_mass_Mpc3 (void)
{ return ncm_c_crit_mass_density_h2 () / ncm_c_mass_solar () * gsl_pow_3 (ncm_c_Mpc ()); }

NCM_INLINE gdouble ncm_c_crit_number_density_p (void)
{ return ncm_c_crit_density_h2 () / ncm_c_rest_energy_p (); }

NCM_INLINE gdouble ncm_c_crit_number_density_n (void)
{ return ncm_c_crit_density_h2 () / ncm_c_rest_energy_n (); }

NCM_INLINE gdouble ncm_c_blackbody_energy_density (void)
{ return 4.0 * (gsl_pow_2 (ncm_c_pi ()) * gsl_pow_4 (ncm_c_kb ()) / (60.0 * gsl_pow_3 (ncm_c_h () / (2.0 * ncm_c_pi ())) * gsl_pow_2 (ncm_c_c ()))) / ncm_c_c (); }

NCM_INLINE gdouble ncm_c_blackbody_per_crit_density_h2 (void)
{ return ncm_c_blackbody_energy_density () / ncm_c_crit_density_h2 (); }

NCM_INLINE gdouble ncm_c_radiation_temp_to_h2omega_r (const gdouble T)
{ return ncm_c_blackbody_per_crit_density_h2 () * gsl_pow_4 (T); }

NCM_INLINE gdouble ncm_c_radiation_h2Omega_r0_to_temp (const gdouble omr)
{ return pow (omr / ncm_c_blackbody_per_crit_density_h2 (), 0.25); }

/*******************************************************************************
 * END: Observational data
 *******************************************************************************/

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_C_INLINE_H_ */
