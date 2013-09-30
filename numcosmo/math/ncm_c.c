/***************************************************************************
 *            ncm_c.c
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

/**
 * SECTION:ncm_c
 * @title: Numerical and Physical Constants
 * @short_description: Numerical constants
 *
 * Mathematical and physical constants and constants manipulation
 * functions.
 * Using 2006 CODATA recommended values, see constants.txt.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcmC, ncm_c, G_TYPE_OBJECT);

static void
ncm_c_init (NcmC *ncm_c)
{
  NCM_UNUSED (ncm_c);
}

static void
ncm_c_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_c_parent_class)->finalize (object);
}

static void
ncm_c_class_init (NcmCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = ncm_c_finalize;
}

/*******************************************************************************
 * Mathematical constants
 *******************************************************************************/

/**
 * ncm_c_sqrt_1_4pi:
 *
 * Returns: sqrt (1 / (4 * pi))
 *
 */
/**
 * ncm_c_sqrt_2pi:
 *
 * Returns: sqrt (2 * pi)
 *
 */
/**
 * ncm_c_sqrt_3_4pi:
 *
 * Returns: sqrt (3 / (4 * pi))
 *
 */
/**
 * ncm_c_lnpi_4:
 *
 * Returns: ln (pi) / 4
 * 1.8378770664093454835606594728112352797227949472755668256343
 */
/**
 * ncm_c_ln2pi:
 *
 * Returns: ln (2pi)
 * 
 */
/**
 * ncm_c_pi:
 *
 * Returns: pi
 *
 */
/**
 * ncm_c_tan_1arcsec:
 *
 * Returns: tan (2 * pi/ (360 * 60 * 60))
 *
 */

/**
 * ncm_c_degree_to_radian:
 * @d: anle in degrees.
 *
 * Returns: d * pi / 180
 *
 */
/**
 * ncm_c_radian_to_degree:
 * @r: angle in radians 
 *
 * Returns: r * 180 / pi
 *
 */
/**
 * ncm_c_radian_0_2pi:
 * @r: angle in radians
 *
 * Returns: the angle in the interval [0, 2pi]
 *
 */
/**
 * ncm_c_sign_sin:
 * @r: angle in radias
 *
 * Returns: the sign of the value of sin(d).
 *
 */

/*******************************************************************************
 * START: 2006 CODATA recommended values (see end of file)
 *******************************************************************************/

/**
 * ncm_c_c:
 *
 * Returns: Speed of light.
 *
 */
/**
 * ncm_c_h:
 *
 * Returns: Planck constant.
 *
 */
/**
 * ncm_c_hbar:
 *
 * Returns: Planck constant over 2 pi.
 *
 */
/**
 * ncm_c_fine_struct:
 *
 * Returns: Fine structure constant.
 *
 */
/**
 * ncm_c_kb:
 *
 * Returns: Boltzmann constant.
 *
 */
/**
 * ncm_c_G:
 *
 * Returns: Newton constant.
 *
 */
/**
 * ncm_c_planck_length:
 *
 * Returns: Planck length.
 *
 */
/**
 * ncm_c_thomson_cs:
 *
 * Returns: Thomson cross section.
 *
 */
/**
 * ncm_c_stefan_boltzmann:
 *
 * Returns: Stefan Boltzmann constant.
 *
 */
/**
 * ncm_c_mass_e:
 *
 * Returns: Electron mass.
 *
 */
/**
 * ncm_c_mass_p:
 *
 * Returns: Proton mass.
 *
 */
/**
 * ncm_c_mass_n:
 *
 * Returns: Neuton mass.
 *
 */
/**
 * ncm_c_mass_ratio_alpha_p:
 *
 * Returns: The proton and alpha particle (helium-4) mass ratio.
 *
 */

/*******************************************************************************
 * END: 2006 CODATA recommended values
 *******************************************************************************/

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

/**
 * ncm_c_hc:
 *
 * Returns: Planck constant times the speed of light.
 *
 */
/**
 * ncm_c_fine_struct_square:
 *
 * Returns: The square of the fine struct constant.
 *
 */
/**
 * ncm_c_kpc:
 *
 * Returns: One kilo parsec.
 *
 */
/**
 * ncm_c_Mpc:
 *
 * Returns: One mega parsec.
 *
 */
/**
 * ncm_c_AR:
 *
 * Returns: Radiation constant AR.
 *
 */
/**
 * ncm_c_c2:
 *
 * Returns: Square of the speed of light.
 *
 */
/**
 * ncm_c_planck_length2:
 *
 * Returns: Square of the Planck length.
 *
 */
/**
 * ncm_c_rest_energy_e:
 *
 * Returns: Electron's rest energy.
 *
 */
/**
 * ncm_c_rest_energy_p:
 *
 * Returns: Proton's rest energy.
 *
 */
/**
 * ncm_c_rest_energy_n:
 *
 * Returns: Neutron's rest energy.
 *
 */

/*******************************************************************************
 * Constants from other places
 *******************************************************************************/

/**
 * ncm_c_decay_H_rate_2s_1s:
 *
 * FIXME: Cite source.
 *
 * Returns: Decay rate of Hydrogen from 2s -> 1s.
 *
 */
/**
 * ncm_c_decay_He_rate_2s_1s:
 *
 * FIXME: Cite source.
 *
 * Returns: Decay rate of Helium from 2s -> 1s.
 *
 */
/**
 * ncm_c_HeI_bind_1s:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI binding energy 1s.
 *
 */
/**
 * ncm_c_HeII_bind_1s:
 *
 * FIXME: Cite source.
 *
 * Returns: HeII binding energy 1s.
 *
 */
/**
 * ncm_c_HeI_Lyman_2s:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI Lyman 2s energy.
 *
 */
/**
 * ncm_c_HeI_Lyman_2p:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI Lyman 2p energy.
 *
 */
/**
 * ncm_c_HeI_Lyman_2s_wl:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI Lyman 2s wave length.
 *
 */
/**
 * ncm_c_HeI_Lyman_2p_wl:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI Lyman 2p wave length.
 *
 */
/**
 * ncm_c_HeI_bind_2s:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI binding energy 2s.
 *
 */
/**
 * ncm_c_HeI_bind_2p:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI binding energy 2p.
 *
 */
/**
 * ncm_c_HeI_2s_m_2p:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI energy difference between states 2s and 2p.
 *
 */
/**
 * ncm_c_HeI_2s_m_2p_kb:
 *
 * FIXME: Cite source.
 *
 * Returns: HeI energy difference between states 2s and 2p divided by Boltzmann constant.
 *
 */
/**
 * ncm_c_HeI_Lyman_2s_wl3_8pi:
 *
 * FIXME: Cite source.
 *
 * Returns: Cubic power of HeI Lyman 2s wave length divided by (8 * pi).
 *
 */
/**
 * ncm_c_HeI_Lyman_2p_wl3_8pi:
 *
 * FIXME: Cite source.
 *
 * Returns: Cubic power of HeI Lyman 2p wave length divided by (8 * pi).
 *
 */
/**
 * ncm_c_H_reduced_mass:
 *
 * FIXME: Cite source.
 *
 * Returns: Hydrogen reduced mass.
 *
 */
/**
 * ncm_c_H_reduced_energy:
 *
 * FIXME: Cite source.
 *
 * Returns: Hydrogen reduced energy.
 *
 */
/**
 * ncm_c_H_bind:
 * @n: FIXME
 * @j: FIXME
 * FIXME: Cite source.
 *
 * Returns: Hydrogen binding energy.
 *
 */
/**
 * ncm_c_H_bind_1s:
 *
 * FIXME: Cite source.
 *
 * Returns: Hydrogen 1s binding energy.
 *
 */
/**
 * ncm_c_H_bind_2s:
 *
 * FIXME: Cite source.
 *
 * Returns: Hydrogen 2s binding energy.
 *
 */
/**
 * ncm_c_H_bind_2p:
 *
 * FIXME: Cite source.
 *
 * Returns: Hydrogen 2p binding energy.
 *
 */
/**
 * ncm_c_H_Lyman_series:
 * @n: FIXME
 * @j: FIXME
 *
 * Energy difference between levels 1s and n,j.
 * FIXME: Cite source.
 *
 * Returns: Hydrogen Lyman series.
 *
 */
/**
 * ncm_c_H_Lyman_2s:
 *
 * FIXME: Cite source.
 *
 * Returns: Energy difference between levels 1s and 2s.
 *
 */
/**
 * ncm_c_H_Lyman_2p:
 *
 * FIXME: Cite source.
 *
 * Returns: Energy difference between levels 1s and 2p.
 *
 */
/**
 * ncm_c_H_Lyman_series_wl:
 * @n: FIXME
 * @j: FIXME
 *
 * FIXME: Cite source.
 *
 * Returns: Wavelenght relative to the energy difference between levels 1s and n,j.
 *
 */
/**
 * ncm_c_H_Lyman_2s_wl:
 *
 * FIXME: Cite source.
 *
 * Returns: Wavelenght relative to the energy difference between levels 1s and 2s.
 *
 */
/**
 * ncm_c_H_Lyman_2p_wl:
 *
 * FIXME: Cite source.
 *
 * Returns: Wavelenght relative to the energy difference between levels 1s and 2s.
 *
 */
/**
 * ncm_c_H_Lyman_2s_wl3_8pi:
 *
 * FIXME: Cite source.
 *
 * Returns: Cubic power of the Wavelenght relative to the energy difference between levels 1s and 2s divided by (8*pi).
 *
 */
/**
 * ncm_c_H_Lyman_2p_wl3_8pi:
 *
 * FIXME: Cite source.
 *
 * Returns: Cubic power of the wavelenght relative to the energy difference between levels 1s and 2p divided by (8*pi).
 *
 */
/**
 * ncm_c_thermal_wl_e:
 *
 * FIXME: Cite source.
 *
 * Returns: Thermal electron wavelenght.
 *
 */
/**
 * ncm_c_thermal_wl_p:
 *
 * FIXME: Cite source.
 *
 * Returns: Thermal proton wavelenght.
 *
 */
/**
 * ncm_c_thermal_wl_n:
 *
 * FIXME: Cite source.
 *
 * Returns: Thermal neutron wavelenght.
 *
 */
/**
 * ncm_c_thermal_wn_e:
 *
 * FIXME: Cite source.
 *
 * Returns: Thermal eletron wavenumber.
 *
 */
/**
 * ncm_c_thermal_wn_p:
 *
 * FIXME: Cite source.
 *
 * Returns: Thermal proton wavenumber.
 *
 */
/**
 * ncm_c_thermal_wn_n:
 *
 * FIXME: Cite source.
 *
 * Returns: Thermal neutron wavenumber.
 *
 */
/**
 * ncm_c_boltzmann_factor_H_1s:
 * @T: temperature.
 *
 * FIXME: Cite source.
 *
 * Returns: Boltzmann factor for Hydrogen 1s level.
 *
 */
/**
 * ncm_c_boltzmann_factor_H_2s:
 * @T: temperature.
 *
 * FIXME: Cite source.
 *
 * Returns: Boltzmann factor for Hydrogen 2s level.
 *
 */
/**
 * ncm_c_boltzmann_factor_H_2p:
 * @T: temperature.
 *
 * FIXME: Cite source.
 *
 * Returns: Boltzmann factor for Hydrogen 2p level.
 *
 */
/**
 * ncm_c_boltzmann_factor_HeI_1s:
 * @T: temperature.
 *
 * FIXME: Cite source.
 *
 * Returns: Boltzmann factor for HeI 1s level.
 *
 */
/**
 * ncm_c_boltzmann_factor_HeI_2s:
 * @T: temperature.
 *
 * FIXME: Cite source.
 *
 * Returns: Boltzmann factor for HeI 2s level.
 *
 */
/**
 * ncm_c_boltzmann_factor_HeI_2p:
 * @T: temperature.
 *
 * FIXME: Cite source.
 *
 * Returns: Boltzmann factor for HeI 2p level.
 *
 */
/**
 * ncm_c_AU:
 *
 * Returns: Astronomical unit (http://ssd.jpl.nasa.gov/?constants).
 *
 */
/**
 * ncm_c_pc:
 *
 * Returns: Parsec unit 1 AU / tan (1 arcsec) - Copied from CAMB/constants.f90 to facilitate comparison.
 *
 */
/**
 * ncm_c_mass_solar:
 *
 * Returns: One solar mass.
 *
 */

/* Statistics */

/**
 * ncm_c_stats_1sigma:
 *
 * The integral of a gaussian distribution with mean mu
 * and standard deviation sigma in (mu - 1 * sigma, mu + 1 * sigma)
 *
 * Returns: P (mu - 1 * sigma, mu + 1 * sigma)
 *
 */
/**
 * ncm_c_stats_2sigma:
 *
 * The integral of a gaussian distribution with mean mu
 * and standard deviation sigma in (mu - 2 * sigma, mu + 2 * sigma)
 *
 * Returns: P (mu - 2 * sigma, mu + 2 * sigma)
 *
 */
/**
 * ncm_c_stats_3sigma:
 *
 * The integral of a gaussian distribution with mean mu
 * and standard deviation sigma in (mu - 3 * sigma, mu + 3 * sigma)
 *
 * Returns: P (mu - 3 * sigma, mu + 3 * sigma)
 *
 */

/*******************************************************************************
 * Observational data
 *******************************************************************************/

/**
 * ncm_c_wmap3_cmb_z:
 *
 * Returns: Wmap3 last scatering redshift.
 *
 */
/**
 * ncm_c_wmap3_cmb_R:
 *
 * Returns: Wmap3 last scatering shift parameter.
 *
 */
/**
 * ncm_c_wmap3_cmb_sigma_R:
 *
 * Returns: Wmap3 last scatering shift parameter standard deviation.
 *
 */
/**
 * ncm_c_wmap5_cmb_z:
 *
 * Returns: Wmap5 last scatering redshift.
 *
 */
/**
 * ncm_c_wmap5_cmb_R:
 *
 * Returns: Wmap5 last scatering shift parameter.
 *
 */
/**
 * ncm_c_wmap5_cmb_sigma_R:
 *
 * Returns: Wmap5 last scatering shift parameter standard deviation.
 *
 */
/**
 * ncm_c_wmap5_coadded_I_K:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_wmap5_coadded_I_Ka:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_wmap5_coadded_I_Q:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_wmap5_coadded_I_V:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_wmap5_coadded_I_W:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_eisenstein_z:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_eisenstein_A:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_eisenstein_sigma_A:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_eisenstein_DV:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_eisenstein_sigma_DV:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_percival_DV_DV:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_bao_percival_sigma_DV_DV:
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_hubble_cte_wmap:
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_hubble_cte_hst:
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_hubble_cte_msa:
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_neutrino_n_eff:
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
/**
 * ncm_c_prim_He_Yp:
 *
 * The primoridial helium mass fraction 
 * $$Y_p = \frac{m_\text{He}n_\text{He}}
 * {m_\text{He}n_\text{He}+m_\text{H}n_\text{H}},$$ where $m_\text{He}$, 
 * n_\text{He}, m_\text{H} and m_\text{H} are respectively helium mass and 
 * number density and hydrogen mass and number density.
 *
 * Returns: The primordial helium mass abundance.
 *
 */
/**
 * ncm_c_prim_H_Yp:
 *
 * The primordial hydrogen mass fraction $$Y_{\text{H}p} = 1 - Y_p,$$
 * where $Y_p$ is the helium mass fraction, see ncm_c_prim_He_Yp ().
 *
 * Returns: The primordial hydrogen mass abundance.
 *
 */
/**
 * ncm_c_prim_XHe:
 * 
 * The primordial helium to hydrogen ratio $$X_\text{He} = 
 * \frac{n_\text{He}}{n_\text{H}} = \frac{m_\text{H}}{m_\text{He}}
 * \frac{Y_p}{Y_{\text{H}p}},$$ see ncm_c_prim_H_Yp () and ncm_c_prim_He_Yp ().
 * 
 * Returns: The primordial helium to hydrogen ratio.
 *
 */
/**
 * ncm_c_hubble_radius:
 *
 * FIXME
 *
 * Returns: Hubble radius
 *
 */
/**
 * ncm_c_hubble_radius_planck:
 *
 * FIXME
 *
 * Returns: Hubble radius
 *
 */
/**
 * ncm_c_crit_density:
 *
 * FIXME
 *
 * Returns: Critical density in ... units.
 *
 */
/**
 * ncm_c_crit_mass_density:
 *
 * FIXME
 *
 * Returns: Critical mass density in ... units.
 *
 */
/**
 * ncm_c_crit_mass_density_solar_Mpc:
 *
 * FIXME
 *
 * Returns: Critical mass density in ... units.
 *
 */
/**
 * ncm_c_crit_number_density_p:
 *
 * FIXME
 *
 * Returns: Critical proton number density in ... units.
 *
 */
/**
 * ncm_c_crit_number_density_n:
 *
 * FIXME
 *
 * Returns: Critical neutron number density in ... units.
 *
 */
/**
 * ncm_c_blackbody_energy_density:
 *
 * FIXME
 *
 * Returns: Blackbody energy density in ... units.
 *
 */
/**
 * ncm_c_radiation_temp_to_h2Omega_r:
 * @T: FIXME
 *
 * FIXME
 *
 * Returns: .
 *
 */
/**
 * ncm_c_radiation_h2Omega_r_to_temp:
 * @omr: FIXME
 *
 * FIXME
 *
 * Returns: .
 *
 */

