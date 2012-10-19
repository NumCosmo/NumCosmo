/***************************************************************************
 *            nc_constants.h
 *
 *  Wed Oct 15 17:31:25 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:nc_constants
 * @short_description: Physical, mathematical and observations constants
 *
 * FIXME
 */

#ifndef _NC_CONSTANTS_H
#define _NC_CONSTANTS_H

#include <glib.h>
#include <gsl/gsl_const_num.h>

/**
 * SECTION:nc_constants
 * @title: Numerical and Physical Constants
 * @short_description: FIXME
 *
 * FIXME
 *
 */

G_BEGIN_DECLS

/* Mathematical constants */
#define NC_C_SQRT_1_4PI   (0.28209479177387814347403972578038630L)
#define NC_C_SQRT_3_4PI   (0.48860251190291992158638462283834700L)
#define NC_C_LNPI_4       (0.28618247146235004353585683783826468L)
#define NC_C_PIL          (3.1415926535897932384626433832795029L)
#define NC_C_TAN_1ARCSECL (4.8481368111333441675396429478852853e-6L)
#define NC_C_TAN_1ARCSEC  (4.8481368111333441675396429478852853e-6)

/* START: 2006 CODATA recommended values (see end of file) */

#define NC_C_c (299792458.0)                       /**< Speed of light */
#define NC_C_h (6.62606896e-34)                    /**< Planck constant */
#define NC_C_hbar (1.054571628e-34)                /**< Planck constant over 2pi */
#define NC_C_FINE_STRUCT (7.2973525376e-3)         /**< Fine structure constant */
#define NC_C_kb (1.3806504e-23)                    /**< Boltzmann constant */
#define NC_C_G (6.67428e-11)                       /**< Newton constant */
#define NC_C_PLANCK_LENGTH (1.616252e-35)          /**< Planck length */
#define NC_C_THOMPSON_CS (0.6652458558e-28)        /**< Thompson cross section */
#define NC_C_STEFAN_BOLTZMANN (5.670400e-8)        /**< Stefan Boltzmann constant */
#define NC_C_MASS_e (9.10938215e-31)               /**< Electron mass */
#define NC_C_MASS_p (1.672621637e-27)              /**< Proton mass */
#define NC_C_MASS_n (1.674927211e-27)              /**< Neuton mass */
#define NC_C_MASS_RATIO_alpha_p (3.97259968951)    /**< The proton and alpha particle (helium-4) mass ratio. */

/* END: 2006 CODATA recommended values */

#define NC_C_AU (1.49597870691e11)                 /**< Astronomical unit (http://ssd.jpl.nasa.gov/?constants) */
#define NC_C_PARSEC (3.085678e16)                  /**< Parsec unit 1 AU / tan (1 arcsec) - Copied from CAMB/constants.f90 to facilitate comparison */
/*#define NC_C_PARSEC (NC_C_AU / NC_C_TAN_1ARCSEC)*/   /**< Parsec unit 1 AU / tan (1 arcsec)*/
#define NC_C_MASS_SOLAR (1.98892e30)               /**< FIXME */

#define NC_C_hc (NC_C_c * NC_C_h)
#define NC_C_FINE_STRUCT_SQUARE (NC_C_FINE_STRUCT * NC_C_FINE_STRUCT)
#define NC_C_KILO_PARSEC (GSL_CONST_NUM_KILO * NC_C_PARSEC)
#define NC_C_MEGA_PARSEC (GSL_CONST_NUM_MEGA * NC_C_PARSEC)
#define NC_C_AR (4.0 * NC_C_STEFAN_BOLTZMANN / NC_C_c)
#define NC_C_c2 (NC_C_c * NC_C_c)
#define NC_C_PLANCK_LENGTH2 (NC_C_PLANCK_LENGTH * NC_C_PLANCK_LENGTH)

#define NC_C_ENERGY_e (NC_C_MASS_e * NC_C_c2)
#define NC_C_ENERGY_p (NC_C_MASS_p * NC_C_c2)
#define NC_C_ENERGY_n (NC_C_MASS_n * NC_C_c2)

/* Statistics */
#define NC_C_STATS_1SIGMA 0.6827
#define NC_C_STATS_2SIGMA 0.9545
#define NC_C_STATS_3SIGMA 0.9973

/* Observational data  */
#define NC_C_WMAP3_REDSHIFT 1089.0
#define NC_C_WMAP3_R 1.70
#define NC_C_WMAP3_SIGMA_R 0.03

#define NC_C_WMAP5_REDSHIFT 1090.0
#define NC_C_WMAP5_R 1.71
#define NC_C_WMAP5_SIGMA_R 0.013

#define NC_C_WMAP5_COADDED_I_K  1.436
#define NC_C_WMAP5_COADDED_I_Ka 1.470
#define NC_C_WMAP5_COADDED_I_Q  2.197
#define NC_C_WMAP5_COADDED_I_V  3.133
#define NC_C_WMAP5_COADDED_I_W  6.538

#define NC_C_BAO_EISENSTEIN_REDSHIFT 0.35
#define NC_C_BAO_EISENSTEIN_A 0.469
#define NC_C_BAO_EISENSTEIN_SIGMA_A 0.017
#define NC_C_BAO_EISENSTEIN_DV 1334.0
#define NC_C_BAO_EISENSTEIN_SIGMA_DV 88.0

#define NC_C_BAO_PERCIVAL_DV_DV 1.812
#define NC_C_BAO_PERCIVAL_SIGMA_DV_DV 0.060

#define NC_C_HUBBLE_CTE_WMAP 73.0
#define NC_C_HUBBLE_CTE_HST  72.0
#define NC_C_HUBBLE_CTE_MSA  68.0

#define NC_C_NEUTRINO_N_EFF 3.04

#define NC_C_PRIM_HE_Y_P 0.24                     /**< The primordial helium abundance.    */
#define NC_C_PRIM_H_FRAC (1.0 - NC_C_PRIM_HE_Y_P) /**< The primordial hydrogen abundance.  */
#define NC_C_PRIM_HE_XHe (NC_C_PRIM_HE_Y_P / (NC_C_MASS_RATIO_alpha_p * (1.0 - NC_C_PRIM_HE_Y_P)))

#define NC_C_DECAY_H_2S_1S 8.22458
#define NC_C_DECAY_He_2S_1S 51.3

#define NC_C_HeI_BINDING_ENERGY_1s (1.98310772e7 * NC_C_hc)
#define NC_C_HeII_BINDING_ENERGY_1s (4.389088863e7 * NC_C_hc)

#define NC_C_HeI_LYMAN_2s (1.66277434e7 * NC_C_hc)
#define NC_C_HeI_LYMAN_2p (1.71134891e7 * NC_C_hc)

#define NC_C_HeI_LYMAN_2s_WL (1.0 / 1.66277434e7)
#define NC_C_HeI_LYMAN_2p_WL (1.0 / 1.71134891e7)

#define NC_C_HeI_BINDING_ENERGY_2s (3.2033338e6 * NC_C_hc)
#define NC_C_HeI_BINDING_ENERGY_2p (2.7175881e6 * NC_C_hc)

#define NC_C_HeI_2s_m_2p ((3.2033338e6 - 2.7175881e6) * NC_C_hc)
#define NC_C_HeI_2s_m_2p_Kb ((3.2033338e6 - 2.7175881e6) * NC_C_hc / NC_C_kb)

#define NC_C_HeI_LYMAN_2s_WL3_8PI \
(NC_HeI_LYMAN_2s_WL*NC_HeI_LYMAN_2s_WL*NC_HeI_LYMAN_2s_WL / (8.0 * M_PI))
#define NC_C_HeI_LYMAN_2p_WL3_8PI \
(NC_C_HeI_LYMAN_2p_WL*NC_C_HeI_LYMAN_2p_WL*NC_C_HeI_LYMAN_2p_WL / (8.0 * M_PI))

#define NC_C_HYDROGEN_REDUCED_MASS (NC_C_MASS_e / (1.0 + NC_C_MASS_e / NC_C_MASS_p))
#define NC_C_HYDROGEN_REDUCED_ENERGY (NC_C_HYDROGEN_REDUCED_MASS * NC_C_c2)

#define NC_C_HYDROGEN_BINDING_ENERGY(n,j) \
(NC_C_HYDROGEN_REDUCED_ENERGY *  \
 (1.0 - 1.0 / sqrt (1.0 + \
                    NC_C_FINE_STRUCT_SQUARE / \
                    pow (n - j - 0.5 + sqrt(pow(j + 0.5, 2.0) - NC_C_FINE_STRUCT_SQUARE), 2.0) \
                    )) \
)

#define NC_C_HYDROGEN_BINDING_ENERGY_1s NC_C_HYDROGEN_BINDING_ENERGY(1.0,0.5)
#define NC_C_HYDROGEN_BINDING_ENERGY_2s NC_C_HYDROGEN_BINDING_ENERGY(2.0,0.5)
#define NC_C_HYDROGEN_BINDING_ENERGY_2p NC_C_HYDROGEN_BINDING_ENERGY(2.0,1.5)

#define NC_C_HYDROGEN_LYMAN_SERIE(n,j) \
(NC_C_HYDROGEN_BINDING_ENERGY(1,0.5) - NC_C_HYDROGEN_BINDING_ENERGY(n,j))

#define NC_C_HYDROGEN_LYMAN_2s NC_C_HYDROGEN_LYMAN_SERIE(2,0.5)
#define NC_C_HYDROGEN_LYMAN_2p NC_C_HYDROGEN_LYMAN_SERIE(2,1.5)

#define NC_C_HYDROGEN_LYMAN_SERIE_WL(n,j) (NC_C_hc /NC_C_HYDROGEN_LYMAN_SERIE(n,j))

#define NC_C_HYDROGEN_LYMAN_2s_WL (NC_C_hc/NC_C_HYDROGEN_LYMAN_SERIE(2,0.5)) /**< 2s -> 1s transition wavelength */
#define NC_C_HYDROGEN_LYMAN_2p_WL (NC_C_hc/NC_C_HYDROGEN_LYMAN_SERIE(2,1.5)) /**< 2p -> 1s transition wavelength */

#define NC_C_HYDROGEN_LYMAN_2s_WL3_8PI \
(NC_C_HYDROGEN_LYMAN_2s_WL*NC_C_HYDROGEN_LYMAN_2s_WL*NC_C_HYDROGEN_LYMAN_2s_WL / (8.0 * M_PI))
#define NC_C_HYDROGEN_LYMAN_2p_WL3_8PI \
(NC_C_HYDROGEN_LYMAN_2p_WL*NC_C_HYDROGEN_LYMAN_2p_WL*NC_C_HYDROGEN_LYMAN_2p_WL / (8.0 * M_PI))

#define NC_C_THERMAL_WAVELENGTH_e (sqrt((2.0 * M_PI * NC_C_hbar * NC_C_hbar ) / (NC_C_MASS_e * NC_C_kb)))
#define NC_C_THERMAL_WAVELENGTH_p (sqrt((2.0 * M_PI * NC_C_hbar * NC_C_hbar ) / (NC_C_MASS_p * NC_C_kb)))
#define NC_C_THERMAL_WAVELENGTH_n (sqrt((2.0 * M_PI * NC_C_hbar * NC_C_hbar ) / (NC_C_MASS_n * NC_C_kb)))

#define NC_C_THERMAL_WAVENUMBER_e (1.0/NC_C_THERMAL_WAVELENGTH_e)
#define NC_C_THERMAL_WAVENUMBER_p (1.0/NC_C_THERMAL_WAVELENGTH_p)
#define NC_C_THERMAL_WAVENUMBER_n (1.0/NC_C_THERMAL_WAVELENGTH_n)

#define NC_C_BOLTZMAN_FACTOR_H_1s(T) \
(pow(NC_C_THERMAL_WAVENUMBER_e, 3.0) * \
 exp(-NC_C_HYDROGEN_BINDING_ENERGY_1s / NC_C_kb / (T)))

#define NC_BOLTZMAN_FACTOR_H_2s(T) \
(pow(NC_C_THERMAL_WAVENUMBER_e, 3.0) * \
 exp(-NC_C_HYDROGEN_BINDING_ENERGY_2s / NC_C_kb / (T)))

#define NC_BOLTZMAN_FACTOR_H_2p(T) \
(pow(NC_C_THERMAL_WAVENUMBER_e, 3.0) * \
 exp(-NC_C_HYDROGEN_BINDING_ENERGY_2p / NC_C_kb / (T)))

#define NC_BOLTZMAN_FACTOR_HeI_1s(T) \
(pow(NC_C_THERMAL_WAVENUMBER_e, 3.0) * \
 exp(-NC_C_HeI_BINDING_ENERGY_1s / NC_C_kb / (T)))

#define NC_BOLTZMAN_FACTOR_HeI_2s(T) \
(pow(NC_C_THERMAL_WAVENUMBER_e, 3.0) * \
 exp(-NC_C_HeI_BINDING_ENERGY_2s / NC_C_kb / (T)))

#define NC_BOLTZMAN_FACTOR_HeI_2p(T) \
(pow(NC_C_THERMAL_WAVENUMBER_e, 3.0) * \
 exp(-NC_HeI_BINDING_ENERGY_2p / NC_C_kb / (T)))

#define NC_C_HUBBLE_RADIUS (NC_C_c / (100.0e3))
#define NC_C_HUBBLE_RADIUS_PLANCK (NC_C_HUBBLE_RADIUS * NC_C_MEGA_PARSEC / NC_C_PLANCK_LENGTH)

#define NC_C_CRIT_DENSITY (3.0 * pow (NC_C_c / (10.0 * NC_C_PARSEC), 2.0) / (8.0 * M_PI * NC_C_G))
#define NC_C_CRIT_MASS_DENSITY (3.0 * pow (1.0 / (10.0 * NC_C_PARSEC), 2.0) / (8.0 * M_PI * NC_C_G))
#define NC_C_CRIT_MASS_DENSITY_SOL_MPC (NC_C_CRIT_MASS_DENSITY / NC_C_MASS_SOLAR * pow (NC_C_MEGA_PARSEC, 3.0))

#define NC_C_CRIT_NUMBER_DENSITY_p (NC_C_CRIT_DENSITY / NC_C_ENERGY_p)
#define NC_C_CRIT_NUMBER_DENSITY_n (NC_C_CRIT_DENSITY / NC_C_ENERGY_n)

#define NC_C_BLACKBODY_ENERGY_DENSITY (4.0 * NC_C_STEFAN_BOLTZMANN / NC_C_c)

#define NC_C_RADIATION_TEMP_TO_h2OMEGA_R(T) (NC_C_BLACKBODY_ENERGY_DENSITY * ((T * T) * (T * T)) / NC_C_CRIT_DENSITY)
#define NC_C_RADIATION_h2OMEGA_R_TO_TEMP(omr) (pow(NC_C_CRIT_DENSITY * omr / NC_C_BLACKBODY_ENERGY_DENSITY, 0.25))

G_END_DECLS

#endif /* _NC_CONSTANTS_H */
