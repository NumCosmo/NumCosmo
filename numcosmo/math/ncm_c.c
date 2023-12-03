/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_c.c
 *
 *  Wed Oct 15 17:31:25 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * @title: NcmC
 * @short_description: Numerical and physical constants.
 * @stability: Stable
 * @include: numcosmo/math/ncm_c.h
 *
 * Mathematical and physical constants and constants manipulation
 * functions.
 *
 * The sources are:
 *
 * - High precision mathematical constants obtained from [MPFR](http://www.mpfr.org/).
 *
 * - Fundamental constants: 2018 [CODATA](http://physics.nist.gov/cuu/Constants/index.html)
 * recommended values, see constants.txt distributed with NumCosmo sources.
 *
 * - The atomic weights: Commission on Isotopic Abundances and Atomic Weights
 * ([CIAAW](http://www.ciaaw.org/atomic-weights.htm)) of the
 * International Union of Pure and Applied Chemistry (IUPAC). See also from the
 * [NIST compilation](http://www.nist.gov/pml/data/comp.cfm).
 *
 * - Astronomical constants: [IAU 2015](https://www.iau.org/administration/resolutions/general_assemblies/)
 * resolutions for the astronomical unit ncm_c_au(), parsec ncm_c_pc() and derived constants.
 * See also [Luzum 2011][XLuzum2011].
 *
 * - Atomic Spectra: National Institute of Standards and Technology (NIST) [Atomic Spectra](http://www.nist.gov/pml/data/asd.cfm)
 * Standard Reference Database 78 - Version 5.7 (October 2018).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

struct _NcmC
{
  GObject parent_instance;
};

G_DEFINE_TYPE (NcmC, ncm_c, G_TYPE_OBJECT)

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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = ncm_c_finalize;
}

/*******************************************************************************
 * Mathematical constants
 *******************************************************************************/

/**
 * ncm_c_sqrt_1_4pi:
 *
 * Returns: $\sqrt{1 / (4 \pi)}$.
 */

/**
 * ncm_c_sqrt_pi:
 *
 * Returns: $\sqrt{\pi}$.
 */

/**
 * ncm_c_sqrt_2pi:
 *
 * Returns: $\sqrt{2 \pi}$.
 */

/**
 * ncm_c_sqrt_pi_2:
 *
 * Returns: $\sqrt{\pi / 2}$.
 */

/**
 * ncm_c_sqrt_3_4pi:
 *
 * Returns: $\sqrt{3 / (4 \pi)}$.
 */

/**
 * ncm_c_ln2:
 *
 * Returns: $\ln(2)$.
 */

/**
 * ncm_c_ln3:
 *
 * Returns: $\ln(3)$.
 */

/**
 * ncm_c_lnpi_4:
 *
 * Returns: $\ln(\pi) / 4$.
 */

/**
 * ncm_c_ln2pi:
 *
 * Returns: $\ln(2\pi)$.
 */

/**
 * ncm_c_lnpi:
 *
 * Returns: $\ln(\pi)$.
 */

/**
 * ncm_c_pi:
 *
 * Returns: $\pi$.
 */

/**
 * ncm_c_2_pi_2:
 *
 * Returns: $2\pi^2$.
 */

/**
 * ncm_c_tan_1arcsec:
 *
 *
 * Returns: $\tan(2 \pi/ (360 \times 60 \times 60))$.
 */

/**
 * ncm_c_deg2_steradian:
 *
 * The convertion factor from degrees squared to steradian.
 *
 * Returns: $\pi^2/(180)^2$.
 */

/**
 * ncm_c_degree_to_radian:
 * @d: angle in degrees
 *
 * Returns: $d \times \pi / 180$.
 */

/**
 * ncm_c_radian_to_degree:
 * @r: angle in radians
 *
 * Returns: $r \times 180 / \pi$.
 */

/**
 * ncm_c_radian_0_2pi:
 * @r: angle in radians
 *
 * Returns: the angle in the interval $[0, 2\pi]$.
 */

/**
 * ncm_c_sign_sin:
 * @r: angle in radias
 *
 * Returns: the sign of the value of $\sin(r)$.
 */

/*******************************************************************************
 * START: 2018 CODATA recommended values (see constants.txt)
 *******************************************************************************/

/**
 * ncm_c_c:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Speed of light $c = 299792458.0 \,\left[\mathrm{m}\mathrm{s}^{-1}\right]$.
 *
 */

/**
 * ncm_c_h:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Planck constant $h = 6.62607015 \times 10^{-34} \,\left[\mathrm{J}\,\mathrm{s}\right]$.
 */

/**
 * ncm_c_hbar:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Planck constant over $2\pi$, $\hbar \equiv h / (2\pi) = 1.054571817 \times 10^{-34} \,\left[\mathrm{J}\,\mathrm{s}\right]$.
 */

/**
 * ncm_c_fine_struct:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Fine structure constant $\alpha = 7.2973525693 \times 10^{-3} $.
 */

/**
 * ncm_c_kb:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Boltzmann constant $k_\mathrm{B} = 1.380649 \times 10^{-23} \,\left[\mathrm{J}\,\mathrm{K}^{-1}\right]$.
 */

/**
 * ncm_c_G:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Newton's (or gravitational) constant $\mathrm{G} = 6.67430 \times 10^{-11} \,\left[\mathrm{m}^3\,\mathrm{kg}^{-1}\,\mathrm{s}^{-2}\right]$.
 */

/**
 * ncm_c_planck_length:
 *
 * Returns: Planck length $l_\mathrm{P} = 1.616255 \times 10^{-35} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_thomson_cs:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Thomson cross section $\sigma_\mathrm{T} = 0.66524587321 \times 10^{-28} \,\left[\mathrm{m}^2\right]$.
 */

/**
 * ncm_c_stefan_boltzmann:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Stefan Boltzmann constant $\sigma_\mathrm{SB} = 5.670374419 \times 10^{-8} \,\left[\mathrm{W}\,\mathrm{m}^{-2}\,\mathrm{K}^{-4}\right]$.
 */

/**
 * ncm_c_mass_atomic:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Atomic mass constant $m_\mathrm{A} = 1.66053906660 \times 10^{-27} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_e:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Electron mass $m_\mathrm{e} = 9.1093837015 \times 10^{-31} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_p:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Proton mass $m_\mathrm{p} = 1.67262192369 \times 10^{-27} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_n:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: Neuton mass $m_\mathrm{n} = 1.67492749804 \times 10^{-27} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_ratio_alpha_p:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: The proton and alpha particle (Helium-4 III) mass ratio $m_\alpha / m_\mathrm{p} = 3.97259969009$.
 */

/**
 * ncm_c_mass_ratio_e_p:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: The electron and proton mass ratio $m_\mathrm{e} / m_\mathrm{p} = 5.44617021487 \times 10^{-4}$.
 */

/**
 * ncm_c_Rinf:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: The Rydberg constant $\mathrm{R}_\infty = 10973731.568160 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_Ry:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: The Rydberg unity of energy $\mathrm{Ry} = hc\mathrm{R}_\infty = 2.1798723611035 \times 10^{-18} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_eV:
 *
 * Using CODATA values, see [description][NcmC.description].
 *
 * Returns: The value of one electron volt $\mathrm{eV} = 1.602176634 \times 10^{-19} \,\left[\mathrm{J}\right]$.
 */

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

/**
 * ncm_c_year:
 *
 * One year ($365.25$ days) in seconds.
 *
 * Returns: $1$ year $365.25 \times 24 \times 60 \times 60 \,\left[\mathrm{s}\right]$.
 */

/**
 * ncm_c_lightyear:
 *
 * One year times the speed of light ncm_c_c() in meters.
 *
 * Returns: $1$ light-year $365.25 \times 24 \times 60 \times 60 \times c \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_lightyear_pc:
 *
 * One light-year in parsecs.
 *
 * Returns: $1$ light-year $365.25 \times 24 \times 60 \times 60 \times c \,\left[\mathrm{pc}\right]$.
 */

/**
 * ncm_c_Glightyear_Mpc:
 *
 * One giga light-year in mega parsecs.
 *
 * Returns: $1$ giga light-year $10^6 \times 365.25 \times 24 \times 60 \times 60 \times c \,\left[\mathrm{Mpc}\right]$.
 */

/**
 * ncm_c_hc:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Planck constant times the speed of light $hc \,\left[\mathrm{kg}\,\mathrm{m}^3\,\mathrm{s}^{-2}\right]$.
 */

/**
 * ncm_c_fine_struct_square:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: The square of the fine struct constant $\alpha^2$.
 */

/**
 * ncm_c_AR:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Radiation constant AR.
 */

/**
 * ncm_c_c2:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Square of the speed of light $c^2 \,\left[\mathrm{m}^2\,\mathrm{s}^{-2}\right]$.
 */

/**
 * ncm_c_planck_length2:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Square of the Planck length $l_\mathrm{P}^2 \,\left[\mathrm{m}^2\right]$.
 */

/**
 * ncm_c_rest_energy_atomic:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Rest energy of one atomic mass $m_\mathrm{A}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_e:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Electron's rest energy $m_\mathrm{e}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_p:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Proton's rest energy $m_\mathrm{p}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_n:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Returns: Neutron's rest energy $m_\mathrm{n}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_H_reduced_mass:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * Reduced mass for the electron in Hydrogen binding energy calculation, i.e.,
 * $m_\mathrm{r} = m_\mathrm{e} / (1 + m_\mathrm{e}/m_\mathrm{p})$
 *
 * Returns: Electron reduced mass $m_\mathrm{r} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_thermal_wl_e:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * The electron termal wavelength is $\lambda_\mathrm{e} = \sqrt{2\pi\hbar^2/(m_\mathrm{e}k_\mathrm{B}T)} \,\left[\mathrm{m}\right]$.
 *
 * Returns: Thermal electron wavelength times the temperature $\lambda_\mathrm{e}\sqrt{T}$.
 */

/**
 * ncm_c_thermal_wl_p:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * The proton termal wavelength is $\lambda_\mathrm{p} = \sqrt{2\pi\hbar^2/(m_\mathrm{p}k_\mathrm{B}T)} \,\left[\mathrm{m}\right]$.
 *
 * Returns: Thermal electron wavelength times the temperature $\lambda_\mathrm{p}\sqrt{T}$.
 */

/**
 * ncm_c_thermal_wl_n:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * The neutron termal wavelength is $\lambda_\mathrm{n} = \sqrt{2\pi\hbar^2/(m_\mathrm{n}k_\mathrm{B}T)} \,\left[\mathrm{m}\right]$.
 *
 * Returns: Thermal electron wavelength times the temperature $\lambda_\mathrm{n}\sqrt{T}$.
 */

/**
 * ncm_c_thermal_wn_e:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * The electron termal wavenumber is $k_\mathrm{e} = 1/\lambda_\mathrm{e}$,
 * see ncm_c_thermal_wl_e().
 *
 * Returns: Thermal eletron wavenumber $k_\mathrm{e}/\sqrt{T}$.
 */

/**
 * ncm_c_thermal_wn_p:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * The proton termal wavenumber is $k_\mathrm{p} = 1/\lambda_\mathrm{p}$,
 * see ncm_c_thermal_wl_p().
 *
 * Returns: Thermal proton wavenumber $k_\mathrm{e}/\sqrt{T}$.
 */

/**
 * ncm_c_thermal_wn_n:
 *
 * Derived from CODATA values, see [description][NcmC.description].
 *
 * The neutron termal wavenumber is $k_\mathrm{n} = 1/\lambda_\mathrm{n}$,
 * see ncm_c_thermal_wl_n().
 *
 * Returns: Thermal neutron wavenumber $k_\mathrm{e}/\sqrt{T}$.
 */

/**
 * ncm_c_H_reduced_energy:
 *
 * Reduced mass times $c^2$, $m_\mathrm{r}c^2$, see ncm_c_H_reduced_mass().
 *
 * Returns: $m_\mathrm{r}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_H_bind:
 * @n: Principal quantum number
 * @j: Total angular momentum
 *
 * Energy difference from unbounded state to state $(n,\,j)$, i.e., minus the
 * binding energy of the state $(n,\,j)$, calculated from
 * \begin{equation}
 * E^\mathrm{H}_{n,j} = m_\mathrm{e}c^2\left[1 - f(n,j)\right],
 * \end{equation}
 * where
 * \begin{align}
 * f(n, j)   &= \left[1+\left(\frac{\alpha}{n - \delta(j)}\right)^2\right]^{-\frac{1}{2}}, \\\\
 * \delta(j) &= j+\frac{1}{2} + \sqrt{\left(j+1/2\right)^2 - \alpha^2}.
 * \end{align}
 *
 * Returns: Hydrogen binding energy $E^\mathrm{H}_{n,j}$.
 */

/*******************************************************************************
 * END: 2018 CODATA recommended values
 *******************************************************************************/

/*******************************************************************************
 * START: IUPAC related constants
 *******************************************************************************/

/**
 * ncm_c_mass_1H_u:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 *
 * Returns: Hydrogen-1's mass over one atomic mass $m_\mathrm{1H}/m_\mathrm{A} = 1.00782503223$.
 */

/**
 * ncm_c_mass_2H_u:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 *
 * Returns: Hydrogen-2's mass over one atomic mass $m_\mathrm{2H}/m_\mathrm{A} = 2.01410177812$.
 */

/**
 * ncm_c_mass_3H_u:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 *
 * Returns: Hydrogen-3's mass over one atomic mass $m_\mathrm{3H}/m_\mathrm{A} = 3.0160492779$.
 */

/**
 * ncm_c_mass_3He_u:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 *
 * Returns: Helium-3's mass over one atomic mass $m_\mathrm{3He}/m_\mathrm{A} = 3.0160293201$.
 */

/**
 * ncm_c_mass_4He_u:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 *
 * Returns: Helium-4's mass over one atomic mass $m_\mathrm{4He}/m_\mathrm{A} = 4.00260325413$.
 */

/**
 * ncm_c_mass_1H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_1H_u() $\times$ ncm_c_mass_atomic().
 *
 * Returns: Hydrogen-1's mass $m_\mathrm{1H} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_2H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_2H_u() $\times$ ncm_c_mass_atomic().
 *
 * Returns: Hydrogen-2's mass $m_\mathrm{2H} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_3H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_3H_u() $\times$ ncm_c_mass_atomic().
 *
 * Returns: Hydrogen-3's mass $m_\mathrm{3H} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_3He:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_3He_u() $\times$ ncm_c_mass_atomic().
 *
 * Returns: Helium-3's mass $m_\mathrm{3He} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_mass_4He:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_4He_u() $\times$ ncm_c_mass_atomic().
 *
 * Returns: Helium-4's mass $m_\mathrm{4He} \,\left[\mathrm{kg}\right]$.
 */

/**
 * ncm_c_rest_energy_1H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_1H_u() $\times$ ncm_c_rest_energy_atomic().
 *
 * Returns: Hydrogen-1's rest energy $m_\mathrm{1H}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_2H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_2H_u() $\times$ ncm_c_rest_energy_atomic().
 *
 * Returns: Hydrogen-2's rest energy $m_\mathrm{2H}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_3H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_3H_u() $\times$ ncm_c_rest_energy_atomic().
 *
 * Returns: Hydrogen-3's rest energy $m_\mathrm{3H}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_3He:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_3He_u() $\times$ ncm_c_rest_energy_atomic().
 *
 * Returns: Helium-3's rest energy $m_\mathrm{3He}c^2 \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_rest_energy_4He:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_4He_u() $\times$ ncm_c_rest_energy_atomic().
 *
 * Returns: Helium-4's rest energy $m_\mathrm{4He}c^2 \,\left[\mathrm{J}\right]$.
 */

/*******************************************************************************
 * Derived constants
 *******************************************************************************/

/**
 * ncm_c_mass_ratio_4He_1H:
 *
 * Obtained from CIAAW commission of IUPAC, see [description][NcmC.description].
 * Calculated using ncm_c_mass_4He_u() / ncm_c_mass_1H_u().
 *
 * Returns: Helium-4 / Hydrogen-1 mass ratio $m_\mathrm{4He} / m_\mathrm{1H}$.
 */

/*******************************************************************************
 * END: IUPAC related constants
 *******************************************************************************/

/*******************************************************************************
 * START: IAU related constants
 *******************************************************************************/

/**
 * ncm_c_au:
 *
 * Using IAU 2015 recommendation see [description][NcmC.description],
 * compatible with [NASA JPL](http://ssd.jpl.nasa.gov/?constants) recommendations
 * (as in 5 January 2016).
 *
 * Returns: One astronomical unit in meters $\mathrm{au} = 1.49597870700 \times 10^{11} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_pc:
 *
 * Using IAU 2015 recommendation see [description][NcmC.description].
 *
 * Returns: One parsec in meters $\mathrm{pc} = 648000 \mathrm{au} / \pi = 3.0856775814913672789139379577965 \times 10^{16} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_kpc:
 *
 * Using IAU 2015 recommendation see [description][NcmC.description].
 *
 * Returns: One kilo parsec $\mathrm{kpc} = 10^3 \mathrm{pc}$.
 */

/**
 * ncm_c_Mpc:
 *
 * Using IAU 2015 recommendation see [description][NcmC.description].
 *
 * Returns: One mega parsec $\mathrm{Mpc} = 10^6 \mathrm{pc}$.
 */

/**
 * ncm_c_G_mass_solar:
 *
 * Using IAU 2015 recommendation see [description][NcmC.description].
 *
 * IAU recomends the use of a fixed value for the gravitational constant
 * times the solar mass.
 *
 * Returns: One solar mass times the gravitational constant $(\mathcal{GM})_\odot = 1.3271244 \times 10^{20} \,\left[\mathrm{m}^3\,\mathrm{s}^{-2}\right]$.
 */

/**
 * ncm_c_mass_solar:
 *
 * Using IAU 2015 recommendation see [description][NcmC.description].
 *
 * As in the recomendation above $\mathrm{M}_\odot = (\mathcal{GM})_\odot / \mathrm{G}$.
 * Here we use the CODATA 2018 value for $G$, see ncm_c_G().
 *
 * Returns: One solar mass $\mathrm{M}_\odot = (\mathcal{GM})_\odot / \mathrm{G} \,\left[\mathrm{kg}\right]$.
 */

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

/**
 * ncm_c_HI_ion_wn_1s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for H-I $1s\,{}^2\\!S_{1/2}$ state, i.e., $k_{1s\,{}^2\\!S_{1/2}}$.
 *
 * Returns: Hydrogen $1s\,{}^2\\!S_{1/2}$ ionization energy wavelength, $k_{1s\,{}^2\\!S_{1/2}} = 1.0967877174307 \times 10^{7} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_ion_wn_2s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for H-I $2s\,{}^2\\!S_{1/2}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2s\,{}^2\\!S_{1/2}} = k_{1s\,{}^2\\!S_{1/2}} - k_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly}$,
 * see ncm_c_HI_Lyman_wn_2s_2S0_5().
 *
 * Returns: Hydrogen $2s\,{}^2\\!S_{1/2}$ ionization energy wavelength, $k_{2s\,{}^2\\!S_{1/2}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_ion_wn_2p_2P0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for H-I $2p\,{}^2\\!P_{1/2}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^2\\!P_{1/2}} = k_{1s\,{}^2\\!S_{1/2}} - k_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly}$,
 * see ncm_c_HI_Lyman_wn_2p_2P0_5().
 *
 * Returns: Hydrogen $2p\,{}^2\\!P_{1/2}$ ionization energy wavelength, $k_{2p\,{}^2\\!P_{1/2}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_ion_wn_2p_2P3_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for H-I $2p\,{}^2\\!P_{3/2}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^2\\!P_{3/2}} = k_{1s\,{}^2\\!S_{3/2}} - k_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly}$,
 * see ncm_c_HI_Lyman_wn_2p_2P3_5().
 *
 * Returns: Hydrogen $2p\,{}^2\\!P_{3/2}$ ionization energy wavelength, $k_{2p\,{}^2\\!P_{3/2}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_ion_wn_2p_2Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * The mean ionization energy wavenumber for H-I $2p\,{}^2\\!P_{1/2}$ and
 * $2p\,{}^2\\!P_{3/2}$ states , i.e., $k_{2p\,{}^2\\!P_\mathrm{mean}} = (k_{2p\,{}^2\\!P_{1/2}} + k_{2p\,{}^2\\!P_{3/2}}) / 2$,
 * see ncm_c_HI_Lyman_wn_2p_2Pmean().
 *
 * Returns: Hydrogen states $2p\,{}^2\\!P_{1/2}$ and $2p\,{}^2\\!P_{3/2}$ mean ionization energy wavelength, $k_{2p\,{}^2\\!P_\mathrm{mean}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/* Ionization energy: E */

/**
 * ncm_c_HI_ion_E_1s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{1s\,{}^2\\!S_{1/2}}$,
 * see ncm_c_HI_ion_wn_1s_2S0_5().
 *
 * Returns: Hydrogen $1s\,{}^2\\!S_{1/2}$ ionization energy, $E_{1s\,{}^2\\!S_{1/2}} = hc\times{}k_{1s\,{}^2\\!S_{1/2}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HI_ion_E_2s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2s\,{}^2\\!S_{1/2}}$,
 * see ncm_c_HI_ion_wn_2s_2S0_5().
 *
 * Returns: Hydrogen $2s\,{}^2\\!S_{1/2}$ ionization energy, $E_{2s\,{}^2\\!S_{1/2}} = hc\times{}k_{2s\,{}^2\\!S_{1/2}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HI_ion_E_2p_2P0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^2\\!P_{1/2}}$,
 * see ncm_c_HI_ion_wn_2p_2P0_5().
 *
 * Returns: Hydrogen $2p\,{}^2\\!P_{1/2}$ ionization energy, $E_{2p\,{}^2\\!P_{1/2}} = hc\times{}k_{2p\,{}^2\\!P_{1/2}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HI_ion_E_2p_2P3_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^2\\!P_{3/2}}$,
 * see ncm_c_HI_ion_wn_2p_2P3_5().
 *
 * Returns: Hydrogen $2p\,{}^2\\!P_{3/2}$ ionization energy, $E_{2p\,{}^2\\!P_{3/2}} = hc\times{}k_{2p\,{}^2\\!P_{3/2}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HI_ion_E_2p_2Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^2\\!P_\mathrm{mean}}$,
 * see ncm_c_HI_ion_wn_2p_2Pmean().
 *
 * Returns: Hydrogen $2p\,{}^2\\!P_\mathrm{mean}$ ionization energy, $E_{2p\,{}^2\\!P_\mathrm{mean}} = hc\times{}k_{2p\,{}^2\\!P_\mathrm{mean}} \,\left[\mathrm{J}\right]$.
 */

/* Lyman series wavenumber: wn */

/**
 * ncm_c_HI_Lyman_wn_2s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2s\,{}^2\\!S_{1/2} \to 1s\,{}^2\\!S_{1/2}$ transition $k_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly} = 8.22589543992821 \times 10^6 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_Lyman_wn_2p_2P0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2p\,{}^2\\!P_{1/2} \to 1s\,{}^2\\!S_{1/2}$ transition $k_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly} = 8.22589191133 \times 10^6 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_Lyman_wn_2p_2P3_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2p\,{}^2\\!P_{3/2} \to 1s\,{}^2\\!S_{1/2}$ transition $k_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly} = 8.22592850014 \times 10^6 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HI_Lyman_wn_2p_2Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Mean Lyman emission wavenumber for the $2p\,{}^2\\!P_{1/2}$ and $2p\,{}^2\\!P_{3/2}$
 * states, $k_{2p\,{}^2\\!P_{mean}^\mathrm{Ly}} = (k_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly} + k_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly}) / 2$.
 *
 * Returns: $k_{2p\,{}^2\\!P_{mean}}^\mathrm{Ly} \,\left[\mathrm{m}^{-1}\right]$.
 */

/* Lyman series wavelength: wl */

/**
 * ncm_c_HI_Lyman_wl_2s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2s\,{}^2\\!S_{1/2} \to 1s\,{}^2\\!S_{1/2}$ transition
 * $\lambda_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly} = \left(k_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HI_Lyman_wn_2s_2S0_5().
 *
 * Returns: Wavelength for the $2s\,{}^2\\!S_{1/2} \to 1s\,{}^2\\!S_{1/2}$ transition, $\lambda_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HI_Lyman_wl_2p_2P0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^2\\!P_{1/2} \to 1s\,{}^2\\!S_{1/2}$ transition
 * $\lambda_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly} = \left(k_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HI_Lyman_wn_2p_2P0_5().
 *
 * Returns: Wavelength for the $2p\,{}^2\\!P_{1/2} \to 1s\,{}^2\\!S_{1/2}$ transition, $\lambda_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HI_Lyman_wl_2p_2P3_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^2\\!P_{3/2} \to 1s\,{}^2\\!S_{1/2}$ transition
 * $\lambda_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly} = \left(k_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HI_Lyman_wn_2p_2P3_5().
 *
 * Returns: Wavelength for the $2p\,{}^2\\!P_{3/2} \to 1s\,{}^2\\!S_{1/2}$ transition, $\lambda_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HI_Lyman_wl_2p_2Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^2\\!P_\mathrm{mean} \to 1s\,{}^2\\!S_{1/2}$ transition
 * $\lambda_{2p\,{}^2\\!P_\mathrm{mean}}^\mathrm{Ly} = \left(k_{2p\,{}^2\\!P_\mathrm{mean}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HI_Lyman_wn_2p_2Pmean().
 *
 * Returns: Wavelength for the $2p\,{}^2\\!P_\mathrm{mean} \to 1s\,{}^2\\!S_{1/2}$ transition, $\lambda_{2p\,{}^2\\!P_\mathrm{mean}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/* Lyman series factor: wl^3 / (8pi) */

/**
 * ncm_c_HI_Lyman_wl3_8pi_2s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2s\,{}^2\\!S_{1/2}} = \left(\lambda_{2s\,{}^2\\!S_{1/2}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HI_Lyman_wl_2s_2S0_5().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2s\,{}^2\\!S_{1/2}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HI_Lyman_wl3_8pi_2p_2P0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^2\\!P_{1/2}} = \left(\lambda_{2p\,{}^2\\!P_{1/2}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HI_Lyman_wl_2p_2P0_5().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^2\\!P_{1/2}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HI_Lyman_wl3_8pi_2p_2P3_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^2\\!P_{3/2}} = \left(\lambda_{2p\,{}^2\\!P_{3/2}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HI_Lyman_wl_2p_2P3_5().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^2\\!P_{3/2}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HI_Lyman_wl3_8pi_2p_2Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^2\\!P_\mathrm{mean}} = \left(\lambda_{2p\,{}^2\\!P_\mathrm{mean}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HI_Lyman_wl_2p_2Pmean().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^2\\!P_\mathrm{mean}} \,\left[\mathrm{m}^3\right]$.
 */

/* Boltzmann factor */

/**
 * ncm_c_boltzmann_factor_HI_1s_2S0_5:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{1s\,{}^2\\!S_{1/2}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{1s\,{}^2\\!S_{1/2}} / (k_\mathrm{B}T)\right]$,
 * for the $1s\,{}^2\\!S_{1/2}$ hydrogen energy level, see
 * ncm_c_HI_ion_E_1s_2S0_5() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{1s\,{}^2\\!S_{1/2}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HI_2s_2S0_5:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2s\,{}^2\\!S_{1/2}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2s\,{}^2\\!S_{1/2}} / (k_\mathrm{B}T)\right]$,
 * for the $2s\,{}^2\\!S_{1/2}$ hydrogen energy level, see
 * ncm_c_HI_ion_E_2s_2S0_5() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2s\,{}^2\\!S_{1/2}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HI_2p_2P0_5:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^2\\!P_{1/2}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^2\\!P_{1/2}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^2\\!P_{1/2}$ hydrogen energy level, see
 * ncm_c_HI_ion_E_2p_2P0_5() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^2\\!P_{1/2}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HI_2p_2P3_5:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^2\\!P_{3/2}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^2\\!P_{3/2}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^2\\!P_{3/2}$ hydrogen energy level, see
 * ncm_c_HI_ion_E_2p_2P3_5() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^2\\!P_{3/2}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HI_2p_2Pmean:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^2\\!P_\mathrm{mean}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^2\\!P_\mathrm{mean}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^2\\!P_\mathrm{mean}$ hydrogen energy level, see
 * ncm_c_HI_ion_E_2p_2Pmean() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^2\\!P_\mathrm{mean}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/*******************************************************************************
 * -- END: Hydrogen I
 *******************************************************************************/

/*******************************************************************************
 * -- START: Helium I
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

/**
 * ncm_c_HeI_ion_wn_1s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $1s\,{}^1\\!S_{0}$ state, i.e., $k_{1s\,{}^1\\!S_{0}}$.
 *
 * Returns: Helium-I $1s\,{}^1\\!S_{0}$ ionization energy wavelength, $k_{1s\,{}^1\\!S_{0}} = 1.9831066637 \times 10^{7} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2s\,{}^1\\!S_{0}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2s\,{}^1\\!S_{0}} = k_{1s\,{}^1\\!S_{0}} - k_{2s\,{}^1\\!S_{0}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2s_1S0().
 *
 * Returns: Helium-I $2s\,{}^1\\!S_{0}$ ionization energy wavelength, $k_{2s\,{}^1\\!S_{0}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2s\,{}^3\\!S_{1}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2s\,{}^3\\!S_{1}} = k_{1s\,{}^1\\!S_{0}} - k_{2s\,{}^3\\!S_{1}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2s_3S1().
 *
 * Returns: Helium-I $2s\,{}^3\\!S_{1}$ ionization energy wavelength, $k_{2s\,{}^3\\!S_{1}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2p_1P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2p\,{}^1\\!P_{1}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^1\\!P_{1}} = k_{1s\,{}^1\\!S_{0}} - k_{2p\,{}^1\\!P_{1}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2p_1P1().
 *
 * Returns: Helium-I $2p\,{}^1\\!P_{1}$ ionization energy wavelength, $k_{2p\,{}^1\\!P_{1}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2p_3P0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2p\,{}^3\\!P_{0}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^3\\!P_{0}} = k_{1s\,{}^1\\!S_{0}} - k_{2p\,{}^3\\!P_{0}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2p_3P0().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_{0}$ ionization energy wavelength, $k_{2p\,{}^3\\!P_{0}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2p_3P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2p\,{}^3\\!P_{1}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^3\\!P_{1}} = k_{1s\,{}^1\\!S_{0}} - k_{2p\,{}^3\\!P_{1}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2p_3P1().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_{1}$ ionization energy wavelength, $k_{2p\,{}^3\\!P_{1}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2p_3P2:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2p\,{}^3\\!P_{2}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^3\\!P_{2}} = k_{1s\,{}^1\\!S_{0}} - k_{2p\,{}^3\\!P_{2}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2p_3P2().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_{2}$ ionization energy wavelength, $k_{2p\,{}^3\\!P_{2}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_ion_wn_2p_3Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-I $2p\,{}^3\\!P_\mathrm{mean}$ state calculated
 * from the difference between the first state and the corresponding Lyman
 * wavenumber, i.e., $k_{2p\,{}^3\\!P_{0}} = k_{1s\,{}^1\\!S_{0}} - k_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly}$,
 * see ncm_c_HeI_Lyman_wn_2p_3Pmean().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_\mathrm{mean}$ ionization energy wavelength, $k_{2p\,{}^3\\!P_\mathrm{mean}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/* Ionization energy: E */

/**
 * ncm_c_HeI_ion_E_1s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{1s\,{}^1\\!S_{0}}$,
 * see ncm_c_HeI_ion_wn_1s_1S0().
 *
 * Returns: Helium-I $1s\,{}^1\\!S_{0}$ ionization energy, $E_{1s\,{}^1\\!S_{0}} = hc\times{}k_{1s\,{}^1\\!S_{0}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2s\,{}^1\\!S_{0}}$,
 * see ncm_c_HeI_ion_wn_2s_1S0().
 *
 * Returns: Helium-I $2s\,{}^1\\!S_{0}$ ionization energy, $E_{2s\,{}^1\\!S_{0}} = hc\times{}k_{2s\,{}^1\\!S_{0}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2s\,{}^3\\!S_{1}}$,
 * see ncm_c_HeI_ion_wn_2s_3S1().
 *
 * Returns: Helium-I $2s\,{}^3\\!S_{1}$ ionization energy, $E_{2s\,{}^3\\!S_{1}} = hc\times{}k_{2s\,{}^3\\!S_{1}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2p_1P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^1\\!P_{1}}$,
 * see ncm_c_HeI_ion_wn_2p_1P1().
 *
 * Returns: Helium-I $2p\,{}^1\\!P_{1}$ ionization energy, $E_{2p\,{}^1\\!P_{1}} = hc\times{}k_{2p\,{}^1\\!P_{1}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2p_3P0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^3\\!P_{0}}$,
 * see ncm_c_HeI_ion_wn_2p_3P0().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_{0}$ ionization energy, $E_{2p\,{}^3\\!P_{0}} = hc\times{}k_{2p\,{}^3\\!P_{0}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2p_3P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^3\\!P_{1}}$,
 * see ncm_c_HeI_ion_wn_2p_3P1().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_{1}$ ionization energy, $E_{2p\,{}^3\\!P_{1}} = hc\times{}k_{2p\,{}^3\\!P_{1}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2p_3P2:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^3\\!P_{2}}$,
 * see ncm_c_HeI_ion_wn_2p_3P2().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_{2}$ ionization energy, $E_{2p\,{}^3\\!P_{2}} = hc\times{}k_{2p\,{}^3\\!P_{2}} \,\left[\mathrm{J}\right]$.
 */

/**
 * ncm_c_HeI_ion_E_2p_3Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy calculated from the wavenumber $k_{2p\,{}^3\\!P_\mathrm{mean}}$,
 * see ncm_c_HeI_ion_wn_2p_3Pmean().
 *
 * Returns: Helium-I $2p\,{}^3\\!P_\mathrm{mean}$ ionization energy, $E_{2p\,{}^3\\!P_\mathrm{mean}} = hc\times{}k_{2p\,{}^3\\!P_\mathrm{m}} \,\left[\mathrm{J}\right]$.
 */

/* Lyman series wavenumber: wn */

/**
 * ncm_c_HeI_Lyman_wn_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2s\,{}^1\\!S_{0} \to 1s\,{}^1\\!S_{0}$ transition $k_{2s\,{}^1\\!S_{0}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2s\,{}^1\\!S_{0}}^\mathrm{Ly} = 1.66277440141 \times 10^7 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wn_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2s\,{}^3\\!S_{1} \to 1s\,{}^1\\!S_{0}$ transition $k_{2s\,{}^3\\!S_{1}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2s\,{}^3\\!S_{1}}^\mathrm{Ly} = 1.598559743297 \times 10^7 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wn_2p_1P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2p\,{}^1\\!P_{1} \to 1s\,{}^1\\!S_{0}$ transition $k_{2p\,{}^1\\!P_{1}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^1\\!P_{1}}^\mathrm{Ly} = 1.71134896946 \times 10^7 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wn_2p_3P0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2p\,{}^3\\!P_{0} \to 1s\,{}^1\\!S_{0}$ transition $k_{2p\,{}^3\\!P_{0}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^3\\!P_{0}}^\mathrm{Ly} = 1.690878308131 \times 10^7 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wn_2p_3P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2p\,{}^3\\!P_{1} \to 1s\,{}^1\\!S_{0}$ transition $k_{2p\,{}^3\\!P_{1}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^3\\!P_{1}}^\mathrm{Ly} = 1.690868428979 \times 10^7 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wn_2p_3P2:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Lyman emission wavenumber for the $2p\,{}^3\\!P_{2} \to 1s\,{}^1\\!S_{0}$ transition $k_{2p\,{}^3\\!P_{2}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^3\\!P_{2}}^\mathrm{Ly} = 1.690867664725 \times 10^7 \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wn_2p_3Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Mean Lyman emission wavenumber for the $2p\,{}^3\\!P_{*}$, i.e.,
 * $k_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly} = \left(k_{2p\,{}^3\\!P_{0}}^\mathrm{Ly} + k_{2p\,{}^3\\!P_{1}}^\mathrm{Ly} + k_{2p\,{}^3\\!P_{2}}^\mathrm{Ly}\right) / 3$.
 * See ncm_c_HeI_Lyman_wn_2p_3P0(), ncm_c_HeI_Lyman_wn_2p_3P1() and ncm_c_HeI_Lyman_wn_2p_3P2().
 *
 * Returns: $k_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly} \,\left[\mathrm{m}^{-1}\right]$.
 */

/* Lyman series wavelength: wl */

/**
 * ncm_c_HeI_Lyman_wl_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2s\,{}^1\\!S_{0} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2s\,{}^1\\!S_{0}}^\mathrm{Ly} = \left(k_{2s\,{}^1\\!S_{0}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2s_1S0().
 *
 * Returns: Wavelength for the $2s\,{}^1\\!S_{0} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2s\,{}^1\\!S_{0}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2s\,{}^3\\!S_{1} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2s\,{}^3\\!S_{1}}^\mathrm{Ly} = \left(k_{2s\,{}^3\\!S_{1}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2s_3S1().
 *
 * Returns: Wavelength for the $2s\,{}^3\\!S_{1} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2s\,{}^3\\!S_{1}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl_2p_1P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^1\\!P_{1} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2p\,{}^1\\!P_{1}}^\mathrm{Ly} = \left(k_{2p\,{}^1\\!P_{1}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2p_1P1().
 *
 * Returns: Wavelength for the $2p\,{}^1\\!P_{1} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2p\,{}^1\\!P_{1}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl_2p_3P0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^3\\!P_{0} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2p\,{}^3\\!P_{0}}^\mathrm{Ly} = \left(k_{2p\,{}^3\\!P_{0}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2p_3P0().
 *
 * Returns: Wavelength for the $2p\,{}^3\\!P_{0} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2p\,{}^3\\!P_{0}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl_2p_3P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^3\\!P_{1} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2p\,{}^3\\!P_{1}}^\mathrm{Ly} = \left(k_{2p\,{}^3\\!P_{1}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2p_3P1().
 *
 * Returns: Wavelength for the $2p\,{}^3\\!P_{1} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2p\,{}^3\\!P_{1}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl_2p_3P2:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^3\\!P_{2} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2p\,{}^3\\!P_{2}}^\mathrm{Ly} = \left(k_{2p\,{}^3\\!P_{2}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2p_3P2().
 *
 * Returns: Wavelength for the $2p\,{}^3\\!P_{2} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2p\,{}^3\\!P_{2}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl_2p_3Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Wavelength for the $2p\,{}^3\\!P_\mathrm{mean} \to 1s\,{}^1\\!S_{0}$ transition
 * $\lambda_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly} = \left(k_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly}\right)^{-1}$,
 * see ncm_c_HeI_Lyman_wn_2p_3Pmean().
 *
 * Returns: Wavelength for the $2p\,{}^3\\!P_\mathrm{mean} \to 1s\,{}^1\\!S_{0}$ transition, $\lambda_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly} \,\left[\mathrm{m}\right]$.
 */

/* Lyman series factor: wl^3 / (8pi) */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2s\,{}^1\\!S_{0}} = \left(\lambda_{2s\,{}^1\\!S_{0}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2s_1S0().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2s\,{}^1\\!S_{0}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2s\,{}^3\\!S_{1}} = \left(\lambda_{2s\,{}^3\\!S_{1}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2s_3S1().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2s\,{}^3\\!S_{1}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2p_1P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^1\\!P_{1}} = \left(\lambda_{2p\,{}^1\\!P_{1}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2p_1P1().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^1\\!P_{1}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2p_3P0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^3\\!P_{0}} = \left(\lambda_{2p\,{}^3\\!P_{0}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2p_3P0().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^3\\!P_{0}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2p_3P1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^3\\!P_{1}} = \left(\lambda_{2p\,{}^3\\!P_{1}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2p_3P1().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^3\\!P_{1}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2p_3P2:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^3\\!P_{2}} = \left(\lambda_{2p\,{}^3\\!P_{2}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2p_3P2().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^3\\!P_{2}} \,\left[\mathrm{m}^3\right]$.
 */

/**
 * ncm_c_HeI_Lyman_wl3_8pi_2p_3Pmean:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Effective volume of the Lyman wavelength $V^\mathrm{Ly}_{2p\,{}^3\\!P_\mathrm{mean}} = \left(\lambda_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly}\right)^{3} / (8\pi)$,
 * see ncm_c_HeI_Lyman_wl_2p_3Pmean().
 *
 * Returns: Effective volume $V^\mathrm{Ly}_{2p\,{}^3\\!P_\mathrm{mean}} \,\left[\mathrm{m}^3\right]$.
 */

/* Boltzmann factor */

/**
 * ncm_c_boltzmann_factor_HeI_1s_1S0:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{1s\,{}^1\\!S_{0}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{1s\,{}^1\\!S_{0}} / (k_\mathrm{B}T)\right]$,
 * for the $1s\,{}^1\\!S_{0}$ helium energy level, see
 * ncm_c_HeI_ion_E_1s_1S0() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{1s\,{}^1\\!S_{0}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2s_1S0:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2s\,{}^1\\!S_{0}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2s\,{}^1\\!S_{0}} / (k_\mathrm{B}T)\right]$,
 * for the $2s\,{}^1\\!S_{0}$ helium energy level, see
 * ncm_c_HeI_ion_E_2s_1S0() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2s\,{}^1\\!S_{0}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2s_3S1:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2s\,{}^3\\!S_{1}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2s\,{}^3\\!S_{1}} / (k_\mathrm{B}T)\right]$,
 * for the $2s\,{}^3\\!S_{1}$ helium energy level, see
 * ncm_c_HeI_ion_E_2s_3S1() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2s\,{}^3\\!S_{1}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2p_1P1:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^1\\!P_{1}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^1\\!P_{1}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^1\\!P_{1}$ helium energy level, see
 * ncm_c_HeI_ion_E_2p_1P1() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^1\\!P_{1}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2p_3P0:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^3\\!P_{0}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^3\\!P_{0}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^3\\!P_{0}$ helium energy level, see
 * ncm_c_HeI_ion_E_2p_3P0() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^3\\!P_{0}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2p_3P1:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^3\\!P_{1}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^3\\!P_{1}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^3\\!P_{1}$ helium energy level, see
 * ncm_c_HeI_ion_E_2p_3P1() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^3\\!P_{1}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2p_3P2:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^3\\!P_{2}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^3\\!P_{2}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^3\\!P_{2}$ helium energy level, see
 * ncm_c_HeI_ion_E_2p_3P2() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^3\\!P_{2}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/**
 * ncm_c_boltzmann_factor_HeI_2p_3Pmean:
 * @T: temperature $T$
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Calculates the Boltzmann factor $B_{2p\,{}^3\\!P_\mathrm{mean}}(T) = k_\mathrm{e}^3 T^{-3/2}\,\exp\left[-E_{2p\,{}^3\\!P_\mathrm{mean}} / (k_\mathrm{B}T)\right]$,
 * for the $2p\,{}^3\\!P_\mathrm{mean}$ helium energy level, see
 * ncm_c_HeI_ion_E_2p_3Pmean() and ncm_c_thermal_wn_e().
 *
 * Returns: Boltzmann factor $B_{2p\,{}^3\\!P_\mathrm{mean}}(T) \,\left[\mathrm{m}^3\,\mathrm{K}^{-3/2}\right]$.
 */

/* Balmer series wavenumber: wn */

/**
 * ncm_c_HeI_Balmer_wn_2p_1P1_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Balmer emission wavenumber for the $2p\,{}^1\\!P_{1} \to 2s\,{}^1\\!S_{0}$ transition $k_{2p\,{}^1\\!P_{1}}^{2s\,{}^1\\!S_{0}}$,
 * calculated from the difference between the Lyman lines $2s\,{}^1\\!S_{0}$ state and the
 * corresponding Lyman wavenumber, i.e.,
 * $k_{2p\,{}^1\\!P_{1}}^{2s\,{}^1\\!S_{0}} = k_{2p\,{}^1\\!P_{1}}^\mathrm{Ly} - k_{2s\,{}^1\\!S_{0}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^1\\!P_{1}}^{2s\,{}^1\\!S_{0}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/**
 * ncm_c_HeI_Balmer_wn_2p_3Pmean_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Balmer emission wavenumber for the $2p\,{}^3\\!P_\mathrm{mean} \to 2s\,{}^3\\!S_{1}$ transition $k_{2p\,{}^3\\!P_\mathrm{mean}}^{2s\,{}^3\\!S_{1}}$,
 * calculated from the difference between the Lyman lines $2s\,{}^1\\!S_{0}$ state and the
 * corresponding Lyman wavenumber, i.e.,
 * $k_{2p\,{}^3\\!P_\mathrm{mean}}^{2s\,{}^3\\!S_{1}} = k_{2p\,{}^3\\!P_\mathrm{mean}}^\mathrm{Ly} - k_{2s\,{}^3\\!S_{1}}^\mathrm{Ly}$.
 *
 * Returns: $k_{2p\,{}^3\\!P_\mathrm{mean}}^{2s\,{}^3\\!S_{1}} \,\left[\mathrm{m}^{-1}\right]$.
 */

/* Balmer series: E / k_B */

/**
 * ncm_c_HeI_Balmer_E_kb_2p_1P1_2s_1S0:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Balmer emission energy $E_{2p\,{}^1\\!P_{1}}^{2s\,{}^1\\!S_{0}} = hc\times{}k_{2p\,{}^1\\!P_{1}}^{2s\,{}^1\\!S_{0}}$
 * over $k_\mathrm{B}$.
 *
 * Returns: $E_{2p\,{}^1\\!P_{1}}^{2s\,{}^1\\!S_{0}} / k_\mathrm{B}$.
 */

/**
 * ncm_c_HeI_Balmer_E_kb_2p_3Pmean_2s_3S1:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Balmer emission energy $E_{2p\,{}^3\\!P_\mathrm{mean}}^{2s\,{}^3\\!S_{1}} = hc\times{}k_{2p\,{}^3\\!P_\mathrm{mean}}^{2s\,{}^3\\!S_{1}}$
 * over $k_\mathrm{B}$.
 *
 * Returns: $E_{2p\,{}^3\\!P_\mathrm{mean}}^{2s\,{}^3\\!S_{1}} / k_\mathrm{B}$.
 */

/*******************************************************************************
 * -- END: Helium I
 *******************************************************************************/

/*******************************************************************************
 * -- START: Helium II
 *******************************************************************************/
/* Ionization energy wavenumber: wn */

/**
 * ncm_c_HeII_ion_wn_1s_2S0_5:
 *
 * NIST compilation of atomic spectra see [description][NcmC.description].
 *
 * Ionization energy wavenumber for He-II $1s\,{}^2\\!S_{1/2}$ state, i.e., $k_{1s\,{}^2\\!S_{1/2}}$.
 *
 * Returns: Helium-II $1s\,{}^2\\!S_{1/2}$ ionization energy wavelength, $k_{1s\,{}^2\\!S_{1/2}} = 1.0967877174307 \times 10^{7} \,\left[\mathrm{m}^{-1}\right]$.
 */

/* Ionization energy: E */

/**
 * ncm_c_HeII_ion_E_1s_2S0_5:
 *
 * Ionization energy for He-II $1s\,{}^2\\!S_{1/2}$ state, i.e., $E_{1s\,{}^2\\!S_{1/2}} = hc \times k_{1s\,{}^2\\!S_{1/2}}$.
 *
 * Returns: Helium-II $1s\,{}^2\\!S_{1/2}$ ionization energy $E_{1s\,{}^2\\!S_{1/2}} \,\left[\mathrm{J}\right]$.
 */

/*******************************************************************************
 * -- END: Helium II
 *******************************************************************************/

/*******************************************************************************
 * END: NIST Atomic Spectra database
 *******************************************************************************/

/*******************************************************************************
 * Constants from other sources
 *******************************************************************************/

/**
 * ncm_c_decay_H_rate_2s_1s:
 *
 * Theoretical value for the two photons decay rate for Hydrogen
 * $2\mathrm{s} \to 1\mathrm{s}$ states [Goldman 1989][XGoldman1989].
 *
 * Returns: Decay rate of Hydrogen from $\Lambda_{2\mathrm{s} \to 1\mathrm{s}} = 8.2245809 \,\left[\mathrm{s}^{-1}\right]$.
 */

/**
 * ncm_c_decay_He_rate_2s_1s:
 *
 * Theoretical value for the two photons decay rate for Helium
 * $2\mathrm{s} \to 1\mathrm{s}$ states [Drake 1969][XDrake1969].
 *
 * Returns: Decay rate of Helium from $\Lambda_{2\mathrm{s} \to 1\mathrm{s}} = 51.3 \,\left[\mathrm{s}^{-1}\right]$.
 */

/*******************************************************************************
 * START: Statistics
 *******************************************************************************/

/**
 * ncm_c_stats_1sigma:
 *
 * The integral of a Gaussian distribution with mean $\mu$
 * and standard deviation $\sigma$ in $(\mu - 1 \sigma, \mu + 1 \sigma)$.
 *
 * Returns: $P (\mu - 1 \sigma, \mu + 1 \sigma)$
 *
 */

/**
 * ncm_c_stats_2sigma:
 *
 * The integral of a Gaussian distribution with mean $\mu$
 * and standard deviation $\sigma$ in $(\mu - 2 \sigma, \mu + 2 \sigma)$.
 *
 * Returns: $P (\mu - 2 \sigma, \mu + 2 \sigma)$
 *
 */

/**
 * ncm_c_stats_3sigma:
 *
 * The integral of a Gaussian distribution with mean $\mu$
 * and standard deviation sigma in $(\mu - 3 \sigma, \mu + 3 \sigma)$.
 *
 * Returns: $P (\mu - 3 \sigma, \mu + 3 \sigma)$
 *
 */

/*******************************************************************************
 * END: Statistics
 *******************************************************************************/

/*******************************************************************************
 * START: Observational data
 *******************************************************************************/

/**
 * ncm_c_wmap5_coadded_I_K:
 *
 * WMAP-5 noise per observation for resolution 9 and 10 $I$ maps, to $\sim$0.1\% uncertainty for differencing assembly (DA) $\text{K1}$. See [Hinshaw (2009)][X2009ApJS..180..225H] [[arXiv](https://arxiv.org/abs/0803.0732)].
 *
 * Returns: $\sigma_0(I) = 1.436 \left[ \text{mK} \right]$.
 *
 */

/**
 * ncm_c_wmap5_coadded_I_Ka:
 *
 * WMAP-5 noise per observation for resolution 9 and 10 $I$ maps, to $\sim$0.1\% uncertainty for differencing assembly (DA) $\text{Ka1}$. See [Hinshaw (2009)][X2009ApJS..180..225H] [[arXiv](https://arxiv.org/abs/0803.0732)].
 *
 * Returns: $\sigma_0(I) = 1.470 \left[ \text{mK} \right]$.
 *
 */

/**
 * ncm_c_wmap5_coadded_I_Q:
 *
 * WMAP-5 noise per observation for resolution 9 and 10 $I$ maps, to $\sim$0.1\% uncertainty for differencing assembly (DA) $\text{Q}$. See [Smidt (2009)][X2009PhRvD..80l3005S] [[arXiv](https://arxiv.org/abs/0907.4051)].
 *
 * Returns: $\sigma_0(I) = 2.197 \left[ \text{mK} \right]$.
 *
 */

/**
 * ncm_c_wmap5_coadded_I_V:
 *
 * WMAP-5 noise per observation for resolution 9 and 10 $I$ maps, to $\sim$0.1\% uncertainty for differencing assembly (DA) $\text{V}$. See [Smidt (2009)][X2009PhRvD..80l3005S] [[arXiv](https://arxiv.org/abs/0907.4051)].
 *
 * Returns: $\sigma_0(I) = 3.133 \left[ \text{mK} \right]$.
 *
 */

/**
 * ncm_c_wmap5_coadded_I_W:
 *
 * WMAP-5 noise per observation for resolution 9 and 10 $I$ maps, to $\sim$0.1\% uncertainty for differencing assembly (DA) $\text{W}$. See [Smidt (2009)][X2009PhRvD..80l3005S] [[arXiv](https://arxiv.org/abs/0907.4051)].
 *
 * Returns: $\sigma_0(I) = 6.538 \left[ \text{mK} \right]$.
 *
 */

/**
 * ncm_c_hubble_cte_planck6_base:
 *
 * Planck VI Hubble constant base-$\Lambda$CDM model TT,TE,EE$+$lowE$+$lensing final result. See [Planck Collaboration (2018)][X2018arXiv180706209P] [[arXiv](https://arxiv.org/abs/1807.06209)].
 *
 * Returns: $H_0 = 67.36 \left[\text{km}\,\text{s}^{-1}\,\text{Mpc}^{-1}\right]$.
 *
 */

/**
 * ncm_c_hubble_cte_hst:
 *
 * HST Hubble constant final result. See [Freedman (2001)][X2001ApJ...553...47F] [[arXiv](https://arxiv.org/abs/astro-ph/0012376)].
 *
 * Returns: $H_0 = 72 \left[\text{km}\,\text{s}^{-1}\,\text{Mpc}^{-1}\right]$.
 *
 */

/**
 * ncm_c_hubble_radius_hm1_Mpc:
 *
 * Hubble radius in units of $\mathsf{h}^{-1} \, \text{Mpc}$ defined as
 * \begin{equation}
 *   R_H  h^{-1} = \frac{c}{100 \mathsf{h} \, \text{km} \, \text{sec}^{-1} \, \text{Mpc}^{-1}} \, ,
 * \end{equation}
 * where $c$ is the speed of light (#ncm_c_c()). Calculated using ncm_c_c() $/$ $10^{5}$.
 *
 * Returns: Hubble radius $R_H \mathsf{h}^{-1} \left[\text{Mpc}\right]$.
 *
 */

/**
 * ncm_c_hubble_radius_hm1_planck:
 *
 * Hubble radius in units of $\mathsf{h}^{-1} \, l_{\text{p}}$. Calculated using ncm_c_hubble_radius_hm1_Mpc () $\times$ ncm_c_Mpc () $/$ ncm_c_planck_length ().
 *
 * Returns: Hubble radius $R_H \mathsf{h}^{-1} \left[l_\text{p}\right]$.
 *
 */

/**
 * ncm_c_crit_density_h2:
 *
 * The critical density is defined as
 * \begin{equation}
 * \rho_{\mathrm{crit}0} = \frac{3 c^2 H_0^2}{8\pi G},
 * \end{equation}
 * where $G$ is the gravitational constant (#ncm_c_G()), $c$ is the speed of light
 * (#ncm_c_c()) and $H_0$ is the Hubble parameter,
 * $$H_0 = 100 \times \mathsf{h} \,\left[\text{km}\,\text{s}^{-1}\,\text{Mpc}^{-1}\right].$$
 *
 * Returns: Critical density over $\mathsf{h}^2$, $$\frac{\rho_{\mathrm{crit}0}}{\mathsf{h}^2} \left[\frac{\text{kg}}{\text{m}^3} \frac{\text{m}^2}{\text{s}^2}\right].$$
 */

/**
 * ncm_c_crit_mass_density_h2:
 *
 * This function computes the critical mass density over $\mathsf{h}^2 \times c^2$.
 *
 * Returns: Critical mass density over $\mathsf{h}^2 \times c^2$, $$\frac{\rho_{\mathrm{crit}0}}{c^2\mathsf{h}^2} \,\left[\frac{\text{kg}}{\text{m}^3}\right].$$
 *
 */

/**
 * ncm_c_crit_mass_density_h2_solar_mass_Mpc3:
 *
 * This function computes the critical mass density in units of solar mass $M_\odot$ and Mpc.
 *
 * Returns: Critical mass density in $M_\odot$ and Mpc units $\frac{\rho_{\mathrm{crit}0}}{\mathsf{h}^2 M_\odot} \left(1 \mathrm{Mpc}\right)^3$.
 */

/**
 * ncm_c_crit_number_density_p:
 *
 * This function computes the proton number density in units of its rest energy. Calculated using ncm_c_crit_density_h2 () $/$ ncm_c_rest_energy_p ().
 *
 * Returns: Critical proton number density in units of its rest mass.
 *
 */

/**
 * ncm_c_crit_number_density_n:
 *
 * This function computes the neutron number density in units of its rest energy. Calculated using ncm_c_crit_density_h2 () $/$ ncm_c_rest_energy_n ().
 *
 * Returns: Critical neutron number density in units of its rest mass.
 *
 */

/**
 * ncm_c_blackbody_energy_density:
 *
 * This functions returns the black body energy density divided by $T^4$. Defined as
 * \begin{equation}
 *   \frac{\rho_{\text{BL}}}{T^4} = \frac{8\pi^2k_{\text{b}}^4}{15 \left( hc \right)^3},
 * \end{equation}
 * where $\rho_{\text{BL}}$ is the black body energy density, $T$ is the temperature, $k_{\text{b}}$ is the Boltzmann constant (#ncm_c_kb()), $h$ is the Planck constant (#ncm_c_h()) and $c$ is the speed of light (#ncm_c_c()).
 *
 * Returns: Blackbody energy density in $\left[ \text{J} \, \text{m}^{-3} \, \text{K}^{-4}  \right]$units.
 *
 */

/**
 * ncm_c_blackbody_per_crit_density_h2:
 *
 * This functions returns ncm_c_blackbody_energy_density () $/$ ncm_c_crit_density_h2 ().
 *
 * Returns: Blackbody energy density over critical density times $\mathsf{h}^2$.
 */

/**
 * ncm_c_radiation_temp_to_h2Omega_r0:
 * @T: temperature $T$
 *
 * Returns: ncm_c_blackbody_per_crit_density_h2 () $\times$ $T^4$.
 *
 */

/*******************************************************************************
 * END: Observational data
 *******************************************************************************/

