#!/usr/bin/env python
#
# test_py_constants.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_constants.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests for the Python bindings for constants."""

import math
from numpy.testing import assert_allclose
import numpy as np
from scipy import constants
from scipy.constants import physical_constants
from astropy import units as u

from numcosmo_py import Ncm

Ncm.cfg_init()


def test_constants_object():
    """Test constants from ncm_c.h"""
    c = Ncm.C()
    assert c is not None


def test_constants_math():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.sqrt_1_4pi(), math.sqrt(1.0 / (4.0 * math.pi)))
    assert_allclose(Ncm.C.sqrt_pi(), math.sqrt(math.pi))
    assert_allclose(Ncm.C.sqrt_2pi(), math.sqrt(2.0 * math.pi))
    assert_allclose(Ncm.C.sqrt_pi_2(), math.sqrt(math.pi / 2.0))
    assert_allclose(Ncm.C.sqrt_3_4pi(), math.sqrt(3.0 / (4.0 * math.pi)))
    assert_allclose(Ncm.C.ln2(), math.log(2.0))
    assert_allclose(Ncm.C.ln3(), math.log(3.0))
    assert_allclose(Ncm.C.lnpi_4(), math.log(math.pi) / 4.0)
    assert_allclose(Ncm.C.ln2pi(), math.log(2.0 * math.pi))
    assert_allclose(Ncm.C.lnpi(), math.log(math.pi))
    assert_allclose(Ncm.C.pi(), math.pi)
    assert_allclose(Ncm.C.two_pi_2(), 2.0 * math.pi**2)

    assert_allclose(Ncm.C.tan_1arcsec(), math.tan(math.pi / (3600.0 * 180.0)))
    assert_allclose(Ncm.C.deg2_steradian(), (math.pi / 180.0) ** 2)

    angles = np.linspace(-10 * math.pi, 10.0 * math.pi, 311, endpoint=False)
    angles_deg = angles * 180.0 / math.pi

    assert_allclose(
        [Ncm.C.radian_0_2pi(angle) for angle in angles], angles % (2.0 * math.pi)
    )
    assert_allclose(
        [Ncm.C.sign_sin(angle) for angle in angles], np.sign(np.sin(angles))
    )

    assert_allclose(
        [Ncm.C.degree_to_radian(angle) for angle in angles_deg],
        angles_deg * math.pi / 180.0,
    )
    assert_allclose(
        [Ncm.C.radian_to_degree(angle) for angle in angles],
        angles * 180.0 / math.pi,
    )


def test_constants_scipy():
    """Test constants from scipy.constants"""

    assert_allclose(Ncm.C.c(), constants.c)
    assert_allclose(Ncm.C.h(), constants.h)
    assert_allclose(Ncm.C.hbar(), constants.hbar)
    assert_allclose(Ncm.C.fine_struct(), constants.alpha)
    assert_allclose(Ncm.C.kb(), constants.k)
    assert_allclose(Ncm.C.G(), constants.G)
    assert_allclose(Ncm.C.hbar(), constants.hbar)
    assert_allclose(Ncm.C.thomson_cs(), physical_constants["Thomson cross section"][0])
    assert_allclose(Ncm.C.magnetic_constant(), constants.mu_0)
    assert_allclose(Ncm.C.mass_p(), constants.m_p)
    assert_allclose(Ncm.C.mass_e(), constants.m_e)
    assert_allclose(Ncm.C.mass_n(), constants.m_n)

    assert_allclose(Ncm.C.mass_ratio_alpha_p(), 3.97259969009)
    assert_allclose(Ncm.C.Rinf(), physical_constants["Rydberg constant"][0])
    assert_allclose(Ncm.C.Ry(), physical_constants["Rydberg constant times hc in J"][0])
    assert_allclose(Ncm.C.eV(), physical_constants["electron volt"][0])


# NCM_INLINE gdouble ncm_c_H_bind (const gdouble n, const gdouble j) G_GNUC_CONST;


def test_constants_derived():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.year(), u.year.to(u.s))
    assert_allclose(Ncm.C.lightyear(), u.lyr.to(u.m))
    assert_allclose(Ncm.C.lightyear_pc(), u.lyr.to(u.pc))
    assert_allclose(Ncm.C.Glightyear_Mpc(), 1.0e9 * u.lyr.to(u.Mpc))
    assert_allclose(Ncm.C.hc(), constants.h * constants.c)
    assert_allclose(Ncm.C.fine_struct_square(), constants.alpha**2)
    assert_allclose(Ncm.C.electric_constant(), constants.epsilon_0)
    assert_allclose(
        Ncm.C.AR(),
        4.0 * physical_constants["Stefan-Boltzmann constant"][0] / constants.c,
    )
    assert_allclose(Ncm.C.c2(), constants.c**2)
    assert_allclose(
        Ncm.C.planck_length2(), constants.hbar * constants.G / constants.c**3
    )
    assert_allclose(Ncm.C.rest_energy_atomic(), constants.atomic_mass * constants.c**2)
    assert_allclose(Ncm.C.rest_energy_e(), constants.m_e * constants.c**2)
    assert_allclose(Ncm.C.rest_energy_p(), constants.m_p * constants.c**2)
    assert_allclose(Ncm.C.rest_energy_n(), constants.m_n * constants.c**2)
    factor = math.sqrt(2.0 * math.pi * constants.hbar**2 / constants.k)
    assert_allclose(Ncm.C.thermal_wl_e(), factor / math.sqrt(constants.m_e))
    assert_allclose(Ncm.C.thermal_wl_p(), factor / math.sqrt(constants.m_p))
    assert_allclose(Ncm.C.thermal_wl_n(), factor / math.sqrt(constants.m_n))
    assert_allclose(Ncm.C.thermal_wn_e(), math.sqrt(constants.m_e) / factor)
    assert_allclose(Ncm.C.thermal_wn_p(), math.sqrt(constants.m_p) / factor)
    assert_allclose(Ncm.C.thermal_wn_n(), math.sqrt(constants.m_n) / factor)

    H_reduced_mass = constants.m_p * constants.m_e / (constants.m_p + constants.m_e)
    assert_allclose(Ncm.C.H_reduced_mass(), H_reduced_mass)
    assert_allclose(Ncm.C.H_reduced_energy(), H_reduced_mass * constants.c**2)
    assert_allclose(Ncm.C.H_bind(1.0, 0.5), 2.178715e-18)
    assert_allclose(Ncm.C.H_bind(2.0, 0.5), 5.446805e-19)


def test_constants_masses_atomic_unit():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.mass_1H_u(), 1.00782503)
    assert_allclose(Ncm.C.mass_2H_u(), 2.0141017778)
    assert_allclose(Ncm.C.mass_3H_u(), 3.0160492777)
    assert_allclose(Ncm.C.mass_3He_u(), 3.0160293191)
    assert_allclose(Ncm.C.mass_4He_u(), 4.00260325413)


def test_constants_masses():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.mass_1H(), 1.00782503 * constants.atomic_mass)
    assert_allclose(Ncm.C.mass_2H(), 2.0141017778 * constants.atomic_mass)
    assert_allclose(Ncm.C.mass_3H(), 3.0160492777 * constants.atomic_mass)
    assert_allclose(Ncm.C.mass_3He(), 3.0160293191 * constants.atomic_mass)
    assert_allclose(Ncm.C.mass_4He(), 4.00260325413 * constants.atomic_mass)


def test_constants_rest_energy():
    """Test constants from ncm_c.h"""

    assert_allclose(
        Ncm.C.rest_energy_1H(), 1.00782503 * constants.atomic_mass * constants.c**2
    )
    assert_allclose(
        Ncm.C.rest_energy_2H(), 2.0141017778 * constants.atomic_mass * constants.c**2
    )
    assert_allclose(
        Ncm.C.rest_energy_3H(), 3.0160492777 * constants.atomic_mass * constants.c**2
    )
    assert_allclose(
        Ncm.C.rest_energy_3He(), 3.0160293191 * constants.atomic_mass * constants.c**2
    )
    assert_allclose(
        Ncm.C.rest_energy_4He(),
        4.00260325413 * constants.atomic_mass * constants.c**2,
    )


def test_constants_mass_ratio():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.mass_ratio_4He_1H(), 4.00260325413 / 1.00782503)


def test_constants_distances():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.au(), constants.au)
    assert_allclose(Ncm.C.pc(), u.pc.to(u.m))
    assert_allclose(Ncm.C.kpc(), u.kpc.to(u.m))
    assert_allclose(Ncm.C.Mpc(), u.Mpc.to(u.m))
    assert_allclose(Ncm.C.G_mass_solar(), constants.G * u.M_sun.to(u.kg))
    assert_allclose(Ncm.C.mass_solar(), u.M_sun.to(u.kg))


def test_constants_HI_ion_wn():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HI_ion_wn_1s_2S0_5(), 10967877.174307)
    assert_allclose(Ncm.C.HI_ion_wn_2s_2S0_5(), 2741981.734379)
    assert_allclose(Ncm.C.HI_ion_wn_2p_2P0_5(), 2741985.262977)
    assert_allclose(Ncm.C.HI_ion_wn_2p_2P3_5(), 2741948.674167)
    assert_allclose(Ncm.C.HI_ion_wn_2p_2Pmean(), 2741966.968572)


def test_constants_HI_ion_E():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HI_ion_E_1s_2S0_5(), 2.1787094174620437e-18)
    assert_allclose(Ncm.C.HI_ion_E_2s_2S0_5(), 5.44679825663478e-19)
    assert_allclose(Ncm.C.HI_ion_E_2p_2P0_5(), 5.446805266004078e-19)
    assert_allclose(Ncm.C.HI_ion_E_2p_2P3_5(), 5.446732584314035e-19)
    assert_allclose(Ncm.C.HI_ion_E_2p_2Pmean(), 5.446768925159056e-19)


def test_constants_HI_Lyman_wn():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HI_Lyman_wn_2s_2S0_5(), 8.22589543992821e6)
    assert_allclose(Ncm.C.HI_Lyman_wn_2p_2P0_5(), 8.22589191133e6)
    assert_allclose(Ncm.C.HI_Lyman_wn_2p_2P3_5(), 8.22592850014e6)
    assert_allclose(Ncm.C.HI_Lyman_wn_2p_2Pmean(), 8.225910205735e6)


def test_constants_HI_Lyman_wl():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HI_Lyman_wl_2s_2S0_5(), 1.0 / 8.22589543992821e6)
    assert_allclose(Ncm.C.HI_Lyman_wl_2p_2P0_5(), 1.0 / 8.22589191133e6)
    assert_allclose(Ncm.C.HI_Lyman_wl_2p_2P3_5(), 1.0 / 8.22592850014e6)
    assert_allclose(Ncm.C.HI_Lyman_wl_2p_2Pmean(), 1.0 / 8.225910205735e6)


def test_constants_HI_Lyman_wl3_8pi():
    """Test constants from ncm_c.h"""

    assert_allclose(
        Ncm.C.HI_Lyman_wl3_8pi_2s_2S0_5(),
        (1.0 / 8.22589543992821e6) ** 3 / (8.0 * math.pi),
    )
    assert_allclose(
        Ncm.C.HI_Lyman_wl3_8pi_2p_2P0_5(),
        (1.0 / 8.22589191133e6) ** 3 / (8.0 * math.pi),
    )
    assert_allclose(
        Ncm.C.HI_Lyman_wl3_8pi_2p_2P3_5(),
        (1.0 / 8.22592850014e6) ** 3 / (8.0 * math.pi),
    )
    assert_allclose(
        Ncm.C.HI_Lyman_wl3_8pi_2p_2Pmean(),
        (1.0 / 8.225910205735e6) ** 3 / (8.0 * math.pi),
    )


def test_constants_boltzmann_factor():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.boltzmann_factor_HI_1s_2S0_5(5.0e3), 47450557.836754)
    assert_allclose(Ncm.C.boltzmann_factor_HI_2s_2S0_5(5.0e3), 9.040447e17)
    assert_allclose(Ncm.C.boltzmann_factor_HI_2p_2P0_5(5.0e3), 9.040355e17)
    assert_allclose(Ncm.C.boltzmann_factor_HI_2p_2P3_5(5.0e3), 9.041307e17)
    assert_allclose(Ncm.C.boltzmann_factor_HI_2p_2Pmean(5.0e3), 9.040831e17)


def test_constants_HeI_ion_wn():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HeI_ion_wn_1s_1S0(), 19831066.637)
    assert_allclose(Ncm.C.HeI_ion_wn_2s_1S0(), 3203322.6229)
    assert_allclose(Ncm.C.HeI_ion_wn_2s_3S1(), 3845469.20403)
    assert_allclose(Ncm.C.HeI_ion_wn_2p_1P1(), 2717576.9424)
    assert_allclose(Ncm.C.HeI_ion_wn_2p_3P0(), 2922283.55569)
    assert_allclose(Ncm.C.HeI_ion_wn_2p_3P1(), 2922382.34721)
    assert_allclose(Ncm.C.HeI_ion_wn_2p_3P2(), 2922389.98975)
    assert_allclose(Ncm.C.HeI_ion_wn_2p_3Pmean(), 2922351.964217)


def test_constants_HeI_Lyman_wn():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HeI_Lyman_wn_2s_1S0(), 16627744.0141)
    assert_allclose(Ncm.C.HeI_Lyman_wn_2s_3S1(), 15985597.43297)
    assert_allclose(Ncm.C.HeI_Lyman_wn_2p_1P1(), 17113489.6946)
    assert_allclose(Ncm.C.HeI_Lyman_wn_2p_3P0(), 16908783.08131)
    assert_allclose(Ncm.C.HeI_Lyman_wn_2p_3P1(), 16908684.28979)
    assert_allclose(Ncm.C.HeI_Lyman_wn_2p_3P2(), 16908676.64725)
    assert_allclose(Ncm.C.HeI_Lyman_wn_2p_3Pmean(), 16908714.672783)


def test_constants_HeI_Lyman_wl():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HeI_Lyman_wl_2s_1S0(), 1.0 / 16627744.0141)
    assert_allclose(Ncm.C.HeI_Lyman_wl_2s_3S1(), 1.0 / 15985597.43297)
    assert_allclose(Ncm.C.HeI_Lyman_wl_2p_1P1(), 1.0 / 17113489.6946)
    assert_allclose(Ncm.C.HeI_Lyman_wl_2p_3P0(), 1.0 / 16908783.08131)
    assert_allclose(Ncm.C.HeI_Lyman_wl_2p_3P1(), 1.0 / 16908684.28979)
    assert_allclose(Ncm.C.HeI_Lyman_wl_2p_3P2(), 1.0 / 16908676.64725)
    assert_allclose(Ncm.C.HeI_Lyman_wl_2p_3Pmean(), 1.0 / 16908714.672783)


def test_constants_HeI_Lyman_wl3_8pi():
    """Test constants from ncm_c.h"""

    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2s_1S0(), (1.0 / 16627744.0141) ** 3 / (8.0 * math.pi)
    )
    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2s_3S1(), (1.0 / 15985597.43297) ** 3 / (8.0 * math.pi)
    )
    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2p_1P1(), (1.0 / 17113489.6946) ** 3 / (8.0 * math.pi)
    )
    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2p_3P0(), (1.0 / 16908783.08131) ** 3 / (8.0 * math.pi)
    )
    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2p_3P1(), (1.0 / 16908684.28979) ** 3 / (8.0 * math.pi)
    )
    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2p_3P2(), (1.0 / 16908676.64725) ** 3 / (8.0 * math.pi)
    )
    assert_allclose(
        Ncm.C.HeI_Lyman_wl3_8pi_2p_3Pmean(),
        (1.0 / 16908714.672783) ** 3 / (8.0 * math.pi),
    )


def test_constants_boltzmann_factor_HeI():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.boltzmann_factor_HeI_1s_1S0(7.0e3), 4794.233313)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2s_1S0(7.0e3), 3.337521e18)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2s_3S1(7.0e3), 8.916898e17)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2p_1P1(7.0e3), 9.057814e18)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2p_3P0(7.0e3), 5.946928e18)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2p_3P1(7.0e3), 5.94572e18)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2p_3P2(7.0e3), 5.945627e18)
    assert_allclose(Ncm.C.boltzmann_factor_HeI_2p_3Pmean(7.0e3), 5.946092e18)


def test_constants_HeI_Balmer_wn():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HeI_Balmer_wn_2p_1P1_2s_1S0(), 485745.6805)
    assert_allclose(Ncm.C.HeI_Balmer_wn_2p_3Pmean_2s_3S1(), 923117.239813)


def test_constants_HeI_Balmer_E_kb():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HeI_Balmer_E_kb_2p_1P1_2s_1S0(), 6988.796535)
    assert_allclose(Ncm.C.HeI_Balmer_E_kb_2p_3Pmean_2s_3S1(), 13281.597399)


def test_constants_HeII_ion():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.HeII_ion_wn_1s_2S0_5(), 43890887.85)
    assert_allclose(Ncm.C.HeII_ion_E_1s_2S0_5(), 8.718687e-18)


def test_constants_decay():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.decay_H_rate_2s_1s(), 8.224581)
    assert_allclose(Ncm.C.decay_He_rate_2s_1s(), 51.3)


def test_constants_stats():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.stats_1sigma(), 0.682689492137)
    assert_allclose(Ncm.C.stats_2sigma(), 0.954499736104)
    assert_allclose(Ncm.C.stats_3sigma(), 0.997300203937)


def test_constants_hubble():
    """Test constants from ncm_c.h"""

    assert_allclose(Ncm.C.hubble_cte_planck6_base(), 67.36)
    assert_allclose(Ncm.C.hubble_cte_hst(), 72.0)
    assert_allclose(Ncm.C.hubble_radius_hm1_Mpc(), 2997.92458)
    assert_allclose(Ncm.C.hubble_radius_hm1_planck(), 5.723496e60)


def test_constants_crit_density():
    """Test constants from ncm_c.h"""

    print(Ncm.C.crit_mass_density_h2())
    assert_allclose(Ncm.C.crit_density_h2(), 1.6881692556555728e-09)
    assert_allclose(Ncm.C.crit_mass_density_h2(), 1.8783416169331677e-26)
    assert_allclose(Ncm.C.crit_mass_density_h2_solar_mass_Mpc3(), 2.775366e11)
    assert_allclose(Ncm.C.crit_number_density_p(), 11.229923)
    assert_allclose(Ncm.C.crit_number_density_n(), 11.214465)
    assert_allclose(Ncm.C.blackbody_energy_density(), 7.565733e-16)
    assert_allclose(Ncm.C.blackbody_per_crit_density_h2(), 4.48162e-07)
    assert_allclose(Ncm.C.radiation_temp_to_h2Omega_r0(200.0), 717.059214)
