#
# test_py_de_cont.py
#
# Sun Apr 21 21:40:35 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_de_cont.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcDECont class."""

from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_de_cont_setget():
    """Test basic functionality of NcDECont."""
    de_cont = Nc.DECont.new(Omegaw=0.3, OmegaL=0.7, cs2=1.0e-5, w=1.0e-5)

    de_cont.set_reltol(1.0e-12)
    de_cont.set_ti(1.0e-8)
    de_cont.set_tf(1.0e10)
    de_cont.set_k(1.0)

    assert_allclose(de_cont.get_reltol(), 1.0e-12)
    assert_allclose(de_cont.get_ti(), 1.0e-8)
    assert_allclose(de_cont.get_tf(), 1.0e10)
    assert_allclose(de_cont.get_k(), 1.0)


def test_de_cont_serialize():
    """Test serialization of NcDECont."""
    de_cont = Nc.DECont.new(Omegaw=0.3, OmegaL=0.7, cs2=1.0e-5, w=1.0e-5)

    de_cont.set_reltol(1.0e-12)
    de_cont.set_ti(1.0e-8)
    de_cont.set_tf(1.0e10)
    de_cont.set_k(1.0)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    de_cont2 = ser.dup_obj(de_cont)
    assert de_cont2 is not None
    assert de_cont2 is not de_cont

    assert_allclose(de_cont.get_reltol(), de_cont2.get_reltol())
    assert_allclose(de_cont.get_ti(), de_cont2.get_ti())
    assert_allclose(de_cont.get_tf(), de_cont2.get_tf())
    assert_allclose(de_cont.get_k(), de_cont2.get_k())


def test_de_cont_prepare_prop():
    """Test propagator preparation of NcDECont."""
    de_cont = Nc.DECont.new(Omegaw=0.3, OmegaL=0.7, cs2=1.0e-5, w=1.0e-5)

    de_cont.set_reltol(1.0e-12)
    de_cont.set_ti(1.0e-8)
    de_cont.set_tf(1.0e10)
    de_cont.set_k(1.0)

    de_cont.prepare_prop(None, 0.0, 1.0e-30, 1.0e1)


def test_de_cont_nonadiab():
    """Test nonadiabatic vacuum approximation of NcDECont."""
    de_cont = Nc.DECont.new(Omegaw=0.3, OmegaL=0.7, cs2=1.0e-5, w=1.0e-5)

    ti = 1.0e-8
    tf = 1.0e-1
    de_cont.set_reltol(1.0e-12)
    de_cont.set_ti(ti)
    de_cont.set_tf(tf)
    de_cont.set_k(1.0)

    de_cont.prepare_prop(None, 0.0, 1.0e-30, 1.0e1)

    state_prop0 = Ncm.CSQ1DState.new()
    state_prop0.set_up(Ncm.CSQ1DFrame.NONADIAB1, 0.0, 0.0, 0.0)

    state_prop = Ncm.CSQ1DState.new()
    state_nonadiab = Ncm.CSQ1DState.new()

    N = 100
    t_a = np.geomspace(ti, tf, N)

    for frame in [
        Ncm.CSQ1DFrame.ORIG,
        Ncm.CSQ1DFrame.NONADIAB1,
        Ncm.CSQ1DFrame.NONADIAB2,
    ]:
        for t in t_a:
            de_cont.evolve_prop_vector(None, state_prop0, frame, t, state_prop)
            de_cont.compute_nonadiab(None, t, state_nonadiab)

            de_cont.change_frame(None, state_nonadiab, frame)

            assert_allclose(
                state_prop.get_phi_Pphi(), state_nonadiab.get_phi_Pphi(), rtol=1.0e-4
            )


def test_de_cont_eval():
    """Test evaluation of NcDECont."""
    de_cont = Nc.DECont.new(Omegaw=0.3, OmegaL=0.7, cs2=1.0e-5, w=1.0e-5)

    ti = 1.0e0
    tf = 1.0e10
    de_cont.set_reltol(1.0e-12)
    de_cont.set_ti(ti)
    de_cont.set_tf(tf)
    de_cont.set_k(1.0)

    N = 100
    t_a = np.geomspace(ti, tf, N)

    for t in t_a:
        assert np.isfinite(de_cont.eval_xi(None, t))
        assert np.isfinite(de_cont.eval_F1(None, t))
        assert np.isfinite(de_cont.eval_F2(None, t))
