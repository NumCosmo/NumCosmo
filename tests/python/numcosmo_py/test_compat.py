#
# test_compat.py
#
# Sat Oct 3 10:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_compat.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on the deprecated NumCosmo 0.27 Python API shims."""

import warnings

import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()


def test_add_submodel_after_construction():
    """A late attach fills the empty slot and warns."""
    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw()
    reion = Nc.HIReionCamb()

    with pytest.warns(DeprecationWarning, match="after construction"):
        cosmo.add_submodel(prim)
    with pytest.warns(DeprecationWarning, match="after construction"):
        cosmo.add_submodel(reion)

    assert cosmo.peek_prim() is prim
    assert cosmo.peek_reion() is reion
    assert prim.peek_host() is cosmo
    assert reion.peek_host() is cosmo


def test_yp_with_default_bbn_is_ignored():
    """With the default NcBBNParthenope, Yp is ignored, as a fixed Yp was in 0.27."""
    cosmo = Nc.HICosmoDEXcdm()
    yp_before = cosmo.Yp_4He()

    with pytest.warns(DeprecationWarning, match="Yp is no longer"):
        cosmo.param_set_by_name("Yp", 0.30)

    assert cosmo.Yp_4He() == yp_before


def test_yp_with_bbn_parametrized_is_routed():
    """With an NcBBNParametrized, Yp is set on the BBN submodel."""
    cosmo = Nc.HICosmoDEXcdm(bbn=Nc.BBNParametrized())

    with pytest.warns(DeprecationWarning, match="Yp is no longer"):
        cosmo.param_set_by_name("Yp", 0.2454)

    assert_allclose(cosmo.Yp_4He(), 0.2454, rtol=1.0e-12)


def test_param_set_by_name_other_names_do_not_warn():
    """Names other than Yp go straight to the original method."""
    cosmo = Nc.HICosmoDEXcdm()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cosmo.param_set_by_name("H0", 70.0)

    assert cosmo.param_get_by_name("H0") == 70.0


def test_set_z_from_tau_old_signature_without_host():
    """The 0.27 call on an unattached reion matches the current call on an attached one."""
    reion_new = Nc.HIReionCamb()
    host = Nc.HICosmoDEXcdm(reion=reion_new)
    reion_new.set_z_from_tau(0.06)

    cosmo = Nc.HICosmoDEXcdm()
    reion_old = Nc.HIReionCamb()
    with pytest.warns(DeprecationWarning, match="set_z_from_tau"):
        reion_old.set_z_from_tau(cosmo, 0.06)

    assert_allclose(
        reion_old.orig_param_get_by_name("z_re"),
        reion_new.orig_param_get_by_name("z_re"),
        rtol=1.0e-12,
    )
    assert reion_new.peek_host() is host


def test_set_z_from_tau_old_signature_with_host():
    """The 0.27 call on an attached reion uses the host."""
    reion = Nc.HIReionCamb()
    cosmo = Nc.HICosmoDEXcdm(reion=reion)

    with pytest.warns(DeprecationWarning, match="set_z_from_tau"):
        reion.set_z_from_tau(cosmo, 0.06)

    assert_allclose(reion.get_tau(cosmo), 0.06, rtol=1.0e-6)


def test_set_z_from_tau_old_signature_with_reparam():
    """The 0.27 call on a tau-reparametrized reion sets tau."""
    reion = Nc.HIReionCamb()
    cosmo = Nc.HICosmoDEXcdm(reion=reion)
    reion.z_to_tau()

    with pytest.warns(DeprecationWarning, match="set_z_from_tau"):
        reion.set_z_from_tau(cosmo, 0.06)

    assert_allclose(reion["tau_reion"], 0.06, rtol=1.0e-10)


def test_set_z_from_tau_old_signature_wrong_cosmo():
    """The 0.27 call with a cosmology other than the host raises."""
    reion = Nc.HIReionCamb()
    host = Nc.HICosmoDEXcdm(reion=reion)
    other = Nc.HICosmoDEXcdm()
    assert reion.peek_host() is host

    with (
        pytest.warns(DeprecationWarning, match="set_z_from_tau"),
        pytest.raises(ValueError, match="not the host"),
    ):
        reion.set_z_from_tau(other, 0.06)


def test_set_z_from_tau_new_signature_does_not_warn():
    """The current call goes straight to the original method."""
    reion = Nc.HIReionCamb()
    cosmo = Nc.HICosmoDEXcdm(reion=reion)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        reion.set_z_from_tau(0.06)

    assert_allclose(reion.get_tau(cosmo), 0.06, rtol=1.0e-6)
