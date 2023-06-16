#!/usr/bin/env python
#
# py_hiprim_example.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# py_hiprim_example.py
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

"""Example implementation of a primordial cosmology model."""

import math
from numcosmo_py import Nc, Ncm, GObject

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()

#
# New ModelBuilder object, defines a new model NcHIPrimExample implementing
# the Nc.HIPrim abstract class.
#
mb = Ncm.ModelBuilder.new(Nc.HIPrim, "NcHIPrimExample", "A example primordial model")

#
# New parameter A_s to describe the spectrum amplitud (it is usually better
# to work with ln(10^10As))
#
mb.add_sparam("A_s", "As", 0.0, 1.0, 0.1, 0.0, 1.0e-9, Ncm.ParamType.FREE)

#
# New parameter n_s to describe the spectral index
#
mb.add_sparam("n_s", "ns", 0.5, 1.5, 0.1, 0.0, 0.96, Ncm.ParamType.FREE)

#
# New parameters a, b and c to describe the spectrum modification
#
mb.add_sparam("a", "a", 0.0, 1.0, 0.01, 0.0, 0.5, Ncm.ParamType.FREE)
mb.add_sparam("b", "b", 0.0, 1.0e4, 0.10, 0.0, 100.0, Ncm.ParamType.FREE)
mb.add_sparam("c", "c", 0.0, 6.0, 0.10, 0.0, 0.0, Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet!
#
GNcHIPrimExample = mb.create()

#
# Workaraound to make the pygobject "see" the new object
#
GObject.new(GNcHIPrimExample)

#
# Gets the Python version of the object (.pytype) and register
# it as a PyGObject object.
#
NcHIPrimExample = GNcHIPrimExample.pytype
GObject.type_register(NcHIPrimExample)


#
# Creating a new class implementing our object NcHIPrimExample
#
class PyHIPrimExample(NcHIPrimExample):  # type: ignore
    """Example implementation of a primordial cosmology model."""

    #
    # Defining some property which is not part of the model paramers.
    # All model parameters must be defined by the ModelBuilder.
    #
    some_property = GObject.Property(type=str)

    #
    # Implementing the adiabatic spectrum function
    #
    def do_lnSA_powspec_lnk(self, lnk):
        """Return the ln of the adiabatic spectrum at a given lnk."""
        lnk0 = self.get_lnk_pivot()
        lnka = lnk - lnk0
        ka = math.exp(lnka)
        As = self.props.As
        ns = self.props.ns
        a = self.props.a
        b = self.props.b
        c = self.props.c
        a2 = a * a

        return (
            (ns - 1.0) * lnka
            + math.log(As)
            + math.log1p(a2 * math.cos(b * ka + c) ** 2)
        )

    def get_lnSA_powspec_lnk0(self):
        """Return the pivot lnk."""
        As = self.props.As
        a = self.props.a
        b = self.props.b
        c = self.props.c
        a2 = a * a

        return math.log(As) + math.log1p(a2 * math.cos(b + c) ** 2)


#
# Register our new Python class PyHIPrimExample
#
GObject.type_register(PyHIPrimExample)
