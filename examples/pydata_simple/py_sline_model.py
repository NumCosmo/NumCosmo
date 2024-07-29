#
# py_sline_mfunc.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# py_sline_mfunc.py
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

"""Example of defining a new Ncm.Model for SLine."""

from numcosmo_py import Ncm, GObject

#
# New ModelBuilder object, defines a new model NcHIPrimExample
# implementing the base Ncm.Model abstract class.
#
mb = Ncm.ModelBuilder.new(Ncm.Model, "NcPySLineModel", "A simple python example model")

#
# New parameter `alpha' to describe the ln-slope
# Allowed interval: [0, 5]; Default scale: 0.1
# Absolute tolerance: 0; Default value: 2
#
mb.add_sparam(r"\alpha", "alpha", 0.0, 5.0, 0.1, 0.0, 2.0, Ncm.ParamType.FREE)

#
# New parameter `a' to describe the amplitude
# Allowed interval: [0.2, 2]; Default scale: 0.1
# Absolute tolerance: 0; Default value: 1
#
mb.add_sparam("a", "a", 0.2, 2.0, 0.1, 0.0, 1.0, Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet! Then register
# the new object in the GObject type system by creating a new
# instance. Finally, gets the Python version of the object (.pytype)
# and register it as a PyGObject object.
#
GNcPySLineModel = mb.create()
GObject.new(GNcPySLineModel)
NcPySLineModel = GNcPySLineModel.pytype
GObject.type_register(NcPySLineModel)


#
# Creating a new class implementing our object NcPySLineModel
#
class PySLineModel(NcPySLineModel):  # type: ignore
    """A simple python example model."""

    #
    # Defining some property which is not part of the model paramers.
    # All model parameters must be defined by the ModelBuilder.
    #
    some_property = GObject.Property(type=str)

    #
    # Calling the father's constructor
    #
    def __init__(self):
        NcPySLineModel.__init__(self)

    #
    # Method to calculate the y(x)
    #
    def f_x(self, x):
        """Method to calculate the y(x)."""
        return self.props.alpha * x + self.props.a


#
# Register our new Python class PyNcPySLineModel
#
GObject.type_register(PySLineModel)
