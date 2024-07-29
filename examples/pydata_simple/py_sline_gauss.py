#
# py_sline_gauss.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# py_sline_gauss.py
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

"""Example of a data class for the SLine model based on Ncm.DataGaussCov."""

from py_sline_model import PySLineModel

from numcosmo_py import Ncm, GObject


#
# Creating a new class implementing our object Ncm.Data
#
class PySLineGauss(Ncm.DataGaussCov):
    """A simple data class for the SLine model based on Ncm.DataGaussCov."""

    #
    # We need one vector property to save the independent variables x
    #
    xv = GObject.Property(type=Ncm.Vector, flags=GObject.PARAM_READWRITE)

    #
    # The contructor assigns some default values and calls the father's constructor.
    #
    def __init__(self, length=600):
        Ncm.DataGaussCov.__init__(self, n_points=length)
        n_points = super().get_size()

        if n_points > 0:
            self.xv = Ncm.Vector.new(n_points)
        else:
            self.xv = None

        self.cov_init = False

        #
        # Initializing to sane values
        #
        if n_points > 0:
            self.peek_cov().set_identity()
            self.xv.set_zero()

    #
    # Implements the virtual method get_length.
    #
    def do_get_length(self) -> int:  # pylint: disable-msg=arguments-differ
        return super().get_size()

    #
    # Implements the virtual method get_dof.
    #
    def do_get_dof(self) -> int:  # pylint: disable-msg=arguments-differ
        return super().get_size()

    #
    # Implements the virtual method `begin'.
    # This method usually do some groundwork in the data
    # before the actual calculations. For example, if the likelihood
    # involves the decomposition of a constant matrix, it can be done
    # during `begin' once and then used afterwards.
    #
    def do_begin(self):  # pylint: disable-msg=arguments-differ
        return

    #
    # Implements the virtual method `prepare'.
    # This method should do all the necessary calculations using mset
    # to be able to calculate the likelihood afterwards.
    #
    def do_prepare(self, _mset):  # pylint: disable-msg=arguments-differ
        return

    #
    # Implements the virtual method `mean_func'.
    # This method should compute the theoretical mean for the gaussian
    # distribution.
    #
    # pylint: disable-next=arguments-differ
    def do_mean_func(self, mset: Ncm.MSet, vp: Ncm.Vector) -> None:
        mid = mset.get_id_by_ns("NcPySLineModel")
        slm = mset.peek(mid)
        assert isinstance(slm, PySLineModel)
        n_points = super().get_size()

        for i in range(n_points):
            x = self.xv.get(i)
            vp.set(i, slm.f_x(x))

    def create_random_cov(self, _slm, rng):
        """Create a random covariance matrix."""

        # ya = [slm.f_x(x) for x in self.xv.dup_array()]
        # yv = Ncm.Vector.new_array(ya)

        # self.cov.fill_rand_cov2(yv, 0.5, 2.0, 15.0, rng)
        self.peek_cov().fill_rand_cov(0.5, 2.0, 15.0, rng)

        # Remove comment to print the covariance matrix
        # self.cov.log_vals ("# COV: ", "% 10.3g")


#
# Register our new Python class PySLineGauss
#
GObject.type_register(PySLineGauss)
