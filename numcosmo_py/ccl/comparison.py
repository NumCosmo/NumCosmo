#
# comparison.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# comparison.py
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

"""Compare CCL and NumCosmo results."""

import numpy as np
import pylab as plt


# pylint: disable-next=dangerous-default-value
def compare_ccl_nc_func(
    x,  # pylint: disable=invalid-name
    y_ccl,
    y_nc,
    x_name="x",
    y_name="func",
    subplots_pars={"figsize": (12, 6)},
    xscale="linear",
    yscale="log",
):
    """Compare CCL and NumCosmo results for a function."""

    ccl_name, nc_name = f"{y_name}^\\mathrm{{ccl}}", f"{y_name}^\\mathrm{{nc}}"

    x = np.array(x)
    y_ccl = np.array(y_ccl)
    y_nc = np.array(y_nc)
    diff = np.zeros_like(y_ccl)

    non_zind = np.where(y_ccl != 0.0)[0]
    zind = np.where(y_ccl == 0.0)[0]
    diff[non_zind] = y_nc[non_zind] / y_ccl[non_zind] - 1.0
    diff[zind] = y_nc[zind] - y_ccl[zind]
    print(
        f"[{y_name:10}]: rel diff min: {min(abs(diff)):.e}\t"
        f"rel diff max: {max(abs(diff)):.e}"
    )

    fig, axs = plt.subplots(2, sharex=True, **subplots_pars)
    fig.subplots_adjust(hspace=0)

    axs[0].plot(x, y_ccl, label="ccl", lw=3)
    axs[0].plot(x, y_nc, label="nc")
    axs[1].plot(x, np.abs(diff), c="r")
    axs[1].set_xscale(xscale)
    axs[1].set_yscale(yscale)

    axs[0].legend()
    axs[0].set_ylabel(f"${y_name}$")
    axs[1].set_xlabel(f"${x_name}$")
    axs[1].set_ylabel(f"${nc_name}/{ccl_name}-1$")

    return fig, axs
