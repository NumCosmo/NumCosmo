#
# tools.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# tools.py
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

"""NumCosmo's plotting tools."""

import math
import numpy as np

from matplotlib.patches import Ellipse
from matplotlib import transforms
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm


def confidence_ellipse(mu, cov, ax, n_std=1.0, facecolor="none", **kwargs):
    """Adds a confidence ellipse to the given axis based on the given
    covariance matrix.
    """
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse(
        (0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs,
    )

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = mu[0]

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = mu[1]

    transf = (
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(scale_x, scale_y)
        .translate(mean_x, mean_y)
    )

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def latex_float(value: float):
    """Convert a float to a string with a fixed number of decimal places."""
    float_str = f"{value:.2g}"
    if "e" in float_str:
        base, exponent = float_str.split("e")
        if base == 1.0:
            if exponent == 0.0:
                return r"1"
            return f"10^{{{int(exponent)}}}"

        if exponent == 0.0:
            return f"{base}"

        return f"{base} \times 10^{{{int(exponent)}}}"

    return float_str


def set_rc_params_article(column_width: float = 246.0, ncol: int = 2, nrows: int = 1):
    """Set matplotlib rcParams for a LaTeX article."""
    fig_width_pt = column_width * ncol  # \showthe\columnwidth
    inches_per_pt = 1.0 / 72.27  # Convert pt to inch
    golden_mean = (math.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean  # height in inches
    fig_size = [fig_width, fig_height * nrows]

    params = {
        "axes.labelsize": 8,
        "font.size": 8,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "text.usetex": True,
        "figure.figsize": fig_size,
    }

    plt.rcParams.update(params)


def plot_m2lnp(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    ax: plt.Axes,
    *,
    plotn: int = 150,
    vmin: float = 1.0e-12,
    vmax: float = 1.0,
):
    """Plot the -2lnp."""

    z = z - np.min(z)
    z = np.exp(-0.5 * z)
    exp_z = z.reshape(plotn, plotn)

    img = ax.imshow(
        exp_z,
        interpolation="bicubic",
        origin="lower",
        cmap=cm.gray_r,  # type: ignore # pylint:disable-msg=no-member
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extent=(x[0], x[-1], y[0], y[-1]),
        aspect="auto",
        rasterized=True,
    )

    return img
