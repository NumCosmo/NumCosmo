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

import shutil
import math
import numpy as np
import numpy.typing as npt

from matplotlib.patches import Ellipse
from matplotlib import transforms
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm


def confidence_ellipse(mu, cov, ax, n_std=1.0, facecolor="none", **kwargs):
    """Add confidence ellipse.

    This function adds an confidence ellipse to the given axes, given the mean
    and covariance of the data.
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


def add_ellipse_from_ellipticity(
    ax: plt.Axes,
    ra_a: npt.NDArray[np.float64],
    dec_a: npt.NDArray[np.float64],
    e1_a: npt.NDArray[np.float64],
    e2_a: npt.NDArray[np.float64],
    ellipse_scale: float = 0.5,
    ellipse_area: float = 1.0,
    edgecolor: str = "black",
    facecolor: str = "none",
) -> None:
    r"""Add an ellipse to a plot based on the ellipticity.

    The ellipticity is given by the e1 and e2 components of the ellipticity vector. The
    semi-major and semi-minor axes are calculated from the ellipticity and the angle of
    the ellipse is calculated from the ellipticity components. The equation for the
    angle is given by:
    $$
    \theta = 0.5 \arctan2(e2, e1),$$ $$q = \frac{1 - \epsilon}{1 + \epsilon},
    $$
    where $q$ is the axis ratio, $a$ is the semi-major axis, and $b$ is the semi-minor
    axis. The area of the ellipse is given by:
    $$
    A = \pi a b.
    $$

    :param ax: The axes to add the ellipse to.
    :param ra_a: The right ascension of the center of the ellipse.
    :param dec_a: The declination of the center of the ellipse.
    :param e1_a: The e1 component of the ellipticity.
    :param e2_a: The e2 component of the ellipticity.
    :param ellipse_scale: The scale factor for the ellipse axes.
    :param ellipse_area: The area of the ellipse.
    :param edgecolor: The color of the ellipse edge.
    :param facecolor: The color of the ellipse face.

    """
    for ra, dec, e1, e2 in zip(ra_a, dec_a, e1_a, e2_a):
        # Calculate the ellipse parameters
        ellip: float = np.hypot(e1, e2)
        theta: float = 0.5 * np.arctan2(e2, e1)

        # Convert ellipticity to semi-major and semi-minor axes
        q = (1.0 - ellip) / (1.0 + ellip)
        a = np.sqrt(ellipse_area / q)
        b = np.sqrt(ellipse_area * q)

        # Create an ellipse
        ellipse = Ellipse(
            xy=(ra, dec),
            width=ellipse_scale * a,
            height=ellipse_scale * b,
            angle=np.degrees(theta),
            edgecolor=edgecolor,
            facecolor=facecolor,
        )

        # Add the ellipse to the plot
        ax.add_patch(ellipse)


def latex_float(value: float, precision: int = 2, convert_g: bool = True) -> str:
    """Convert a float to a string with a fixed number of decimal places."""
    if convert_g:
        float_str = f"{value:.{precision}g}"
    else:
        float_str = f"{value:.{precision}e}"
    if "e" in float_str:
        base, exponent = float_str.split("e")
        if math.isclose(float(base), 1.0):
            if exponent == 0.0:
                return r"1"
            return f"10^{{{int(exponent)}}}"

        if exponent == 0.0:
            return f"{base}"

        return rf"{base} \times 10^{{{int(exponent)}}}"

    return float_str


def format_time(value: float) -> str:
    """Format a time value in seconds with appropriate units."""
    units = [
        (1.0e-9, r"$\mu$s"),  # Microseconds
        (1.0e-6, "ms"),  # Milliseconds
        (1.0e-3, "s"),  # Seconds
        (1.0, "s"),  # Keep in seconds for larger values
    ]

    # Choose the most appropriate unit
    for factor, unit in units:
        if value < factor * 1000:  # Keep within 3 significant digits
            return f"{value / factor:.2f} {unit}"

    return f"{value:.2f} s"  # Fallback (large values)


def set_rc_params_article(
    column_width: float = 246.0,
    ncol: int = 2,
    nrows: int = 1,
    fontsize: int = 8,
    use_tex: bool | None = None,
) -> None:
    """Set matplotlib rcParams for a LaTeX article."""
    fig_width_pt = column_width * ncol  # \showthe\columnwidth
    inches_per_pt = 1.0 / 72.27  # Convert pt to inch
    golden_mean = (math.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean  # height in inches
    fig_size = [fig_width, fig_height * nrows]
    if use_tex is None:
        use_tex = bool(shutil.which("latex"))

    params = {
        "axes.labelsize": fontsize,
        "font.size": fontsize,
        "legend.fontsize": fontsize,
        "xtick.labelsize": fontsize,
        "ytick.labelsize": fontsize,
        "text.usetex": use_tex,
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
