#
# kernels.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# kernels.py
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

"""The library kernel of a vocabulary window.

Every parameter a window carries in :mod:`numcosmo_py.certified.xcor.windows` is a
property of its analytic kernel, under the same name, and every kernel carries
``bessel-deriv`` and ``scale-dependence``. One keyword construction therefore covers
every window, its derivative weights and its scale-dependent variant.
"""

import typing

from numcosmo_py import Nc, Ncm
from numcosmo_py.certified.xcor import windows as W

FORM_CLASSES: typing.Final[dict[str, typing.Any]] = {
    "gauss": Nc.XcorKernelAnalyticGauss,
    "tophat": Nc.XcorKernelAnalyticTophat,
    "tophat_smooth": Nc.XcorKernelAnalyticTophatSmooth,
    "student_t": Nc.XcorKernelAnalyticStudentT,
    "power_exp": Nc.XcorKernelAnalyticPowerExp,
    "lensing": Nc.XcorKernelAnalyticLensing,
    "multi": Nc.XcorKernelAnalyticMulti,
}


def make_cosmo_bits() -> tuple[Nc.HICosmo, Nc.Distance, Ncm.Powspec]:
    """Return a cosmology, a prepared distance, and a closed-form power spectrum.

    ``NcmPowspecAnalytic`` rather than a transfer-function spline: sqrt(P) sits
    *inside* the function the closure fits, so a splined P would put its own
    interpolation error under everything measured, and it has no closed form for a
    reference to reproduce. BBKS is elementary, hence exactly reproducible in
    extended precision, and keeps the ln^2 k tail that a rational form loses. The
    analytic windows are closed forms in chi alone, so neither the cosmology nor the
    distance reaches a radial integral; the power spectrum enters C_ell only.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    return cosmo, dist, ps


def build_kernel(
    window: str | W.Window,
    dist: Nc.Distance,
    powspec: Ncm.Powspec,
    integrator: Ncm.SBesselIntegrator | None = None,
) -> Nc.XcorKernelRadial:
    """Construct the analytic kernel of a window, not yet prepared.

    The caller sets what it measures at (``l_limber``, ``reltol``, ``peak-epsilon``,
    the closure) and prepares it.

    :param window: a vocabulary name, current or historical, or a
        :class:`~numcosmo_py.certified.xcor.windows.Window`.
    :param dist: the distance object the kernel is constructed with.
    :param powspec: the power spectrum the kernel is constructed with.
    :param integrator: the radial integrator; a new ``NcmSBesselIntegratorLevin``
        for multipoles 0 to 8 when omitted. Pass one to share it across the
        kernels of a block, which is what makes the solver's factorization reuse
        possible; the ODE solver is not reentrant, so a shared integrator must not
        be used concurrently.
    :return: the kernel.
    """
    spec = W.get(window) if isinstance(window, str) else window
    props: dict[str, typing.Any] = dict(zip(W.FORM_PARAMS[spec.form], spec.params))

    if spec.form == "multi":
        props = {
            key: (
                Ncm.Vector.new_array(list(value)) if isinstance(value, tuple) else value
            )
            for key, value in props.items()
        }

    if spec.kdep:
        props["scale_dependence"] = Nc.XcorKernelRadialKDepGrowth.new(*spec.kdep)

    if integrator is None:
        integrator = Ncm.SBesselIntegratorLevin.new(0, 8)

    return FORM_CLASSES[spec.form](
        dist=dist,
        powspec=powspec,
        integrator=integrator,
        bessel_deriv=spec.deriv,
        **props,
    )
