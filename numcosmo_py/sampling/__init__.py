#
# __init__.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# __init__.py
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

from typing import Optional, Union, Type
from .. import Ncm, GEnum


class FitRunner(GEnum):
    """Fit algorithm for Ncm.Fit."""

    GSL_LS = Ncm.FitType.GSL_LS
    GSL_MM = Ncm.FitType.GSL_MM
    GSL_MMS = Ncm.FitType.GSL_MMS
    LEVMAR = Ncm.FitType.LEVMAR
    NLOPT = Ncm.FitType.NLOPT


class FitRunMessages(GEnum):
    """Fit messages for Ncm.Fit."""

    NONE = Ncm.FitRunMsgs.NONE
    SIMPLE = Ncm.FitRunMsgs.SIMPLE
    FULL = Ncm.FitRunMsgs.FULL


class FitGradType(GEnum):
    """Fit gradient type for Ncm.Fit."""

    NUMDIFF_FORWARD = Ncm.FitGradType.NUMDIFF_FORWARD
    NUMDIFF_CENTRAL = Ncm.FitGradType.NUMDIFF_CENTRAL
    NUMDIFF_ACCURATE = Ncm.FitGradType.NUMDIFF_ACCURATE


def get_algorithms(
    runner: FitRunner,
) -> Optional[
    Union[
        Type[Ncm.FitNloptAlgorithm],
        Type[Ncm.FitLevmarAlgos],
        Type[Ncm.FitGSLMMSAlgos],
        Type[Ncm.FitGSLMMAlgos],
    ]
]:
    """Get algorithms for a given runner."""

    if runner == FitRunner.NLOPT:
        return Ncm.FitNloptAlgorithm
    if runner == FitRunner.LEVMAR:
        return Ncm.FitLevmarAlgos
    if runner == FitRunner.GSL_MMS:
        return Ncm.FitGSLMMSAlgos
    if runner == FitRunner.GSL_MM:
        return Ncm.FitGSLMMAlgos
    if runner == FitRunner.GSL_LS:
        return None
    raise RuntimeError(f"Runner {runner} not found.")


def check_runner_algorithm(runner: FitRunner, algorithm: str):
    """Check if algorithm is valid."""

    if algorithm is not None:
        algorithms = get_algorithms(runner)
        if algorithms is None:
            raise RuntimeError(f"Runner {runner} do not support algorithms.")
        if Ncm.cfg_get_enum_by_id_name_nick(algorithms, algorithm) is None:
            Ncm.cfg_enum_print_all(algorithms, "Allowed algorithms")
            raise RuntimeError(f"Algorithm {algorithm} not found for runner {runner}.")
