#
# planck_native_release.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck_native_release.py
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

"""Curated native Planck 2018 likelihoods as downloadable serialized objects.

Mirrors the weak-lensing catalog release mechanism (``NcGalaxyWLObsCatalogId`` /
``nc_galaxy_wl_obs_new_from_catalog_id``): each native Planck likelihood is
pre-built once from the public ``plc_3.0`` clik data (the converters are the
data-reduction step) and serialized to a self-contained ``.gvar`` object. These
lightweight value-added products are hosted on a NumCosmo GitHub release;
``load_planck_release`` downloads (and caches) one by id and deserializes it, so
users can assemble native Planck experiments with **no clik data and no PLC
library**.

NumCosmo provides independent reimplementations of the Planck likelihoods,
including support for alternative resampled representations. The underlying
likelihood methodology and original data products are due to the Planck
Collaboration. Analyses using these likelihoods should cite the corresponding
Planck publications in addition to NumCosmo; see
``planck_native_provenance.md`` for the source data and the reference list.

Maintainers rebuild the release artifacts with :func:`build_release` and upload
them to the ``RELEASE_TAG`` GitHub release.
"""

from typing import Callable
import os
import urllib.request
from enum import StrEnum

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lite import (
    find_baseline_file,
    PLIK_LITE_TT_RELPATH,
    PLIK_LITE_TTTEEE_RELPATH,
    build_plik_lite,
    build_plik_lite_tt,
)
from numcosmo_py.experiments.planck_simall import (
    SIMALL_EE_RELPATH,
    SIMALL_BB_RELPATH,
    SIMALL_EEBB_RELPATH,
    build_simall,
)
from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH, build_commander
from numcosmo_py.experiments.planck_smica import (
    PLIK_TT_RELPATH,
    PLIK_TTTEEE_RELPATH,
    build_smica_tt,
    build_smica_ttteee,
)
from numcosmo_py.experiments.planck_lensing import (
    LENSING_FULL_RELPATH,
    LENSING_MARGED_RELPATH,
    build_lensing,
)

# The shared NumCosmo data-asset release, the same one the SNIa covariances, the
# weak-lensing catalogs and the Planck clik baseline tarball are served from
# (see nc_data_snia_cov.c, nc_galaxy_wl_obs.c and nc_data_planck_lkl.c).
RELEASE_TAG = "datafile-release-v1.0.0"
RELEASE_URL = f"https://github.com/NumCosmo/NumCosmo/releases/download/{RELEASE_TAG}"


class PlanckReleaseId(StrEnum):
    """Identifier of a curated native Planck likelihood object.

    The value carries the Planck data release the object was reduced from
    (``pr3`` = the 2018 ``plc_3.0`` package, the ``R3.00`` tarball on the same
    GitHub release), so a set built from another release lives alongside it
    instead of colliding with it.
    """

    PR3_COMMANDER = "pr3_commander"
    PR3_SIMALL_EE = "pr3_simall_ee"
    PR3_SIMALL_BB = "pr3_simall_bb"
    PR3_SIMALL_EEBB = "pr3_simall_eebb"
    PR3_PLIK_TT = "pr3_plik_tt"
    PR3_PLIK_TTTEEE = "pr3_plik_ttteee"
    PR3_PLIK_LITE_TT = "pr3_plik_lite_tt"
    PR3_PLIK_LITE_TTTEEE = "pr3_plik_lite_ttteee"
    PR3_LENSING = "pr3_lensing"
    PR3_LENSING_MARGED = "pr3_lensing_marged"


# id -> (source clik relpath, builder(clik_path, pb) -> NcmData). The builder is
# always called with pb=None here so the serialized object is data-only.
_REGISTRY: dict[PlanckReleaseId, tuple[str, Callable]] = {
    PlanckReleaseId.PR3_COMMANDER: (COMMANDER_RELPATH, build_commander),
    PlanckReleaseId.PR3_SIMALL_EE: (SIMALL_EE_RELPATH, build_simall),
    PlanckReleaseId.PR3_SIMALL_BB: (SIMALL_BB_RELPATH, build_simall),
    PlanckReleaseId.PR3_SIMALL_EEBB: (SIMALL_EEBB_RELPATH, build_simall),
    PlanckReleaseId.PR3_PLIK_TT: (PLIK_TT_RELPATH, build_smica_tt),
    PlanckReleaseId.PR3_PLIK_TTTEEE: (PLIK_TTTEEE_RELPATH, build_smica_ttteee),
    PlanckReleaseId.PR3_PLIK_LITE_TT: (PLIK_LITE_TT_RELPATH, build_plik_lite_tt),
    PlanckReleaseId.PR3_PLIK_LITE_TTTEEE: (PLIK_LITE_TTTEEE_RELPATH, build_plik_lite),
    PlanckReleaseId.PR3_LENSING: (LENSING_FULL_RELPATH, build_lensing),
    PlanckReleaseId.PR3_LENSING_MARGED: (LENSING_MARGED_RELPATH, build_lensing),
}


def release_filename(rid: PlanckReleaseId) -> str:
    """Release/cache filename for a given id.

    The id already carries the Planck data release (``pr3_...``), so the name is
    unique across releases: ``planck_native_pr3_commander.gvar``.
    """
    return f"planck_native_{rid.value}.gvar"


def _cache_path(rid: PlanckReleaseId, cache_dir: str | None) -> str:
    base = cache_dir if cache_dir is not None else Ncm.cfg_get_fullpath_base()
    return os.path.join(base, release_filename(rid))


def _ensure_downloaded(rid: PlanckReleaseId, cache_dir: str | None = None) -> str:
    """Return the local path for @rid, downloading it from the release if absent."""
    path = _cache_path(rid, cache_dir)
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        url = f"{RELEASE_URL}/{release_filename(rid)}"
        # Fetched under a private name and published by rename, so a reader
        # never sees a partial file and an interrupted fetch leaves nothing
        # that a later run would take for a complete download. The name carries
        # the pid: several processes may fetch the same id at once.
        tmp = f"{path}.{os.getpid()}.part"
        try:
            print(f"# Downloading file [{url}]...", flush=True)
            urllib.request.urlretrieve(url, tmp)  # nosec B310 - fixed https host
            os.replace(tmp, path)
        except Exception as exc:  # pylint: disable=broad-except
            if os.path.exists(tmp):
                os.remove(tmp)
            raise RuntimeError(
                f"cannot download {url}: {exc}. Download it manually to {path}, "
                f"or rebuild the release with build_release()."
            ) from exc
    return path


def load_planck_release(
    rid: PlanckReleaseId,
    pb: Nc.HIPertBoltzmann | None = None,
    cache_dir: str | None = None,
) -> Ncm.Data:
    """Download (and cache) a curated native Planck likelihood and deserialize it.

    With @pb set, the shared Boltzmann source is attached; the native object
    self-configures its targets/lmax in prepare(), so no clik data or PLC library
    is needed.
    """
    path = _ensure_downloaded(rid, cache_dir)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    data = ser.from_binfile(path)
    assert isinstance(
        data,
        Nc.DataPlanckCommander
        | Nc.DataPlanckLensing
        | Nc.DataPlanckLKL
        | Nc.DataPlanckPlikLite
        | Nc.DataPlanckSimall
        | Nc.DataPlanckSmica,
    )

    if pb is not None:
        data.set_hipert_boltzmann(pb)

    return data


def build_release(out_dir: str | None = None) -> list[str]:
    """Rebuild the release artifacts from the local ``plc_3.0`` clik data.

    Each curated likelihood is built data-only (no Boltzmann) and serialized to
    ``out_dir`` (default: the NumCosmo cache dir, so the files double as the local
    cache). Ids whose source clik data is missing are skipped. Returns the list of
    written paths. Upload the results to the ``RELEASE_TAG`` GitHub release.
    """
    if out_dir is None:
        out_dir = Ncm.cfg_get_fullpath_base()
    os.makedirs(out_dir, exist_ok=True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    written: list[str] = []

    for rid, (relpath, builder) in _REGISTRY.items():
        clik = find_baseline_file(relpath)
        if clik is None:
            continue
        data = builder(clik, None)  # data-only
        out = os.path.join(out_dir, release_filename(rid))
        ser.to_binfile(data, out)
        written.append(out)

    return written
