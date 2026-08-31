#!/usr/bin/env python3
#
# prefetch_data.py
#
# Sun August 30 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# prefetch_data.py
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

"""Download every data file NumCosmo fetches on demand into its data directory.

CI caches that directory between runs (see .github/actions/data-cache). The
cache is written under an exact key, so it is never refreshed while the key
holds: whatever the writing job downloaded is what every other job gets. A job
only downloads what its own tests touch, so a cache written from an ordinary
test run would leave the rest of the assets re-downloaded on every run,
forever. This fills the directory with all of them first, so the one job that
writes the cache writes a complete one.

The asset lists come from the enums themselves, so a catalog added to any of
them is picked up here without editing this script.
"""

import sys

from numcosmo_py import Nc, Ncm
from numcosmo_py.experiments import planck_native_release as pnr


def _enum_members(enum_cls, prefix: str = "") -> list[int]:
    """Enum values of @enum_cls whose name is uppercase and starts with @prefix.

    Introspected enums are plain ints carrying the members as class attributes,
    alongside int's own lowercase methods -- hence the isupper() filter.
    """
    return [
        getattr(enum_cls, name)
        for name in dir(enum_cls)
        if name.isupper() and name.startswith(prefix)
    ]


def main() -> int:
    """Fetch every on-demand data file, reporting each one."""
    Ncm.cfg_init()
    base = Ncm.cfg_get_fullpath_base()
    print(f"# Prefetching NumCosmo data files into {base}", flush=True)

    # Planck baseline (plc_3.0): a tarball unpacked into <base>/baseline.
    Nc.DataPlanckLKL.download_baseline(base)

    # SNIa covariance catalogs. Only the COV_* ids are downloaded; the SIMPLE_*
    # ones are built in.
    for cat_id in _enum_members(Nc.DataSNIAId, "COV_"):
        print(f"# {Nc.DataSNIACov.get_catalog_by_id(cat_id)}", flush=True)

    # Curated weak-lensing catalogs.
    for cat_id in _enum_members(Nc.GalaxyWLObsCatalogId):
        print(f"# {Nc.galaxy_wl_obs_catalog_id_get_filename(cat_id)}", flush=True)

    # Curated native Planck likelihoods. Loading, not just downloading: a
    # truncated .gvar then fails here instead of inside a test.
    for release_id in pnr.PlanckReleaseId:
        pnr.load_planck_release(release_id)
        print(f"# {pnr.release_filename(release_id)}", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
