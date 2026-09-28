#!/usr/bin/env python
#
# make_sphere_healpy_truth_table.py
#
# Sun Sep 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sphere_healpy_truth_table.py
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

"""Generate the healpy truth tables of NcmSphereMap.

Writes, in ``data/truth_tables/sphere/``:

``healpy_pixels.bin``
    For nside 1, 2, 4 and 8, one row per RING index ``i``: ``ring2nest(i)``,
    ``nest2ring(i)`` (``i`` read as a NESTED index), and the centre ``theta, phi`` of
    RING pixel ``i``. For nside 128 and 256, 256 seeded pixels, with the index itself
    as a first column.
``healpy_ang2pix.bin``
    For nside 1, 16 and 256: ``theta, phi, ring pixel, nest pixel`` for 150 seeded
    directions and 50 chosen ones (the poles, the edge of the polar caps at
    ``z = 2/3``, ``phi`` at 0, just below 2 pi, negative and beyond 2 pi).
``healpy_transforms.bin``
    A seeded nside-8 map (``map``) and, for lmax 23 (3 nside - 1) and 32 (4 nside),
    healpy's ``map2alm`` at iter 0 and 3 (``alm_lmax<L>_iter<I>``, rows ``re, im`` in
    healpy's index order, m-major), ``anafast`` (``cl_lmax<L>_iter<I>``) and the
    ``alm2map`` of the iter-3 coefficients (``synth_lmax<L>``), and with a second seeded
    map (``map2``) the cross spectrum ``anafast(map, map2, iter=3)``
    (``cross_lmax<L>``). All transforms use healpy's defaults except
    ``use_weights=False``, the quadrature NcmSphereMap uses.
``healpy_map_nest_nside16.fits``
    A NESTED nside-16 celestial map written by healpy with 1024 pixels per row, its
    value at NESTED pixel ``i`` being ``(i % 97) / 7 - 3``, exact in double precision.
``fits_noorder.fits``, ``fits_explicit.fits``, ``fits_car.fits``, ``fits_short.fits``
    nside-4 tables written with astropy for the header cases of
    ``ncm_sphere_map_load_fits()``: no ORDERING (loads as RING), a partial-sky map
    (INDXSCHM = EXPLICIT), another pixelization (PIXTYPE = CAR) and a column of 100 values.

The binfiles are NcmObjDictStr of NcmMatrix/NcmVector, read with
``ncm_serialize_dict_str_from_binfile()``. Integers are stored as doubles, exact far
beyond these sizes. The consumer is ``tests/c/ncm/sphere/test_ncm_sphere_map.c``; it
needs neither healpy nor this script, which is re-run only when a case changes.
"""

import argparse
from pathlib import Path

import healpy
import numpy as np

from numcosmo_py import Ncm

SEED = 20260927


def _matrix(rows: np.ndarray) -> Ncm.Matrix:
    rows = np.ascontiguousarray(rows, dtype=np.float64)
    return Ncm.Matrix.new_array(rows.ravel().tolist(), rows.shape[1])


def _vector(values: np.ndarray) -> Ncm.Vector:
    return Ncm.Vector.new_array(np.asarray(values, dtype=np.float64).tolist())


def _pixels(rng: np.random.Generator) -> Ncm.ObjDictStr:
    table = Ncm.ObjDictStr.new()
    for nside in (1, 2, 4, 8):
        idx = np.arange(healpy.nside2npix(nside))
        theta, phi = healpy.pix2ang(nside, idx)
        rows = np.column_stack(
            [healpy.ring2nest(nside, idx), healpy.nest2ring(nside, idx), theta, phi]
        )
        table.add(f"nside{nside}", _matrix(rows))
    for nside in (128, 256):
        idx = np.sort(rng.choice(healpy.nside2npix(nside), 256, replace=False))
        theta, phi = healpy.pix2ang(nside, idx)
        rows = np.column_stack(
            [
                idx,
                healpy.ring2nest(nside, idx),
                healpy.nest2ring(nside, idx),
                theta,
                phi,
            ]
        )
        table.add(f"nside{nside}", _matrix(rows))
    return table


def _directions(rng: np.random.Generator) -> np.ndarray:
    theta = np.arccos(rng.uniform(-1.0, 1.0, 150))
    phi = rng.uniform(0.0, 2.0 * np.pi, 150)
    cap = np.arccos(2.0 / 3.0)
    edge_theta = [0.0, 1.0e-10, np.pi, np.pi - 1.0e-10, cap, np.pi - cap]
    edge_theta += [cap - 1.0e-12, cap + 1.0e-12, np.pi / 2.0, 1.0e-3, np.pi - 1.0e-3]
    edge_phi = [0.0, np.nextafter(2.0 * np.pi, 0.0), -np.pi / 4.0, 7.0 * np.pi, 0.3]
    pairs = [(t, p) for t in edge_theta for p in edge_phi][:50]
    return np.vstack([np.column_stack([theta, phi]), np.array(pairs, dtype=np.float64)])


def _ang2pix(rng: np.random.Generator) -> Ncm.ObjDictStr:
    table = Ncm.ObjDictStr.new()
    directions = _directions(rng)
    theta, phi = directions[:, 0], directions[:, 1]
    for nside in (1, 16, 256):
        rows = np.column_stack(
            [
                theta,
                phi,
                healpy.ang2pix(nside, theta, phi, nest=False),
                healpy.ang2pix(nside, theta, phi, nest=True),
            ]
        )
        table.add(f"nside{nside}", _matrix(rows))
    return table


def _transforms(rng: np.random.Generator) -> Ncm.ObjDictStr:
    nside = 8
    table = Ncm.ObjDictStr.new()
    skymap = rng.standard_normal(healpy.nside2npix(nside))
    skymap2 = rng.standard_normal(healpy.nside2npix(nside))
    table.add("map", _vector(skymap))
    table.add("map2", _vector(skymap2))
    for lmax in (3 * nside - 1, 4 * nside):
        for it in (0, 3):
            alm = healpy.map2alm(skymap, lmax=lmax, iter=it, use_weights=False)
            cl = healpy.anafast(skymap, lmax=lmax, iter=it, use_weights=False)
            table.add(
                f"alm_lmax{lmax}_iter{it}",
                _matrix(np.column_stack([alm.real, alm.imag])),
            )
            table.add(f"cl_lmax{lmax}_iter{it}", _vector(cl))
        synth = healpy.alm2map(alm, nside, lmax=lmax)
        table.add(f"synth_lmax{lmax}", _vector(synth))
        cross = healpy.anafast(skymap, skymap2, lmax=lmax, iter=3, use_weights=False)
        table.add(f"cross_lmax{lmax}", _vector(cross))
    return table


def _header_fixtures(output_dir: Path) -> None:
    from astropy.io import fits  # pylint: disable=import-outside-toplevel

    cases = {
        "fits_noorder.fits": (192, {}),
        "fits_explicit.fits": (192, {"ORDERING": "RING", "INDXSCHM": "EXPLICIT"}),
        "fits_car.fits": (192, {"ORDERING": "RING", "PIXTYPE": "CAR"}),
        "fits_short.fits": (100, {"ORDERING": "RING"}),
    }
    for name, (npix, keys) in cases.items():
        values = np.arange(float(npix))
        table = fits.BinTableHDU.from_columns(
            fits.ColDefs([fits.Column(name="SIGNAL", format="D", array=values)])
        )
        table.header["NSIDE"] = 4
        table.header["COORDSYS"] = "C"
        for key, value in keys.items():
            table.header[key] = value
        fits.HDUList([fits.PrimaryHDU(), table]).writeto(
            output_dir / name, overwrite=True
        )


def main() -> None:
    """Write the tables."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "data"
        / "truth_tables"
        / "sphere",
    )
    args = parser.parse_args()

    Ncm.cfg_init()
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    rng = np.random.default_rng(SEED)

    ser.dict_str_to_binfile(_pixels(rng), str(args.output_dir / "healpy_pixels.bin"))
    ser.dict_str_to_binfile(_ang2pix(rng), str(args.output_dir / "healpy_ang2pix.bin"))
    ser.dict_str_to_binfile(
        _transforms(rng), str(args.output_dir / "healpy_transforms.bin")
    )

    nest_map = (np.arange(healpy.nside2npix(16)) % 97) / 7.0 - 3.0
    healpy.write_map(
        str(args.output_dir / "healpy_map_nest_nside16.fits"),
        nest_map,
        nest=True,
        coord="C",
        dtype=np.float64,
        overwrite=True,
    )
    _header_fixtures(args.output_dir)


if __name__ == "__main__":
    main()
