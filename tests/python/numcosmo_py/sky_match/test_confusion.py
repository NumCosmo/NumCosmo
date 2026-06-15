#!/usr/bin/env python
#
# test_confusion.py
#
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
#

"""Tests for the catalog confusion-matrix metrics."""

import pytest
import numpy as np

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table

from numcosmo_py.catalog import (
    calculate_catalog_metrics,
    calculate_split_metrics,
    get_ratios,
)


def _object_column(values):
    """Build an astropy-compatible object column of variable-length lists."""
    col = np.empty(len(values), dtype=object)
    col[:] = values
    return col


def _split_halo_table() -> Table:
    """Halo match table with one of each TP/FP/TN/FN by construction."""
    table = Table()
    table["halo_id"] = [1, 2, 3, 4]
    #               TP    FN    FP    TN
    table["is_detected"] = [1, 1, 0, 0]
    table["RA_matched"] = _object_column([[10.0], [], [11.0], []])
    return table


def _split_cluster_table() -> Table:
    """Cluster match table with one of each TP/FP/TN/FN by construction."""
    table = Table()
    table["cluster_id"] = [10, 11, 12, 13]
    # parent_id != 0 -> real (positive); 0 -> fake (negative)
    table["parent_id"] = [100, 200, 0, 0]
    table["RA_matched"] = _object_column([[1.0], [], [2.0], []])
    return table


def test_split_metrics_halo() -> None:
    """One TP/FP/TN/FN each gives 0.5 ratios across the board."""
    result = calculate_split_metrics(_split_halo_table(), "halo")
    assert isinstance(result, dict)
    assert result["counts"] == {"TP": 1, "FP": 1, "TN": 1, "FN": 1}
    assert result["ratios"] == {"TP": 0.5, "FP": 0.5, "TN": 0.5, "FN": 0.5}
    assert result["sums"] == {"positive": 1.0, "negative": 1.0}


def test_split_metrics_cluster() -> None:
    """Cluster ground truth uses parent_id (!= 0 is a real object)."""
    result = calculate_split_metrics(_split_cluster_table(), "cluster")
    assert isinstance(result, dict)
    assert result["counts"] == {"TP": 1, "FP": 1, "TN": 1, "FN": 1}
    assert result["ratios"] == {"TP": 0.5, "FP": 0.5, "TN": 0.5, "FN": 0.5}


def test_split_metrics_empty_split_returns_error() -> None:
    """When every object is matched there is no undetected set; an error string is returned."""
    table = Table()
    table["halo_id"] = [1, 2]
    table["is_detected"] = [1, 0]
    table["RA_matched"] = _object_column([[10.0], [11.0]])

    result = calculate_split_metrics(table, "halo")
    assert isinstance(result, str)
    assert "empty" in result


def test_split_metrics_invalid_cat_type() -> None:
    """An unknown category type raises ValueError."""
    with pytest.raises(ValueError, match="Invalid cat_type"):
        calculate_split_metrics(_split_halo_table(), "galaxy")


def test_catalog_metrics_halo() -> None:
    """Unique-detection catalog scored against the ground-truth mock."""
    mock = Table()
    mock["halo_id"] = [1, 2, 3, 4]
    mock["is_detected"] = [1, 1, 0, 0]  # detected: {1, 2}; undetected: {3, 4}

    catalog = Table()
    catalog["halo_id"] = [1, 5]
    catalog["parent_id"] = [1, 0]  # row 1 points to its parent (TP); row 2 is a ghost (FP)

    result = calculate_catalog_metrics(catalog, mock, "halo")
    assert result["counts"] == {"TP": 1, "FP": 1, "TN": 2, "FN": 1}
    assert result["ratios"] == {"TP": 0.5, "FP": 0.5, "TN": 1.0, "FN": 0.5}
    # n_total_detections / len(mock)
    assert result["overall_ratio"] == pytest.approx(0.5)


def test_catalog_metrics_invalid_cat_type() -> None:
    """An unknown category type raises ValueError."""
    mock = Table()
    mock["halo_id"] = [1]
    mock["is_detected"] = [1]
    with pytest.raises(ValueError, match="Invalid cat_type"):
        calculate_catalog_metrics(mock, mock, "galaxy")


def test_get_ratios_passthrough() -> None:
    """get_ratios returns the ratios sub-dict when present."""
    result = calculate_split_metrics(_split_halo_table(), "halo")
    assert get_ratios(result) == {"TP": 0.5, "FP": 0.5, "TN": 0.5, "FN": 0.5}


def test_get_ratios_handles_error_string() -> None:
    """get_ratios degrades to NaNs for error strings or malformed input."""
    nan_ratios = get_ratios("Error: something went wrong")
    assert set(nan_ratios) == {"TP", "FP", "TN", "FN"}
    assert all(np.isnan(v) for v in nan_ratios.values())

    assert all(np.isnan(v) for v in get_ratios({"counts": {}}).values())
