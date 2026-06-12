"""Confusion-matrix metrics for sky-match catalogs.

Given a catalog matched with :mod:`numcosmo_py.catalog.sky_match` and the
ground-truth mock it was drawn from, these helpers score the match as a binary
classification problem (true/false positives and negatives).

Two complementary views are provided:

- :func:`calculate_split_metrics` scores a *complete* match table (one row per
  query object, with a possibly empty list of matched candidates) by splitting
  it into objects the algorithm matched and objects it did not.
- :func:`calculate_catalog_metrics` scores a *unique* detection catalog against
  the full ground-truth mock catalog.

Both return a nested ``dict`` with ``counts``, ``ratios`` and ``sums`` entries;
:func:`get_ratios` safely extracts the ``ratios`` sub-dict (returning NaNs when a
metric could not be computed).
"""

from typing import Literal

import numpy as np
from astropy.table import Table

CatalogType = Literal["halo", "cluster"]


def get_ratios(metrics_result: dict | str) -> dict[str, float]:
    """Extract the TP/FP/TN/FN ratios from a metrics result.

    :param metrics_result: the return value of :func:`calculate_split_metrics` or
        :func:`calculate_catalog_metrics`. When a metric could not be computed the
        functions return an error string instead of a dict.
    :return: the ``ratios`` sub-dict, or a dict of NaNs if it is unavailable.
    """
    if isinstance(metrics_result, dict) and "ratios" in metrics_result:
        return metrics_result["ratios"]
    # Returns NaN if metrics failed or returned a string, so callers don't crash.
    return {"TP": np.nan, "FP": np.nan, "TN": np.nan, "FN": np.nan}


def calculate_split_metrics(table: Table, cat_type: CatalogType = "halo") -> dict | str:
    """Confusion-matrix metrics from a complete (split) match table.

    Splits ``table`` into detected (rows with at least one matched candidate) and
    undetected (rows with an empty ``RA_matched``) sets, then compares the object
    IDs against the ground truth (``is_detected`` for halos, ``parent_id`` for
    clusters).

    :param table: complete match table, with an ``RA_matched`` column of
        candidate lists and the relevant ID/ground-truth columns.
    :param cat_type: ``"halo"`` or ``"cluster"``.
    :return: a dict with ``counts``, ``ratios`` and ``sums`` entries, or an error
        string when one of the split tables is empty.
    """
    # 1. Configuration based on category type
    if cat_type == "halo":
        id_col = "halo_id"
        gt_col = "is_detected"
    elif cat_type == "cluster":
        id_col = "cluster_id"
        gt_col = "parent_id"
    else:
        raise ValueError(f"Invalid cat_type: '{cat_type}'. Use 'halo' or 'cluster'.")

    # 2. Split the table by algorithm performance (Did we find a match?)
    # mask_has_data: rows where the algorithm matched a candidate
    mask_has_data = np.array([len(d) > 0 for d in table["RA_matched"]])
    detected_table = table[mask_has_data]
    undetected_table = table[~mask_has_data]

    n_detected = len(detected_table)
    n_undetected = len(undetected_table)

    if n_detected == 0 or n_undetected == 0:
        return "Error: One of the split tables is empty. check your distance data."

    # 3. Define Ground Truth subsets
    if cat_type == "halo":
        gt_positives = table[table[gt_col] == 1]
        gt_negatives = table[table[gt_col] == 0]
    else:
        gt_positives = table[table[gt_col] != 0]
        gt_negatives = table[table[gt_col] == 0]

    # 4. Confusion Matrix Logic

    # TP: Real objects found in the detected table
    tp_count = np.sum(np.isin(gt_positives[id_col], detected_table[id_col]))

    # FP: "Ghost" objects (GT says 0) found in the detected table
    fp_count = np.sum(np.isin(gt_negatives[id_col], detected_table[id_col]))

    # TN: "Ghost" objects correctly left in the undetected table
    tn_count = np.sum(np.isin(gt_negatives[id_col], undetected_table[id_col]))

    # FN: Real objects incorrectly left in the undetected table
    fn_count = np.sum(np.isin(gt_positives[id_col], undetected_table[id_col]))

    # 5. Ratios
    tp_ratio = tp_count / n_detected
    fp_ratio = fp_count / n_detected
    tn_ratio = tn_count / n_undetected
    fn_ratio = fn_count / n_undetected

    return {
        "counts": {"TP": tp_count, "FP": fp_count, "TN": tn_count, "FN": fn_count},
        "ratios": {"TP": tp_ratio, "FP": fp_ratio, "TN": tn_ratio, "FN": fn_ratio},
        "sums": {"positive": tp_ratio + fp_ratio, "negative": tn_ratio + fn_ratio},
    }


def calculate_catalog_metrics(
    catalog_table: Table, mock_table: Table, cat_type: CatalogType = "halo"
) -> dict:
    """Confusion-matrix metrics for a unique detection catalog.

    Compares a catalog of unique detections against the ground-truth mock it was
    drawn from. True/false positives are decided by whether each detection points
    back to its true parent (``halo_id == parent_id``); true/false negatives are
    decided by the ground-truth objects that were, or were not, recovered.

    :param catalog_table: table of unique detections (e.g. ``unique_halos`` or
        ``unique_detections``).
    :param mock_table: the ground-truth catalog (e.g. ``halos_mock`` or
        ``clusters_mock``).
    :param cat_type: ``"halo"`` or ``"cluster"``.
    :return: a dict with ``counts``, ``ratios``, ``sums`` and ``overall_ratio``.
    """
    # 1. Handle cat_type specific logic
    if cat_type == "halo":
        # Halos use 'is_detected' 0 or 1
        undetected_mask = mock_table["is_detected"] == 0
        detected_mask = mock_table["is_detected"] == 1
        id_col = "halo_id"
    elif cat_type == "cluster":
        # Clusters use 'parent_id' 0 or not 0
        undetected_mask = mock_table["parent_id"] == 0
        detected_mask = mock_table["parent_id"] != 0
        id_col = "cluster_id"
    else:
        raise ValueError(
            f"Invalid cat_type: '{cat_type}'. Expected 'halo' or 'cluster'."
        )

    # Filter mock tables based on the logic above
    mock_undetected_table = mock_table[undetected_mask]
    mock_detected_table = mock_table[detected_mask]

    # 2. Total Counts for Denominators
    n_total_detections = len(catalog_table[id_col])
    # n_non_detections is the difference between found matches and the ground truth catalog
    n_non_detections = abs(len(catalog_table[id_col]) - len(mock_table[id_col]))

    # 3. True Positives (TP) and False Positives (FP)
    # Check if detection matches its original source
    tp_count = np.sum(catalog_table["halo_id"] == catalog_table["parent_id"])
    fp_count = n_total_detections - tp_count

    tp_ratio = tp_count / n_total_detections if n_total_detections > 0 else 0
    fp_ratio = fp_count / n_total_detections if n_total_detections > 0 else 0

    # 4. True Negatives (TN)
    # Correctly missed: in mock_undetected but NOT in catalog_table
    tn_count = np.sum(~np.isin(mock_undetected_table[id_col], catalog_table[id_col]))
    tn_ratio = tn_count / n_non_detections if n_non_detections > 0 else 0

    # 5. False Negatives (FN)
    # Incorrectly missed: in mock_detected but NOT in catalog_table
    fn_count = np.sum(~np.isin(mock_detected_table[id_col], catalog_table[id_col]))
    fn_ratio = fn_count / n_non_detections if n_non_detections > 0 else 0

    # 6. Additional Stats
    purity_or_completeness = (
        n_total_detections / len(mock_table) if len(mock_table) > 0 else 0
    )

    return {
        "counts": {"TP": tp_count, "FP": fp_count, "TN": tn_count, "FN": fn_count},
        "ratios": {"TP": tp_ratio, "FP": fp_ratio, "TN": tn_ratio, "FN": fn_ratio},
        "sums": {"positive": tp_ratio + fp_ratio, "negative": tn_ratio + fn_ratio},
        "overall_ratio": purity_or_completeness,
    }
