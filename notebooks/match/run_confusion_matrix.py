from mock_generator import MockGenerator
from numcosmo_py import Nc, Ncm, sky_match
import numpy as np
from astropy.table import unique, Table

def calculate_split_metrics(table, cat_type="halo"):
    """
    Calculates confusion matrix metrics by splitting a table into 
    detected (has distance data) and undetected (empty distance data) sets.
    """
    # 1. Configuration based on category type
    if cat_type == "halo":
        id_col = 'halo_id'
        # Ground truth: 1 is a real object that should be found, 0 is not
        true_pos_condition = 1 
        true_neg_condition = 0
        gt_col = 'is_detected'
    elif cat_type == "cluster":
        id_col = 'cluster_id'
        # Ground truth: non-zero parent_id should be found, 0 should not
        true_pos_condition = 0 # used for '!= 0' logic below
        gt_col = 'parent_id'
    else:
        raise ValueError(f"Invalid cat_type: '{cat_type}'. Use 'halo' or 'cluster'.")

    # 2. Split the table by algorithm performance (Did we find a match?)
    # mask_has_data: rows where the algorithm matched a candidate
    mask_has_data = np.array([len(d) > 0 for d in table['RA_matched']])
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
    tp_mask = np.isin(gt_positives[id_col], detected_table[id_col])
    tp_count = np.sum(tp_mask)
    
    # FP: "Ghost" objects (GT says 0) found in the detected table
    fp_mask = np.isin(gt_negatives[id_col], detected_table[id_col])
    fp_count = np.sum(fp_mask)

    # TN: "Ghost" objects correctly left in the undetected table
    tn_mask = np.isin(gt_negatives[id_col], undetected_table[id_col])
    tn_count = np.sum(tn_mask)

    # FN: Real objects incorrectly left in the undetected table
    fn_mask = np.isin(gt_positives[id_col], undetected_table[id_col])
    fn_count = np.sum(fn_mask)

    # 5. Ratios
    tp_ratio = tp_count / n_detected
    fp_ratio = fp_count / n_detected
    tn_ratio = tn_count / n_undetected
    fn_ratio = fn_count / n_undetected

    return {
        "counts": {"TP": tp_count, "FP": fp_count, "TN": tn_count, "FN": fn_count},
        "ratios": {"TP": tp_ratio, "FP": fp_ratio, "TN": tn_ratio, "FN": fn_ratio},
        "sums": {"positive": tp_ratio + fp_ratio, "negative": tn_ratio + fn_ratio}
    }



def calculate_catalog_metrics(catalog_table, mock_table, cat_type="halo"):
    """
    Calculates TP, FP, TN, and FN for either 'halo' or 'cluster' catalogs.
    
    Parameters:
    catalog_table: The table of unique detections (unique_halos or unique_detections).
    mock_table: The ground truth catalog (halos_mock or clusters_mock).
    cat_type: String, either 'halo' or 'cluster'.
    """
    
    # 1. Handle cat_type specific logic
    if cat_type == "halo":
        # Halos use 'is_detected' 0 or 1
        undetected_mask = mock_table['is_detected'] == 0
        detected_mask = mock_table['is_detected'] == 1
        id_col = 'halo_id'
    elif cat_type == "cluster":
        # Clusters use 'parent_id' 0 or not 0
        undetected_mask = mock_table['parent_id'] == 0
        detected_mask = mock_table['parent_id'] != 0
        id_col = 'cluster_id'
    else:
        raise ValueError(f"Invalid cat_type: '{cat_type}'. Expected 'halo' or 'cluster'.")

    # Filter mock tables based on the logic above
    mock_undetected_table = mock_table[undetected_mask]
    mock_detected_table = mock_table[detected_mask]

    # 2. Total Counts for Denominators
    n_total_detections = len(catalog_table[id_col])
    # n_non_detections is the difference between found matches and the ground truth catalog
    n_non_detections = abs(len(catalog_table[id_col]) - len(mock_table[id_col]))

    # 3. True Positives (TP) and False Positives (FP)
    # Check if detection matches its original source
    tp_mask = catalog_table['halo_id'] == catalog_table['parent_id']
    tp_count = np.sum(tp_mask)
    fp_count = n_total_detections - tp_count
    
    tp_ratio = tp_count / n_total_detections if n_total_detections > 0 else 0
    fp_ratio = fp_count / n_total_detections if n_total_detections > 0 else 0

    # 4. True Negatives (TN)
    # Correctly missed: in mock_undetected but NOT in catalog_table
    tn_mask = ~np.isin(mock_undetected_table[id_col], catalog_table[id_col])
    tn_count = np.sum(tn_mask)
    tn_ratio = tn_count / n_non_detections if n_non_detections > 0 else 0

    # 5. False Negatives (FN)
    # Incorrectly missed: in mock_detected but NOT in catalog_table
    fn_mask = ~np.isin(mock_detected_table[id_col], catalog_table[id_col])
    fn_count = np.sum(fn_mask)
    fn_ratio = fn_count / n_non_detections if n_non_detections > 0 else 0

    # 6. Additional Stats
    purity_or_completeness = n_total_detections / len(mock_table) if len(mock_table) > 0 else 0

    return {
        "counts": {"TP": tp_count, "FP": fp_count, "TN": tn_count, "FN": fn_count},
        "ratios": {"TP": tp_ratio, "FP": fp_ratio, "TN": tn_ratio, "FN": fn_ratio},
        "sums": {"positive": tp_ratio + fp_ratio, "negative": tn_ratio + fn_ratio},
        "overall_ratio": purity_or_completeness
    }