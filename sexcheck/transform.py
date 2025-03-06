import os
import json
import numpy as np
from typing import Tuple
from snputils.snp import SNPObject

def load_reference_data() -> Tuple[np.ndarray, np.ndarray]:
    """Loads reference SNP data from a JSON file."""
    REF_PATH = os.path.join(os.path.dirname(__file__), "data", "reference.json")
    with open(REF_PATH, "r") as f:
        data = json.load(f)
    return np.array(data["ids"]), np.array(data["refs"])

def preprocess(snpobj: SNPObject) -> np.ndarray:
    X = snpobj.calldata_gt
    X = np.where((X < 0) | (X > 1), np.nan, X)
    X = X.sum(axis=2).T
    return X

def compute_null_percentage(data: np.ndarray) -> float:
    """Computes the percentage of null values in a given dataset."""
    return np.isnan(data).sum() / data.size * 100

def transform_input(snpobj: SNPObject) -> np.ndarray:
    """
    Maps input SNP data to a reference set, ensuring consistent ordering.
    
    - Builds a dictionary mapping SNP IDs to reference indices.
    - Iterates over input SNPs, aligning those present in the reference.
    """
    ref_ids, ref_refs = load_reference_data()
    
    X = preprocess(snpobj)  # (n_muestras, n_snps_input)
    
    input_ids = snpobj.variants_id   # (n_snps_input,)
    input_refs = snpobj.variants_ref  # (n_snps_input,)
    n_samples = X.shape[0]
    n_ref = len(ref_ids)
    
    # Build dictionary: id -> (idx in reference, ref)
    ref_dict = {rid: (i, rref) for i, (rid, rref) in enumerate(zip(ref_ids, ref_refs))}
    
    snps_used = 0
    X_mapped = np.full((n_samples, n_ref), np.nan, dtype=X.dtype)
    for idx_input, (snp_id, snp_ref) in enumerate(zip(input_ids, input_refs)):
        if snp_id in ref_dict:
            ref_index, ref_allele = ref_dict[snp_id]
            if snp_ref == ref_allele:
                X_mapped[:, ref_index] = X[:, idx_input]
                snps_used += 1
    
    initial_null = compute_null_percentage(X)
    final_null = compute_null_percentage(X_mapped)
    overlap_percentage = (snps_used / len(input_ids)) * 100
    
    stats = {
        "num_samples": n_samples,
        "num_snps_input": len(input_ids),
        "num_snps_used": snps_used,
        "overlap_percentage": overlap_percentage,
        "initial_null_percentage": initial_null,
        "final_null_percentage": final_null
    }
    
    return X_mapped, stats
