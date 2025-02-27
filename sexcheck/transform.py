import os
import json
import numpy as np
from snputils.snp import SNPObject

def load_reference_data():
    REF_PATH = os.path.join(os.path.dirname(__file__), "data", "reference.json")
    with open(REF_PATH, "r") as f:
        data = json.load(f)
    return np.array(data["ids"]), np.array(data["refs"])

def preprocess(snpobj: SNPObject) -> np.ndarray:
    X = snpobj.calldata_gt
    X = np.where(X == -9, np.nan, X)
    X = X.sum(axis=2).T
    return X

def transform_input(snpobj: SNPObject) -> np.ndarray:
    """
    Versión optimizada usando un diccionario de la referencia.
    Se construye un diccionario (id -> (índice, ref)) a partir de la referencia y se itera
    sobre los SNPs del input (que se asume son muchos menos que en la referencia).
    """
    ref_ids, ref_refs = load_reference_data()
    
    X = preprocess(snpobj)  # (n_muestras, n_snps_input)
    
    input_ids = snpobj.variants_id   # (n_snps_input,)
    input_refs = snpobj.variants_ref  # (n_snps_input,)
    
    n_samples = X.shape[0]
    n_ref = len(ref_ids)
    
    # Construir diccionario: id -> (índice en referencia, ref)
    ref_dict = {rid: (i, rref) for i, (rid, rref) in enumerate(zip(ref_ids, ref_refs))}
    
    X_mapped = np.full((n_samples, n_ref), np.nan, dtype=X.dtype)
    
    # Iterar sobre el input (menos SNPs)
    for idx_input, (snp_id, snp_ref) in enumerate(zip(input_ids, input_refs)):
        if snp_id in ref_dict:
            ref_index, ref_allele = ref_dict[snp_id]
            if snp_ref == ref_allele:
                X_mapped[:, ref_index] = X[:, idx_input]
    return X_mapped
