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

def transform_input1(snpobj: SNPObject) -> np.ndarray:
    """
    Versión original no optimizada.
    Para cada SNP de referencia, se busca en el snpobj y se asigna la columna correspondiente si el REF coincide.
    """
    ref_ids, ref_refs = load_reference_data()
    
    # (n_muestras, n_snps_input)
    X = preprocess(snpobj)
    
    input_ids = snpobj.variants_id  # (n_snps_input,)
    input_refs = snpobj.variants_ref  # (n_snps_input,)
    
    n_samples = X.shape[0]
    n_ref_snps = len(ref_ids)
    
    X_mapped = np.full((n_samples, n_ref_snps), np.nan)
    
    for j, (ref_id, ref_ref) in enumerate(zip(ref_ids, ref_refs)):
        indices = np.where(input_ids == ref_id)[0]
        if indices.size > 0:
            idx = indices[0]
            if input_refs[idx] == ref_ref:
                X_mapped[:, j] = X[:, idx]
            else:
                X_mapped[:, j] = np.nan
    return X_mapped

def transform_input2(snpobj: SNPObject) -> np.ndarray:
    """
    Versión optimizada usando np.argsort y np.searchsorted.
    Se ordena el array de IDs de entrada para buscar vectorizadamente los SNPs de referencia.
    """
    ref_ids, ref_refs = load_reference_data()
    
    X = preprocess(snpobj)  # (n_muestras, n_snps_input)
    
    input_ids = snpobj.variants_id   # (n_snps_input,)
    input_refs = snpobj.variants_ref  # (n_snps_input,)
    
    # Ordenar input_ids para usar searchsorted
    sort_order = np.argsort(input_ids)
    sorted_ids = input_ids[sort_order]
    
    # Buscar cada ref_id en el array ordenado
    pos = np.searchsorted(sorted_ids, ref_ids)
    # Verificar que se encuentre la coincidencia exacta
    found = (pos < len(sorted_ids)) & (sorted_ids[pos] == ref_ids)
    
    mapped_idx = np.full(ref_ids.shape, -1, dtype=int)
    mapped_idx[found] = sort_order[pos[found]]
    
    # Verificar que, para los SNP encontrados, el alelo REF coincida
    allele_match = (mapped_idx != -1) & (input_refs[mapped_idx] == ref_refs)
    
    n_samples = X.shape[0]
    n_ref_snps = len(ref_ids)
    X_mapped = np.full((n_samples, n_ref_snps), np.nan, dtype=X.dtype)
    
    # Asignar vectorizadamente las columnas donde el alelo coincide
    X_mapped[:, allele_match] = X[:, mapped_idx[allele_match]]
    
    return X_mapped

def transform_input3(snpobj: SNPObject) -> np.ndarray:
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
