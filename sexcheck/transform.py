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

def transform_input(snpobj:SNPObject) -> np.ndarray:
    """
    Transforma el objeto SNPObject:
      1. Preprocesa la matriz de genotipos para obtener X de forma (n_muestras, n_snps_input).
      2. Carga la referencia con la lista de SNP IDs y alelos REF (usada en entrenamiento).
      3. Crea una nueva matriz X_mapped de forma (n_muestras, n_snps_referencia) en la que:
           - Para cada SNP de referencia, se busca en el snpobj (por su id).
           - Si se encuentra y el alelo REF coincide, se copia la fila (la columna correspondiente de X).
           - Si no se encuentra o el alelo no coincide, se deja np.nan.
    """
    ref_ids, ref_refs = load_reference_data()
    
    # (n_samples, n_snps_input)
    X = preprocess(snpobj)
    
    input_ids = snpobj.variants_id  # Se asume array de forma (n_snps_input,)
    input_refs = snpobj.variants_ref  # Se asume array de forma (n_snps_input,)
    
    n_samples = X.shape[0]
    n_ref_snps = len(ref_ids)
    
    # Inicializar la matriz final con np.nan
    X_mapped = np.full((n_samples, n_ref_snps), np.nan)
    
    # Mapear cada SNP de referencia a la posición correspondiente en X
    for j, (ref_id, ref_ref) in enumerate(zip(ref_ids, ref_refs)):
        # Buscar el índice del SNP de entrada con el mismo ID
        indices = np.where(input_ids == ref_id)[0]
        if indices.size > 0:
            idx = indices[0]
            # Si el alelo REF coincide, copiar los genotipos
            if input_refs[idx] == ref_ref:
                X_mapped[:, j] = X[:, idx]
            else:
                # Si no coincide, dejar np.nan para este SNP en todas las muestras
                X_mapped[:, j] = np.nan
        # Si el SNP de referencia no se encuentra, se mantiene np.nan
    return X_mapped
