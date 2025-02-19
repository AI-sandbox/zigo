import os
import json
import numpy as np

def load_reference_data():
    REF_PATH = os.path.join(os.path.dirname(__file__), "data", "reference.json")
    with open(REF_PATH, "r") as f:
        data = json.load(f)
    return np.array(data["ids"]), np.array(data["refs"])

def transform_input(qsnpobj):
    """
    Dada la instancia qsnpobj leída con snputils.read_snp,
    genera una matriz ordenada de genotipos correspondiente únicamente a los SNPs
    presentes en la referencia y valida la coincidencia del alelo REF.
    """
    ref_ids, ref_refs = load_reference_data()
    n_samples = qsnpobj.n_samples

    transformed = np.full((len(ref_ids), n_samples), np.nan)

    for i, (rid, rref) in enumerate(zip(ref_ids, ref_refs)):
        indices = np.where(qsnpobj.variants_id == rid)[0]
        if indices.size > 0:
            idx = indices[0]
            if qsnpobj.variants_ref[idx] == rref:
                transformed[i, :] = qsnpobj.calldata_gt[idx, :]
    return transformed
