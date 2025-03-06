import onnxruntime as ort
import numpy as np
from typing import List, Optional, Tuple

def format_sexcheck_report(
    outputs: Tuple[np.ndarray], 
    sample_ids: Optional[List[str]] = None, 
    pedsex: Optional[List[str]] = None
) -> str:
    """
    Formats inference outputs into a PLINK2-style .sexcheck report.
    
    Outputs:
      - outputs[0]: Predicted sex codes (0: male, 1: female) converted to (1: male, 2: female).
      - outputs[1] (optional): Probability data, not used, assigned "NA".
    
    Columns in the report:
      - FID, IID, SID: Sample identifiers.
      - PEDSEX: Input sex code (1, 2, or NA); assigned "NA" if not provided.
      - SNPSEX: Imputed sex code (1: male, 2: female).
      - STATUS: "OK" if PEDSEX matches SNPSEX, "PROBLEM" otherwise.
      - F, YCOUNT, YRATE, YOBS: Assigned "NA" since only chrX is used.
    
    If sample_ids are not provided, they default to "sample1", "sample2", etc.
    """
    sex_pred = outputs[0]
    n_samples = sex_pred.shape[0]
    
    if sample_ids is None:
        sample_ids = [f"sample{i+1}" for i in range(n_samples)]
    pedsex_list = pedsex if pedsex != None else ["NA"] * n_samples
    
    snpsex_list = []
    status_list = []
    for i in range(n_samples):
        try:
            sp = str(int(round(sex_pred[i])) + 1)
        except Exception:
            sp = "NA"
        snpsex_list.append(sp)
        status_list.append("OK" if pedsex_list[i] != "NA" and pedsex_list[i] == sp else "PROBLEM")
    
    header = "\t".join(["FID", "IID", "SID", "PEDSEX", "SNPSEX", "STATUS", "F", "YCOUNT", "YRATE", "YOBS"])
    lines = [header] + [
        "\t".join([
            sample_ids[i], sample_ids[i], sample_ids[i], pedsex_list[i],
            snpsex_list[i], status_list[i], "NA", "NA", "NA", "NA"
        ])
        for i in range(n_samples)
    ]
    
    return "\n".join(lines)

def run_inference(
    data: np.ndarray, 
    model_path: str, 
    sample_ids: Optional[List[str]] = None, 
    pedsex: Optional[List[str]] = None
) -> str:
    """
    Runs inference using ONNX Runtime and formats the results into a .sexcheck report.
    
    Parameters:
      - data: Input data matrix.
      - model_path: Path to the ONNX model.
      - sample_ids: Optional list of sample identifiers.
      - pedsex: Optional list of input sex codes (1/2/NA).
    
    Returns:
      - A formatted string report.
    """
    session = ort.InferenceSession(model_path)
    input_name = session.get_inputs()[0].name
    data = data.astype(np.float32)
    outputs = session.run(None, {input_name: data})
    
    return format_sexcheck_report(outputs, sample_ids=sample_ids, pedsex=pedsex)
