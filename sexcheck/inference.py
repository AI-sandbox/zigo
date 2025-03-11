import onnxruntime as ort
import numpy as np
from typing import List, Optional, Tuple

def format_sexcheck_report(
    outputs: Tuple[np.ndarray], 
    sample_ids: Optional[List[str]] = None, 
    pedsex: Optional[dict] = None
) -> Tuple[str, str]:
    """
    Formats inference outputs into a PLINK2-style .sexcheck and .nosex report.
    """
    sex_pred = outputs[0]
    probabilities = outputs[1]
    n_samples = sex_pred.shape[0]
    
    if sample_ids is None:
        sample_ids = [f"sample{i+1}" for i in range(n_samples)]
    
    snpsex_list = []
    status_list = []
    pedsex_list = []
    nosex_list  = []

    for i in range(n_samples):
        sp = str(int(round(sex_pred[i])) + 1)
        snpsex_list.append(sp)
        
        ped_sex = pedsex.get(sample_ids[i], "NA")
        pedsex_list.append(ped_sex)
        
        if ped_sex == "NA" or ped_sex not in ["1", "2"]:
            nosex_list.append(f"{sample_ids[i]}\t{sample_ids[i]}")
        
        status_list.append("OK" if ped_sex == sp else "PROBLEM")
    
    header_sexcheck = "\t".join(["FID", "IID", "SID", "PEDSEX", "SNPSEX", "STATUS", "F"])
    lines_sexcheck = [header_sexcheck] + [
        "\t".join([
            sample_ids[i], sample_ids[i], sample_ids[i], pedsex_list[i],
            snpsex_list[i], status_list[i], str(probabilities[i][int(snpsex_list[i]) - 1])
        ])
        for i in range(n_samples)
    ]
    
    return "\n".join(lines_sexcheck), "\n".join(nosex_list)

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
