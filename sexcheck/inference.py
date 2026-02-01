import numpy as np
import re
import json
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
        sp = str(int(round(sex_pred[i])) + 1) # 0 -> 1 = male, 1 -> 2 = female
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
            snpsex_list[i], status_list[i], str(probabilities[i][0])
        ])
        for i in range(n_samples)
    ]
    
    return "\n".join(lines_sexcheck), "\n".join(nosex_list)

def _clip01(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    return np.clip(p, eps, 1.0 - eps)

def _sigmoid(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-z))

_POLY_TOKEN_RE = re.compile(r"x(\d+)(?:\^(\d+))?")

def _eval_poly_feature(X: np.ndarray, name: str) -> np.ndarray:
    """
    Evalúa una feature tipo:
        "1", "x0", "x1", "x2", "x0^2", "x0 x1", "x1^2", "x2^2", ...
    en modo vectorizado sobre X (N,3).
    """
    name = name.strip()
    n = X.shape[0]
    if name == "1":
        return np.ones(n, dtype=np.float64)

    phi = np.ones(n, dtype=np.float64)
    for tok in name.split():
        m = _POLY_TOKEN_RE.fullmatch(tok)
        if m is None:
            raise ValueError(f"Token de feature no reconocido: {tok}")
        idx = int(m.group(1))
        pow_str = m.group(2)
        power = int(pow_str) if pow_str is not None else 1
        phi *= np.power(X[:, idx], power, dtype=np.float64)
    return phi

def _poly_logit_from_model(X: np.ndarray, model: dict) -> np.ndarray:
    """
    X: (N,3)
    model: dict con keys: intercept, coefs, feature_names
    devuelve z = intercept + sum_j coefs_j * phi_j(X)
    """
    intercept = float(model["intercept"])
    coefs = np.asarray(model["coefs"], dtype=np.float64)
    names = model["feature_names"]
    z = np.full(X.shape[0], intercept, dtype=np.float64)
    for c, name in zip(coefs, names):
        if abs(c) < 1e-12:
            continue
        phi = _eval_poly_feature(X, name)
        z += c * phi
    return z

def parse_poly_json(model_path: str) -> dict:
    """Carga modelo polinómico desde el JSON de distill_poly.py."""
    with open(model_path, "r") as f:
        d = json.load(f)
    required = {"intercept", "coefs", "feature_names"}
    if not required.issubset(d.keys()):
        raise ValueError("JSON no parece un modelo polinómico válido.")
    return {
        "intercept": d["intercept"],
        "coefs": d["coefs"],
        "feature_names": d["feature_names"],
    }

def run_poly_inference(data: np.ndarray, poly_model: dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Ejecuta inferencia con el modelo polinómico:
        z = intercept + sum coefs_j * phi_j(x)
        p_male = sigmoid(z)   # z viene de ajustar sobre p(male) del CatBoost
    Queremos que:
        probabilities[:,0] = P(male)
        probabilities[:,1] = P(female)
    para que la columna F siga siendo P(male), igual que en ONNX/symbolic.
    """
    X = data.astype(np.float64, copy=False)
    z = _poly_logit_from_model(X, poly_model)

    p_male = _sigmoid(z)
    p_male = _clip01(p_male)
    p_female = 1.0 - p_male

    # Haploids are always male
    is_haploid = (X[:, 0] == 0) & (X[:, 1] == 0)
    p_male[is_haploid] = 1.0
    p_female[is_haploid] = 0.0

    probabilities = np.column_stack([p_male, p_female])
    predictions = (p_female >= 0.5).astype(int)  # 0 = male (1 in PLINK), 1 = female (2 in PLINK)

    return predictions, probabilities

def run_inference(
    data: np.ndarray,
    model_path: str,
    sample_ids: Optional[List[str]] = None,
    pedsex: Optional[List[str]] = None,
) -> str:
    """
    Ejecuta inferencia (ONNX, simbólica o polinómica)
    """
    poly_model = parse_poly_json(model_path)
    x = data.astype(np.float64, copy=False)
    predictions, probabilities = run_poly_inference(x, poly_model)
    outputs = (predictions, probabilities)

    return format_sexcheck_report(outputs, sample_ids=sample_ids, pedsex=pedsex)
