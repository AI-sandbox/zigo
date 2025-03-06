import onnxruntime as ort
import numpy as np

def format_sexcheck_report(outputs, sample_ids=None, pedsex=None):
    """
    Formatea los outputs de la inferencia en un reporte al estilo .sexcheck de PLINK2.
    
    Se esperan:
      - outputs[0]: Array con el código imputado de sexo para cada muestra (0: macho, 1: hembra),
                    que se transformará a (1: macho, 2: hembra).
      - outputs[1] (opcional): Se espera que contenga información de probabilidad, pero en este caso
                    no se utiliza para calcular F y se asigna "NA".
    
    Para cada muestra se generan las siguientes columnas:
      FID: se asigna igual que IID.
      IID: identificador individual (sample id).
      SID: se asigna igual que IID.
      PEDSEX: código de sexo en el input (1, 2 o NA); si no se provee se asigna "NA".
      SNPSEX: código de sexo imputado (1 para macho y 2 para hembra, derivado de outputs[0]).
      STATUS: "OK" si PEDSEX es distinto de "NA" y coincide con SNPSEX; "PROBLEM" en otro caso.
      F: se asigna "NA" ya que no se dispone de este valor.
      YCOUNT, YRATE, YOBS: se asignan "NA" ya que solo se usa chrX.
    
    Si no se pasan sample_ids se generan nombres de muestra como "sample1", "sample2", ...
    """
    # Obtener el array de predicciones de sexo
    sex_pred = outputs[0]
    n_samples = sex_pred.shape[0]
    
    # Si no se proporcionan sample_ids, generar unos por defecto
    if sample_ids is None:
        sample_ids = [f"sample{i+1}" for i in range(n_samples)]
    
    # Si no se provee PEDSEX, asignar "NA" a todas
    if pedsex is None:
        pedsex_list = ["NA"] * n_samples
    else:
        pedsex_list = pedsex  # se asume que es una lista de longitud n_samples
    
    # En este ejemplo, no disponemos de un valor xf (inbreeding coefficient), se asigna "NA"
    xf_arr = np.full((n_samples,), "NA")
    
    # Las columnas de chrY serán todas "NA"
    ycount = ["NA"] * n_samples
    yrate  = ["NA"] * n_samples
    yobs   = ["NA"] * n_samples
    
    # Formatear SNPSEX: se transforman los valores 0->1 y 1->2
    snpsex_list = []
    status_list = []
    for i in range(n_samples):
        try:
            # Se obtiene la predicción (0 o 1) y se transforma: 0 -> 1 (macho), 1 -> 2 (hembra)
            sp = str(int(round(sex_pred[i])) + 1)
        except Exception:
            sp = "NA"
        snpsex_list.append(sp)
        if pedsex_list[i] != "NA" and pedsex_list[i] == sp:
            status_list.append("OK")
        elif pedsex_list[i] == "NA":
            status_list.append("NA")
        else:
            status_list.append("PROBLEM")
    
    # Asignar FID y SID (se usa el mismo sample_id)
    fid_list = sample_ids
    sid_list = sample_ids
    
    # Construir la cabecera y cada línea del reporte (columnas separadas por tabulador)
    header = "\t".join(["FID", "IID", "SID", "PEDSEX", "SNPSEX", "STATUS", "F", "YCOUNT", "YRATE", "YOBS"])
    lines = [header]
    for i in range(n_samples):
        line = "\t".join([
            fid_list[i],
            sample_ids[i],
            sid_list[i],
            pedsex_list[i],
            snpsex_list[i],
            status_list[i],
            str(xf_arr[i]),
            ycount[i],
            yrate[i],
            yobs[i]
        ])
        lines.append(line)
    report = "\n".join(lines)
    return report

def run_inference(data, model_path, sample_ids=None, pedsex=None):
    """
    Ejecuta la inferencia usando ONNX Runtime y, a continuación,
    formatea los outputs en un reporte al estilo .sexcheck.
    
    Parámetros:
      - data: matriz de datos de entrada.
      - model_path: ruta al modelo ONNX.
      - sample_ids (opcional): lista de identificadores de muestra.
      - pedsex (opcional): lista de códigos de sexo provenientes del input (1/2/NA).
      
    Retorna:
      - Un string con el reporte formateado.
    """
    session = ort.InferenceSession(model_path)
    
    input_name = session.get_inputs()[0].name
    data = data.astype(np.float32)
    
    outputs = session.run(None, {input_name: data})
    
    report = format_sexcheck_report(outputs, sample_ids=sample_ids, pedsex=pedsex)
    return report
