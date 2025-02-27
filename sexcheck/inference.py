# inference.py
import onnxruntime as ort
import numpy as np

def run_inference(data, model_path):
    session = ort.InferenceSession(model_path)
    
    input_name = session.get_inputs()[0].name
    data = data.astype(np.float32)
    
    outputs = session.run(None, {input_name: data})
    return outputs