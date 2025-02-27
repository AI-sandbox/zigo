import argparse
import os
import sys
import time
import numpy as np
import snputils as su

from .transform import transform_input1, transform_input2, transform_input3
from .inference import run_inference

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Ruta al archivo de entrada genómico (.vcf, .vcf.gz, .pgen, .bed, etc.)"
    )
    parser.add_argument(
        "-o", "--output",
        default="results.txt",
        help="Ruta del archivo de salida para guardar los resultados"
    )
    args = parser.parse_args()

    print(f"Leyendo archivo de entrada: {args.input}")
    try:
        qsnpobj = su.read_snp(args.input)
    except Exception as e:
        print(f"[ERROR] No se pudo leer el archivo {args.input}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Se han leído {qsnpobj.n_snps} SNPs y {qsnpobj.n_samples} muestras.")
    print("Ejemplos de IDs de SNP:", qsnpobj.variants_id[:5])
    print("Ejemplos de alelos REF:", qsnpobj.variants_ref[:5])
    print("Forma de la matriz de genotipos (calldata_gt):", qsnpobj.calldata_gt.shape)
    
    # Ejecutar transform_input1 y medir tiempo
    t0 = time.time()
    output1 = transform_input1(qsnpobj)
    t1 = time.time() - t0
    print(f"transform_input1: {t1:.4f} segundos. Resultado: {output1.shape}")
    
    # Ejecutar transform_input2 y medir tiempo
    t0 = time.time()
    output2 = transform_input2(qsnpobj)
    t2 = time.time() - t0
    print(f"transform_input2: {t2:.4f} segundos. Resultado: {output2.shape}")
    
    # Ejecutar transform_input3 y medir tiempo
    t0 = time.time()
    output3 = transform_input3(qsnpobj)
    t3 = time.time() - t0
    print(f"transform_input3: {t3:.4f} segundos. Resultado: {output3.shape}")
    
    # Comparar los resultados: deben ser iguales (o muy cercanos)
    if np.allclose(output1, output2, equal_nan=True):
        print("Output1 y Output2 son iguales.")
    else:
        print("¡Atención! Output1 y Output2 difieren.")
    
    if np.allclose(output1, output3, equal_nan=True):
        print("Output1 y Output3 son iguales.")
    else:
        print("¡Atención! Output1 y Output3 difieren.")
    
    # Utilizar, por ejemplo, output1 para la inferencia (deben ser equivalentes)
    MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "model.onnx")
    predictions = run_inference(output1, MODEL_PATH)
    
    with open(args.output, "w") as f:
        f.write(str(predictions))
    print(f"Resultados guardados en: {args.output}")

if __name__ == '__main__':
    main()
