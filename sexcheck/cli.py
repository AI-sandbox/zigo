import argparse
import os
import sys
import time
import numpy as np
import snputils as su

from .transform import transform_input
from .inference import run_inference

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Ruta al archivo de entrada genómico (.vcf, .vcf.gz, .pgen, .bed, etc.)"
    )
    parser.add_argument( # path
        "-o", "--output",
        default="results.txt",
        help="Ruta del archivo de salida para guardar los resultados"
    )
    args = parser.parse_args()

    print(f"Leyendo archivo de entrada: {args.input}")
    t0 = time.time()
    try:
        qsnpobj = su.read_snp(args.input)
    except Exception as e:
        print(f"[ERROR] No se pudo leer el archivo {args.input}: {e}", file=sys.stderr)
        sys.exit(1)
    t1 = time.time() - t0

    ## Warning / Sanity Check de ids rsXXX...

    print(f"Se han leído {qsnpobj.n_snps} SNPs y {qsnpobj.n_samples} muestras.")
    print("Ejemplos de IDs de SNP:", qsnpobj.variants_id[:5])
    print("Ejemplos de alelos REF:", qsnpobj.variants_ref[:5])
    print("Forma de la matriz de genotipos (calldata_gt):", qsnpobj.calldata_gt.shape)
    print(f"read_snp: {t1:.4f} segundos.")
    
    # Ejecutar transform_input1 y medir tiempo
    t0 = time.time()
    output1 = transform_input(qsnpobj)
    t1 = time.time() - t0
    print(f"transform_input: {t1:.4f} segundos. Resultado: {output1.shape}")
    
    # Utilizar, por ejemplo, output1 para la inferencia (deben ser equivalentes)
    MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "model.onnx")
    t0 = time.time()
    predictions = run_inference(output1, MODEL_PATH)
    t1 = time.time() - t0
    print(f"run_inference: {t1:.4f} segundos. Resultado: {len(predictions)}")
    
    with open(args.output, "w") as f:
        f.write(str(predictions))
    print(f"Resultados guardados en: {args.output}")

    # .log ## llibreria logging ### info warning debugging error
    # dicc json with stats (% overlap, nsnps, nan percentage, ...)
    # .sexcheck

if __name__ == '__main__':
    main()
