import argparse
import os
import sys
import numpy as np
import snputils as su

from .transform import transform_input
from .inference import run_inference

def main():
    '''
    qsnpobj: SNPObject instance of the query.
    '''
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
    
    transformed_data = transform_input(qsnpobj)
    print("Transformación completada. Matriz transformada:", transformed_data.shape)

    MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "model.onnx")
    predictions = run_inference(transformed_data, MODEL_PATH)
    
    with open(args.output, "w") as f:
        f.write(str(predictions))
    print(f"Resultados guardados en: {args.output}")

if __name__ == '__main__':
    main()