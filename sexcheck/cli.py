import argparse
import os
import sys
import time
import numpy as np
import logging
import snputils as su

from .transform import transform_input
from .inference import run_inference

def setup_logging(output_dir):
    """Configura el sistema de logging."""
    log_file = os.path.join(output_dir, "sexcheck.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def check_rs_format(variants_id):
    """Verifica cuántos IDs de SNPs tienen formato rsXXXX."""
    if variants_id is None or len(variants_id) == 0:
        return 0.0

    count_rs = sum(1 for vid in variants_id if isinstance(vid, str) and vid.startswith("rs") and vid[2:].isdigit())
    percentage = (count_rs / len(variants_id)) * 100
    return percentage

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Ruta al archivo de entrada genómico (.vcf, .vcf.gz, .pgen, .bed, etc.)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Directorio donde se guardarán los resultados"
    )
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)
    logger = setup_logging(output_dir)

    # Load
    logger.info(f"Leyendo archivo de entrada: {args.input}")
    t0 = time.time()
    try:
        qsnpobj = su.read_snp(args.input)
    except Exception as e:
        logger.error(f"No se pudo leer el archivo {args.input}: {e}")
        sys.exit(1)
    t1 = time.time() - t0
    logger.info(f"Archivo leído en {t1:.4f} segundos.")
    logger.info(f"Se han leído {qsnpobj.n_snps} SNPs y {qsnpobj.n_samples} muestras.")

    # rsID: Sanity Check
    rs_percentage = check_rs_format(qsnpobj.variants_id)
    logger.info(f"Se han detectado un {rs_percentage:.2f}% de IDs en formato rs.")
    if rs_percentage < 50:
        logger.warning("Menos del 50% de los SNPs tienen un identificador en formato rsXXXX. "
                       "Revisar la calidad de los datos de entrada.")

    # Transform
    t0 = time.time()
    output1 = transform_input(qsnpobj)
    t1 = time.time() - t0
    logger.info(f"Transformación de input completada en {t1:.4f} segundos. Resultado: {output1.shape}")

    # Inference
    MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "model.onnx")
    t0 = time.time()
    predictions = run_inference(output1, MODEL_PATH, sample_ids=qsnpobj.samples)
    t1 = time.time() - t0
    logger.info(f"Inferencia completada en {t1:.4f} segundos. Resultado: {len(predictions)}")

    # Save
    output_file = os.path.join(output_dir, "results.sexcheck")
    with open(output_file, "w") as f:
        f.write(str(predictions))
    logger.info(f"Resultados guardados en: {output_file}")

if __name__ == '__main__':
    main()
