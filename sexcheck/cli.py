import argparse
import sys
import numpy as np
import snputils as su

from .transform import *
from .inference import *

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
    

if __name__ == '__main__':
    main()