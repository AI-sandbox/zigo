#!/usr/bin/env bash
# -------------------------------------------------
# Script: plink-bench.sh
# Uso:    ./plink-bench.sh <vcf_directory> <output_directory> <path_to_ped_update>
#
# Ejemplo de uso:
#   ./plink-bench.sh /home/oscar/chr/data/eval-test/ \
#                    /home/oscar/chr/sex-check/benchmarks/plink1.9/ \
#                    /home/oscar/chr/data/raw/labels/sex_update.txt
#
# -------------------------------------------------

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <vcf_directory> <output_directory> <path_to_ped_update>"
  exit 1
fi

VCF_DIR=$1
OUTPUT_DIR=$2
PED_UPDATE=$3

if [ ! -d "$VCF_DIR" ]; then
  echo "Error: The VCF directory '$VCF_DIR' does not exist."
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Procesa cada archivo test-*.vcf en VCF_DIR
for vcf_file in "$VCF_DIR"/test-*.vcf; do
    if [ ! -e "$vcf_file" ]; then
        echo "No .vcf files found in '$VCF_DIR'."
        exit 1
    fi

    filename=$(basename "$vcf_file")
    # Extrae el id quitando el prefijo "test-" y la extensión ".vcf"
    id=${filename#test-}
    id=${id%.vcf}

    echo "Processing file: $vcf_file (ID: $id)"
    # Crea el subdirectorio para este test
    mkdir -p "$OUTPUT_DIR/$id"
    
    # Ejecuta PLINK 1.9 usando el archivo VCF y el archivo de actualización de sexo.
    /home/oscar/plink1.9/plink --vcf "$vcf_file" \
        --update-sex "$PED_UPDATE" \
        --check-sex \
        --out "$OUTPUT_DIR/$id/results"
done

echo "All files processed. Results stored in '$OUTPUT_DIR'."

# Activa el entorno virtual (si es necesario, según la configuración de tu sistema)
source ~/.bashrc
activatevenv sex

# Ejecuta el script Python para evaluar los resultados. Este script buscará
# cualquier archivo .sexcheck en los subdirectorios (0,1,2,...) y generará results.csv.
python3 analyze-results.py "$OUTPUT_DIR"
