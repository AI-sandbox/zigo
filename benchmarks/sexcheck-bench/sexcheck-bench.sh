#!/usr/bin/env bash

# -------------------------------------------------
# Script: sexcheck-bench.sh
# Use:    ./sexcheck-bench.sh <vcf_directory> <output_directory> <path_to_ped>
# -------------------------------------------------

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <vcf_directory> <output_directory> <path_to_ped>"
  exit 1
fi

VCF_DIR=$1
OUTPUT_DIR=$2
PED_PATH=$3

if [ ! -d "$VCF_DIR" ]; then
  echo "Error: The VCF directory '$VCF_DIR' does not exist."
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

for vcf_file in "$VCF_DIR"/*.vcf
do
  if [ ! -e "$vcf_file" ]; then
    echo "No .vcf files found in '$VCF_DIR'."
    exit 1
  fi

  filename=$(basename "$vcf_file")
  id=${filename#val-}
  id=${id%.vcf}

  echo "Processing file: $vcf_file"
  sex-check -i "$vcf_file" -o "$OUTPUT_DIR/$id" --ped "$PED_PATH"
done

echo "Processing finished. Results can be found in '$OUTPUT_DIR'."

# alias made by oscar in his own .bashrc to activate a python virtual environment
source ~/.bashrc
activatevenv sex

python3 analyze-results.py "$OUTPUT_DIR"
