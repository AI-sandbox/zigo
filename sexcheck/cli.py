import argparse
import os
import sys
import time
import logging
import json
import pandas as pd

from sexcheck.reader import read_vcf
from sexcheck.inference import run_inference

MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "zigo.json")

def setup_logging(output_dir: str) -> logging.Logger:
    """
    Configures the logging system.
    
    Parameters:
        output_dir: Directory where log file will be saved
    
    Returns:
        Configured logger instance
    """
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

def process_ped(ped_path: str, logger: logging.Logger) -> dict:
    """
    Processes the PED file, ensuring required columns exist.
    
    Parameters:
        ped_path: Path to the PED file
        logger: Logger instance for logging messages
    
    Returns:
        Dictionary mapping sample IDs to sex labels (1 or 2)
    """
    logger.info(f"Reading PED file: {ped_path}")
    try:
        ped_df = pd.read_csv(ped_path, sep="\t")
        required_columns = ["Individual ID", "Gender"]
        
        if not all(col in ped_df.columns for col in required_columns):
            logger.error("The PED file must contain 'Individual ID' and 'Gender' columns.")
            sys.exit(1)
        
        ped_dict = ped_df.set_index("Individual ID")["Gender"].to_dict()
        
        return {str(k): str(v) for k, v in ped_dict.items() if v in [1, 2]}
    
    except Exception as e:
        logger.error(f"Error processing PED file: {e}")
        sys.exit(1)

def main() -> None:
    """
    Main entry point for the sex-check command-line tool.
    """
    parser = argparse.ArgumentParser(
        description="Sex inference from genomic data using zygosity distributions"
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the genomic input file (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Directory for saving results"
    )
    parser.add_argument(
        "--ped", 
        required=False, 
        help="Path to the PED file containing 'Individual ID' and 'Gender'"
    )
    
    args = parser.parse_args()

    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)
    logger = setup_logging(output_dir)

    # Use fast C-based VCF processor with normalization
    input_file = args.input
    start_time = time.time()
    
    try:
        transformed_data, vcf_data, stats = read_vcf(input_file)
    except Exception as e:
        logger.error(f"Failed to process VCF: {e}")
        sys.exit(1)
    
    logger.info(f"File processed in {time.time() - start_time:.4f} seconds.")
    logger.info(f"Loaded {stats['num_snps_input']} SNPs and {stats['num_samples']} samples.")
    
    sample_ids = vcf_data.samples
    ped_sex = process_ped(args.ped, logger) if args.ped else None

    start_time = time.time()

    sexcheck, nosex = run_inference(data=transformed_data, model_path=MODEL_PATH, sample_ids=sample_ids, pedsex=ped_sex)
    logger.info(f"Inference completed in {time.time() - start_time:.4f} seconds.")

    sexcheck_file = os.path.join(output_dir, "results.sexcheck")
    with open(sexcheck_file, "w") as f:
        f.write(str(sexcheck))
    logger.info(f"Results saved to: {sexcheck_file}")

    nosex_file = os.path.join(output_dir, "results.nosex")
    with open(nosex_file, "w") as f:
        f.write(nosex)
    logger.info(f"No sex information saved to: {nosex_file}")

if __name__ == '__main__':
    main()
