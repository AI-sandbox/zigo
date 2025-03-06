import argparse
import os
import sys
import time
import logging
import json
from typing import List, Optional
import pandas as pd
import snputils as su

from .transform import transform_input
from .inference import run_inference

def setup_logging(output_dir: str) -> logging.Logger:
    """Configures the logging system."""
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

def process_ped(ped_path: str, samples: List[str], logger: logging.Logger) -> List[str]:
    """Processes the PED file, ensuring required columns exist and reordering samples."""
    logger.info(f"Reading PED file: {ped_path}")
    try:
        ped_df = pd.read_csv(ped_path, sep="\t")
        required_columns = ["Individual ID", "Gender"]
        
        if not all(col in ped_df.columns for col in required_columns):
            logger.error("The PED file must contain 'Individual ID' and 'Gender' columns.")
            sys.exit(1)
        
        ped_df = ped_df.set_index("Individual ID").reindex(samples)
        
        if ped_df.isnull().any().any():
            logger.warning("Some individuals in the PED file were not found in the input samples.")
        
        return ped_df["Gender"].tolist()
    except Exception as e:
        logger.error(f"Error processing PED file: {e}")
        sys.exit(1)

def check_rs_format(variants_id: Optional[List[str]]) -> float:
    """Calculates the percentage of SNP IDs that follow the rsXXXX format."""
    if variants_id is None or len(variants_id) == 0:
        return 0.0
    
    count_rs = sum(vid.startswith("rs") and vid[2:].isdigit() for vid in variants_id if isinstance(vid, str))
    return (count_rs / len(variants_id)) * 100

def main() -> None:
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the genomic input file (.vcf, .pgen, etc.)"
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

    # logger.info(f"Reading input file: {args.input}")
    start_time = time.time()
    try:
        snp_obj = su.read_snp(args.input)
    except Exception as e:
        logger.error(f"Failed to read {args.input}: {e}")
        sys.exit(1)
    logger.info(f"File loaded in {time.time() - start_time:.4f} seconds.")
    logger.info(f"Loaded {snp_obj.n_snps} SNPs and {snp_obj.n_samples} samples.")

    ped_sex = process_ped(args.ped, snp_obj.samples, logger) if args.ped else None

    rs_percentage = check_rs_format(snp_obj.variants_id)
    logger.info(f"{rs_percentage:.2f}% of SNP IDs follow the rsXXXX format.")
    if rs_percentage < 50:
        logger.warning("Less than 50% of SNPs have rsXXXX identifiers. Check data quality.")

    start_time = time.time()
    transformed_data, stats = transform_input(snp_obj)
    logger.info(f"Input transformation completed in {time.time() - start_time:.4f} seconds.")
    logger.info(f"Using {stats['num_snps_used']} out of {stats['num_snps_input']} SNPs ({stats['overlap_percentage']:.2f}% overlap).")
    logger.info(f"Initial null percentage: {stats['initial_null_percentage']:.2f}%, Final null percentage: {stats['final_null_percentage']:.2f}%.")

    
    stats_path = os.path.join(output_dir, "info.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=4)
    logger.info(f"Data statistics saved to: {stats_path}")

    model_path = os.path.join(os.path.dirname(__file__), "models", "model.onnx")
    start_time = time.time()
    predictions = run_inference(data=transformed_data, model_path=model_path, sample_ids=snp_obj.samples, pedsex=ped_sex)
    logger.info(f"Inference completed in {time.time() - start_time:.4f} seconds.")

    output_file = os.path.join(output_dir, "results.sexcheck")
    with open(output_file, "w") as f:
        f.write(str(predictions))
    logger.info(f"Results saved to: {output_file}")

if __name__ == '__main__':
    main()
