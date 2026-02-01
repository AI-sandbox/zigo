# Zigo: Sex Checking by Zigosity Distributions

A command-line tool for sex inference from genomic data. This tool uses a distilled polynomial equation to predict genetic sex from SNP zygosity distributions and can compare predictions with provided PED files.

## Installation

Installation instructions will be added soon

## Key Features

- **Genetic Sex Prediction**: Accurately predicts genetic sex from genomic data using a machine learning model
- **PED File Integration**: Optional comparison with provided sex information from PED files (requires "Individual ID" and "Gender" columns)
- **Comprehensive Logging**: Detailed logs of the analysis process

## Usage

```bash:sex-check/README.md
sex-check -i INPUT_FILE -o OUTPUT_DIR [--ped PED_FILE]
```

### Arguments

- `-i, --input`: Path to the genomic input file (.vcf, .pgen, etc.)
- `-o, --output`: Directory for saving results
- `--ped`: (Optional) Path to the PED file containing 'Individual ID' and 'Gender' columns

### Input Format Requirements

- **Genomic Data**: `.vcf` or `.vcf.gz` 
- **PED File**: Tab-separated file with at least two columns:
  - 'Individual ID': Sample identifiers matching those in the genomic data
  - 'Gender': Sex information coded as 1 (male) or 2 (female)

### Output Files

- `results.sexcheck`: Main results file with sex predictions, containing the following columns:
  - `FID`: Family ID (same as IID in this implementation)
  - `IID`: Individual ID
  - `SID`: Sample ID (same as IID in this implementation)
  - `PEDSEX`: Sex as provided in the PED file (1=male, 2=female, NA=not available)
  - `SNPSEX`: Predicted sex from SNP data (1=male, 2=female)
  - `STATUS`: Comparison result ("OK" if PEDSEX matches SNPSEX, "PROBLEM" if they differ)
  - `F`: Probability score for the predicted sex label
  
  Columns related to Y chromosome have been skipped. 
- `results.nosex`: List of samples with no sex information
- `sexcheck.log`: Detailed log of the analysis process

## Example

```bash
sex-check -i samples.vcf.gz -o results --ped samples.ped
```

## Dependencies

- **pandas**: Data manipulation and file I/O
- **numpy**: Numerical computations
- **GCC compiler**: Required to compile the C-based VCF processor

## License

[License information]