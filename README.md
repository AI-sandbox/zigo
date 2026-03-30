# Zigo: Sex Checking by Zigosity Distributions

A command-line tool for sex inference from genomic data. This tool uses a distilled polynomial equation to predict genetic sex from SNP zygosity distributions and can compare predictions with provided PED files.

## Installation

We are working on publishing this package to PyPI. In the meantime, you can install the latest version directly from the repository using `pip`:

```bash
pip install git+https://github.com/AI-sandbox/zigo.git
```

## Key Features

- **Genetic Sex Prediction**: Accurately predicts genetic sex from genomic data using a machine learning model
- **PED File Integration**: Optional comparison with provided sex information from PED files (requires "Individual ID" and "Gender" columns)
- **Comprehensive Logging**: Detailed logs of the analysis process

## Usage

```bash:sex-check/README.md
zigo -i INPUT_FILE -o OUTPUT_DIR [--ped PED_FILE]
```

### Arguments

- `-i, --input`: Path to the genomic input file (.vcf, .vcf.gz, .bed, .pgen)
- `-o, --output`: Directory for saving results
- `--ped`: (Optional) Path to the PED file containing 'Individual ID' and 'Gender' columns

### Input Format Requirements

- **Genomic Data**: Supports `.bed`, `.pgen`, and includes a high-performance C-based reader for `.vcf` and `.vcf.gz`.
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

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Cite

If you use **Zigo** in your research, please cite [our paper](https://www.biorxiv.org/content/10.64898/2026.03.15.711924v1).

**BibTeX:**

```bibtex
@article {zigo2026,
	author = {Molina-Sedano, Oscar and Mas Montserrat, Daniel and Ioannidis, Alexander G.},
	title = {Sex checking by zygosity distributions},
	year = {2026},
	doi = {10.64898/2026.03.15.711924},
	URL = {https://www.biorxiv.org/content/early/2026/03/18/2026.03.15.711924},
	journal = {bioRxiv}
}
