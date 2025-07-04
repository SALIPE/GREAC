# GREAC - Genomic Region Extraction and Classifier

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.6%2B-blue.svg)](https://julialang.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://doi.org/pending)

## Overview

GREAC (Genomic Region Extraction and Classifier) is a novel computational methodology for discriminative genomic region extraction and classification. This data-driven, non-parametric approach employs exclusive k-mers through strategically defined sliding windows to identify genomic regions with high concentrations of variant-specific signatures.

### Key Features

- **K-mer based analysis**: Utilizes exclusive k-mers for discriminative feature extraction
- **Sliding window approach**: Strategic window-based analysis for region identification
- **Non-parametric methodology**: Data-driven approach without strong distributional assumptions
- **High-performance**: Optimized Julia implementation for computational efficiency
- **Flexible classification**: Multiple distance metrics and evaluation tools
- **Benchmarking suite**: Comprehensive performance evaluation capabilities

## Usage

### Pre-compiled Binary

Download the latest pre-compiled binary from the [Releases](https://github.com/SALIPE/genomic-extractor/releases) page:

### External Dependencies (for complete workflow)

For the full workflow using `scripts/local/doall.sh`, install these external tools:

### [GRAMEP](https://github.com/omatheuspimenta/GRAMEP) - For Exclusive K-mers Search

### [FastasSplitter](https://github.com/SALIPE/Fasta-splitter) - For Dataset Balancing

## Quick Start

```bash
# Basic feature extraction
GREAC extract-features --group-name denv --window 0.1 --train-dir /path/to/training/data

# Run benchmark with classification
GREAC benchmark --group-name denv --window 0.1 --train-dir /path/to/train --test-dir /path/to/test

# Brute Force Parameters search
GREAC fit-parameters --group-name denv --window 0.05 --train-dir /path/to/train --test-dir /path/to/test
```

## Usage

### Command Structure

```bash
GREAC [GLOBAL OPTIONS] COMMAND [COMMAND OPTIONS]
```

### Global Options

| Option | Description | Required | Type |
|--------|-------------|----------|------|
| `--group-name` | Process group name for organizing results | Yes | String |
| `-w, --window` | Sliding window percent size (0.0-1.0) | Yes | Float |
| `--no-cache` | Remove cached files before processing | No | Flag |

### Commands

#### 1. Extract Features

Extract k-mer features from genomic sequences for downstream analysis.

```bash
GREAC --group-name denv --window 0.1 extract-features  --train-dir /data/training/
```

**Options:**
- `--train-dir`: Path to training dataset directory (required)

#### 2. Benchmark

Perform classification benchmark with confusion matrix generation.

```bash
GREAC --group-name denv --window 0.1 benchmark \
  --train-dir /data/training/ --test-dir /data/testing/ \
  --metric manhattan --threshold 0.05 --output-directory /results/
```

**Options:**
- `--train-dir`: Training dataset path (required)
- `--test-dir`: Test dataset path (required)
- `-m, --metric`: Distance metric (`manhattan`, `euclidian`, `chisquared`, `mahalanobis`, `kld`)
- `--threshold`: Window threshold for consideration (Float16)
- `-o, --output-directory`: Output directory for results

#### 3. Fit Parameters

Optimize model parameters using training and test datasets.

```bash
GREAC --group-name denv --window 0.05 fit-parameters \
  --train-dir /data/training/ --test-dir /data/testing/
```

**Options:**
- `--train-dir`: Training dataset path (required)
- `--test-dir`: Test dataset path (required)

#### 4. Export reduced FASTA files

Create the datasets FASTA files with cutted sequences on the extract regions.

```bash
GREAC --group-name denv --window 0.05 fasta-regions --input /data/training/
```

**Options:**
- `--input`: Training dataset path (required)


## Data Format

### Input Requirements

- **FASTA/FASTQ files**: Genomic sequences in standard format
- **Directory structure**: Organize sequences by class/group in separate directories
- **File naming**: Consistent naming convention for reproducibility

### Example Directory Structure

The exclusive k-mers files are the output from GRAMEP, and the *.fasta files have all the sequences from training/extraction.

```
training_data/
├── class_A/
│   ├── class_A_ExclusiveKmers.sav
│   ├── class_A_ExclusiveKmers.txt
│   └── class_A.fasta
├── class_B/
│   ├── class_B_ExclusiveKmers.sav
│   ├── class_B_ExclusiveKmers.txt
│   └── class_B.fasta
└── class_C/
│   ├── class_C_ExclusiveKmers.sav
│   ├── class_C_ExclusiveKmers.txt
    └── class_C.fasta
```

## Examples

### Complete Workflow Example Script

The `doall.sh` script includes:
- Data preprocessing with FastasSplitter
- K-mer extraction using GRAMEP
- Complete GREAC workflow execution
- Result organization and summary generation


## Performance Considerations

- **Memory usage**: Scales with k-mer size and sequence length
- **Processing time**: Linear with dataset size and window overlap
- **Disk space**: Caching improves speed but requires storage
- **Parallelization**: Multi-threading support for large datasets

## Troubleshooting

### Common Issues

**Memory errors with large datasets:**
```bash
# Use smaller window sizes or process in batches
GREAC benchmark --window 0.05 --no-cache [other options]
```

**File format errors:**
- Ensure FASTA files are properly formatted
- Check file permissions and paths
- Verify directory structure matches expected format

## Citation

If you use GREAC in your research, please cite:


## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request


### ⚠️ Important Disclaimers

#### Running Example Scripts
The `scripts/` directory contains example shell scripts demonstrating various GREAC workflows. These scripts are provided as **templates and examples only**. Users should:

- **Review and modify** all paths, parameters, and configurations before execution
- **Test with small datasets** before running on production data
- **Understand each command** and its implications for your specific use case
- **Backup your data** before running any automated workflows

#### Complete Workflow Script (`doall.sh`)
The `scripts/doall.sh` script provides a **comprehensive end-to-end example** of the entire GREAC workflow, including data preprocessing, feature extraction, parameter optimization, and benchmarking. 

**⚠️ Prerequisites Warning**: This script requires two external tools that must be installed separately:

1. **[GRAMEP](https://github.com/omatheuspimenta/GRAMEP)** - Tool for exclusive k-mers search
   - Used for k-mer extraction and analysis
   - Must be installed and accessible in your PATH
   - Follow GRAMEP installation instructions before using `doall.sh`

2. **[FastasSplitter](https://github.com/SALIPE/Fasta-splitter)** - Dataset balancing and structuring tool
   - Used to balance datasets and export structured files for GREAC execution
   - Required for proper data organization and preprocessing
   - Install according to FastasSplitter documentation

**The `doall.sh` script will NOT work without these dependencies properly installed and configured.**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Acknowledgments

- Julia community for the excellent computational framework
---

**Note**: This tool is under active development. Please check the [releases page](https://github.com/SALIPE/genomic-extractor/releases) for the latest stable version.