# Examples - Shell Scripts Usage

This `scripts/` directory contains example shell scripts demonstrating various GREAC workflows. These scripts are provided as **templates and examples only**. Users should:

- **Review and modify** all paths, parameters, and configurations before execution
- **Test with small datasets** before running on production data
- **Understand each command** and its implications for your specific use case
- **Backup your data** before running any automated workflows

## Scripts 

- [`doall.sh`](./local/doall.sh): provides a **comprehensive end-to-end example** of the entire GREAC workflow, including data preprocessing, feature extraction, parameter optimization, and benchmarking. 

### ⚠️ Important Disclaimers

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

- [`benchmark.sh`](./local/benchmark.sh): provides a **complete execution** of GREAC, benchmarking datasets. 
- [`extract-features.sh`](./local/extract-features.sh): extraction of features and model fitting (regions and frequency behavior). 
- [`fasta-regions.sh`](./local/fasta-regions.sh): create reduced FASTA files with the extract regions. 
- [`file-classification.sh`](./local/file-classification.sh): classify FASTA sequences files based on pre-trained model by `extracted-features`. 