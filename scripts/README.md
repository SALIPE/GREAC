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

- [`clean_dataset.py`](./clean_dataset.py): **dataset cleaning** with [SeqKit2](https://bioinf.shenwei.me/seqkit/). Takes the directory holding the FASTA files, walks it recursively and writes a `<dataset>_clean` copy with the very same tree and file names, so downstream GREAC runs need no other change. Applies `rmdup` (duplicated sequences inside each file), `common` (sequences shared between the files of a group, compared pairwise), `seq` (the shared sequences are extracted to a reference file under `.clean_meta/`) and `grep -v` (the filtered files, without the shared sequences, to avoid over-fitting).

  ```bash
  # one group per directory (default)
  python scripts/clean_dataset.py ~/Desktop/datasets/virus_dataset

  # every class in its own subdirectory: compare the whole tree, by sequence
  python scripts/clean_dataset.py ~/Desktop/datasets/virus_dataset --group all --by seq -j 8
  ```

  **Requires `seqkit` (SeqKit2) in the PATH** (or `--seqkit /path/to/seqkit`). Runs on Python 3.6+ (standard library only). See `--help` for `--by {id,name,seq}`, `--group {dir,all}`, `-j/--threads`, `-o/--output-dir` and `-f/--force`.

- [`benchmark.sh`](./local/benchmark.sh): provides a **complete execution** of GREAC, benchmarking datasets. 
- [`extract-features.sh`](./local/extract-features.sh): extraction of features and model fitting (regions and frequency behavior). 
- [`fasta-regions.sh`](./local/fasta-regions.sh): create reduced FASTA files with the extract regions. 
- [`file-classification.sh`](./local/file-classification.sh): classify FASTA sequences files based on pre-trained model by `extracted-features`. 