# Code

This directory contains the active simulation and real-data code for
"Dimension-Free Principal Component Analysis with High Complexity".

- `simulation/main.R`: unified Monte Carlo, summary, and plotting entry point.
- `simulation/pca_function/`: experiment implementations and shared inference utilities.
- `real_data/main.R`: unified entry point for EEGMMIDB and CAPTURE-24.
- `real_data/`: dataset-specific download/preprocessing code and shared analysis helpers.
- `install_packages.R`: installs the additional R packages used by real-data preprocessing.

Run commands from the repository root.

```bash
# Lightweight code check
Rscript code/simulation/main.R mode=smoke

# Full simulation sequence
Rscript code/simulation/main.R experiment=sequence ncores=20

# Analyze the processed real-data caches included in the repository
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10

# Complete real-data rebuild from public source data
Rscript code/real_data/main.R dataset=eeg action=all ncores=10
Rscript code/real_data/main.R dataset=capture24 action=all ncores=10
```

The public repository intentionally omits raw public datasets and Monte Carlo
replication-level output. The corresponding directories are generated locally
when the download or simulation pipelines are run. Compact processed data,
summary CSV files, and manuscript figures are retained in the repository.

Both simulation and real-data analyses call
`simulation/pca_function/inference_utils.R`; there is no separate real-data
implementation of the estimators. CAPTURE-24 uses Python only for streaming and
compressing the raw accelerometry archive into the deterministic five-second
analysis cache.
