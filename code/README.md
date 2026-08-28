# Code

This directory contains all active code used for the numerical supplement.

- `simulation/main.R` is the unified Monte Carlo, summary, and plotting entry point.
- `simulation/pca_function/` contains the experiment functions and shared utilities.
- `real_data/main.R` is the unified entry point for all real-data applications.
- `real_data/` also contains preprocessing and shared inference helpers.
- `install_packages.R` installs the packages required by real-data preprocessing.

Run both entry points from the project root:

```bash
Rscript code/simulation/main.R experiment=sequence ncores=20
Rscript code/real_data/main.R dataset=eeg action=all ncores=10
Rscript code/real_data/main.R dataset=capture24 action=all ncores=10
Rscript code/real_data/main.R dataset=sequence action=all ncores=10
```

Both pipelines call `simulation/pca_function/inference_utils.R`; there is no
separate Python or real-data reimplementation of the estimators. CAPTURE-24
uses a Python helper only to stream and compress its roughly one billion raw
accelerometer rows into a compact five-second matrix.
