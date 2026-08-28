# Reproducibility Materials for Dimension-Free PCA Inference

This repository contains the code and compact analysis outputs for the paper
**“Dimension-Free Principal Component Analysis with High Complexity.”** The
implementation covers the simulation studies, the CAPTURE-24 accelerometry
application, and the EEG Motor Movement/Imagery Database application.

The public repository is intentionally distributed without the large source
data files and without Monte Carlo replication-level output. Both are
reproducible from the source code. Deterministic processed real-data caches,
manuscript-facing summary CSV files, and publication figures are retained.

## Repository contents

```text
code/                  simulation and real-data source code
data/                  processed real-data caches and documentation
manuscript/            manuscript/supplement figure files
output/                compact manuscript-facing numerical summaries
software_versions.txt  recorded computation-server environment
```

Not versioned:

```text
data/eegmmidb/raw/        downloaded PhysioNet EDF files
data/capture24/source/    downloaded CAPTURE-24 archive
output/data/replicate/    Monte Carlo replication-level results
output/data/truth/        generated truth/configuration snapshots
output/data/manifest/     generated run registries
output/table/             generated LaTeX tables
output/log/               local scheduler/log files
output/real_data/manifest/ real-data run registries
```

The `.gitignore` records these exclusions. Text files use LF line endings;
binary data and figure formats are marked in `.gitattributes`.

## 1. Quick start

Run commands from the repository root. Install non-base packages used by the
real-data workflow with:

```bash
Rscript code/install_packages.R
```

A lightweight simulation check is:

```bash
Rscript code/simulation/main.R mode=smoke
```

A complete Monte Carlo sequence is:

```bash
Rscript code/simulation/main.R experiment=sequence ncores=20
```

Analyze the included processed real-data caches without downloading the raw
data:

```bash
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10
```

Use `action=all` only when a complete download and preprocessing rebuild is
required. CAPTURE-24 raw preprocessing is substantially more expensive than
analysis of the included cache.

## 2. Simulation interface

`code/simulation/main.R` is the simulation entry point. Available experiment
selections are:

```text
1a, 1b, 2a, 2b, 3a, 3b, 3c, 3d, k0, universality,
1, 2, 3, robustness, all, sequence
```

The main families are:

- **1A:** panel/longitudinal phase-transition examples.
- **1B:** controlled phase-margin and conditioning experiment.
- **2A:** functional inference across covariance and score designs.
- **2B:** functional eigenvalue inference, including classical FPCA versus
  Proposed coverage and confidence-interval length.
- **3A:** large-domain simple-direction recovery.
- **3B:** large-domain repeated-eigenspace recovery.
- **3C:** dense-grid/discretization saturation.
- **3D:** multiplicity-two eigengap inference under local-to-tie alternatives.
- **k0 / universality:** K0 sensitivity and non-Gaussian robustness.

Production defaults use master seed `20260810`, deterministic population-design
seeds, replication-specific seeds derived from the experiment and replication
indices, covariance divisor `n`, and Gram-matrix PCA. Changing `ncores` does not
change the Monte Carlo draws on Unix-like systems.

## 3. Core Proposed inference

For ordered nonzero sample eigenvalues, the first `K0` components are excluded
from empirical bulk transforms

```text
mhat_s(z) = (1/n) * sum_{l > K0} (hatlambda[l] - z)^(-s).
```

The common Proposed estimators are

```text
hatalpha = -1 / mhat_1
hatDelta = mhat_1^2 / mhat_2
hatr2    = -mhat_1 / (hatlambda * mhat_2).
```

`code/simulation/pca_function/inference_utils.R` is the shared implementation
used by simulations and real-data analyses. It also contains the universal
standard errors, Wald intervals, one-sided phase lower bound, and the
multiplicity-two eigengap procedures.

## 4. Classical FPCA versus Proposed inference

The two labels denote independently justified inferential procedures.

For a simple classical FPCA eigenvalue, the fixed-covariance first-order
expansion estimates the asymptotic variance by the empirical variance of squared
FPC scores divided by `n`. This implementation follows the classical Hilbert-
space PCA/FPCA asymptotic framework of Dauxois, Pousse and Romain (1982) and
Hall and Hosseini-Nasab (2006).

Equality of two classical FPCA roots is nonregular. The primary general
classical benchmark estimates the covariance of the two-dimensional score
perturbation block and simulates its Gaussian repeated-root limit. Anderson's
(1963) chi-square equality-of-roots statistic is saved as the Gaussian
benchmark. Classical FPCA inference never calls the Proposed `Delta` or `r2`
quantities.

The Proposed procedure instead uses empirical-spectrum inversion. For a
multiplicity-two cluster it evaluates pooled bulk transforms at the mean sample
root and uses the high-complexity Gaussian-matrix limit. The formal adjacent-
pair target is the symmetric relative eigengap

```text
g(a,b) = 2*(a-b)/(a+b).
```

Under the universal condition, the two-root corrected statistic has the
Rayleigh/Rice, equivalently central/noncentral chi-square with two degrees of
freedom, representation used by the code. Older asymmetric-gap quantities are
legacy diagnostics only.

## 5. Experiment 3D: local-to-tie eigengap inference

Experiment 3D uses

```text
alpha1 = 6 * (1 + delta/(2*sqrt(n)))
alpha2 = 6 * (1 - delta/(2*sqrt(n)))
delta  = 0, 0.5, 1, 2, 3.
```

The true symmetric population gap is exactly `delta/sqrt(n)`. Because the roots
approach one another on the `n^(-1/2)` scale, the ordinary classical
simple-root FPCA eigengap Wald interval is nonregular throughout Experiment 3D.
It is retained only as QA metadata and is **not** used for manuscript-facing
coverage. The primary FPCA comparison uses the fixed-covariance repeated-root
equality test and Anderson benchmark. Proposed reports its equality test and
local-gap magnitude inference.

The compact final outputs are the `eigengap_*_symmetric.csv` files under
`output/data/summary/`, with corresponding publication figures under
`manuscript/simulation/`.

## 6. Real-data applications

### EEGMMIDB

The EEG analysis uses R04 for 105 prespecified usable subjects, with 64 channels
and 8--30 Hz right-minus-left log-power trajectories. The primary analysis uses
`K0=10` and reports PC1--PC6 population-scale eigenvalue, phase-margin and
reliability estimates; K0 sensitivity uses `8,10,12`. PC3--PC4 is the only
prespecified candidate joint eigenspace. The PC2 lateralization analysis also
contains matched 500-subject-bootstrap FPCA and Proposed percentile intervals.

### CAPTURE-24

The CAPTURE-24 analysis uses 150 wake-aligned participants after the
prespecified wake-anchor exclusion. Five-second log-ENMO curves are analyzed
over nested 6, 12, 18 and 24 hour domains. The primary analysis uses `K0=8`,
with sensitivity values `6,8,10`, and reports PC1--PC6 population-scale
eigenvalue, phase-margin and reliability estimates. PC3--PC4 and PC5--PC6 are
the prespecified non-overlapping candidate joint eigenspaces.

### Interpretation rules

Phase separation, sample--population alignment, and adjacent-root separation are
kept as distinct diagnostics. Pairwise nonrejection means only insufficient
evidence to distinguish the roots; it does not automatically create a cluster.
For classical FPCA near-tie pairs, the equality p-value remains available while
the ordinary distinct-root eigengap Wald interval is suppressed and marked
`not valid under near-tie`. FPCA decisions do not gate the validity of Proposed
phase or component estimates. The historical 0.05 relative-gap threshold is
retained only in legacy/run-specific diagnostics.

## 7. Outputs and figures

Compact public outputs are retained under `output/data/summary/` and
`output/real_data/<dataset>/`. Replication-level simulation files, run
registries, truth snapshots, generated tables/logs, and run-ID-prefixed
real-data files are generated locally and excluded from version control.

Every manuscript PDF is an independent title-free panel with a corresponding
plot-data CSV. `output/figure_output_summary.csv` links publication figures to
their setting, quantity, comparison, plot-data file, and source run.

## 8. Software environment

The recorded computation-server environment is listed in
`software_versions.txt`. The install helper installs missing current CRAN
packages and is not a lockfile. Simulation uses base/recommended R facilities;
EEG preprocessing additionally uses `edfReader` and `signal`, and CAPTURE-24
preprocessing uses `jsonlite` plus Python with `numpy` and `pandas`.
