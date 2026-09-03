# Reproducibility Materials for Dimension-Free PCA Inference

This repository contains the code and compact analysis outputs for the paper
**“Dimension-Free Principal Component Analysis with High Complexity.”** The
implementation covers the simulation studies, the CAPTURE-24 accelerometry
application, and the EEG Motor Movement/Imagery Database application.

The repository is intentionally distributed without the large public raw data
files and without Monte Carlo replication-level output. Both are reproducible
from the source code. Compact processed real-data caches, manuscript-facing
summary CSV files, and publication figures are included so that the numerical
results can be inspected without committing multi-gigabyte data or very large
simulation output.

## Repository contents

```text
code/          simulation and real-data source code
data/          processed real-data caches and data documentation
manuscript/    manuscript/supplement figure files
output/        compact manuscript-facing numerical summaries
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
```

The `.gitignore` records these exclusions. The `.gitattributes` file keeps
source and text files in LF format across platforms; PDF, PNG, RDS, BIN, ZIP,
and EDF files are treated as binary.

## 1. Quick start

Run commands from the repository root. Install the additional packages used by
the real-data pipeline with

```bash
Rscript code/install_packages.R
```

A lightweight simulation check is

```bash
Rscript code/simulation/main.R mode=smoke
```

The complete simulation sequence is

```bash
Rscript code/simulation/main.R experiment=sequence ncores=20
```

The real-data analyses can either use the included processed caches or be
rebuilt from the public source data. To analyze the included caches without
redownloading the raw data, run

```bash
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10
```

For complete reproduction from the public raw datasets, use `action=all`; the
pipeline downloads, preprocesses, analyzes, summarizes, and plots the selected
dataset:

```bash
Rscript code/real_data/main.R dataset=eeg action=all ncores=10
Rscript code/real_data/main.R dataset=capture24 action=all ncores=10
```

CAPTURE-24 is a multi-gigabyte download and preprocessing is substantially more
expensive than analysis of the included processed cache.

## 2. Simulation interface

`code/simulation/main.R` is the entry point for simulation, summary, and
plotting tasks. The available experiment selections are

```text
1a, 1b, 2a, 2b, 3a, 3b, 3c, k0, universality,
1, 2, 3, robustness, all, sequence
```

Common examples are

```bash
Rscript code/simulation/main.R experiment=1a ncores=20
Rscript code/simulation/main.R experiment=1b n=200 p=200 R=500 Delta='c(0.2,0.4,0.6,0.8)'
Rscript code/simulation/main.R experiment=3a n=300 R=500 ncores=10
Rscript code/simulation/main.R experiment=sequence ncores=20
```

The default action is `all`, which runs the selected Monte Carlo experiment and
then writes summaries, tables, and figures. Other actions are `simulate`,
`summarize`, and `plot`. Because replication-level files are intentionally not
committed to this public repository, a fresh clone must run the corresponding
simulation before using `action=summarize` or `action=plot` for that run.

Every production run receives a configuration fingerprint automatically. The
fingerprint is used only for generated local artifacts; users do not enter it
on the command line.

Common configuration keys include `ncores`, `master_seed`, `confidence_level`,
`replications`, `M`, `K0`, and `certification_threshold`. CLI aliases include
`R=replications`, `cores=ncores`, `seed=master_seed`, `m=m_values`,
`Delta=Delta_grid`, and `K0s=K0_values`.

## 3. Reproducibility rules

The production defaults use the following deterministic rules:

- master seed `20260810`;
- fixed population-design seeds within each simulation scenario;
- replication seed `MASTER_SEED + 100000 * experiment_id + replication_id`;
- population covariance, Haar matrices, and spike directions held fixed across
  replications within each scenario;
- confidence level `0.95`;
- sample-size normalization by `n` throughout the Gram-matrix calculations.

Parallel simulation uses `parallel::mclapply` on Unix-like systems. Each worker
sets its assigned deterministic seed, so changing `ncores` does not alter the
simulated populations or Monte Carlo draws. On Windows the helper falls back to
serial execution.

## 4. Core inference implementation

`code/simulation/pca_function/inference_utils.R` contains the common
implementation of the empirical bulk transforms, population-spike inversion,
phase-margin and reliability estimators, derivative factors, universal standard
errors, Wald intervals, and phase lower bound. The simulation and real-data
pipelines call the same implementation.

PCA is computed from the `n x n` Gram matrix, allowing vector, matrix,
coefficient-space functional, and discretized functional observations to share
the same inference API without forming a potentially very large ambient sample
covariance matrix.

## 5. Experiment 1A: Panel Phase Transition

Default hyperparameters:

```text
n = 300, p = 20, T = 30, p*T = 600
replications = 500, M = 3, K0 = 5
score distribution = Gaussian
```

`P1_independent` uses fixed diagonal cross-sectional and temporal covariance
factors with entries drawn from `Uniform[0.5,1.5]`. Three random Frobenius
directions are orthonormalized once, held fixed, and assigned signal variances
`(6,4,2)`.

`P2_block_AR1` uses cross-sectional covariance
`A[k,l] = 0.6^abs(k-l)` and three compound-symmetry temporal blocks
`0.55 * 11' + 0.45 * I`. The first two bulk eigenvalues use multipliers `(5,4)`.

In both P1 and P2, the final third population eigenvalue is calibrated after
constructing the covariance so that `Delta3=-0.10`. Consequently, both panels
contain two supercritical components and one genuinely subcritical component,
while `alpha3` remains above the largest population bulk eigenvalue. The implied
spectral adjustment is saved in the truth CSV and `target_Delta3` is configurable.

Regular outlier inference is computed only for spikes 1 and 2. For spike 3, the
locally generated replication-level CSV saves `hat_lambda`, population truth, and phase status, sets
`regular_inference_valid=FALSE`, and leaves `hat_alpha`, `hat_Delta`, `hat_r2`,
SE, CI, and coverage fields as `NA`.

The locally generated primary CSV contains one row per Monte Carlo replication and spike. The auxiliary bulk
spectrum CSV contains sample ranks 4 through `n` in long form for density plots.
The corresponding spectral figure reconstructs the replication-level sample bulk edge as
`hat_lambda3 - sample_gap_to_lambda4` and overlays its density as a dashed line.

## 6. Experiment 1B: Controlled Phase Margin

Default hyperparameters:

```text
n = p = 300, M = K0 = 1
main Delta grid = 0.40, 0.50, ..., 0.90
boundary Delta grid = 0.10, 0.20, 0.30
replications per Delta = 2000
C = diag(alpha,1,...,1), Gaussian scores
```

For `c = (p-1)/n`, each spike is set by
`alpha = 1 + sqrt(c/(1-Delta))`. Saved diagnostics include naive and corrected
relative bias, scaled empirical variance, all three interval lengths, coverage,
and the indicator that the phase lower bound exceeds `0.20`. Every saved row is
labelled `phase_region=main` or `phase_region=boundary`. Figure 2(c,d) uses only
the regular-interior main grid; the boundary rows remain available for a
separate finite-sample stress analysis. Setting `K0=M` avoids extra leading-component removal in this baseline conditioning experiment.

## 7. Experiment 2A: Functional Inference Across Designs

Default hyperparameters:

```text
n = 300, p = 100, N = 200, decay = 1.1
mu[k] = k^(-1.1), M = 3, K0 = 5
regular inference is evaluated for spikes 1 and 2
```

The primary experiment uses exact orthonormal coefficient-space inner products,
so discretization error is not mixed into inference calibration.

Covariance scenarios:

- `F1_diagonal`: `diag(6,4,3,1,...,1)`.
- `F2_block`: ten compound-symmetry correlation blocks
  `0.55 * 11' + 0.45 * I`, with their three largest eigenvalues multiplied
  by `(5,4,3)`.
- `F3_Haar`: one fixed Haar rotation of `diag(6,4,3,1,...,1)`.

Score distributions:

- `Gaussian`: standard normal.
- `Uniform`: `Uniform[-sqrt(3),sqrt(3)]`.
- `t12`: `sqrt(5/6) * t_12`.

The `t12` design is intentionally retained as a robustness experiment. It is not
sub-Gaussian, so its empirical behavior should be described separately from a
literal verification of assumptions that require sub-Gaussian scores.

All `3 x 3` combinations are run. The five combinations emphasized in the
specification default to 2000 replications; the remaining four default to 1000.
Both counts are configurable. Each result row also stores the actual squared
sample/population eigenvector alignment.

Experiment 2A remains the broad `3 x 3` robustness comparison and has its own
coverage diagnostic figure. The representative theorem-calibration Figure 3
is generated from Experiment 2B, not selected retrospectively from this grid.

## 8. Experiment 2B: Functional Inference Asymptotics

Experiment 2B is a separate triangular-array experiment for convergence in
bias, RMSE, SE calibration, and coverage. It does not replace Experiment 2A.

```text
n = 150, 300, 600
p = nearest multiple of 10 to n/3 = 50, 100, 200
N = 200, decay = 1.1
replications = 2000, 2000, 2000
M = K0 = 3
scenario = F3_Haar + standardized Uniform
evaluated spike = 2
```

The population covariance, `true_alpha`, `true_psi`, `true_Delta`, `true_r2`,
and all population derivative factors are rebuilt separately for every `(n,p)`.
The block-compatible rounding satisfies `p/n -> 1/3` and is exact for the
default sample sizes. For example, `n=1000` gives `p=330` rather than raising
an integrality error. Fixed `N=200` isolates statistical asymptotics from
functional truncation changes.

The middle spike is specified in advance as the representative interior regular
component, away from both the phase boundary and the extreme strong-signal
regime. Figure 3 uses its `n=300` studentized errors. The summary CSV reports
relative bias and relative RMSE for `alpha`; bias and RMSE for `Delta` and
`r2`; empirical SD; average estimated SE; SD/mean-SE; coverage; `sqrt(n)`
times RMSE; and `sqrt(n)` times mean SE. Its dedicated LaTeX writer produces
three panels: population spike, phase margin, and reliability.
The final Figure 3 panel plots `Z_alpha` against `Z_Delta` and `Z_r2`, with
`Z_T=(hat_T-T)/se(hat_T)` and a `y=x` reference, so it directly displays the
joint rank-one first-order fluctuation without dividing by a small derivative.

## 9. Experiment 3A: Simple Direction

Default asymptotic grid:

```text
n = 150, 300, 600
T = n/3, N = 2n
replications = 2000, 2000, 2000
M = 3, K0 = 5, Gaussian scores
```

The spectrum is `(6,4,3)` followed by repeated decay levels
`(floor((j-4)/T)+1)^(-1.1)`. The fixed unit target is

```text
f = 0.6 phi1 + 0.3 phi2 + sqrt(0.55) phi13,
theta_true = 0.6.
```

The sample eigenvector is oriented by `phi1`. Conventional FPCA estimates the
target by the raw projection; Proposed divides that same projection by the
estimated reliability. The publication summary labels these methods `FPCA` and
`Proposed` and reports mean, bias, and RMSE for every `n`.

## 10. Experiment 3B: Repeated Eigenspace

This experiment uses the same grid, decay, and 2000 replications per sample size as Experiment 3A but changes the
leading spectrum to `(6,6,3)`. The first cluster has multiplicity two and

```text
Q_true(f) = 0.6^2 + 0.3^2 = 0.45.
```

The FPCA projector functional is divided by the mean estimated reliability of
sample PCs 1 and 2. The squared individual alignment between sample PC 1 and
population direction 1 is saved only as a non-identifiability diagnostic; it is
not treated as an estimator.

The independent comparison figure plots Monte Carlo mean against `n` for FPCA
and Proposed, with the population truth as a horizontal reference. Its CSV also
retains bias and RMSE over the full `n=150,300,600` convergence grid.

## 11. Experiment 3C: Dense-Grid Approximation

Default hyperparameters:

```text
n = 300, T = 100, N = 200
m = 500, 1000, 1500, 2000, 3000
replications = 1000, M = 3, K0 = 5
```

Each replication generates one latent coefficient matrix. The oracle PCA and
all discretized PCAs reuse those identical latent curves. For every `m`, the CSV
stores `sqrt(n)` times the absolute oracle/discrete discrepancy for corrected
`alpha`, `Delta`, `r2`, and `theta`.

## 12. Robustness Experiments

`robustness_k0` uses `(n,p)=(150,50),(300,100),(600,200)`, `N=200`,
`decay=1.1`, `M=3`, and 2000 replications at each sample size.
It uses `F3_Haar + Uniform`, evaluates the middle spike, and applies
`K0=3,5,8` to the same sample in each replication. The locally generated replication-level CSV saves paired differences between each
`K0>M` estimate and its `K0=M` baseline. The summary reports bias, RMSE,
`sqrt(n)`-scaled RMSE, and `n`-scaled RMSE, so sensitivity to the number of
leading components removed from the bulk estimate can be assessed directly
along the asymptotic sequence. The dedicated figure uses `n` times paired RMSE
to check the predicted first-order stability. The LaTeX writer presents these
quantities in separate panels for `alpha`, `Delta`, and `r2`.

`robustness_universality` uses aligned diagonal functional covariance with
`n=300`, `p=100`, `N=200`, Gaussian, `t12`, and Uniform scores,
evaluates spikes 1 and 2, and uses `K0=5`.
The summary compares empirical SD with the observable universal theoretical SD
and with the full fourth-cumulant theoretical SD. Because all three estimators
share one rank-one first-order fluctuation, their full SD is obtained by
multiplying the corresponding universal SE by the common factor
`sqrt(1 + kappa * true_Delta / 2)`. The figure reports both theoretical-to-
empirical SD ratios; coverage remains available as a supplementary diagnostic.
The fourth cumulants used are `0`, `0.75`, and `-1.2`, respectively.

## 13. Output layout and public-repository policy

A full local simulation run creates run-specific replication, truth,
configuration, manifest, summary, and table files. These generated artifacts
use complete configuration fingerprints so that distinct runs cannot overwrite
one another. The replication-level and run-specific files are intentionally not
committed here because they are large and some generated names are unsuitable
for portable Windows paths.

The repository retains the compact manuscript-facing outputs:

```text
output/data/summary/       compact simulation plot/table data
output/real_data/eegmmidb/ compact EEG analysis outputs
output/real_data/capture24/ compact CAPTURE-24 analysis outputs
manuscript/simulation/     publication simulation figures
manuscript/real_data/      publication real-data figures
```

The compact CSV files are sufficient to inspect the values displayed in the
paper. To recreate the full run registry, truth/configuration snapshots,
replication-level files, or LaTeX tables, rerun the corresponding analysis.

### FPCA comparison outputs

The conventional baseline is the ordinary empirical FPCA decomposition of the
same centered or weighted data. For population-spike comparison it uses the
sample eigenvalue as the usual plug-in quantity. For a fixed directional target
it uses the raw empirical projection, and for a repeated eigenspace target it
uses the raw empirical projector functional. The proposed procedure estimates
the population spike and corrects the corresponding fixed functional by the
estimated reliability. It does not define a corrected full eigenfunction.

The simulation comparison files report Monte Carlo means, bias or relative
bias, empirical SD where relevant, and RMSE against known population targets.
The real-data comparisons are descriptive because population truth is unknown.

## 14. Rebuilding summaries and figures

After a simulation has been run locally, its summaries and figures can be
regenerated without repeating the Monte Carlo draws:

```bash
Rscript code/simulation/main.R experiment=all action=summarize
Rscript code/simulation/main.R experiment=all action=plot
```

A fresh clone does not contain `output/data/replicate/`, so these commands
require the relevant simulation to have been run first. The committed compact
CSV and figure files record the final manuscript-facing results from the
production runs used for the paper.

For real data, `action=analyze` recreates analysis outputs from the included
processed caches. `action=all` performs a complete rebuild from the public raw
data. `action=plot` and `action=summarize` operate on locally generated run
manifests and therefore should be used after `analyze` or `all` on a fresh clone.

## 15. Real-data applications

### EEG Motor Movement/Imagery Database

The analysis uses run R04 from PhysioNet EEGMMIDB v1.0.0 for 105 prespecified
subjects. The public EDF files are not committed. Download them with

```bash
Rscript code/real_data/main.R dataset=eeg action=download ncores=10
```

The repository includes the deterministic prepared cache
`data/eegmmidb/prepared/eegmmidb_r04__0003bb07.rds`, so the analysis can also be
run directly with `action=analyze`. The preprocessing forms 64-channel
right-minus-left log-power contrasts over 8--30 Hz and outputs 80 time points
per channel. The primary analysis uses `K0=10` with sensitivity values
`8, 10, 12`.

### CAPTURE-24

The public CAPTURE-24 archive is not committed. It can be downloaded and
verified with

```bash
Rscript code/real_data/main.R dataset=capture24 action=download
```

The repository includes the deterministic wake-aligned five-second processed
cache under `data/capture24/processed/cfg0003e622/`, including the log-ENMO
matrix, coverage mask, subject QC, metadata, and preprocessing manifest. The
primary analysis uses 150 participants and nested 6-, 12-, 18-, and 24-hour
wake-aligned domains, with `K0=8` and sensitivity values `6, 8, 10`.

The real-data analyses are deterministic sample analyses. CAPTURE-24 reports
how spectral complexity and leading-component reliability change with domain
length and treats closely spaced component pairs through their joint
eigenspaces. EEG reports the population-scale correction for the prespecified
PC2 lateralization functional and treats the closely spaced PC3--PC4 pair as an
eigenspace target.

## 16. Data and storage notes

The raw datasets are omitted solely to keep the repository portable:

- EEGMMIDB is downloaded directly from PhysioNet by the R pipeline.
- CAPTURE-24 is downloaded from the Oxford University Research Archive by the R
  pipeline and is streamed during preprocessing rather than permanently
  extracted as a second copy.

The included processed caches are substantially smaller than the source data
and allow the manuscript analyses to be rerun without a fresh download. See
`data/README.md`, `data/eegmmidb/README.md`, and `data/capture24/README.md` for
more detail.

## 17. Requirements and computational notes

Simulation code uses base R and the recommended `parallel` package. Real-data
preprocessing additionally uses `edfReader`, `signal`, and `jsonlite`.
CAPTURE-24 preprocessing uses Python 3 with `numpy` and `pandas`, plus the
system commands `curl`, `unzip`, and `sha256sum` for download/integrity checks.
The statistical analysis itself remains in R.

The larger functional experiments are computationally intensive. Memory use
grows with both the observation size and the number of active workers. In
particular, the `n=600` large-domain experiments and Experiment 2B at
`(n,p)=(600,200)` require repeated eigendecompositions of `600 x 600` Gram
matrices. Increase `ncores` only after checking available RAM.

## GitHub / Windows note

The repository includes `.gitattributes` so text files are normalized to LF.
Messages such as `LF will be replaced by CRLF` are line-ending warnings rather
than analysis errors; they should disappear when this cleaned repository is
added from a fresh checkout. The generated long run-specific output files have
also been excluded, so the repository does not require Git's Windows
`core.longpaths` option.
