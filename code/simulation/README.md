# Simulation Code

`pca_function/` contains one R file per experiment family plus shared model,
inference, summary, and plotting utilities. `main.R` provides experiment
selection, parameter overrides, sequential execution, saved-result summary/plot
actions, run IDs, and elapsed-time logging.

The production sequence is:

```text
1a -> 1b -> 2a -> 2b -> 3a -> 3b -> 3c -> 3d -> k0 -> universality
```

Run from the project root:

```bash
Rscript code/simulation/main.R experiment=sequence ncores=20
Rscript code/simulation/main.R mode=smoke
```

Replication-level Monte Carlo files are intentionally not committed to the
public repository. Consequently, `action=summarize` or `action=plot` on a fresh
clone requires the corresponding local simulation run first. The compact final
summaries and publication figures from the manuscript runs are retained.

Experiment 3D studies a local-to-tie multiplicity-two design using the symmetric
relative eigengap `2*(a-b)/(a+b)`. The classical FPCA comparison uses
fixed-covariance repeated-root equality inference: a Dauxois-type Gaussian
perturbation-block test and Anderson's Gaussian equality-of-roots benchmark.
The ordinary classical simple-root FPCA gap Wald interval is explicitly marked
`nonregular_local_to_tie` in this experiment and is not used for primary
coverage. Proposed inference uses deflated spectrum inversion and the
high-complexity Rayleigh/Rice/noncentral-chi-square limit.

Classical FPCA functions do not use the Proposed phase-margin or reliability
quantities. See the root `README.md` for the experiment overview.
