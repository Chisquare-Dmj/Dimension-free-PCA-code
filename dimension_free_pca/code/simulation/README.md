# Simulation Code

`pca_function/` contains the experiment families and shared model, inference,
summary, and plotting utilities. `main.R` provides experiment selection,
parameter overrides, deterministic run IDs, parallel Monte Carlo execution, and
saved-result summary/plot actions.

The default production sequence is

```text
1a -> 1b -> 2a -> 2b -> 3a -> 3b -> 3c -> k0 -> universality
```

Run from the repository root:

```bash
Rscript code/simulation/main.R experiment=sequence ncores=20
Rscript code/simulation/main.R mode=smoke
```

Individual branches can be run with, for example,

```bash
Rscript code/simulation/main.R experiment=1b ncores=20
Rscript code/simulation/main.R experiment=2b ncores=20
Rscript code/simulation/main.R experiment=3a ncores=20
```

A full run writes replication-level results, truth/configuration snapshots,
manifests, summaries, tables, and figures. The public repository intentionally
does **not** commit `output/data/replicate/`, `output/data/truth/`, or the
run-specific manifests/tables. These files are recreated automatically by the
simulation pipeline. Compact manuscript-facing summary CSV files and figures
are retained.

Consequently, on a fresh clone the commands

```bash
Rscript code/simulation/main.R experiment=all action=summarize
Rscript code/simulation/main.R experiment=all action=plot
```

should be used only after the corresponding simulations have been run locally.
See the root `README.md` for scientific purpose, default hyperparameters, and
resource requirements for each experiment.
