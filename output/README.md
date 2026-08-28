# Output

This public directory retains compact manuscript-facing numerical summaries and
real-data outputs. Large or run-specific generated artifacts are reproducible
from the code and are intentionally not versioned.

Included:

- `data/summary/`: compact simulation summaries and figure-data CSV files.
- `real_data/eegmmidb/`: named EEG summaries used by the manuscript.
- `real_data/capture24/`: named CAPTURE-24 summaries used by the manuscript.
- `figure_output_summary.csv`: index linking publication figures to plot-data
  CSV files and analysis settings.
- `consistency_report.csv`: checks against locked manuscript-facing numerical
  results.

Not versioned:

- `data/replicate/`: Monte Carlo replication-level files.
- `data/truth/`: generated truth/configuration snapshots.
- `data/manifest/`: run registries and latest-run pointers.
- `table/`: generated LaTeX tables.
- `log/`: scheduler and batch logs.
- `real_data/manifest/`: run registries.
- run-ID-prefixed `real_*` files under the real-data output directories.

The revised eigengap analysis is represented by the `*_symmetric.csv` Experiment
3D summaries. Older asymmetric-gap files are not part of the public release.
