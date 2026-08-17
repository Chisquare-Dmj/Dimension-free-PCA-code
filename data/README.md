# Data

Large public raw data files are intentionally **not included** in this GitHub
repository. The download code and deterministic processed caches needed for the
manuscript analyses are included.

Included:

- `eegmmidb/prepared/eegmmidb_r04__0003bb07.rds`: prepared 64-channel
  right-minus-left log-power contrasts for the 105-subject EEG analysis.
- `capture24/processed/cfg0003e622/`: wake-aligned five-second CAPTURE-24
  log-ENMO matrix, coverage mask, subject QC, metadata, activity-label
  dictionary, and preprocessing manifest.

Not included:

- `eegmmidb/raw/`: PhysioNet R04 EDF files.
- `capture24/source/`: the official CAPTURE-24 archive and download manifest.

Download the public source data from the repository root with

```bash
Rscript code/real_data/main.R dataset=eeg action=download ncores=10
Rscript code/real_data/main.R dataset=capture24 action=download
```

To reproduce preprocessing as well as the analysis, use `action=all`. To rerun
the manuscript analysis directly from the included processed caches, use
`action=analyze`.

Processed cache names contain a preprocessing-configuration fingerprint and are
reused only when the settings match.
