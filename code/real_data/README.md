# Real-Data Code

- `common_real_data.R`: centered Gram PCA, separation diagnostics, K0 sensitivity,
  run manifests, tables, and figure output.
- `eegmmidb.R`: PhysioNet R04 download, EDF+ event parsing, 8--30 Hz log-power
  preprocessing, subject contrasts, lateralization direction, and analysis.
- `capture24.R`: CAPTURE-24 archive verification, wake-aligned five-second
  log-ENMO processing, nested-domain complexity diagnostics, K0 sensitivity,
  preprocessing sensitivities, and deterministic analysis output.
- `capture24_preprocess.py`: memory-bounded conversion of the 100 Hz CAPTURE-24
  participant files. It does not implement PCA or inference.

The public repository does not include the large raw datasets. It does include
the processed EEG and CAPTURE-24 caches used by the manuscript analyses.
Therefore two workflows are available.

Analyze the included caches:

```bash
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10
```

Rebuild from the public source data:

```bash
Rscript code/real_data/main.R dataset=eeg action=all ncores=10
Rscript code/real_data/main.R dataset=capture24 action=all ncores=10
```

Separate `download` and `preprocess` actions are also available. Run-specific
manifests are generated locally and are not committed; compact named outputs in
`output/real_data/` and publication figures in `manuscript/real_data/` are
retained in the repository.

The real-data analyses are deterministic sample analyses. CAPTURE-24 uses K0
values 6, 8, and 10 for sensitivity checks and reports the 30-minute
missingness and three-hour main-sleep sensitivities. EEG uses K0 values 8, 10,
and 12 and reports the conventional and corrected PC2 lateralization point
estimates. Because population truth is unavailable in real data, these
comparisons are descriptive rather than accuracy or bias assessments.
