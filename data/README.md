# Data

This public repository contains deterministic processed caches needed to rerun
the real-data analyses without committing the multi-gigabyte source archives.

Included:

- `eegmmidb/prepared/`: cached 64-channel right-minus-left log-power contrasts.
- `capture24/processed/`: wake-aligned five-second log-ENMO matrix, coverage
  mask, participant QC, metadata, label dictionary, and preprocessing manifest.

Not versioned:

- `eegmmidb/raw/`: downloaded PhysioNet EDF files.
- `capture24/source/`: the downloaded CAPTURE-24 archive.

The raw sources can be reconstructed with the real-data pipeline. For example:

```bash
Rscript code/real_data/main.R dataset=eeg action=download ncores=10
Rscript code/real_data/main.R dataset=capture24 action=download
```

The prepared cache names contain preprocessing-configuration fingerprints and
are reused only when the corresponding settings match.
