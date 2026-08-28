# PhysioNet EEGMMIDB R04

Source: <https://physionet.org/content/eegmmidb/1.0.0/>

The real-data pipeline uses run R04 (imagined left/right fist) for 105
prespecified usable subjects. The public repository does not version the EDF
files; they are downloaded from PhysioNet when `action=download` or `action=all`
is requested.

`prepared/eegmmidb_r04__0003bb07.rds` is the deterministic cache of 8--30 Hz
right-minus-left log-power contrasts used by the manuscript analysis. A fresh
analysis can therefore run directly from the cache with:

```bash
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
```

Delete the prepared RDS file, or disable cache reuse in the pipeline, only when
the preprocessing itself needs to be rebuilt.
