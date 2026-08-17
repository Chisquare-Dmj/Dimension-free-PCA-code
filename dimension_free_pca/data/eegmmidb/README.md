# PhysioNet EEG Motor Movement/Imagery Database (EEGMMIDB), run R04

Source: <https://physionet.org/content/eegmmidb/1.0.0/>

The raw EDF files are public but are not committed to this repository. The
pipeline downloads run R04 (imagined left/right fist) for the 105 prespecified
usable subjects and writes `raw/download_manifest.csv` with source URLs, byte
counts, and MD5 checksums.

Download the raw files with

```bash
Rscript code/real_data/main.R dataset=eeg action=download ncores=10
```

`prepared/eegmmidb_r04__0003bb07.rds` is included in the repository. It is the
deterministic cache of the 8--30 Hz right-minus-left log-power contrasts used
for the manuscript analysis. Run the analysis directly from this cache with

```bash
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
```

Delete the prepared RDS file, set `use_cache=FALSE`, or run `action=all` to
force preprocessing from the downloaded EDF files after changing a
preprocessing parameter.
