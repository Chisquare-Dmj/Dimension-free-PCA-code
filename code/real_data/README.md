# Real-Data Code

- `common_real_data.R`: centered Gram PCA, classical FPCA and Proposed
  adjacent-eigengap inference, prespecified candidate-cluster output, K0
  sensitivity, tables, and publication-figure helpers.
- `eegmmidb.R`: PhysioNet R04 download, EDF+ event parsing, 8--30 Hz log-power
  preprocessing, subject contrasts, lateralization direction, and analysis.
- `capture24.R`: CAPTURE-24 archive handling, wake-aligned five-second log-ENMO
  construction, nested-domain diagnostics, sensitivity analyses, and inference.
- `capture24_preprocess.py`: memory-bounded deterministic conversion of raw
  accelerometry to the processed cache; it does not implement PCA inference.

Normal usage is through `main.R`:

```bash
Rscript code/real_data/main.R dataset=eeg action=analyze ncores=10
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10
```

The public repository includes processed caches, so `action=analyze` does not
redownload the raw data. Use `action=all` only for a complete source-data rebuild.

All PC1--PC6 estimates are saved independently of adjacent equality-test
outcomes. A failure to reject equality means only that the data provide
insufficient evidence to distinguish the two population roots; it does not
automatically create a cluster. The only prespecified non-overlapping candidate
clusters are PC3--PC4 and PC5--PC6 for CAPTURE-24, and PC3--PC4 for EEG.

For classical FPCA, near-tie pairs retain the fixed-covariance equality p-value
but the regular distinct-root Wald eigengap interval is suppressed in the
manuscript-facing output and marked `not valid under near-tie`. Proposed
component and phase diagnostics are not gated by FPCA decisions. The older 0.05
descriptive relative-gap rule survives only in legacy/run-specific columns.

EEG additionally reports matched 500-draw percentile intervals for the raw FPCA
and reliability-corrected Proposed PC2 lateralization projections. These
bootstrap intervals are an empirical comparison and are not consequences of the
eigengap limit theory.
