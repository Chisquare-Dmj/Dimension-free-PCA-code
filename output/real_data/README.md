# Real-Data Numerical Output

This directory contains the compact, manuscript-facing numerical outputs retained
in the public repository. Full run-ID-prefixed artifacts and real-data manifests
are generated locally by the analysis pipeline and are intentionally not
committed.

The retained files include:

- `eegmmidb/eegmmidb_spectral_inference.csv`: leading-component estimates and
  separation/phase diagnostics for the EEG analysis;
- `eegmmidb/eegmmidb_k0_sensitivity.csv`: EEG estimates across the requested
  `K0` values;
- `eegmmidb/eeg_fpca_vs_proposed.csv` and
  `eegmmidb/eeg_pc2_lateralization_fpca_vs_proposed.csv`: conventional-FPCA and
  corrected summaries used in the manuscript comparison;
- `capture24/capture24_spectral_inference.csv`: leading-component inference for
  the 24-hour CAPTURE-24 analysis;
- `capture24/capture24_k0_sensitivity.csv`: CAPTURE-24 `K0` sensitivity;
- `capture24/capture24_domain_analysis.csv`, `capture24_effective_rank.csv`,
  `capture24_pc1_pve.csv`, and `capture24_pc_reliability.csv`: the retained
  large-domain summaries;
- `capture24/capture24_missingness_sensitivity.csv` and
  `capture24/capture24_sleep_duration_sensitivity.csv`: preprocessing
  sensitivity analyses; and
- `capture24/capture24_fpca_vs_proposed.csv`: the descriptive comparison of
  conventional sample FPCA quantities with the corrected spike/reliability
  analysis.

`figure_output_manifest.csv` indexes the retained real-data figure outputs and
plot-data files. The repository-level `output/figure_output_summary.csv` combines
these entries with the simulation figure index.

Because population truth is unavailable for the two real datasets, the retained
FPCA comparisons are descriptive and are not presented as empirical bias or
accuracy assessments.
