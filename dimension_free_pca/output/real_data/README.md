# Real-Data Numerical Output

`eegmmidb/` and `capture24/` contain complete run-ID-prefixed deterministic
artifacts. `manifest/` maps each run to its inference, spectrum, sensitivity,
QC, metadata, configuration, table, and figure inputs.

Common files include:

- `__spectral_inference.csv`: point estimates, theoretical uncertainty, and
  separation/phase diagnostics.
- `__k0_sensitivity.csv`: observable estimates under all requested K0 values.
- `__config.csv` and `__metadata.json`: exact computational and data settings.
- `__spectral_inference.tex`: manuscript-ready spectral inference table.

CAPTURE-24 additionally saves domain-specific spectra, effective-rank results,
participant QC, principal functions, the 30-minute missingness sensitivity, the
three-hour main-sleep sensitivity, and near-cluster diagnostics.
`capture24_fpca_vs_proposed.csv` reports PC1--PC2 sample FPCA eigenvalues/PVE
and Proposed spike, phase-margin, reliability, gap, and correction magnitudes.

`eegmmidb/eeg_fpca_vs_proposed.csv` and
`eeg_pc2_lateralization_fpca_vs_proposed.csv` contain the two deterministic PC2
lateralization point estimates. Real-data files do not report empirical bias or
relative accuracy because population truth is unknown.
