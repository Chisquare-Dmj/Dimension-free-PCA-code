# Manuscript Figures

- `simulation/` contains the simulation PDF/PNG figure files.
- `real_data/eegmmidb/` contains the EEG spectral and inference figures.
- `real_data/capture24/` contains the CAPTURE-24 spectral, effective-rank, PVE,
  reliability, eigenvalue-comparison, and functional-curve figures.

Each PDF is a title-free statistical panel. Multi-panel layout, panel labels,
and captions are handled in the manuscript source. PNG copies are retained for
quick inspection; PDF files remain vector graphics.

The corresponding compact plot-data CSV files are kept under
`output/data/summary/` or `output/real_data/<dataset>/`. Full run-specific
artifacts and Monte Carlo replication-level files are regenerated locally and
are not committed to the public repository.
