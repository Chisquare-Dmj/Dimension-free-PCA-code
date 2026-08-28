# Manuscript Figures

- `simulation/` contains independent simulation PDF/PNG figures.
- `real_data/eegmmidb/` contains independent EEG spectral and inference figures.
- `real_data/capture24/` contains independent CAPTURE spectral, effective-rank,
  PVE, reliability, FPCA-versus-Proposed, and functional-curve figures.

Every PDF is one title-free statistical panel. Panel labels, captions, and
multi-panel layout belong in LaTeX. Every figure has a matching CSV under
`output/data/summary/` or `output/real_data/<dataset>/` and is indexed by
`output/figure_output_summary.csv`.

PNG figures are exported at 200 dpi on a compact single-panel canvas. Matching
PDF files remain vector graphics.
