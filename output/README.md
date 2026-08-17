# Output

This public repository keeps compact manuscript-facing numerical outputs while
omitting large or machine-specific generated artifacts.

Included:

- `data/summary/`: compact simulation plot/table data used for the manuscript;
- `real_data/eegmmidb/`: compact EEG analysis outputs;
- `real_data/capture24/`: compact CAPTURE-24 analysis outputs;
- `figure_output_summary.csv`: index of manuscript figures and their plot-data files;
- `consistency_report.csv`: numerical consistency checks for the retained outputs.

Generated locally and not versioned:

- `data/replicate/`: Monte Carlo replication-level results;
- `data/truth/`: population truth and configuration snapshots;
- `data/manifest/`: run registry and latest-run pointers;
- `table/`: generated LaTeX tables;
- `log/`: scheduler or batch logs;
- `real_data/manifest/`: real-data run registry;
- run-specific `real_*` real-data artifacts with full configuration fingerprints.

The excluded files are recreated automatically by the simulation or real-data
pipelines. Their omission keeps the public repository small and avoids very long
generated paths on Windows. Compact final summaries and figures are retained so
that manuscript values can be inspected without rerunning the full Monte Carlo
program.
