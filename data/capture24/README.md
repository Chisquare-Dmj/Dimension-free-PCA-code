# CAPTURE-24 Wake-Aligned ENMO Application

The public CAPTURE-24 source archive is not committed to this repository because
of its size. The analysis code downloads the official archive from the Oxford
University Research Archive, verifies it, and streams the compressed
participant files during preprocessing.

Download the archive with

```bash
Rscript code/real_data/main.R dataset=capture24 action=download
```

The deterministic processed cache used by the manuscript is included under
`processed/cfg0003e622/`. It contains

- the 151-by-17280 five-second wake-aligned `log1p` ENMO matrix;
- the corresponding coverage mask;
- participant QC and metadata;
- the released activity-label dictionary; and
- the preprocessing manifest.

The preprocessing computes ENMO from the released triaxial acceleration,
averages in non-overlapping five-second epochs, applies `log1p`, identifies the
main sleep episode from the released sleep labels, and circularly aligns time
zero at the end of that episode. Missing epochs are retained through an
explicit coverage mask and filled with their pointwise cohort mean before the
across-participant centering used in PCA.

One participant does not satisfy the prespecified one-hour main-sleep rule, so
the primary analysis uses `n=150`. The manuscript considers nested 6-, 12-,
18-, and 24-hour wake-aligned domains, with primary `K0=8`, K0 sensitivity
`6, 8, 10`, a missingness sensitivity restricted to participants with at most
30 minutes missing, and a main-sleep sensitivity requiring at least three
hours.

Analyze the included processed cache with

```bash
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10
```

Use `action=all` for a complete download, preprocessing, analysis, summary, and
plotting rebuild from the public source archive.
