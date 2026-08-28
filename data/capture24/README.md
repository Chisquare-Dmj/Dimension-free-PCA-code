# CAPTURE-24 Wake-Aligned ENMO Application

The CAPTURE-24 source archive is obtained from the Oxford University Research
Archive. The large source ZIP is intentionally not versioned in this public
repository; `code/real_data/main.R` downloads it when a raw-data rebuild is
requested.

The active construction is prespecified:

- compute ENMO from released triaxial acceleration;
- average in non-overlapping five-second epochs and apply `log1p`;
- infer the main sleep episode from released labels, joining interruptions no
  longer than 20 minutes;
- circularly align time zero at the end of that episode while preserving
  acquisition gaps as missing bins;
- retain an explicit coverage mask and require at least one hour in the main
  sleep episode for the primary wake anchor;
- fill missing epochs with their pointwise cohort mean so that the missing
  contribution is zero after across-participant centering;
- do not center, scale, or residualize individual participant curves before
  the common sample centering used by PCA.

`processed/cfg0003e622/` contains the deterministic public cache used by the
manuscript: the wake-aligned five-second log-ENMO matrix, coverage mask, subject
QC, metadata, label dictionary, and preprocessing manifest. The primary cohort
contains 150 participants after the prespecified wake-anchor exclusion.

Analyze the included cache without redownloading the raw archive with:

```bash
Rscript code/real_data/main.R dataset=capture24 action=analyze ncores=10
```

A complete raw-data rebuild is available with `action=all`, but it is
substantially more expensive because the source archive is several gigabytes.
