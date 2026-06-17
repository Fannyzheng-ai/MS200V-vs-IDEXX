# MS200V vs IDEXX Catalyst One — Method Comparison Data

This repository contains the de-identified data and reproducibility code for:

> **Analytical Performance of a Disk-Based Microfluidic Veterinary Biochemistry Analyzer: A Multi-Site Comparison Study with the IDEXX Catalyst One**
>
> Zheng Y, Kim YJ, Chavarria-Marron A, Clark D, Lee L, Rakovski C, Yaghmaei E, Zheng J, Hao J.


All patient identifiers have been replaced with anonymous study IDs of the form `SITE-NNN` (e.g., `XIN-001`, `WES-014`). No personally identifiable information about clients or animals is retained.

---

## Contents

```
.
├── README.md                       — this file
├── LICENSE                         — CC BY 4.0 (data) / MIT (code)
├── data/
│   ├── paired_results.csv          — master long-format CSV (1,374 paired MS200V/IDEXX measurements)
│   ├── ALB_paired.csv … UREA_paired.csv — per-analyte subsets (convenience)
│   ├── patient_demographics.csv    — anonymized patient roster
│   ├── imprecision_qc.csv          — QC imprecision (n = 18 replicates × 2 levels × 11 analytes)
│   ├── expected_outputs.csv        — reference output to verify reproduction
│   └── data_dictionary.md          — column definitions and units
└── code/
    ├── analysis_reproducibility.py — Python script to reproduce all main tables
    ├── analysis_reproducibility.R  — equivalent R script
    └── generate_figures.py         — regenerates Figures 1–12 (PDF/PNG/SVG/TIFF)
```

## Quick start

### Python

```bash
pip install pandas numpy scipy matplotlib pillow
python code/analysis_reproducibility.py
python code/generate_figures.py
```

### R

```r
install.packages(c("dplyr", "readr", "tibble"))
source("code/analysis_reproducibility.R")
```

## Data summary

- **11 analytes**: ALB, ALP, AMY, Ca, CHOL, CREA, GLU, PHOS, TBIL, TP, UREA
- **1,374** paired MS200V vs IDEXX Catalyst One measurements
- 4 clinical sites in 2 countries

## Citation

If you use this data, please cite the manuscript above and this repository.

## License

- **Data** — Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Code** — MIT License

See `LICENSE` for full terms.

## Contact

Corresponding author: Jijun Hao, jhao@westernu.edu
