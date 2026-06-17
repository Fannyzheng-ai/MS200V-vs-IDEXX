# Data Dictionary

## `paired_results.csv` — paired MS200V vs IDEXX measurements

One row per (patient × analyte) pair. 1,374 rows total.

| Column | Type | Description |
|---|---|---|
| `PatientID` | string | Anonymous patient ID in the form `P###` (e.g., `P001`, `P014`). |
| `Species` | string | `Canine` or `Feline`. |
| `Age_years` | float (or empty) | Patient age at sample collection, in years. Empty when not recorded. |
| `Analyte` | string | One of: ALB, ALP, AMY, Ca, CHOL, CREA, GLU, PHOS, TBIL, TP, UREA. |
| `Unit` | string | Unit of measurement (see *Units* below). |
| `IDEXX` | float | Result from the IDEXX Catalyst One reference analyzer. |
| `MS200V` | float | Result from the MS200V. For UREA, values have been converted from raw urea (mmol/L) to BUN-equivalent mmol/L by dividing by 2.14. |
| `AbsoluteDiff` | float | `MS200V − IDEXX`. |
| `PercentBias` | float (or empty) | `100 × (MS200V − IDEXX) / IDEXX`. Empty when `IDEXX = 0`. |

### Units

| Analyte | Full name | Unit |
|---|---|---|
| ALB  | Albumin | g/L |
| ALP  | Alkaline phosphatase | U/L |
| AMY  | Amylase | U/L |
| Ca   | Calcium (total) | mmol/L |
| CHOL | Cholesterol | mmol/L |
| CREA | Creatinine | µmol/L |
| GLU  | Glucose | mmol/L |
| PHOS | Inorganic phosphate | mmol/L |
| TBIL | Total bilirubin | µmol/L |
| TP   | Total protein | g/L |
| UREA | Blood urea nitrogen | mmol/L (BUN-equivalent) |

### Notes

- UREA values are in BUN-equivalent mmol/L (MS200V's native urea reading divided by 2.14 prior to deposit).

---

## `patient_demographics.csv` — anonymized patient roster

| Column | Type | Description |
|---|---|---|
| `PatientID` | string | Joins with `paired_results.csv`. |
| `Species` | string | Canine or Feline. |
| `Age_years` | float (or empty) | Patient age. |

---

## `imprecision_qc.csv` — MS200V QC imprecision

| Column | Type | Description |
|---|---|---|
| `Analyte` | string | Analyte abbreviation. |
| `Level` | string | QC concentration level: `L1` (within reference interval) or `L2` (outside). |
| `N` | int | Number of replicate measurements (typically 18). |
| `Target` | float | Manufacturer-assigned target concentration for the QC material. |
| `Observed_Mean` | float | Observed mean across replicates. |
| `Observed_SD` | float | Observed standard deviation. |
| `CV_pct` | float | Observed coefficient of variation (%). |
| `ASVCP_TEa_pct` | float | ASVCP-recommended total allowable error (%). |
| `Allowable_CV_pct` | float | Allowable CV = ASVCP_TEa_pct / 1.65. |
| `PassFail` | string | `Pass` if observed CV ≤ allowable CV. |

---

## `expected_outputs.csv` — reference output

Reference values produced by `code/analysis_reproducibility.py`. Running the script on a fresh checkout should reproduce these values to within numerical rounding.
