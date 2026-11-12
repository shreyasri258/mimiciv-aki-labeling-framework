# mimiciv-aki-labeling-framework
This is a reproducible KDIGO acute kidney injury (AKI) labeling for MIMIC-IV v3.1 with a 24-hour landmark design. It includes label-level sensitivity analysis for 23 implementation variants that are compared against the mimiciv_derived implementation.

## 1. Overview

The KDIGO definition of AKI leaves several implementation details such as the baseline creatinine, the body weight used for urine output, the handling of gaps in urine charting, and the treatment of stays that cannot be evaluated open. 

This repository:
specifies one complete labeling rule, V5, for adult first ICU stays in MIMIC-IV v3.1.
- compares V5 with the mimiciv_derived KDIGO tables. It is treated as a comparison implementation and not as ground truth.
- evaluates 23 variants that each change a few  implementation choices, reporting label flips on common stays and the paired change in Cohen's kappa with stay-level bootstrap confidence intervals.
- reports evaluability and the kind of evidence behind each label.


**Design in brief**

| Item | Definition |
|---|---|
| Source population | Adult first ICU stays with ICU length of stay of at least 24 h |
| Time zero | ICU admission |
| Early window | [0, 24) h. Stays with KDIGO-defined AKI here are excluded |
| Label window | [24, 48) h. A stay is AKI-positive if the rule is met in this window |
| Windows | Half-open throughout |

---
## 2. Repository structure

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── src/
│   
├── code/
│   ├── Cohort Construction.ipynb# cohort construction
│   ├── AKI Labeling and Outcome Prevalance.ipynb# V5 labeling and comparison with mimiciv_derived
│   ├── AKI Label Sensitivity Analysis.ipynb # 24 labelings, paired bootstrap, registry
|   |── kdigo_pipeline.py # labeling functions shared by all notebooks
│   └── tables_figures.ipynb #generated figures and tables
results
├── results/
│   ├── tables/# aggregate tables (CSV)
│   └── figures/# final figures (PNG, SVG)

```
*Patient-level files (labels, cohort, caches) are produced locally and are not part of this repository.*

---
## 3. Data access
This work uses **MIMIC-IV version 3.1** (PhysioNet, released October 11, 2024; data 2008–2022).

1. Complete the PhysioNet credentialing process and sign the data use agreement for MIMIC-IV.
2. Download MIMIC-IV v3.1 from PhysioNet.
3. Build the `mimiciv_derived` concepts with the [MIMIC code repository](https://github.com/MIT-LCP/mimic-code) into a DuckDB database (`mimic-iv.duckdb`).
4. Export the derived tables used here as CSV into `Mimic/derived/`.

**Expected local layout** (paths are relative to the project root. They should not be renamed):

```
Mimic/
├── hosp/        # labevents.csv, admissions.csv, patients.csv, diagnoses_icd.csv, ...
├── icu/         # outputevents.csv, chartevents.csv, d_items.csv, icustays.csv, ...
└── derived/     # kdigo_stages.csv, kdigo_creatinine.csv, kdigo_uo.csv,
                 # weight_durations.csv, crrt.csv
mimic-iv.duckdb  # only needed to export mimiciv_derived.crrt if crrt.csv is missing
processed_outputs/   #outputs created by the notebooks
```
**Data are not distributed with this repository.** To reproduce the analysism data must be obtained from MIMIC-IV through PhysioNet.

---

## 4. Environment

- Python **3.13.9**
- Main packages: `pandas`, `numpy`, `scikit-learn`, `duckdb`, `matplotlib`, `pyarrow`, `openpyxl`
- Versions are mentioned in `requirements.txt`

## 5. Steps to Run

Run from the project root so that relative paths resolve.

| Step | Code / script | Produces |
|---|---|---|
| 1 | `code/Cohort Construction.ipynb` | `processed_outputs/cohort.parquet` |
| 2 | `code/AKI Labeling and Outcome Prevalance.ipynb` | V5 labels (`labels_v5.parquet`), cohort flow, comparison with `mimiciv_derived`, `run_log_primary.json` |
| 3 | `code/AKI Label Sensitivity Analysis .ipynb` | variant registry, summary and paired kappa tables, prognostic summary |
| 4 | `code/tables_figures.py` | manuscript tables and figures in `processed_outputs/Tables_and_Figures/` |

**Notes**

- AKI Label Sensitivity Analysis .ipynb uses a cached results file. With `RERUN = False` (default) it loads the cache and **recomputes nothing from MIMIC**. Set `RERUN = True` to recompute all 24 labelings (this is slow and needs the full raw files).
- AKI Labeling and Outcome Prevalance.ipynb freezes the V5 labels the first time it runs. On later runs it checks the rebuilt labels against the frozen file and stops if they differ.
- tables_figures.py reads uses saved results.

---

## 6. Expected results

A successful reproduction gives these values for the V5 reference run.

**Cohort flow**

| Step | Stays |
|---|---|
| ICU stays | 94,458 |
| First ICU stays (one per patient) | 65,366 |
| Adult first stays, length of stay ≥ 24 h | 51,838 |
| Early AKI in [0, 24) h (excluded) | 22,829 (44.04%) |
| After early exclusion | 29,009 |
| No baseline creatinine and urine not evaluable (excluded) | 32 |
| **Analysis cohort** | **28,977** |

**V5 labels in [24, 48) h**

| Quantity | Value |
|---|---|
| AKI-positive | 7,624 (26.31%) |
| AKI-negative | 21,353 |
| Urine-evaluable | 17,045 |
| Not urine-evaluable | 11,932 |
| Stays with a reference label | 26,610 (91.8%) |

**Agreement of V5 with `mimiciv_derived` (n = 26,610)**

| Criterion | Cohen's κ |
|---|---|
| Combined | 0.6234 |
| Creatinine | 0.8548 |
| Urine output | 0.5779 |

Confusion matrix against the comparison implementation: TP 6,405, FP 1,082, FN 3,343, TN 15,780 (sensitivity 0.657, specificity 0.936, PPV 0.855, NPV 0.825).

**Frozen file checksums (reference run)**

| File | SHA-256 |
|---|---|
| `labels_v5.parquet` | `8cac4a06d959973981a1f6e905830ba78361db0d464374c06a02efe5e6265c07` |
| `kdigo_pipeline.py` | `50886662207fe2f4a3bd1a29ca8404a4b7100798164f6be9653943d8b928d077` |

A parquet file hash can differ across `pandas` and `pyarrow` versions even when the data are identical.

---

## 7. Reproducibility

- Fixed seeds with bootstrap seed 42, 1,000 resamples is used.
- `RERUN = False` uses the saved results. The cache is not distributed because it contains patient-level labels.
- The V5 labels are written once, with a checksum, and checked on later runs.
- Tables and figures use data from `results/tables/` without patient-level data.

---