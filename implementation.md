# Phase 26 — Selectivity Ratio Calculator

**Version:** 1.1 (final as-built)
**Author:** Kerwyn Medrano
**Date:** 2026-03-26
**Track:** Track 1 — Cheminformatics Core
**Tier:** Micro (2–3 hrs)
**API Cost:** $0.00 — pure pandas + matplotlib + seaborn

---

## 1. Project Overview

### Goal

Given a multi-target pIC50 table (compounds × targets), compute
selectivity ratios — how much more potent each compound is on the
primary target vs each off-target — and visualize the selectivity
landscape as a ranked heatmap and waterfall.

```bash
python main.py --input data/activity.csv --primary CETP
```

Outputs:
- `output/selectivity_matrix.csv` — per-compound SI and ΔpIC50 for each off-target
- `output/selectivity_heatmap.png` — heatmap: compounds (rows, ranked by primary SI) × targets (cols); values = ΔpIC50
- `output/selectivity_waterfall.png` — waterfall bar chart: ΔpIC50 vs off-target for top-N compounds

### Distinction from Phase 14

Phase 14 (Selectivity Profiling) built radar charts and a broad SI analysis
on a single fixed dataset. Phase 26 is a reusable **calculator**: any CSV,
any primary target, all off-targets auto-detected, log-scale ratios, ranked
output table. It is a utility module, not an exploratory analysis.

### What This Phase Teaches

| Concept | Detail |
|---|---|
| Selectivity ratio (SI) | `SI = 10^(pIC50_primary - pIC50_off)` = IC50_off / IC50_primary |
| ΔpIC50 | Additive form of SI on log scale; directly interpretable |
| Selectivity window | Range of ΔpIC50 values across all off-targets for a compound |
| Wide-to-long pivot | `pd.melt` for off-target unpacking |
| Multi-target ranking | Sort compounds by mean ΔpIC50 across all off-targets |

---

## 2. Architecture

```
selectivity-ratio/
├── main.py
├── requirements.txt
├── README.md
├── .gitignore
├── data/
│   └── activity.csv          — compounds × targets pIC50 table
└── output/
    ├── selectivity_matrix.csv
    ├── selectivity_heatmap.png
    └── selectivity_waterfall.png
```

---

## 3. Input Format

### `data/activity.csv`

```csv
compound_name,CETP,LPL,HL,EL
benz_001_F,7.25,5.10,4.80,5.50
benz_002_CN,7.40,5.30,4.90,5.60
...
```

- First column: `compound_name`
- Remaining columns: one per target (pIC50 values)
- `--primary` flag selects which column is the on-target
- All other columns are treated as off-targets
- Missing values (NaN) skipped per compound-target pair

### Synthetic dataset: 45 compounds, 4 targets

Reuse Phase 25 compound names. Generate 3 off-target pIC50 columns
(LPL, HL, EL) synthetically: each = primary pIC50 - ΔpIC50, where
ΔpIC50 drawn from N(1.5, 0.5) for typical selectivity in CETP series.
A few compounds (bzim family) should show poor selectivity (ΔpIC50 < 0.5).

---

## 4. Calculations

### ΔpIC50 (primary metric)
```python
delta = df[primary] - df[off_target]
```
Positive = selective (primary more potent than off-target).

### Selectivity Index (SI)
```python
SI = 10 ** delta
```
SI > 10 (ΔpIC50 > 1.0) = meaningfully selective.
SI > 100 (ΔpIC50 > 2.0) = highly selective.

### Selectivity window
```python
window = df_delta.min(axis=1)  # worst-case off-target selectivity
```

### Ranking
Sort compounds by `selectivity_window` descending.

---

## 5. Module Specification

### `load_activity(path, primary)` → pd.DataFrame
- Read CSV; validate `primary` column exists
- Identify off-target columns = all columns except `compound_name` and `primary`
- Return DataFrame with compound_name, primary pIC50, and off-target columns

### `compute_selectivity(df, primary)` → pd.DataFrame
- For each off-target: compute ΔpIC50 and SI
- Add `selectivity_window` = min ΔpIC50 across all off-targets
- Add `mean_delta` = mean ΔpIC50 across all off-targets
- Sort by `selectivity_window` descending

### `plot_selectivity_heatmap(sel_df, primary, output_path)`
- seaborn heatmap; rows = compounds ranked by selectivity_window
- Columns = off-targets; cell values = ΔpIC50
- Diverging colormap (RdYlGn): green = selective, red = non-selective
- Center at 0 (= equipotent on primary and off-target)
- Annotate cells with ΔpIC50 (1 decimal)
- Title: "Selectivity Profile: ΔpIC50 vs Off-Targets"

### `plot_selectivity_waterfall(sel_df, primary, output_path, top_n=20)`
- Horizontal grouped bar for top-N compounds by selectivity_window
- One bar per off-target per compound (grouped by compound)
- Color = off-target; bar width = ΔpIC50
- Reference line at ΔpIC50 = 1.0 (SI = 10×)
- Reference line at ΔpIC50 = 2.0 (SI = 100×)
- Title: "Top Compounds: Selectivity Window"

### `main()`
- `--input` (required): activity CSV
- `--primary` (required): primary target column name
- `--output-dir` (default: output)
- `--top-n` (default: 20): compounds to show in waterfall
- Print: top-10 compounds by selectivity_window; counts at SI>10 and SI>100

---

## 6. Expected Results

For the synthetic dataset (ΔpIC50 ~ N(1.5, 0.5)):
- ~30/45 compounds with selectivity_window > 1.0 (SI > 10)
- ~5–10 compounds with selectivity_window > 2.0 (SI > 100)
- bzim family: lowest selectivity (bzim compounds have PosIonizable feature,
  hinting at off-target activity)

---

## 7. Verification Checklist

```bash
python main.py --input data/activity.csv --primary CETP

# Expected:
# - selectivity_matrix.csv: 45 rows, columns: compound_name, CETP,
#   delta_LPL, SI_LPL, delta_HL, SI_HL, delta_EL, SI_EL,
#   selectivity_window, mean_delta
# - selectivity_heatmap.png: green dominant (good selectivity); red patches
#   in bzim rows
# - selectivity_waterfall.png: top 20 compounds; reference lines at 1.0 and 2.0
# - Console: top-10 table; "N compounds with SI > 10×" count
```

---

## 8. Risks / Assumptions / Next Step

**Risks:**
- If a compound has NaN for a target, skip that target in its window/mean calc
- Diverging colormap center=0: use `sns.heatmap(..., center=0, vmin=-2, vmax=4)`

**Assumptions:**
- pIC50 units; ΔpIC50 > 0 = compound is more potent on primary
- Synthetic off-target data drawn from N(primary - 1.5, 0.5) per compound

---

## 9. Actual Results (v1.1)

### Run
```
Loaded 45 compounds; primary: CETP; off-targets: LPL, HL, EL
```

### Top 10 compounds by selectivity window

| compound_name | selectivity_window | mean_delta | CETP |
|---|---|---|---|
| benz_002_Cl | 2.14 | 2.26 | 7.46 |
| benz_008_OMe | 2.02 | 2.20 | 7.01 |
| naph_002_F | 2.00 | 2.08 | 7.44 |
| ind_001_H | 1.84 | 2.00 | 7.62 |
| pyr_003_F | 1.74 | 2.00 | 6.80 |
| naph_004_Me | 1.71 | 1.90 | 6.77 |
| benz_001_F | 1.66 | 2.00 | 7.65 |
| quin_004_Me | 1.65 | 1.89 | 7.01 |
| benz_003_Br | 1.56 | 1.88 | 7.65 |
| benz_006_NO2 | 1.53 | 2.11 | 7.47 |

### Counts
- 39 / 45 compounds with selectivity_window >= 1.0 (SI > 10×)
- 3 / 45 compounds with selectivity_window >= 2.0 (SI > 100×)

### Key insights
- benz/naph compounds dominate the top — substituent effects (Cl, F, Br) boost both CETP potency and selectivity window
- bzim family clusters at the bottom as designed (poor selectivity, window ~0.2–0.5)
- 86% of compounds (39/45) achieve meaningful selectivity (SI > 10×) — consistent with the synthetic delta distribution

### Deviations from plan
- None — pipeline ran cleanly on first attempt

**Next step:** Phase 27 — Matched Molecular Series (MMS) Extractor
