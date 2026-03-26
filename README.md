# selectivity-ratio

**Phase 26 — Selectivity Ratio Calculator**
Track 1 — Cheminformatics Core

Computes per-compound selectivity ratios (SI, ΔpIC50) from a multi-target
pIC50 table and visualizes the selectivity landscape as a ranked heatmap
and grouped waterfall chart.

## Selectivity metrics

| Metric | Formula | Interpretation |
|---|---|---|
| ΔpIC50 | pIC50(primary) − pIC50(off-target) | Log-scale selectivity; positive = selective |
| SI | 10^(ΔpIC50) = IC50_off / IC50_primary | Fold-selectivity ratio |
| selectivity_window | min ΔpIC50 across all off-targets | Worst-case selectivity for that compound |
| mean_delta | mean ΔpIC50 across all off-targets | Average selectivity profile |

**Thresholds:** ΔpIC50 ≥ 1.0 (SI ≥ 10×) = meaningful selectivity;
ΔpIC50 ≥ 2.0 (SI ≥ 100×) = high selectivity.

## Quickstart

```bash
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt

PYTHONUTF8=1 .venv\Scripts\python main.py --input data/activity.csv --primary CETP
```

## CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input` | required | CSV with compound_name + target pIC50 columns |
| `--primary` | required | Column name for the on-target |
| `--output-dir` | `output` | Output directory (created at runtime) |
| `--top-n` | `20` | Compounds to show in the waterfall chart |

## Off-target auto-detection

Off-targets = all columns except `compound_name` and `--primary`.
No configuration needed; works with any number of off-target columns.

## Outputs

| File | Description |
|---|---|
| `output/selectivity_matrix.csv` | Per-compound delta/SI for each off-target + window + mean |
| `output/selectivity_heatmap.png` | Heatmap ranked by selectivity_window; RdYlGn diverging at 0 |
| `output/selectivity_waterfall.png` | Grouped bar for top-N; reference lines at ΔpIC50 = 1.0 and 2.0 |

## Plot interpretation

- **Heatmap**: green cells = selective (ΔpIC50 > 0); red = non-selective or inverse.
  Compounds sorted top-to-bottom by selectivity_window (best first).
- **Waterfall**: each bar = ΔpIC50 against one off-target. Dashed line at 1.0 = SI 10×
  threshold; dotted line at 2.0 = SI 100× threshold.

## Notes

- Project-local `.venv` — do not share across projects.
- `output/` created at runtime and gitignored.
- Mandatory import order: Windows UTF-8 fix → `matplotlib.use("Agg")` → pyplot.
