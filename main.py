import sys

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_activity(path: str, primary: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "compound_name" not in df.columns:
        raise ValueError("CSV must contain a 'compound_name' column.")
    if primary not in df.columns:
        raise ValueError(f"Primary target column '{primary}' not found in CSV.")
    off_targets = [c for c in df.columns if c not in ("compound_name", primary)]
    if not off_targets:
        raise ValueError("CSV must contain at least one off-target column.")
    for col in [primary] + off_targets:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df[["compound_name", primary] + off_targets].copy()


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------

def compute_selectivity(df: pd.DataFrame, primary: str) -> pd.DataFrame:
    off_targets = [c for c in df.columns if c not in ("compound_name", primary)]
    result = df.copy()
    delta_cols = []
    for off in off_targets:
        dcol = f"delta_{off}"
        result[dcol] = result[primary] - result[off]
        result[f"SI_{off}"] = 10 ** result[dcol]
        delta_cols.append(dcol)
    result["selectivity_window"] = result[delta_cols].min(axis=1)
    result["mean_delta"] = result[delta_cols].mean(axis=1)
    result = result.sort_values(
        ["selectivity_window", "mean_delta", "compound_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return result


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_selectivity_heatmap(sel_df: pd.DataFrame, primary: str, output_path: str) -> None:
    off_targets = [
        c.replace("delta_", "")
        for c in sel_df.columns
        if c.startswith("delta_")
    ]
    heat_data = sel_df[[f"delta_{o}" for o in off_targets]].copy()
    heat_data.columns = off_targets
    heat_data.index = sel_df["compound_name"]

    fig, ax = plt.subplots(figsize=(8, max(8, len(sel_df) * 0.26)))
    sns.heatmap(
        heat_data,
        ax=ax,
        cmap="RdYlGn",
        center=0,
        vmin=-2,
        vmax=4,
        annot=True,
        fmt=".1f",
        linewidths=0.3,
        linecolor="white",
    )
    ax.set_xlabel("Off-Target", fontsize=11)
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title("Selectivity Profile: \u0394pIC50 vs Off-Targets",
                 fontsize=13, fontweight="bold", pad=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_selectivity_waterfall(
    sel_df: pd.DataFrame,
    primary: str,
    output_path: str,
    top_n: int = 20,
) -> None:
    off_targets = [
        c.replace("delta_", "")
        for c in sel_df.columns
        if c.startswith("delta_")
    ]
    n_off = len(off_targets)
    top = sel_df.head(top_n).reset_index(drop=True)
    n = len(top)

    colors = plt.cm.tab10.colors
    bar_h = 0.8 / n_off

    fig, ax = plt.subplots(figsize=(9, max(6, n * 0.35)))

    for j, off in enumerate(off_targets):
        dcol = f"delta_{off}"
        color = colors[j % len(colors)]
        for i, row in top.iterrows():
            y = i + (j - n_off / 2 + 0.5) * bar_h
            ax.barh(y, row[dcol], height=bar_h * 0.9,
                    color=color, edgecolor="white",
                    label=off if i == 0 else "")

    ax.axvline(1.0, color="gray", linestyle="--", linewidth=1.0)
    ax.axvline(2.0, color="gray", linestyle=":",  linewidth=1.0)
    ax.text(1.02, n - 0.5, "SI=10\u00d7", fontsize=8, color="gray", va="bottom")
    ax.text(2.02, n - 0.5, "SI=100\u00d7", fontsize=8, color="gray", va="bottom")

    ax.set_yticks(range(n))
    ax.set_yticklabels(top["compound_name"].tolist(), fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("\u0394pIC50 (primary \u2212 off-target)", fontsize=11)
    ax.set_title("Top Compounds: Selectivity Window",
                 fontsize=13, fontweight="bold")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="lower right", fontsize=9, framealpha=0.85)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute selectivity ratios (SI, \u0394pIC50) across multiple targets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input",      required=True,   help="Activity CSV (compound_name + target columns)")
    parser.add_argument("--primary",    required=True,   help="Primary target column name")
    parser.add_argument("--output-dir", default="output", help="Output directory")
    parser.add_argument("--top-n",      default=20, type=int, help="Compounds to show in waterfall")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_activity(args.input, args.primary)
    off_targets = [c for c in df.columns if c not in ("compound_name", args.primary)]
    print(f"  Loaded {len(df)} compounds; primary: {args.primary}; off-targets: {', '.join(off_targets)}")

    sel_df = compute_selectivity(df, args.primary)

    # --- CSV ---
    csv_path = os.path.join(args.output_dir, "selectivity_matrix.csv")
    sel_df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"\nSaved: {csv_path}")

    # --- Heatmap ---
    heatmap_path = os.path.join(args.output_dir, "selectivity_heatmap.png")
    plot_selectivity_heatmap(sel_df, args.primary, heatmap_path)
    print(f"Saved: {heatmap_path}")

    # --- Waterfall ---
    waterfall_path = os.path.join(args.output_dir, "selectivity_waterfall.png")
    plot_selectivity_waterfall(sel_df, args.primary, waterfall_path, top_n=args.top_n)
    print(f"Saved: {waterfall_path}")

    # --- Console summary ---
    print(f"\n--- Top 10 compounds by selectivity window ---")
    top10_cols = ["compound_name", "selectivity_window", "mean_delta", args.primary]
    print(sel_df[top10_cols].head(10).to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    n_10  = (sel_df["selectivity_window"] >= 1.0).sum()
    n_100 = (sel_df["selectivity_window"] >= 2.0).sum()
    print(f"\n  {n_10}  compounds with selectivity_window >= 1.0 (SI > 10x against worst off-target)")
    print(f"  {n_100}  compounds with selectivity_window >= 2.0 (SI > 100x against worst off-target)")
    print("\nDone.")


if __name__ == "__main__":
    main()
