#!/usr/bin/env python3
"""
Describe one sample (or several, pooled): per-cell summaries, summary tables,
redundant features, and per-cell distribution plots.

Outputs (in --out_dir):
- per_cell.csv          one row per cell: object count and the median of each feature
- summary_objects.csv   median and IQR of each feature over all objects
- summary_cells.csv     median and IQR of the per-cell medians, and the
                        coefficient of variation between cells
- redundant_pairs.csv   feature pairs with |Spearman rho| above --corr_threshold
- per_cell_<feature>.svg/.png   every object, grouped by cell, with cell medians

Example:
    python describe_sample.py --objects objects_filtered.csv --out_dir results/description \
        --plot_features volume_um3 sphericity
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_FEATURES = ["volume_um3", "surface_area_um2", "sphericity", "elongation",
                    "flatness", "max_diameter_um"]

# Stanford palette
COLOR = "#8C1515"   # Cardinal Red
INK = "#2E2D29"     # Black

LABELS = {
    "volume_um3": "Object volume (µm³)",
    "surface_area_um2": "Surface area (µm²)",
    "max_diameter_um": "Maximum 3D diameter (µm)",
    "sphericity": "Sphericity",
    "elongation": "Elongation",
    "flatness": "Flatness",
    "aspect_ratio": "Aspect ratio (major / minor axis)",
}


def is_size_feature(name: str) -> bool:
    """Size features (volume, area, length) are right-skewed: plot on a log axis."""
    return name.endswith(("_um3", "_um2", "_um")) or "Volume" in name or "Area" in name \
        or "Diameter" in name or "AxisLength" in name


def set_log_axis(ax):
    """Log y axis with plain-number labels; more labels when the range is narrow."""
    ax.set_yscale("log")
    ax.autoscale_view()
    lo, hi = ax.get_ylim()
    decades = np.log10(hi / lo) if lo > 0 else 2
    subs = (1, 2, 5) if decades > 1.5 else (1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9)
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=subs))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: np.format_float_positional(float(f"{v:.3g}"), trim="-")))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def iqr(x):
    return x.quantile(0.75) - x.quantile(0.25)


def plot_per_cell(objects, cell_medians, feat, out_dir):
    order = cell_medians[feat].sort_values().index
    fig, ax = plt.subplots(figsize=(max(4.0, 0.35 * len(order) + 1.2), 3.2), layout="constrained")
    rng = np.random.default_rng(0)
    for i, cell in enumerate(order):
        v = objects.loc[objects["cell_uid"] == cell, feat].dropna()
        ax.scatter(i + rng.uniform(-0.25, 0.25, len(v)), v, s=8, color=COLOR,
                   alpha=0.4, linewidths=0)
        ax.hlines(v.median(), i - 0.35, i + 0.35, color=INK, lw=2, zorder=3)

    ax.set_xticks(range(len(order)), [str(c) for c in order], rotation=90 if len(order) > 12 else 0,
                  fontsize=7 if len(order) > 12 else None)
    ax.set_xlabel("Cell (sorted by median)", color=INK)
    if is_size_feature(feat):
        set_log_axis(ax)
    ax.set_ylabel(LABELS.get(feat, feat), color=INK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=INK)
    fig.supxlabel("Dots: objects · line: median per cell", fontsize=7, color=INK)
    base = os.path.join(out_dir, f"per_cell_{feat}")
    fig.savefig(base + ".svg")
    fig.savefig(base + ".png", dpi=600)
    plt.close(fig)
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--objects", default="objects_filtered.csv", help="Object table (from filter_objects.py)")
    ap.add_argument("--out_dir", default="results/description", help="Output folder")
    ap.add_argument("--features", nargs="+", default=DEFAULT_FEATURES,
                    help="Features to summarize (column names in the object table)")
    ap.add_argument("--plot_features", nargs="*", default=["volume_um3"],
                    help="Features to plot per cell (default: volume_um3)")
    ap.add_argument("--corr_threshold", type=float, default=0.9,
                    help="Report feature pairs with |Spearman rho| above this (default: 0.9)")
    args = ap.parse_args()

    objects = pd.read_csv(args.objects)
    missing = [f for f in args.features + args.plot_features if f not in objects.columns]
    if missing:
        raise SystemExit(f"[error] features not in {args.objects}: {missing}")
    os.makedirs(args.out_dir, exist_ok=True)
    feats = args.features

    # Per cell: object count and median of each feature
    cell_medians = objects.groupby("cell_uid")[feats].median()
    per_cell = (objects.groupby(["volume_id", "cell_uid"]).size().rename("n_objects")
                .reset_index().set_index("cell_uid")
                .join(cell_medians.add_prefix("median_")).reset_index())
    per_cell.to_csv(os.path.join(args.out_dir, "per_cell.csv"), index=False)

    # Summary tables: all objects pooled, and across cells
    object_summary = objects[feats].agg(["count", "median", iqr]).T
    cell_summary = cell_medians.agg(["count", "median", iqr]).T
    cell_summary["cv_between_cells"] = cell_medians.std() / cell_medians.mean()
    object_summary.to_csv(os.path.join(args.out_dir, "summary_objects.csv"), index_label="feature")
    cell_summary.to_csv(os.path.join(args.out_dir, "summary_cells.csv"), index_label="feature")

    # Redundant features
    rho = objects[feats].corr(method="spearman")
    pairs = (rho.where(np.triu(np.ones(rho.shape, dtype=bool), k=1))
             .stack().rename("spearman_rho").reset_index()
             .rename(columns={"level_0": "feature_1", "level_1": "feature_2"}))
    redundant = pairs[pairs["spearman_rho"].abs() > args.corr_threshold]
    redundant.to_csv(os.path.join(args.out_dir, "redundant_pairs.csv"), index=False)

    print(f"[info] {len(objects)} objects in {objects['cell_uid'].nunique()} cells "
          f"from {objects['volume_id'].nunique()} volume(s)")
    print("\nObjects (pooled):")
    print(object_summary.to_string(float_format=lambda v: f"{v:.4g}"))
    print("\nCells (one median per cell):")
    print(cell_summary.to_string(float_format=lambda v: f"{v:.4g}"))
    if len(redundant):
        print(f"\nRedundant pairs (|rho| > {args.corr_threshold}); report one feature from each:")
        print(redundant.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    for feat in args.plot_features:
        base = plot_per_cell(objects, cell_medians if feat in feats
                             else objects.groupby("cell_uid")[[feat]].median(), feat, args.out_dir)
        print(f"[io] wrote {base}.svg/.png")

    print(f"[done] outputs in {args.out_dir}")


if __name__ == "__main__":
    main()
