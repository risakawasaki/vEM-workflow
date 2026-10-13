#!/usr/bin/env python3
"""
Apply quality filters to the object table and record how many objects each
filter removed.

Filters are applied in this order:
1. unassigned objects (cell_id == 0)            on by default; --keep_unassigned to skip
2. objects touching the volume border            on by default; --keep_border to skip
3. objects smaller than --min_voxels voxels      off by default
4. downsampled objects (did_downsample == True)  off by default; --exclude_downsampled
5. listed cells (--exclude_cells)                off by default

Decide the filters before looking at the results, and use the same filters for
every volume.

Example:
    python filter_objects.py --objects objects.csv --out objects_filtered.csv \
        --min_voxels 200 --exclude_cells vol01_12 vol01_15
"""

import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--objects", default="objects.csv", help="Input object table (from combine_results.py)")
    ap.add_argument("--out", default="objects_filtered.csv", help="Filtered object table")
    ap.add_argument("--log", default=None,
                    help="CSV recording objects removed per filter (default: <out stem>_filter_log.csv)")
    ap.add_argument("--keep_unassigned", action="store_true", help="Keep objects with cell_id == 0")
    ap.add_argument("--keep_border", action="store_true", help="Keep objects that touch the volume border")
    ap.add_argument("--min_voxels", type=int, default=0, help="Remove objects with fewer voxels (default: 0, off)")
    ap.add_argument("--exclude_downsampled", action="store_true",
                    help="Remove objects that were downsampled before feature extraction")
    ap.add_argument("--exclude_cells", nargs="*", default=[],
                    help="cell_uid values to remove, e.g. incompletely imaged cells (vol01_12)")
    args = ap.parse_args()

    df = pd.read_csv(args.objects)
    log = [("input", len(df), 0)]

    def apply(name, keep_mask):
        nonlocal df
        before = len(df)
        df = df[keep_mask.fillna(True)]
        log.append((name, len(df), before - len(df)))

    if not args.keep_unassigned:
        apply("unassigned (cell_id == 0)", df["cell_id"] != 0)
    if not args.keep_border:
        if "touches_border" not in df.columns:
            raise SystemExit("[error] no touches_border column; run combine_results.py first")
        apply("touches volume border", ~df["touches_border"].astype("boolean"))
    if args.min_voxels > 0:
        apply(f"fewer than {args.min_voxels} voxels", df["voxels"] >= args.min_voxels)
    if args.exclude_downsampled:
        apply("downsampled", ~df["did_downsample"].astype("boolean"))
    if args.exclude_cells:
        missing = sorted(set(args.exclude_cells) - set(df["cell_uid"]))
        if missing:
            print(f"[warn] --exclude_cells not found (already removed or misspelled): {missing}")
        apply("excluded cells", ~df["cell_uid"].isin(args.exclude_cells))

    log = pd.DataFrame(log, columns=["filter", "objects_remaining", "objects_removed"])
    print(log.to_string(index=False))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    log_path = args.log or os.path.splitext(args.out)[0] + "_filter_log.csv"
    log.to_csv(log_path, index=False)

    if "did_downsample" in df.columns and not args.exclude_downsampled:
        n_ds = int(df["did_downsample"].astype("boolean").fillna(False).sum())
        if n_ds:
            print(f"[warn] {n_ds} remaining objects were downsampled before feature extraction; "
                  "check that results hold with --exclude_downsampled")

    print(f"[done] wrote {args.out} ({len(df)} objects, {df['cell_uid'].nunique()} cells) and {log_path}")


if __name__ == "__main__":
    main()
