#!/usr/bin/env python3
"""
Combine per-volume PyRadiomics results into one object table.

For each volume:
- reads the CSV written by radiomics_per_mito.py (full-volume or per-cell mode)
- keeps rows with status == "ok" and reports how many failed
- flags objects that touch the edge of the full volume (truncated objects),
  using the full-volume instance labels
- adds volume_id and a cell ID that is unique across volumes (cell_uid)

Writes one CSV with one row per object, with short names for the most-used
shape features.

Example:
    python combine_results.py \
        --volumes vol01 vol02 \
        --radiomics_csv "results/{volume}_radiomics.csv" \
        --instances_mrc "data/{volume}_instances.mrc" \
        --out objects.csv
"""

import argparse
import os

import mrcfile
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

# Short names for the most-used PyRadiomics shape features
SHORT_NAMES = {
    "original_shape_MeshVolume": "volume_um3",
    "original_shape_SurfaceArea": "surface_area_um2",
    "original_shape_Sphericity": "sphericity",
    "original_shape_Elongation": "elongation",
    "original_shape_Flatness": "flatness",
    "original_shape_Maximum3DDiameter": "max_diameter_um",
}


def border_flags(labels: np.ndarray) -> pd.DataFrame:
    """True for objects that touch any face of the volume (truncated)."""
    rows = []
    for iid, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        touches = any(s.start == 0 or s.stop == n for s, n in zip(sl, labels.shape))
        rows.append((iid, touches))
    return pd.DataFrame(rows, columns=["instance_id", "touches_border"])


def read_labels(path: str) -> np.ndarray:
    with mrcfile.open(path, permissive=True) as f:
        # Labels may be stored as float32 when IDs exceed 16 bits
        return np.asarray(f.data).astype(np.int32)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--volumes", nargs="+", required=True,
                    help="Volume IDs, e.g. vol01 vol02. Substituted for {volume} in the paths.")
    ap.add_argument("--radiomics_csv", required=True,
                    help='Path pattern for radiomics results, e.g. "results/{volume}_radiomics.csv"')
    ap.add_argument("--instances_mrc", required=True,
                    help='Path pattern for full-volume instance labels, e.g. "data/{volume}_instances.mrc"')
    ap.add_argument("--out", default="objects.csv", help="Output CSV (default: objects.csv)")
    args = ap.parse_args()

    tables = []
    for vol in args.volumes:
        csv_path = args.radiomics_csv.format(volume=vol)
        mrc_path = args.instances_mrc.format(volume=vol)

        df = pd.read_csv(csv_path)
        n_total = len(df)
        n_err = int((df["status"] != "ok").sum())
        df = df[df["status"] == "ok"].copy()

        # Per-cell exports: use the original full-volume ID
        if "global_instance_id" in df.columns:
            if (df["global_instance_id"] == 0).any():
                print(f"[warn] {vol}: some objects have no global_instance_id "
                      "(export run without --export_relabel?); border flags will be missing for them")
            df = df.drop(columns=[c for c in ["instance_id"] if c in df.columns])
            df = df.rename(columns={"global_instance_id": "instance_id"})

        # Border flags always come from the full volume
        labels = read_labels(mrc_path)
        df = df.merge(border_flags(labels), on="instance_id", how="left")

        df["volume_id"] = vol
        df["cell_uid"] = vol + "_" + df["cell_id"].astype(int).astype(str)
        tables.append(df)

        print(f"[io] {vol}: {n_total} rows, {n_err} failed, {len(df)} kept, "
              f"{int(df['touches_border'].fillna(False).sum())} touch the volume border")
        if n_err:
            print(f"[warn] {vol}: {n_err} objects failed feature extraction; see the 'error' column in {csv_path}")

    objects = pd.concat(tables, ignore_index=True)
    objects = objects.rename(columns=SHORT_NAMES)
    objects["aspect_ratio"] = (objects["original_shape_MajorAxisLength"]
                               / objects["original_shape_MinorAxisLength"].replace(0, np.nan))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    objects.to_csv(args.out, index=False)
    print(f"[done] wrote {args.out}: {len(objects)} objects from {len(args.volumes)} volume(s)")


if __name__ == "__main__":
    main()
