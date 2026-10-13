#!/usr/bin/env python3
"""
Fraction of each cell's volume occupied by objects (e.g. mitochondrial volume
fraction), from a cell label volume and the instance labels.

The cell label IDs must match the cell_id values used elsewhere (e.g. the
--cell volume given to the proofreader).

Only meaningful for cells that are fully inside the imaged volume, or when
equivalent regions (e.g. the same cell compartment) are compared.

Example:
    python volume_fraction.py \
        --instances_mrc data/vol01_instances.mrc \
        --cells_mrc     data/vol01_cells.mrc \
        --volume_id     vol01 \
        --out           results/vol01_volume_fraction.csv
"""

import argparse
import os

import mrcfile
import numpy as np
import pandas as pd


def read_labels(path: str) -> np.ndarray:
    with mrcfile.open(path, permissive=True) as f:
        return np.asarray(f.data).astype(np.int64)


def volume_fraction(instance_labels: np.ndarray, cell_labels: np.ndarray) -> pd.DataFrame:
    """Voxels per cell, object voxels inside each cell, and their ratio."""
    n_cell = np.bincount(cell_labels.ravel())
    n_obj = np.bincount(cell_labels[instance_labels > 0].ravel(), minlength=len(n_cell))
    vf = pd.DataFrame({"cell_id": np.arange(len(n_cell)),
                       "cell_voxels": n_cell, "object_voxels": n_obj})
    vf = vf[(vf["cell_id"] > 0) & (vf["cell_voxels"] > 0)].copy()
    vf["volume_fraction"] = vf["object_voxels"] / vf["cell_voxels"]
    return vf


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instances_mrc", required=True, help="Instance labels (0 = background)")
    ap.add_argument("--cells_mrc", required=True, help="Cell labels, same shape (0 = outside cells)")
    ap.add_argument("--volume_id", required=True, help="Volume ID, used to build cell_uid")
    ap.add_argument("--out", required=True, help="Output CSV")
    ap.add_argument("--spacing_nm", nargs=3, type=float, default=None, metavar=("X", "Y", "Z"),
                    help="Optional voxel size (x y z, nm) to also report volumes in µm³")
    args = ap.parse_args()

    inst = read_labels(args.instances_mrc)
    cells = read_labels(args.cells_mrc)
    if inst.shape != cells.shape:
        raise SystemExit(f"[error] shape mismatch: instances {inst.shape} vs cells {cells.shape}")

    vf = volume_fraction(inst, cells)
    vf.insert(0, "volume_id", args.volume_id)
    vf.insert(1, "cell_uid", args.volume_id + "_" + vf["cell_id"].astype(str))
    if args.spacing_nm:
        voxel_um3 = float(np.prod(args.spacing_nm)) / 1e9
        vf["cell_volume_um3"] = vf["cell_voxels"] * voxel_um3
        vf["object_volume_um3"] = vf["object_voxels"] * voxel_um3

    n_outside = int(((inst > 0) & (cells == 0)).sum())
    if n_outside:
        print(f"[warn] {n_outside} object voxels lie outside any cell and are not counted")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    vf.to_csv(args.out, index=False)
    print(vf.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"[done] wrote {args.out} ({len(vf)} cells)")


if __name__ == "__main__":
    main()
