import os
import json
import argparse
import numpy as np
import pandas as pd
import mrcfile
from scipy import ndimage as ndi


def read_mrc(path):
    with mrcfile.open(path, permissive=True) as mrc:
        data = np.asarray(mrc.data)
        try:
            vx = float(mrc.voxel_size.x)
            vy = float(mrc.voxel_size.y)
            vz = float(mrc.voxel_size.z)
            spacing = (vx / 10.0, vy / 10.0, vz / 10.0)  # Å → nm
        except Exception:
            spacing = (np.nan, np.nan, np.nan)
    return data, spacing


def write_mrc(path, data, spacing_nm_xyz):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = np.ascontiguousarray(data)

    with mrcfile.new(path, overwrite=True) as mrc:
        mrc.set_data(data.astype(np.uint16, copy=False))
        sx, sy, sz = spacing_nm_xyz
        mrc.voxel_size = (sx * 10, sy * 10, sz * 10)
        mrc.update_header_from_data()
        mrc.update_header_stats()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instances_mrc", required=True)
    ap.add_argument("--mapping_csv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--raw_mrc", default=None, help="Optional raw vEM MRC to crop per cell (same shape as instances).")
    ap.add_argument("--export_raw", action="store_true", help="If set (and --raw_mrc provided), export cropped raw per cell.")
    ap.add_argument("--export_relabel", action="store_true", help="If set, export relabelled per-cell instance labels (1..N) and include mapping in meta.")
    args = ap.parse_args()

    inst, spacing = read_mrc(args.instances_mrc)
    mapping = pd.read_csv(args.mapping_csv)

    raw = None
    raw_spacing = None
    if args.raw_mrc is not None:
        raw, raw_spacing = read_mrc(args.raw_mrc)
        raw = np.asarray(raw)
        if raw.shape != inst.shape:
            raise ValueError(f"raw_mrc shape {raw.shape} does not match instances shape {inst.shape}")
        # Prefer raw spacing if it exists, else fall back to instance spacing
        if raw_spacing is None or (isinstance(raw_spacing, tuple) and any(np.isnan(x) for x in raw_spacing)):
            raw_spacing = spacing

    inst = np.asarray(inst)
    inst_int = inst.astype(np.int32, copy=False)
    print("[prep] Precomputing bounding boxes for all instance labels...")
    bboxes = ndi.find_objects(inst_int)
    print(f"[prep] bbox cache ready for labels 1..{len(bboxes)}")
    mapping["instance_id"] = mapping["instance_id"].astype(int)
    mapping["cell_id"] = mapping["cell_id"].astype(int)

    for cell_id, grp in mapping.groupby("cell_id"):
        inst_ids = grp["instance_id"].tolist()

        # Compute union bbox from per-instance bboxes (fast, avoids full-volume masks)
        z0 = y0 = x0 = None
        z1 = y1 = x1 = None
        present_ids = []

        for iid in inst_ids:
            iid = int(iid)
            if iid <= 0 or iid > len(bboxes):
                continue
            slc = bboxes[iid - 1]
            if slc is None:
                continue
            present_ids.append(iid)
            if z0 is None:
                z0, y0, x0 = slc[0].start, slc[1].start, slc[2].start
                z1, y1, x1 = slc[0].stop,  slc[1].stop,  slc[2].stop
            else:
                z0 = min(z0, slc[0].start); y0 = min(y0, slc[1].start); x0 = min(x0, slc[2].start)
                z1 = max(z1, slc[0].stop);  y1 = max(y1, slc[1].stop);  x1 = max(x1, slc[2].stop)

        if z0 is None:
            print(f"[skip] cell {cell_id}: no voxels")
            continue

        crop = inst[z0:z1, y0:y1, x0:x1].copy()
        # Zero out anything not belonging to this cell's instance ids (operate within crop only)
        crop[~np.isin(crop, present_ids)] = 0

        cell_dir = os.path.join(args.out_dir, f"cell_{cell_id:03d}")
        os.makedirs(cell_dir, exist_ok=True)

        mrc_path = os.path.join(cell_dir, "mito_labels.mrc")
        meta_path = os.path.join(cell_dir, "meta.json")

        write_mrc(mrc_path, crop, spacing)

        if args.export_raw and raw is not None:
            raw_crop = raw[z0:z1, y0:y1, x0:x1]
            raw_path = os.path.join(cell_dir, "raw_vem.mrc")
            # Write raw with its original dtype if possible
            with mrcfile.new(raw_path, overwrite=True) as mrc:
                mrc.set_data(np.ascontiguousarray(raw_crop))
                sx, sy, sz = (raw_spacing if raw_spacing is not None else spacing)
                mrc.voxel_size = (sx * 10, sy * 10, sz * 10)
                mrc.update_header_from_data()
                mrc.update_header_stats()
            print(f"[ok] cell {cell_id}: raw crop saved shape={raw_crop.shape}")

        relabel_map = None
        if args.export_relabel:
            # Relabel the present instance ids to 1..N within this cell crop
            uniq = sorted(int(x) for x in np.unique(crop) if int(x) != 0)
            old_to_new = {old: (i + 1) for i, old in enumerate(uniq)}
            relabel = np.zeros_like(crop, dtype=np.uint16)
            for old, new in old_to_new.items():
                relabel[crop == old] = np.uint16(new)
            relabel_path = os.path.join(cell_dir, "mito_labels_relabel.mrc")
            write_mrc(relabel_path, relabel, spacing)
            relabel_map = [{"old_instance_id": int(o), "new_instance_id": int(n)} for o, n in old_to_new.items()]
            print(f"[ok] cell {cell_id}: relabel saved instances={len(uniq)}")

        meta = {
            "cell_id": int(cell_id),
            "full_shape_zyx": [int(inst.shape[0]), int(inst.shape[1]), int(inst.shape[2])],
            "bbox_zyx": [int(z0), int(y0), int(x0), int(z1), int(y1), int(x1)],
            "voxel_spacing_nm_xyz": list(spacing),
            "instance_ids": [int(x) for x in present_ids],
            "has_raw": bool(args.export_raw and raw is not None),
            "has_relabel": bool(args.export_relabel),
        }
        if relabel_map is not None:
            meta["relabel_map"] = relabel_map

        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        print(f"[ok] cell {cell_id}: crop_shape={crop.shape} instances={len(present_ids)} bbox_zyx=({z0},{y0},{x0})-({z1},{y1},{x1})")

    print("[done] Cell-wise exports complete.")


if __name__ == "__main__":
    main()
