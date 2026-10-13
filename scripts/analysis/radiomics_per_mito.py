#!/usr/bin/env python3
"""
Compute PyRadiomics shape features per mitochondria instance from:
- instance labels volume (MRC; each mitochondrion has unique integer id)
- instance->cell_id mapping CSV (columns: instance_id, cell_id)

Only shape features are computed, from the label masks alone. The raw image
volume is not used, so no intensity or texture features are produced.

Outputs a CSV with one row per instance, including cell_id.
"""

import argparse
import os
import sys
import json
import glob
import re
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import mrcfile

import SimpleITK as sitk
from radiomics import featureextractor


def read_mrc(path: str) -> np.ndarray:
    with mrcfile.open(path, permissive=True) as mrc:
        return np.asarray(mrc.data)


def load_cellmap_csv(path: str) -> Dict[int, int]:
    df = pd.read_csv(path)
    if not {"instance_id", "cell_id"}.issubset(df.columns):
        raise ValueError("cellmap CSV must contain columns: instance_id, cell_id")
    out: Dict[int, int] = {}
    for iid, cid in zip(df["instance_id"].to_numpy(), df["cell_id"].to_numpy()):
        try:
            iid_i = int(iid)
            cid_i = int(cid)
        except Exception:
            continue
        out[iid_i] = cid_i
    return out


def load_meta_mapping(meta_path: str) -> Dict[int, int]:
    """
    Best-effort loader for a per-cell export meta.json mapping from local relabeled instance ids
    (1..K in mito_labels_relabel.mrc) back to original/global instance ids.

    Because meta.json formats can vary, we try several common key names.
    Returns dict: local_id -> global_id. Empty dict if not found.
    """
    try:
        with open(meta_path, "r") as f:
            meta = json.load(f)
    except Exception:
        return {}

    # Try several likely keys (keep this flexible)
    candidate_keys = [
        "relabel_to_global",
        "local_to_global",
        "instance_local_to_global",
        "instance_relabel_map",
        "relabel_map",
        "instance_id_map",
        "mito_relabel_to_global",
    ]
    for k in candidate_keys:
        if k in meta and isinstance(meta[k], dict):
            out = {}
            for kk, vv in meta[k].items():
                try:
                    out[int(kk)] = int(vv)
                except Exception:
                    continue
            if out:
                return out

    # Some formats store relabel_map as a list of dicts:
    # [{"old_instance_id": 447, "new_instance_id": 1}, ...]
    if "relabel_map" in meta and isinstance(meta["relabel_map"], list):
        out = {}
        for item in meta["relabel_map"]:
            if not isinstance(item, dict):
                continue
            if "old_instance_id" in item and "new_instance_id" in item:
                try:
                    global_id = int(item["old_instance_id"])
                    local_id = int(item["new_instance_id"])
                    if local_id > 0 and global_id > 0:
                        out[local_id] = global_id
                except Exception:
                    continue
        if out:
            return out

    # Some formats store as list of pairs
    candidate_list_keys = [
        "relabel_pairs",
        "local_global_pairs",
        "instance_pairs",
    ]
    for k in candidate_list_keys:
        if k in meta and isinstance(meta[k], list):
            out = {}
            for item in meta[k]:
                # list/tuple pair: [local, global]
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    try:
                        out[int(item[0])] = int(item[1])
                    except Exception:
                        pass
                    continue

                # dict style: {"local": 1, "global": 447} or similar
                if isinstance(item, dict):
                    for lk, gk in [("local", "global"), ("new", "old"), ("new_instance_id", "old_instance_id")]:
                        if lk in item and gk in item:
                            try:
                                out[int(item[lk])] = int(item[gk])
                            except Exception:
                                pass
                            break
                    continue
            if out:
                return out

    return {}


def compute_bboxes_streaming(lbl_zyx: np.ndarray, verbose: bool = False) -> Dict[int, Tuple[slice, slice, slice]]:
    max_id = int(lbl_zyx.max())
    INF = 10**18
    zmin = np.full(max_id + 1, INF, dtype=np.int64)
    zmax = np.full(max_id + 1, -1, dtype=np.int64)
    ymin = np.full(max_id + 1, INF, dtype=np.int64)
    ymax = np.full(max_id + 1, -1, dtype=np.int64)
    xmin = np.full(max_id + 1, INF, dtype=np.int64)
    xmax = np.full(max_id + 1, -1, dtype=np.int64)

    Z = lbl_zyx.shape[0]
    for z in range(Z):
        sl = np.asarray(lbl_zyx[z])  # reads one z-slice
        ids = np.unique(sl)
        ids = ids[ids != 0]
        for iid in ids.tolist():
            yy, xx = np.where(sl == iid)
            if yy.size == 0:
                continue
            zmin[iid] = min(zmin[iid], z)
            zmax[iid] = max(zmax[iid], z)
            ymin[iid] = min(ymin[iid], int(yy.min()))
            ymax[iid] = max(ymax[iid], int(yy.max()))
            xmin[iid] = min(xmin[iid], int(xx.min()))
            xmax[iid] = max(xmax[iid], int(xx.max()))
        if verbose and (z % 25 == 0):
            print(f"[bbox] z={z}/{Z}", flush=True)

    bboxes: Dict[int, Tuple[slice, slice, slice]] = {}
    for iid in range(1, max_id + 1):
        if zmax[iid] >= 0:
            bboxes[iid] = (
                slice(int(zmin[iid]), int(zmax[iid]) + 1),
                slice(int(ymin[iid]), int(ymax[iid]) + 1),
                slice(int(xmin[iid]), int(xmax[iid]) + 1),
            )
    return bboxes


def to_sitk(img_zyx: np.ndarray, spacing_um_zyx: Tuple[float, float, float]) -> sitk.Image:
    """Create SimpleITK image with spacing in micrometers, in (z,y,x) numpy order."""
    itk = sitk.GetImageFromArray(img_zyx)  # assumes z,y,x
    # SimpleITK spacing expects (x,y,z)
    sz_um, sy_um, sx_um = spacing_um_zyx
    itk.SetSpacing((sx_um, sy_um, sz_um))
    return itk

def resample_itk_pair(itk_img: sitk.Image, itk_msk: sitk.Image, ds_zyx):
    dz, dy, dx = [int(v) for v in ds_zyx]
    if dz <= 1 and dy <= 1 and dx <= 1:
        return itk_img, itk_msk

    # SimpleITK is (x,y,z)
    sx, sy, sz = itk_img.GetSpacing()
    ox, oy, oz = itk_img.GetOrigin()
    direction = itk_img.GetDirection()
    size_x, size_y, size_z = itk_img.GetSize()

    new_size = (
        max(1, int(np.ceil(size_x / dx))),
        max(1, int(np.ceil(size_y / dy))),
        max(1, int(np.ceil(size_z / dz))),
    )
    new_spacing = (sx * dx, sy * dy, sz * dz)

    tx = sitk.Transform()

    img_rs = sitk.Resample(
        itk_img, new_size, tx, sitk.sitkLinear, (ox, oy, oz), new_spacing, direction, 0.0, sitk.sitkFloat32
    )
    msk_rs = sitk.Resample(
        itk_msk, new_size, tx, sitk.sitkNearestNeighbor, (ox, oy, oz), new_spacing, direction, 0, sitk.sitkUInt8
    )
    return img_rs, msk_rs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inst_mrc", required=False, default=None, help="Instance label MRC (0=bg, 1..N mitochondria).")
    ap.add_argument("--cellmap_csv", required=False, default=None, help="CSV with columns instance_id, cell_id.")
    ap.add_argument("--cell_exports_dir", required=False, default=None,
                    help="Directory containing per-cell exports (cell_###/mito_labels_relabel.mrc + meta.json). "
                         "If provided, features are computed per cell crop and concatenated; cell_id is taken from folder name.")
    ap.add_argument("--out_csv", required=True, help="Output CSV path.")
    ap.add_argument("--spacing_nm", nargs=3, type=float, default=[5.0, 5.0, 20.0], help="Voxel spacing (x y z) in nm.")
    ap.add_argument("--min_voxels", type=int, default=50, help="Skip instances smaller than this many voxels.")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument(
        "--only_cells",
        nargs="*",
        default=None,
        help="In --cell_exports_dir mode, process only these cell ids (e.g. 8 12 25) or folder names (e.g. cell_008).",
    )
    ap.add_argument("--max_voxels_pyrad", type=int, default=500000,
                help="If an instance has more voxels than this, downsample before PyRadiomics.")
    ap.add_argument("--ds_large", nargs=3, type=int, default=[1, 2, 2],
                help="Downsample factors (z y x) used when voxels > max_voxels_pyrad. Recommended: 1 2 2 or 1 4 4.")
    
    args = ap.parse_args()

    # Expand ~ and resolve relative paths for all file/dir inputs
    def _expand(p):
        return os.path.abspath(os.path.expanduser(p.strip())) if p else p

    args.inst_mrc = _expand(args.inst_mrc)
    args.cellmap_csv = _expand(args.cellmap_csv)
    args.cell_exports_dir = _expand(args.cell_exports_dir)
    args.out_csv = _expand(args.out_csv)

    if args.cell_exports_dir:
        if not os.path.isdir(args.cell_exports_dir):
            raise ValueError(f"--cell_exports_dir not found: {args.cell_exports_dir}")
    else:
        # Full-volume mode requires inst/cellmap
        if not (args.inst_mrc and args.cellmap_csv):
            raise ValueError("Provide either --cell_exports_dir OR both of: --inst_mrc --cellmap_csv")
        for label, path in [("--inst_mrc", args.inst_mrc), ("--cellmap_csv", args.cellmap_csv)]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"{label} not found: {path!r}")

    def run_one_volume(lbl_path: str, cellmap: Dict[int, int], *,
                       cell_id_override: int | None = None,
                       global_id_map: Dict[int, int] | None = None,
                       label_id_kind: str = "instance_id") -> List[Dict]:
        """
        Compute shape features for one label volume.
        - cellmap: maps label ids -> cell_id (used in full-volume mode)
        - cell_id_override: if provided, stamps this cell_id for all rows (used in per-cell exports)
        - global_id_map: optional local_id -> global_instance_id mapping (used in per-cell exports)
        - label_id_kind: column name for the label id (default 'instance_id'; per-cell uses 'local_instance_id')
        Returns list of row dicts.
        """
        # Keep the mmapped MRC file open for the full run.
        mlbl = mrcfile.mmap(lbl_path, permissive=True)
        try:
            lbl = mlbl.data

            # Convert spacing nm(x,y,z) -> um(z,y,x) for SimpleITK
            sx_nm, sy_nm, sz_nm = [float(x) for x in args.spacing_nm]
            spacing_um_zyx = (sz_nm / 1000.0, sy_nm / 1000.0, sx_nm / 1000.0)

            if args.verbose:
                print(f"[io] lbl={lbl_path} shape={lbl.shape} dtype={lbl.dtype} max={int(np.max(lbl))}", flush=True)

            bboxes = compute_bboxes_streaming(lbl, verbose=args.verbose)
            instance_ids = sorted(bboxes.keys())

            if args.verbose:
                print(f"[prep] instances with bboxes: {len(instance_ids)}", flush=True)

            # Configure extractor: shape features only (computed from the mask alone)
            settings = {
                "interpolator": sitk.sitkNearestNeighbor,
                "resampledPixelSpacing": None,
                "correctMask": True,
            }
            extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
            extractor.disableAllFeatures()
            extractor.enableFeatureClassByName("shape")

            rows: List[Dict] = []
            t0 = None
            if args.verbose:
                import time
                t0 = time.time()

            for i, iid in enumerate(instance_ids, start=1):
                slc = bboxes[iid]
                sub_lbl = lbl[slc]
                mask = (sub_lbl == iid)
                vox = int(np.count_nonzero(mask))
                if vox < int(args.min_voxels):
                    continue

                # Full-resolution physical metrics (always based on the original mask + original spacing)
                zsl, ysl, xsl = slc
                bbox_sz_vox_zyx = (int(zsl.stop - zsl.start), int(ysl.stop - ysl.start), int(xsl.stop - xsl.start))
                voxel_vol_um3 = float(spacing_um_zyx[0] * spacing_um_zyx[1] * spacing_um_zyx[2])
                volume_um3 = float(vox) * voxel_vol_um3
                bbox_sz_um_zyx = (
                    float(bbox_sz_vox_zyx[0]) * float(spacing_um_zyx[0]),
                    float(bbox_sz_vox_zyx[1]) * float(spacing_um_zyx[1]),
                    float(bbox_sz_vox_zyx[2]) * float(spacing_um_zyx[2]),
                )

                itk_msk = to_sitk(mask.astype(np.uint8, copy=False), spacing_um_zyx)
                # PyRadiomics requires an image alongside the mask; shape features
                # ignore its values, so an all-zero placeholder is used.
                dummy = np.zeros(mask.shape, dtype=np.float32)
                itk_img = to_sitk(dummy, spacing_um_zyx)

                # Downsample only for PyRadiomics if very large ROI (avoid heavy mesh ops / segfaults)
                did_downsample = False
                ds_zyx_used = (1, 1, 1)
                if int(args.max_voxels_pyrad) > 0 and vox > int(args.max_voxels_pyrad):
                    did_downsample = True
                    ds_zyx_used = tuple(int(v) for v in args.ds_large)
                    if args.verbose:
                        zsl, ysl, xsl = slc
                        print(
                            f"[ds] iid={iid} vox={vox} > {int(args.max_voxels_pyrad)} ds_zyx={ds_zyx_used} "
                            f"bbox_zyx=({zsl.start}:{zsl.stop},{ysl.start}:{ysl.stop},{xsl.start}:{xsl.stop})",
                            flush=True,
                        )
                    itk_img, itk_msk = resample_itk_pair(itk_img, itk_msk, ds_zyx_used)

                # Record the effective spacing actually used for PyRadiomics (after any downsampling).
                # SimpleITK spacing is (x,y,z); we store (z,y,x) to match numpy.
                sx_eff, sy_eff, sz_eff = itk_img.GetSpacing()
                spacing_used_um_zyx = (float(sz_eff), float(sy_eff), float(sx_eff))

                try:
                    feats = extractor.execute(itk_img, itk_msk)
                except Exception as e:
                    rows.append({
                        label_id_kind: int(iid),
                        "cell_id": int(cell_id_override) if cell_id_override is not None else int(cellmap.get(iid, 0)),
                        "voxels": vox,
                        "status": "error",
                        "error": str(e),
                        "did_downsample": bool(did_downsample),
                        "ds_zyx": str(ds_zyx_used),
                        "spacing_used_um_z": float(spacing_used_um_zyx[0]),
                        "spacing_used_um_y": float(spacing_used_um_zyx[1]),
                        "spacing_used_um_x": float(spacing_used_um_zyx[2]),
                        "voxel_volume_um3": float(voxel_vol_um3),
                        "volume_um3_voxels": float(volume_um3),
                        "bbox_z_vox": int(bbox_sz_vox_zyx[0]),
                        "bbox_y_vox": int(bbox_sz_vox_zyx[1]),
                        "bbox_x_vox": int(bbox_sz_vox_zyx[2]),
                        "bbox_z_um": float(bbox_sz_um_zyx[0]),
                        "bbox_y_um": float(bbox_sz_um_zyx[1]),
                        "bbox_x_um": float(bbox_sz_um_zyx[2]),
                    })
                    continue

                out = {
                    label_id_kind: int(iid),
                    "cell_id": int(cell_id_override) if cell_id_override is not None else int(cellmap.get(iid, 0)),
                    "voxels": vox,
                    "status": "ok",
                    "did_downsample": bool(did_downsample),
                    "ds_zyx": str(ds_zyx_used),
                    "spacing_used_um_z": float(spacing_used_um_zyx[0]),
                    "spacing_used_um_y": float(spacing_used_um_zyx[1]),
                    "spacing_used_um_x": float(spacing_used_um_zyx[2]),
                    "voxel_volume_um3": float(voxel_vol_um3),
                    "volume_um3_voxels": float(volume_um3),
                    "bbox_z_vox": int(bbox_sz_vox_zyx[0]),
                    "bbox_y_vox": int(bbox_sz_vox_zyx[1]),
                    "bbox_x_vox": int(bbox_sz_vox_zyx[2]),
                    "bbox_z_um": float(bbox_sz_um_zyx[0]),
                    "bbox_y_um": float(bbox_sz_um_zyx[1]),
                    "bbox_x_um": float(bbox_sz_um_zyx[2]),
                }
                if global_id_map is not None:
                    out["global_instance_id"] = int(global_id_map.get(int(iid), 0))

                for k, v in feats.items():
                    if k.startswith("diagnostics_"):
                        continue
                    if isinstance(v, (np.generic,)):
                        v = v.item()
                    out[k] = v

                rows.append(out)

                if args.verbose and (i % 25 == 0):
                    import time
                    print(f"[run] {i}/{len(instance_ids)} done. elapsed={time.time()-t0:.1f}s", flush=True)

            if args.verbose:
                import time
                print(f"[run] finished volume. kept_rows={len(rows)} elapsed={time.time()-t0:.1f}s", flush=True)

            return rows
        finally:
            try:
                mlbl.close()
            except Exception:
                pass

    rows: List[Dict] = []

    if args.cell_exports_dir:
        # Per-cell exports mode: iterate cell_### subfolders.
        cell_dirs = sorted(glob.glob(os.path.join(args.cell_exports_dir, "cell_*")))
        if args.only_cells:
            wanted: set[int] = set()
            for tok in args.only_cells:
                s = str(tok).strip()
                m = re.match(r"^cell_(\d+)$", s)
                if m:
                    wanted.add(int(m.group(1)))
                else:
                    # allow '008' or '8'
                    try:
                        wanted.add(int(s))
                    except Exception:
                        pass
            cell_dirs = [
                d for d in cell_dirs
                if (re.match(r"cell_(\d+)$", os.path.basename(d)) and int(re.match(r"cell_(\d+)$", os.path.basename(d)).group(1)) in wanted)
            ]
        if not cell_dirs:
            raise ValueError(f"No cell_* subfolders found in: {args.cell_exports_dir}")

        if args.verbose:
            print(f"[mode] cell_exports_dir={args.cell_exports_dir} cells={len(cell_dirs)}", flush=True)

        for d in cell_dirs:
            base = os.path.basename(d)
            m = re.match(r"cell_(\d+)$", base)
            if not m:
                continue
            cell_id = int(m.group(1))
            lbl_path = os.path.join(d, "mito_labels_relabel.mrc")
            if not os.path.exists(lbl_path):
                lbl_path = os.path.join(d, "mito_labels.mrc")
            meta_path = os.path.join(d, "meta.json")

            if not os.path.exists(lbl_path):
                if args.verbose:
                    print(f"[skip] missing labels in {d}", flush=True)
                continue

            global_map = load_meta_mapping(meta_path) if os.path.exists(meta_path) else {}
            if args.verbose:
                print(f"[cell] {base} cell_id={cell_id} global_map_entries={len(global_map)}", flush=True)

            if args.verbose:
                print(f"[cell_start] {base} -> running radiomics...", flush=True)

            try:
                cell_rows = run_one_volume(
                    lbl_path,
                    cellmap={},
                    cell_id_override=cell_id,
                    global_id_map=global_map if global_map else None,
                    label_id_kind="local_instance_id",
                )
                rows.extend(cell_rows)
                if args.verbose:
                    ok_n = sum(1 for r in cell_rows if r.get("status") == "ok")
                    err_n = sum(1 for r in cell_rows if r.get("status") == "error")
                    print(f"[cell_done] {base} rows={len(cell_rows)} ok={ok_n} error={err_n} total_rows={len(rows)}", flush=True)
            except Exception as e:
                # Do not abort the entire run if one cell fails; record and continue.
                if args.verbose:
                    print(f"[cell_error] {base} failed: {e}", flush=True)
                rows.append({
                    "cell_id": int(cell_id),
                    "status": "error",
                    "error": f"cell_failed: {e}",
                })
                continue

    else:
        # Full-volume mode (original)
        cellmap = load_cellmap_csv(args.cellmap_csv)
        rows = run_one_volume(
            args.inst_mrc,
            cellmap=cellmap,
            cell_id_override=None,
            global_id_map=None,
            label_id_kind="instance_id",
        )

    df = pd.DataFrame(rows)

    # Helpful ordering
    id_col = "instance_id" if "instance_id" in df.columns else ("local_instance_id" if "local_instance_id" in df.columns else None)
    front = [c for c in [id_col, "global_instance_id", "cell_id", "voxels", "status"] if c and c in df.columns]
    cols = front + [c for c in df.columns if c not in front]
    df = df[cols]

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    print(f"[done] wrote {args.out_csv} rows={len(df)} ok={ok}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[fatal] {e}", file=sys.stderr)
        raise