import numpy as np
import pandas as pd
import mrcfile
from scipy import ndimage as ndi

def border_flags(labels):
    """True for objects that touch any face of the volume (truncated)."""
    rows = []
    for iid, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        touches = any(s.start == 0 or s.stop == n for s, n in zip(sl, labels.shape))
        rows.append((iid, touches))
    return pd.DataFrame(rows, columns=["instance_id", "touches_border"])

volumes = ["vol01"]                     # TODO: our volume IDs

tables = []
for vol in volumes:
    df = pd.read_csv(f"results/{vol}_radiomics.csv")
    n_err = int((df["status"] != "ok").sum())
    if n_err:
        print(f"[warn] {vol}: {n_err} objects failed feature extraction")
    df = df[df["status"] == "ok"].copy()
    if "global_instance_id" in df.columns:          # per-cell exports
        df = df.rename(columns={"global_instance_id": "instance_id"})

    # Border flags always come from the full volume
    with mrcfile.open(f"data/{vol}_instances.mrc", permissive=True) as f:
        labels = np.asarray(f.data).astype(np.int32)
    df = df.merge(border_flags(labels), on="instance_id", how="left")

    df["volume_id"] = vol
    df["cell_uid"] = vol + "_" + df["cell_id"].astype(str)   # unique across volumes
    tables.append(df)

objects = pd.concat(tables, ignore_index=True)

# Short names for the features used most often
objects = objects.rename(columns={
    "original_shape_MeshVolume": "volume_um3",
    "original_shape_SurfaceArea": "surface_area_um2",
    "original_shape_Sphericity": "sphericity",
    "original_shape_Elongation": "elongation",
    "original_shape_Flatness": "flatness",
    "original_shape_Maximum3DDiameter": "max_diameter_um",
})
objects["aspect_ratio"] = (objects["original_shape_MajorAxisLength"]
                           / objects["original_shape_MinorAxisLength"])
objects.to_csv("objects.csv", index=False)
