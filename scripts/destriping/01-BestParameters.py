"""Generate a parameter sweep for destriping a single image.

Applies the wavelet-Fourier filter for every combination of wavelet,
decomposition level, and sigma, and saves each result as a 16-bit PNG
(named <wavelet>-L<level>-S<sigma>.png) for visual comparison.

Parameter guide:
    Level 2-6   -> best for small-to-medium-sized artifacts
    Level 7-10  -> effective for large-scale artifacts, may smooth fine details
    Sigma 1-2   -> preserves fine textures, removes mild stripes
    Sigma 3-5   -> works well for most cases, balancing removal and detail
    Sigma > 7   -> best for irregular or strong stripes, may blur details

Usage:
    python best_parameters.py image.tif \
        --wavelets db2 db15 db30 db42 --levels 2 4 6 8 --sigmas 1 2 4 6 8
"""
import argparse
import itertools
from pathlib import Path

from PIL import Image
from tqdm import tqdm

from destripe import pick_path, read_image, remove_stripes, to_uint16


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("image", nargs="?", help="Input .tif (opens a dialog if omitted)")
    p.add_argument("--wavelets", nargs="+", default=["db2", "db15", "db30", "db42"])
    p.add_argument("--levels", nargs="+", type=int, default=[2, 4, 6, 8])
    p.add_argument("--sigmas", nargs="+", type=float, default=[1, 2, 4, 6, 8])
    p.add_argument("--horizontal", action="store_true", help="Remove horizontal stripes")
    args = p.parse_args()

    path = Path(args.image or pick_path("file"))
    orig = read_image(path)

    suffix = "_parameters_h" if args.horizontal else "_parameters_v"
    save_folder = path.parent.parent / (path.name + suffix)
    save_folder.mkdir(parents=True, exist_ok=True)

    combos = list(itertools.product(args.wavelets, args.levels, args.sigmas))
    for wavelet, level, sigma in tqdm(combos, desc="Parameter sweep"):
        out = remove_stripes(orig, level, wavelet, sigma, vertical=not args.horizontal)
        name = f"{wavelet}-L{level}-S{sigma:g}.png"
        Image.fromarray(to_uint16(out)).save(save_folder / name, dpi=(600, 600))

    print(f"Saved {len(combos)} images to {save_folder}")


if __name__ == "__main__":
    main()
