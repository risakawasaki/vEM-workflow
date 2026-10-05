"""Apply wavelet-Fourier destriping to every image in a folder.

Processes all .tif/.tiff files in the input folder and saves the filtered
stack (16-bit, same filenames) to a sibling folder named
<folder>_wfft_<v|h>_<wavelet>_L<level>_S<sigma>.

Usage:
    python wavelet_fft.py /path/to/stack --wavelet db41 --level 6 --sigma 6
"""
import argparse
import warnings
from pathlib import Path

import tifffile
from tqdm import tqdm

from destripe import im2double, pick_path, read_image, remove_stripes, to_uint16


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("folder", nargs="?", help="Image folder (opens a dialog if omitted)")
    p.add_argument("--wavelet", default="db41", help="Wavelet (default db41)")
    p.add_argument("--level", type=int, default=6, help="Decomposition level (default 6)")
    p.add_argument("--sigma", type=float, default=6, help="Gaussian damping (default 6)")
    p.add_argument("--horizontal", action="store_true", help="Remove horizontal stripes")
    args = p.parse_args()

    folder = Path(args.folder or pick_path("folder")).resolve()
    direction = "h" if args.horizontal else "v"
    out_dir = folder.parent / (f"{folder.name}_wfft_{direction}_{args.wavelet}"
                               f"_L{args.level}_S{args.sigma:g}")
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(f for f in folder.iterdir()
                   if f.is_file() and f.suffix.lower() in (".tif", ".tiff")
                   and not f.name.startswith("."))
    if not files:
        raise SystemExit("No valid .tif or .tiff images were found in the selected folder.")

    for f in tqdm(files, desc="Processing images"):
        try:
            img = im2double(read_image(f))
        except Exception as e:
            warnings.warn(f"Skipping unreadable file: {f}\nReason: {e}")
            continue
        out = remove_stripes(img, args.level, args.wavelet, args.sigma,
                             vertical=not args.horizontal)
        tifffile.imwrite(out_dir / f.name, to_uint16(out))

    print(f"Done. Output: {out_dir}")


if __name__ == "__main__":
    main()
