"""Select the best Daubechies wavelet for destriping a single image.

Applies the wavelet-Fourier filter with each wavelet (db1-db43 by default)
at a fixed decomposition level and sigma, scores the result against the
original with SSIM, relative energy loss, and PSNR, and picks the wavelet
with the best weighted score:

    score = 0.5 * SSIM + 0.3 * (1 - energy loss) + 0.2 * PSNR
    (each metric min-max normalized across wavelets)

Usage:
    python best_wavelet.py image.tif --level 6 --sigma 6 [--horizontal]
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from tqdm import tqdm

from destripe import im2double, pick_path, read_image, remove_stripes


def ssim_matlab(img, ref):
    """SSIM with MATLAB ssim() defaults (Gaussian window, sigma 1.5)."""
    return structural_similarity(img, ref, data_range=1.0, gaussian_weights=True,
                                 sigma=1.5, use_sample_covariance=False)


def normalize_range(x):
    return (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else np.zeros_like(x)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("image", nargs="?", help="Input .tif (opens a dialog if omitted)")
    p.add_argument("--level", type=int, default=6, help="Decomposition level (default 6)")
    p.add_argument("--sigma", type=float, default=6, help="Gaussian damping (default 6)")
    p.add_argument("--max-db", type=int, default=43, help="Test db1..dbN (default 43)")
    p.add_argument("--horizontal", action="store_true", help="Remove horizontal stripes")
    p.add_argument("--csv", help="Save the results table to this CSV file")
    args = p.parse_args()

    path = args.image or pick_path("file")
    orig = im2double(read_image(path))
    print("original image type:", read_image(path).dtype)

    wavelets = [f"db{n}" for n in range(1, args.max_db + 1)]
    ssim_s, loss_s, psnr_s = [], [], []

    for name in tqdm(wavelets, desc="Testing wavelets"):
        out = remove_stripes(orig, args.level, name, args.sigma,
                             vertical=not args.horizontal)
        ssim_s.append(ssim_matlab(out, orig))
        loss_s.append(np.sum((orig - out) ** 2) / np.sum(orig ** 2))
        psnr_s.append(peak_signal_noise_ratio(orig, out, data_range=1.0))

    ssim_s, loss_s, psnr_s = map(np.array, (ssim_s, loss_s, psnr_s))
    score = (0.5 * normalize_range(ssim_s)
             + 0.3 * (1 - normalize_range(loss_s))
             + 0.2 * normalize_range(psnr_s))
    best = int(np.argmax(score))
    best_name = wavelets[best]

    # Results table
    print(f"\n{'Wavelet':>8} {'SSIM':>10} {'Energy_Loss':>12} {'PSNR':>9} {'Final_Score':>12}")
    for i, name in enumerate(wavelets):
        print(f"{name:>8} {ssim_s[i]:10.5f} {loss_s[i]:12.6f} {psnr_s[i]:9.3f} {score[i]:12.4f}")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Wavelet", "SSIM", "Energy_Loss", "PSNR", "Final_Score"])
            for i, name in enumerate(wavelets):
                w.writerow([name, ssim_s[i], loss_s[i], psnr_s[i], score[i]])

    # Plots
    idx = np.arange(1, len(wavelets) + 1)
    ticks = idx[::5]
    fig, axes = plt.subplots(4, 1, figsize=(8, 10), sharex=True, layout="constrained")
    panels = [
        (loss_s * 100, "Energy Loss", "Loss [%]"),
        (ssim_s, "SSIM", "SSIM"),
        (psnr_s, "PSNR", "PSNR [dB]"),
        (score, "Composite Score (min-max normalized metrics)", "Final score"),
    ]
    for ax, (vals, title, ylabel) in zip(axes, panels):
        ax.plot(idx, vals, "-o", ms=4)
        ax.axvline(best + 1, ls="--", color="k")
        ax.plot(best + 1, vals[best], "ks", ms=7)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True)
    axes[0].text(best + 0.5, axes[0].get_ylim()[1], f"Selected: {best_name} ", va="top", ha="right")
    axes[-1].set_xticks(ticks, [wavelets[i - 1] for i in ticks], rotation=45)
    axes[-1].set_xlabel(f"Wavelet type (Daubechies db1-db{args.max_db})")
    axes[-1].text(0.02, 0.95,
                  f"Selected wavelet: {best_name}\nSSIM: {ssim_s[best]:.5f}\n"
                  f"Energy loss: {loss_s[best] * 100:.4f} %\nPSNR: {psnr_s[best]:.3f} dB",
                  transform=axes[-1].transAxes, va="top",
                  bbox=dict(facecolor="w", edgecolor="k"))

    fig_path = Path(path).with_name(Path(path).stem + "_wavelet_selection.png")
    fig.savefig(fig_path, dpi=200)
    print(f"\nBest wavelet based on combined score: {best_name}")
    print(f"Plot saved to {fig_path}")
    plt.show()


if __name__ == "__main__":
    main()
