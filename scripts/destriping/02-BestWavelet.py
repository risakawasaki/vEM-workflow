"""Select the best Daubechies wavelet for destriping a single image.

Applies the wavelet-Fourier filter with each wavelet (db1-db43 by default)
at a fixed decomposition level and sigma, and scores each result with four
metrics:

    Stripe reduction  - how much of the stripe signal was removed (higher is better)
    SSIM              - structural similarity to the original (higher is better)
    Energy loss       - relative signal energy removed (lower is better)
    PSNR              - peak signal-to-noise ratio vs. the original (higher is better)

Stripe reduction is measured on the column-mean profile (row-mean profile for
horizontal stripes). Averaging along the stripe direction suppresses tissue
structure and keeps the stripes. The profile is high-pass filtered to remove
slow illumination gradients, and its variance is the stripe strength:

    stripe reduction = 1 - stripe strength(filtered) / stripe strength(original)

Each metric is min-max normalized across wavelets and combined:

    score = w_stripe * stripe reduction + w_ssim * SSIM
            + w_loss * (1 - energy loss) + w_psnr * PSNR

Usage:
    python best_wavelet.py image.tif --level 6 --sigma 6 [--horizontal]
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from tqdm import tqdm

from destripe import im2double, pick_path, read_image, remove_stripes


def ssim_matlab(img, ref):
    """SSIM with MATLAB ssim() defaults (Gaussian window, sigma 1.5)."""
    return structural_similarity(img, ref, data_range=1.0, gaussian_weights=True,
                                 sigma=1.5, use_sample_covariance=False)


def stripe_strength(img, vertical=True, scale=50):
    """Variance of the high-pass filtered mean profile across the stripes.

    scale : Gaussian sigma (pixels) of the background removed from the
        profile. Stripes much wider than this count as background.
    """
    profile = img.mean(axis=0 if vertical else 1)
    highpass = profile - gaussian_filter1d(profile, scale, mode="reflect")
    return np.var(highpass)


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
    p.add_argument("--stripe-scale", type=float, default=50,
                   help="Background scale in pixels for the stripe metric (default 50)")
    p.add_argument("--weights", type=float, nargs=4, default=[0.4, 0.3, 0.2, 0.1],
                   metavar=("STRIPE", "SSIM", "LOSS", "PSNR"),
                   help="Score weights (default 0.4 0.3 0.2 0.1)")
    p.add_argument("--csv", help="Save the results table to this CSV file")
    args = p.parse_args()

    path = args.image or pick_path("file")
    raw = read_image(path)
    print("original image type:", raw.dtype)
    orig = im2double(raw)
    vertical = not args.horizontal

    s_orig = stripe_strength(orig, vertical, args.stripe_scale)
    wavelets = [f"db{n}" for n in range(1, args.max_db + 1)]
    stripe_s, ssim_s, loss_s, psnr_s = [], [], [], []

    for name in tqdm(wavelets, desc="Testing wavelets"):
        out = remove_stripes(orig, args.level, name, args.sigma, vertical=vertical)
        stripe_s.append(1 - stripe_strength(out, vertical, args.stripe_scale) / s_orig)
        ssim_s.append(ssim_matlab(out, orig))
        loss_s.append(np.sum((orig - out) ** 2) / np.sum(orig ** 2))
        psnr_s.append(peak_signal_noise_ratio(orig, out, data_range=1.0))

    stripe_s, ssim_s, loss_s, psnr_s = map(np.array, (stripe_s, ssim_s, loss_s, psnr_s))
    w_stripe, w_ssim, w_loss, w_psnr = args.weights
    score = (w_stripe * normalize_range(stripe_s)
             + w_ssim * normalize_range(ssim_s)
             + w_loss * (1 - normalize_range(loss_s))
             + w_psnr * normalize_range(psnr_s))
    best = int(np.argmax(score))
    best_name = wavelets[best]

    # Results table
    header = ["Wavelet", "Stripe_Reduction", "SSIM", "Energy_Loss", "PSNR", "Final_Score"]
    print(f"\n{header[0]:>8} {header[1]:>17} {header[2]:>9} {header[3]:>12} "
          f"{header[4]:>8} {header[5]:>12}")
    for i, name in enumerate(wavelets):
        print(f"{name:>8} {stripe_s[i]:17.4f} {ssim_s[i]:9.5f} {loss_s[i]:12.6f} "
              f"{psnr_s[i]:8.3f} {score[i]:12.4f}")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            for i, name in enumerate(wavelets):
                w.writerow([name, stripe_s[i], ssim_s[i], loss_s[i], psnr_s[i], score[i]])

    # Plots
    idx = np.arange(1, len(wavelets) + 1)
    ticks = idx[::5]
    fig, axes = plt.subplots(5, 1, figsize=(8, 12), sharex=True, layout="constrained")
    panels = [
        (stripe_s * 100, "Stripe Reduction", "Reduction [%]"),
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
    axes[0].text(best + 0.5, axes[0].get_ylim()[1], f"Selected: {best_name} ",
                 va="top", ha="right")
    axes[-1].set_xticks(ticks, [wavelets[i - 1] for i in ticks], rotation=45)
    axes[-1].set_xlabel(f"Wavelet type (Daubechies db1-db{args.max_db})")
    axes[-1].text(0.02, 0.95,
                  f"Selected wavelet: {best_name}\n"
                  f"Stripe reduction: {stripe_s[best] * 100:.2f} %\n"
                  f"SSIM: {ssim_s[best]:.5f}\n"
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
