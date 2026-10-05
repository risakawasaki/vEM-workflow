"""Core functions for wavelet-Fourier destriping (Muench et al., 2009).
"""
from functools import lru_cache

import numpy as np
import pywt


# --------------------------------------------------------------------------
# Wavelets
# --------------------------------------------------------------------------
@lru_cache(maxsize=None)
def get_wavelet(name):
    """Return a pywt.Wavelet. Daubechies wavelets above db38 (not built into
    PyWavelets) are computed numerically so db1-db45 all work, as in MATLAB."""
    if name in pywt.wavelist(kind="discrete"):
        return pywt.Wavelet(name)
    if name.startswith("db") and name[2:].isdigit():
        rec_lo = _daubechies_filter(int(name[2:]))
        return pywt.Wavelet(name, filter_bank=pywt.orthogonal_filter_bank(rec_lo))
    raise ValueError(f"Unknown wavelet: {name}")


def _daubechies_filter(n, dps=100):
    """Daubechies dbN reconstruction low-pass filter via spectral
    factorization in high precision (same construction as MATLAB dbaux)."""
    import mpmath as mp

    with mp.workdps(dps):
        # P(y) = sum_{k=0}^{N-1} C(N-1+k, k) y^k
        coeffs = [mp.binomial(n - 1 + k, k) for k in range(n)]
        y_roots = mp.polyroots(coeffs[::-1], maxsteps=500, extraprec=4 * dps)

        # y = (2 - z - 1/z) / 4  ->  z^2 - (2 - 4y) z + 1 = 0; keep |z| < 1
        z_roots = []
        for y in y_roots:
            b = 2 - 4 * y
            d = mp.sqrt(b * b - 4)
            z1, z2 = (b + d) / 2, (b - d) / 2
            z_roots.append(z1 if abs(z1) < 1 else z2)

        poly = [mp.mpf(1)]  # ascending powers of z
        factors = [[1, 1]] * n + [[-z, 1] for z in z_roots]
        for f in factors:
            out = [mp.mpf(0)] * (len(poly) + 1)
            for i, p in enumerate(poly):
                out[i] += p * f[0]
                out[i + 1] += p * f[1]
            poly = out

        h = [mp.re(c) for c in poly]
        s = sum(h)
        h = [c * mp.sqrt(2) / s for c in h]
        return [float(c) for c in h[::-1]]


# --------------------------------------------------------------------------
# Destriping
# --------------------------------------------------------------------------
def remove_stripes_vertical(image, dec_num, wavelet, sigma):
    """Remove vertical stripes with combined wavelet-Fourier filtering.

    Parameters
    ----------
    image : 2D array
    dec_num : int
        Wavelet decomposition level.
    wavelet : str
        Wavelet name, e.g. 'db41'.
    sigma : float
        Gaussian damping parameter in Fourier space.
    """
    w = get_wavelet(wavelet)
    shape = image.shape
    approx = np.asarray(image, dtype=np.float64)

    # Wavelet decomposition
    details = []
    for _ in range(dec_num):
        approx, (ch, cv, cd) = pywt.dwt2(approx, w, mode="symmetric")
        details.append([ch, cv, cd])

    # Damp stripe information in the vertical detail bands
    for d in details:
        cv = d[1]
        my = cv.shape[0]
        k = np.arange(-(my // 2), -(my // 2) + my)
        damp = 1 - np.exp(-(k ** 2) / (2 * sigma ** 2))
        f = np.fft.fftshift(np.fft.fft(cv, axis=0), axes=0)
        f *= damp[:, None]
        d[1] = np.real(np.fft.ifft(np.fft.ifftshift(f, axes=0), axis=0))

    # Wavelet reconstruction
    out = approx
    for ch, cv, cd in reversed(details):
        out = out[: ch.shape[0], : ch.shape[1]]
        out = pywt.idwt2((out, (ch, cv, cd)), w, mode="symmetric")

    return out[: shape[0], : shape[1]]


def remove_stripes(image, dec_num, wavelet, sigma, vertical=True):
    """Remove vertical (default) or horizontal stripes."""
    if vertical:
        return remove_stripes_vertical(image, dec_num, wavelet, sigma)
    return remove_stripes_vertical(image.T, dec_num, wavelet, sigma).T


# --------------------------------------------------------------------------
# MATLAB-equivalent helpers
# --------------------------------------------------------------------------
def im2double(img):
    """Equivalent of MATLAB im2double."""
    img = np.asarray(img)
    if img.dtype == np.uint8:
        return img / 255.0
    if img.dtype == np.uint16:
        return img / 65535.0
    if img.dtype == np.int16:
        return (img.astype(np.float64) + 32768) / 65535.0
    return img.astype(np.float64)


def to_uint16(img):
    """Equivalent of MATLAB im2uint16(mat2gray(img)): min-max rescale to
    [0, 1], then convert to uint16."""
    img = np.asarray(img, dtype=np.float64)
    lo, hi = img.min(), img.max()
    scaled = (img - lo) / (hi - lo) if hi > lo else np.zeros_like(img)
    return np.round(scaled * 65535).astype(np.uint16)


def read_image(path):
    """Read a 2D grayscale TIFF."""
    import tifffile

    img = tifffile.imread(path)
    img = np.squeeze(img)
    if img.ndim != 2:
        raise ValueError(f"Expected a 2D grayscale image, got shape {img.shape}")
    return img


def pick_path(kind="file"):
    """Fallback file/folder picker (like uigetfile/uigetdir) when no path is
    given on the command line."""
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    if kind == "file":
        path = filedialog.askopenfilename(filetypes=[("TIFF", "*.tif *.tiff")])
    else:
        path = filedialog.askdirectory(title="Select the folder containing the images")
    root.destroy()
    if not path:
        raise SystemExit("No path selected.")
    return path
