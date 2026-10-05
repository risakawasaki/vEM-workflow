# 01. Artifact Correction

## Goal
Remove acquisition artifacts from raw cryo-FIB/SEM image stacks, mainly slice-to-slice drift, stripes, and noise. The output is a clean, contrast-enhanced volume.

## Inputs
| Item | Format | Notes
|---|---|---|
| Raw image stack | `.tif` series or multipage `.tif` | One image per FIB slice, output of the imaging step |
| Acquisition metadata | voxel size (x, y, z) nm | Physical size of each pixel, slice thickness |

## Requirements
- Python 3.10+
- Packages: `numpy`, `scipy`, `scikit-image`, `PyWavelets`, `tifffile`
- Optional: ORS Dragonfly 3D World for visual inspection for slice alignment and filtering

```bash
pip install numpy scipy scikit-image PyWavelets tifffile
```

## Overview
 ```mermaid
flowchart LR
  A[Raw volume stack] --> B[1.1 Slice alignment]
  B --> C[1.2 Destriping]
  C --> D[1.3 Denoising]
  D --> E[1.4 Contrast enhancement]
  E --> F[Corrected volume]
```

## Procedure

### 1.1 Slice alignment
Correct lateral drift between consecutive slices.

1. Load the raw volume stack.
2. Register each slice to the previous one using SSD (sum-of-squared differences).
3. Crop to the region shared by all slices.

### 1.2 Destriping

Stripes appear along the milling direction and from electron accumulation. Reduce them with the combined wavelet-Fourier filter (Münch et al., 2009).
The wavelet-FFT filtering is done with:
1. Decomposition of each slice with a 2D Daubechies wavelet transform.
2. Damping of stripe in the detail coefficients with a Gaussian filter in Fourier space.
3. Reconstruction of slice.

The scripts are in [`scripts/destriping/`](../scripts/destriping/):

| Script | Purpose |
|---|---|
| `01-bestParameters.py
