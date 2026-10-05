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
- Recommended: ORS Dragonfly 3D World for visual inspection for slice alignment and filtering

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

Can be done with ORS Dragonfly 3D World. Refer to their documentation [here](https://www.theobjects.com/dragonfly/dfhelp/2021-1/Content/Viewing%20and%20Processing%20Images/Image%20Registration/Registering%20Images%20Automatically.htm) for more information.


### 1.2 Destriping

Stripes appear along the milling direction and from electron accumulation. Reduce them with the combined wavelet-Fourier filter (Münch et al., 2009).
The wavelet-FFT filtering is done with:
1. Decomposition of each slice with a 2D Daubechies wavelet transform.
2. Damping of stripe in the detail coefficients with a Gaussian filter in Fourier space.
3. Reconstruction of slice.

The scripts are in [`scripts/destriping/`](../scripts/destriping/):

| Script | Purpose |
|---|---|
| `01-BestParameters.py` | Saves destriped images for a grid of wavelet, level, and sigma values for visual comparison |
| `02-BestWavelet.py` | Scores Daubechies wavelets `db1` to `db43` on one image and picks the best |
| `03-WaveletFFT.py` | Applies chosen parameters to the whole image stack |
| `HelperFunctions.py` | Core filter and helper functions, imported by the other scripts |

```bash
cd scripts/destriping
pip install -r requirements.txt
```

> **Stripe direction.** All scripts are for reducing vertical stripes by default. If your goal is to reduce horizontal stripes, add `--horizontal`.

#### Parameters
| Parameter | Flag | Effect |
|---|---|---|
| Wavelet | `--wavelet` | Daubechies wavelet `dbN`. Higher orders are smoother and usually preserve structure better. |
| Decomposition level | `--level` | Levels 2 to 6 suits small-to-medium artifacts. Levels 7 to 10 handles large-scale artifacts but may smooth fine details. |
| Sigma | `--sigma` | Gaussian damping width. Sigma values 1 to 2 preserves fine textures and reduces mild stripes. Values 3 to 5 balances reduction and detail for most data. Values above 7 handles strong or irregular stripes but may blur details. |

#### Step 1: Choose a test slice
Pick one representative slice from the stack with clearly visible striping, and use it for Steps 1 and 2.

#### Step 2: Tune level and sigma
```bash
python 01-BestParameters.py /path/to/slice.tif \
 --wavelets db2 db15 db30 db42 \
 --levels 2 4 6 8 \
 --sigmas 1 2 4 6 8
```

The script destripes the slice with every combination of a few representative wavelets, levels, and sigmas, and saves one 16-bit PNG per combination, names `<wavelet>-L<level>-S<sigma>.png`, to a folder called `<slice>.tif_parameters_v` next to the slice's folder. Spreading the wavelets across low and high orders shows how level and sigma behave regardless of the wavelet.

Open the images side by side and choose the level and sigma that remove the stripes without blurring fine boundary structures like membranes. Check for:
- Remaining stripes, especially wide, low contrast bands
- Loss of fine detail such as membranes
- New artifacts

`TODO`:add a comparison grid -`images/parameter_sweep.png`

#### Step 3: Select the wavelet
Use the level and sigma chosen in Step 2:
```bash
python 02-BestWavelet.py /path/to/slice.tif --level <level chosen> --sigma <sigma chosen> --csv wavelet_scores.csv
```

The script destripes the slice with each wavelet from `db1` to `db43` and compares each result with the original image using three metrics:
- **SSIM** (structural similarity): higher means structure is better preserved
- **Relative energy loss**: lower means less signal is removed
- **PSNR** (peak signal-to-noise ratio): higher means the result is closer to the original

Each metric is min-max normalized across wavelets and combined into one score:
```
score = 0.5*SSIM + 0.3*(1 - energy loss) + 0.2*PSNR
```

The wavelet with the highest score is printed and marked on a four-panel plot saved as `<slice>_wavelet_selection.png`.
>

