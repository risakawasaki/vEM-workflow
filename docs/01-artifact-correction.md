# 01. Artifact Correction

## Goal
Remove acquisition artifacts from raw cryo-FIB/SEM image stacks, mainly slice-to-slice drift, stripes, and noise. The output is a destiriped, noise reduced, contrast-enhanced volume.

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

Can be done with [ORS Dragonfly 3D World](https://dragonfly.comet.tech/).<br>
Refer to their documentation [here](https://www.theobjects.com/dragonfly/dfhelp/2021-1/Content/Viewing%20and%20Processing%20Images/Image%20Registration/Registering%20Images%20Automatically.htm) for more information on slice alignment using the software.


### 1.2 Destriping

Stripes appear along the milling direction and from electron accumulation. Reduce them with the combined wavelet-Fourier filter (Münch et al., 2009). Check out the paper [here](https://opg.optica.org/oe/fulltext.cfm?uri=oe-17-10-8567).
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

The script destripes the slice with each wavelet from `db1` to `db43` and scores each result with four metrics:

| Metric | What it measures | Better |
|---|---|---|
| **Stripe reduction** | How much of the stripe signal was removed | Higher |
| **SSIM** (structural similarity) | How well structure is preserved compared with the original image | Higher |
| **Relative energy loss** | How much signal was removed overall | Lower |
| **PSNR** (peak signal-to-noise ratio) | How close the result is to the original | Higher |

**How stripe reduction is measured** The script averages the image along the stripe direction. Tissue and cell structure varies along the stripe direction and averages out, but the stripes don't, so this profile mostly contains stripe signal. The profile is high-pass filtered to remove slow image intensity gradients, and its variance is the stripe strength:
```
stripe reduction = 1 - stripe strength (destriped) / stripe strength (original)
```

The other three metrics compare the result with the original, striped image. They penalize wavelets that distort the data, while stripe reduction rewards removing the stripes. The two pull in opposite directions, so the score balances them. Each metric is min-max normalized across wavelets and combined:
```
score = 0.4 ✕ stripe reduction + 0.3 ✕ SSIM + 0.2 ✕ (1 - energy loss) + 0.1 ✕ PSNR
```

The wavelet with the highest score is printed and marked on a five-panel plot saved as `<slice>_wavelet_selection.png`.

| Option | Default | Description |
|---|---|---|
| `--weights` | `0.4 0.3 0.2 0.1` | Weights for stripe reduction, SSIM, energy loss, and PSNR |
| `--stripe-scale` | `50` | Width in pixels of the background removed from the image. Stripe much wider than this are treated as background, so raise it if your stripes are wide. |
| `--csv` | none | Save the full results table |

> **Interpreting the scores.** Min-max normalization stretches small differences. For example, wavelets with stripe reduction of 99.6% and 97% can end up far apart in the score. Look at the raw values in the table or plot as well as the final score, and confirm the selected wavelet visually on the test slice before processing the full stack. The stripe metric also assumes stripes run the full height of the slice; partial-height stripe is measured, but diluted.

`TODO`: add the plot from your data → `images/wavelet-selection.png`
`TODO`: add a before/after comparison → `images/destriping-before-after.png`

#### Step 4: Apply to the full stack
```bash
python 03-WaveletFFT.py /path/to/aligned_stack --wavelet <wavelet chosen> --level <level chosen> --sigma <sigma chosen>
```

The script processes every `.tif` file in the folder and saves the results with the same filenames to a sibling folder named `<stack>_wfft_v_<wavelet>_L<level>_S<sigma>`, so the settings used are recorded in the folder name. Unreadable files are skipped with a warning.

> **Output intensity.** Each slice is rescaled to the full 16-bit range on its own (min-max). Absolute intensities are not comparable between slices afterwards, so apply any intensity normalization across the stack after this step.
>  **Wavelets above db38.** PyWavelets only includes Daubechies wavelets up to db38. `HelperFunctions.py` computes higher orders (such as db41) itself, using the same construction as MATLAB, so results are consistent with the original MATLAB implementation introduced by [Münch et al., 2009](https://opg.optica.org/oe/fulltext.cfm?uri=oe-17-10-8567).

### 1.3 Denoising and contrast enhancement
Apply a median filter to suppress background salt-and-pepper noise while keeping membrane edges. Then, apply unsharp masking to sharpen membranes boundaries. For easy comparison of different parameters for each filter, and exploring other filters of interest, [ORS Dragonfly 3D World](https://dragonfly.comet.tech/) is highly recommended.

If ORS Dragonfly 3D World is unavailable, you can use python APIs such as `scipy` and `scikit-image`:
```python
from scipy.ndimage import median_filter
from skimage.filters import unsharp_mask

denoised_img = median_filter(img, size = <size>)
contrast_enhanced_img = unsharp_mask(denoised_img, radius = <radius>, amount = <amount>)
```
Refer to [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.median_filter.html) for `median_filter`.<br>
Refer to [here](https://scikit-image.org/docs/0.25.x/auto_examples/filters/plot_unsharp_mask.html) for `unsharp_mask`.

### 1.4 Save the filtered volume.
ORS Dragonfly 3D World has export options for saving filtered volumes. <br>

Otherwise, after filtering, you can save the filtered volume like below:
```python
import tifffile
tifffile.imwrite("filtered_volume.tif", volume.astype(np.uint16))
```

## References

- Münch B, Trtik P, Marone F, Stampanoni M (2009). Stripe and ring artifact removal with combined wavelet–Fourier filtering. *Optics Express* 17(10):8567–8591.
- Daubechies I (1992). *Ten Lectures on Wavelets*. SIAM.

---
**Next step:** [2. Segmentation](02-segmentation.md)
