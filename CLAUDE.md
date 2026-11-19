# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Bachelor's thesis work in optics: reconstructing corneal/lens surface topography from Hartmann-Shack spot
patterns (photographs of a grid of light spots projected onto a cornea/lens, distorted by its shape). There
is no application code, package, or CLI — the repository is a collection of Jupyter notebooks, each
representing one experiment or analysis session, plus their input images and output data files.

## Commands

There is no build, lint, or test tooling — no `requirements.txt`/`environment.yml`, no test suite, no CI.
Work happens by running notebooks directly in Jupyter (`jupyter lab` / `jupyter notebook`) or via VS Code's
notebook UI, executing cells top to bottom.

Libraries used across notebooks (install manually, no pinned versions exist anywhere):
`numpy`, `pandas`, `matplotlib`, `scipy`, `opencv-python` (`cv2`), `scikit-learn` (`sklearn.cluster.KMeans`).

Several notebooks set `rc('text', usetex=True)` for plot labels (Times New Roman via LaTeX). This requires a
working LaTeX installation on the machine; without it, any cell that plots will raise a LaTeX/`latex.exe not
found` error unrelated to the actual analysis code.

## Architecture: the experiment pipeline

Each dated top-level folder (`11-11-17x17Analysis`, `22-11-AlgorithmTest`, `23-11-KeratoconusMC`,
`23-10-ToricLens`, `29-10-tests`, `17-11-FourierAnalysis`, ...) is one self-contained experiment, named
`DD-MM-<description>`. The `codefiles/` subfolder holds the notebooks, `data/` holds centroid text files,
`images/` holds source `.tif`/spot-pattern images. Loose notebooks and `.txt` files at the repo root
(`Campbell_Implementation.ipynb`, `imageProcessing.ipynb`, `codefileToroidalSurf.ipynb`, `output*.txt`, etc.)
are earlier/scratch iterations of the same pipeline, superseded by the versions inside the dated folders —
don't treat them as the canonical implementation.

Within an experiment, work flows through two notebook types, always run in this order:

1. **`imgprocessing*.ipynb`** (image processing): loads a grayscale `.tif` spot-pattern image, applies an
   intensity threshold (`cv.threshold`) to isolate spots, then clusters pixels into blobs (`KMeans`, roughly
   one cluster per expected spot — e.g. 252 spots for a 17x17/19x19 pattern) to find each spot's centroid.
   Writes results as a whitespace-delimited `.txt` file with columns `pos, x, y` (pixel coordinates) into the
   experiment's `data/` folder. This step is done twice per experiment: once for the **reference** pattern
   (an undistorted/flat surface) and once for the **test** pattern (the actual eye/lens surface).

2. **`*Analysis.ipynb` / `codefileRevSurf*.ipynb`** (reconstruction): loads both centroid `.txt` files with
   pandas, converts pixel coordinates to physical distance via a fixed `delta` scale factor (mm/pixel,
   derived from the known spacing of the projector grid). The flat list of centroids is then hand-mapped by
   index into a square grid matrix (typically 19x19) representing the spots' spatial layout — the mapping is
   literal per-index assignment (`mat_R_X[row][col] = df["x_d"][i]`) reflecting the physical (often
   circular/diamond-clipped) aperture of the instrument, not a generic reshape. Zero entries in the grid mean
   "no spot here" (outside the aperture).

   From the per-node displacement between the test grid and the reference grid, the notebooks compute local
   wavefront slopes and then reconstruct the surface (elevation/curvature) using and comparing several
   independent zonal integration methods implemented from scratch: **Southwell**, **Campbell**, and
   **Malacara**. Downstream results computed the same way in every experiment notebook: astigmatism,
   max/min curvature directions, focus error, curvature maps, and an elevation map per method — intended to
   be cross-checked against each other for the same input pattern.

`CampbellAlgorithm/` and `Campbell_Implementation.ipynb` hold a standalone/isolated implementation of the
Campbell method used to validate that piece independently of a full experiment notebook.

## Gotchas when editing notebooks

- Many notebooks read input files via **hardcoded absolute Windows paths**
  (e.g. `C:/Users/PC/Desktop/UN/OpticsWork/codefiles/thesisreproduction/...`) rather than relative paths, even
  though relative paths (`./../data/...`) are used elsewhere in the same notebook. Check which style a given
  notebook uses before assuming `data/`/`images/` paths are portable.
- The grid-index mapping (spot `pos` N → grid cell `[row][col]`) is duplicated by hand across every analysis
  notebook for both the X and Y matrices, and separately for the test pattern vs. the reference pattern. When
  fixing a bug in this mapping (as in the Southwell/Campbell/Malacara fix commits), check whether the same
  hand-written block needs the identical fix in the sibling notebooks/cells — it is not shared code.
