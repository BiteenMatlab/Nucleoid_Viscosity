# Macrodomain Mask Workflow (`Macrodomain_mask`)

This README explains how to run the scripts in `Macrodomain_mask` in sequence, including data loading, parameter setup, and expected outputs.

## 1. Prerequisites

- Tested on MATLAB R2024a
- Scripts use helper functions from `Supporting_functions`, add it to path.
- [200 colormap](https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap) from MATLAB Central File Exchange was used for plotting

## 2. Recommended Run Sequence

1. `Locus_filter.mlx`
2. `Heatmap_plot.mlx`
3. `ExtractRegion_HeatMap.mlx`

## 3. Script Details

## 3.1 `Locus_filter.mlx`

Purpose:
- Load raw `PhaseMask` and `tracks` from the locus-tracking dataset.
  - Example dataset can be found at DOI [10.5281/zenodo.18598888](https://zenodo.org/records/18598888)
- Compute cell morphology and apply morphology filter.
- Compute loci information and apply locus-number,  pole-distance and circularity filters.
- Save filtered trajectories.

What to load:
- Select `Folder` containing:
  - `*PhaseMask.mat` (`PhaseMask`)
  - `*fits.mat` (`tracks`)

Parameters to set:
- `pixel_size` (default `49` nm/pixel)
- Morphology filter:
  - `Bndry_condtion` (`1` long axis, `2` short axis, `3` aspect ratio, `4` area)
  - `low_bd`, `up_bd`
- Cluster threshold:
  - `Thres` (minimum localizations per cluster, default `50`)
- Locus filters:
  - `Thres_LocNum` (default `[1,4]`)
  - `Thres_dis2pol` (default `[250,1800]` nm)
  - `Thres_circle` (default `[1,3]`)

Outputs (saved under `Folder\Results_Filter`):
- `Cell Morphology.mat/.png`
- `Length_Filter.mat`
- `Locus_num.mat`
- `Locus_info.mat/.png`
- `newTracks.mat`
- `Locus_Filter.mat`
- `Tracks_passFilter.mat`

---

## 3.2 `Heatmap_plot.mlx`

Purpose:
- Build localization heatmaps from filtered tracks.
- Create non-symmetrized and symmetrized maps.

Important variable dependency:
- This script expects `Folder` from the previous script to still exist in workspace as `...\Results_Filter`.
- It derives:
  - `Filefolder = Folder(1:end-15)` (the dataset root)

What it loads:
- `*fits.mat` and `*PhaseMask.mat` from `Filefolder`
- `Tracks_passFilter.mat` and `Length_Filter.mat` from `Filefolder\Results_Filter`

Parameters to set:
- `aspt_ratio` (based on the aspect ratio from nanocage-tracking dataset at the same condition)
- `gap_y` ($\frac{1}{num\_pix\_y}$, default $\frac{1}{20}$)
- `pixel_size` (pixel size)

Outputs (saved under `Filefolder\Heatmaps`):
- `Localization_HeatMap_noSym.png/.svg` and `Result_cat.mat` (non-symmetrized)
- `Localization_HeatMap_filtering.png/.svg` and `Result_cat_sym.mat` (symmetrized)

Workspace variables used by next step:
- `Result_cat_sym`, `gap_y`, `aspt_ratio`, `Filefolder`, `low_bd`, `up_bd`

---

## 3.3 `ExtractRegion_HeatMap.mlx`

Purpose:
- Extract boundary points of a macrodomain.


Parameters to set:
- `maskArea_ratio` (default `0.10`, meaning the bins with the highest localization density that are added to 0.1 are used as proxy for macrodomain region)
- `scale_rate` (default `10`; scaling factor for Kernel Density Estimation)
- DBSCAN settings used for boundary grouping:
  - `eps = 0.01`, `minpts = 3`


Outputs:
- Coordinates of grids used as macrodomain mask
  - Change filename per domain run (for example `Ori_boundary.mat`, `Ter_boundary.mat`, `Left_boundary.mat`, `Right_boundary.mat`).
- Heatmap overlaid with macrodomain boundary  (`Heatmap_ori.svg` ).

---


## 4. Practical Notes

- Some paths are hard-coded to the author environment; update these before running:
  - `uigetdir(...)` default locations
  - save/load folders in `ExtractRegion_HeatMap.mlx`
- Repeat `ExtractRegion_HeatMap.mlx` for each locus dataset and save with domain-specific filenames.
- Boundary files from this folder are used by `Acc_Vis_pipeline/Vis_Nuc_Macrodomain.mlx`.
