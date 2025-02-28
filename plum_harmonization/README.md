# PLUM_harmonization

This repository contains code written to harmonize PLUM output data (land use areas) with historical inputs used in LPJ-GUESS.

It is based on the [code used in the Land Use Harmonization dataset version 1 (LUH1)](http://luh.umd.edu/code.shtml). Extensive changes have been made to generalize the code to harmonize more than just cropland and pasture, with the aim of harmonizing the area of *each crop* instead of just overall cropland area.

Here's how this works. We begin with the values in year0 (here, 2010) from one dataset (here, LUH1), then apply the year0 to year1 changes (deltas) from the second dataset (here, PLUM). Grid cells can reach limits: The deltas might specify loss of cropland when the grid cell is already 0% cropland, or likewise the deltas might specify expansion of cropland when it's already 100% cropland. In such cases, the algorithm looks for space to apply the remaining "unmet" deltas in the grid cells bordering the cell in question. It expands the radius of this search until all deltas are satisfied.

This code has another feature not present in the LUH1 harmonization, which is to harmonize fertilizer and irrigation (the latter in arbitrary units of intensity\*area). The way it does this is analogous to how it treats land use area. Limits for irrigation in a cell are 0 and 1. The lower limit for fertilizer for any given crop is also zero, but the upper limit varies. It is either the maximum seen for the crop in any gridcell in any PLUM output thus far (i.e., if we're working on deltas for 2020--2021, consider PLUM outputs from 2021 and before), or in LUH2 during or before the base year (here, 2010), or in any harmonized output thus far (although, because of the first two rules, this should never come into play).

## Changes from original LUH1 code

#### Added features

+ Generalized to harmonize any number of arbitrary land use types.
+ Harmonizes not just area but also fertilizer and irrigation.
+ Reads from and saves to LPJ-GUESS-compatible formats.
+ Can specify minimum "natural" area allowed in each grid cell, which allows for PLUM's `MIN_NATURAL_RATE` as well as protected area maps.
+ Added optional debugging calculations and output, globally or for gridcells.

#### Fixes (see below for detailed explanation)

1. Unmet area resulting from excess amounts of a land use is now positive.
2. `j` resets to zero before each land use's "ring redistribution" loop.
3. Now allows cells to absorb as much of their own unmet as possible before going to a ring.
4. In distribution to half-degree cells, available land is updated after every land use.

## Other notes

+ The `archive/` directory contains code from previous versions of this work. It probably will not function (correctly) without edits.

## Detailed explanation of "fixes"

1. Original code's `IAM_preprocessing/IMAGE/new_grids2.m`, lines 127–128. In line 127, the `1–GLMcrop_2deg` will result in `unmetcrop1` being negative wherever it's not zero (i.e., wherever `GLMcrop_2deg > 1`). My understanding from the rest of the code is that it should be positive. A positive `unmetcrop1` would result in `GLMcrop_2deg` being corrected to 1 in line 129, whereas a negative `unmetcrop1` just increases `GLMcrop_2deg` further above 1. A positive `unmetcrop1` would also be consistent with the fact that on line 142, `unmetcrop` is positive. (Line 128 has the same issue as Line 127, but for pasture instead of cropland---same for the other lines I mention here.)
2. Original code's `IAM_preprocessing/IMAGE/new_grids2.m`, line 206: `j` does not reset to 1 before pasture's `while` loop.
3. LUH1 code always started a ring if a cell had any unmet, even if all the unmet could be met by the cell itself.
4. Can be disabled by setting `update_avail_land = false` (not tested).

## Code overview
Whereas the first part of this README provided a high-level overview of how the code works, this section is intended to walk through the code step-by-step at a medium level. Error checking and debugging steps are excluded here, as are steps related to MATLAB array manipulation etc.

When looking through the code, you'll see lots of arrays with the suffix `_YXv`. This indicates that the dimensions of the array are (latitude, longitude, variable). "Variable" here usually refers to different land use types. In `agri_YXv` arrays, "variable" comprises the crop types plus one pasture type. in `nfert_YXv` and `irrig_YXv` arrays, "variable" just refers to the crop types.

Arrays also often have `y0` or `y1` in their names. Those refer to years `N-1` and `N`, respectively, in the following description.

### Initial setup
1. Read options.
2. Import reference data (lat/lon map, land area, area ineligible for land use).
3. Import baseline land use areas and management levels---i.e., the LUH2 or HILDA+ maps from the "baseline" year (call it year `B`) against which PLUM will be harmonized. (These maps have already been processed for use in LPJ-GUESS.)
4. Ensure that all imported maps have the same land/ocean mask.
5. Define dictionaries to map between LPJ-GUESS's crop names and PLUM's crop names.

### Loop through all PLUM directories
For each, e.g., SSP that PLUM processed...

### Loop through all PLUM years
Starting with the year after the baseline (`B+1`), for each year `N`:
1. (If `N == B+1`: Import all PLUM outputs for year `N-1`.)
2. Import PLUM area outputs for year `N`.
3. For each land use type, calculate the area changes ("deltas") in the PLUM outputs between years `N-1` and `N`.
4. Get limits: (1) Area available for land use given harmonized (or, if `N == B+1`, LUH2/HILDA+ baseline) maps for year `N-1`. (2) Max allowable N fertilization rate.
5. Coarsen all maps from 0.5° to 2°.
6. Apply area deltas to harmonized 2° map from year `N-1`.
7. Handle invalid area changes: Loop through every 2° gridcell. If a gridcell has "unmet" crop and/or pasture area change (i.e., PLUM says to expand area beyond what is actually available or to reduce area by more than what the gridcell actually has), look for available space in its neighboring gridcells. Keep expanding that search ring as needed until all unmet area change has been displaced to new 2° cells.
8. Apply management deltas to harmonized 2° map from year `N-1`.
9. Handle "unmet" management changes (N fertilization, irrigation intensity) using the same basic algorithm as step (7).
10. Loop through all 2° gridcells and distribute the new crop and pasture area changes to the 0.5° gridcells within: If the crop or pasture change is a decrease, apply the same PERCENTAGE decrease to all 0.5° crop or pasture gridcells within the 2° cell. If the crop or pasture change is an increase, apply this increase to all 0.5° gridcells within the 2° cell proportionally according to available land area. Make sure that total area within a 0.5° gridcell does not exceed 1 or 0.
11. Write output files with harmonized data for year `N`.
