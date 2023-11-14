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