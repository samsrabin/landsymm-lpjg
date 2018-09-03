# PLUM_harmonization

This repository contains code written to harmonize PLUM output data (land use areas) with historical inputs used in LPJ-GUESS.

It is based on the [code used in the Land Use Harmonization dataset version 1 (LUH1)](http://luh.umd.edu/code.shtml). Extensive changes have been made to generalize the code to harmonize more than just cropland and pasture, with the aim of harmonizing the area of *each crop* instead of just overall cropland area.

## Changes from original LUH1 code

#### Added features

+ Reads from and saves to LPJ-GUESS-compatible formats.
+ Generalized to harmonize any number of arbitrary land use types.
+ Can specify minimum "natural" area allowed in each grid cell, which allows for PLUM's `MIN_NATURAL_RATE` as well as protected area maps.
+ Debugging calculations and output can be turned on (for 2-degree gridcells, at least for now) by specifying the gridcell of interest in `debugIJ_2deg`. Setting `debugIJ_2deg = []` skips debugging calculations, making the processing much faster.

#### Fixes (see below for detailed explanation)

1. Unmet area resulting from excess amounts of a land use is now positive.
2. `j` resets to zero before each land use's "ring redistribution" loop.

## Other notes

+ The archive/ directory contains code from previous versions of this work. It probably will not function (correctly) without edits.

## Detailed explanation of "fixes"

1. IAM_preprocessing/IMAGE/new_grids2.m, lines 127–128. In line 127, the `1–GLMcrop_2deg` will result in `unmetcrop1` being negative wherever it's not zero (i.e., wherever `GLMcrop_2deg > 1`). My understanding from the rest of the code is that it should be positive. A positive `unmetcrop1` would result in `GLMcrop_2deg` being corrected to 1 in line 129, whereas a negative `unmetcrop1` just increases `GLMcrop_2deg` further above 1. A positive `unmetcrop1` would also be consistent with the fact that on line 142, `unmetcrop` is positive. (Line 128 has the same issue as Line 127, but for pasture instead of cropland---same for the other lines I mention here.)
2. IAM_preprocessing/IMAGE/new_grids2.m, line 206: `j` does not reset to 1 before pasture's `while` loop.