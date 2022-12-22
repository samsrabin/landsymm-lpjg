# Processing land-use input for use in GGCMI-to-PLUM LPJ-GUESS runs



## The steps

### Remapping land use

This takes place in the `remap_*.m` script (henceforth, "remap script") .

For traceability and reproducibility, it's important that you create a new remap script for each new set of inputs and/or processing method you want to use. A short version name is first specified as `remapVer` at the top of the remap script. As a best practice, you should name that script to reflect the version. So, e.g., version `8b2oil` is specified and first processed in script `remap_v8b2oil.m`.

What we're doing here is translating land use areas (from LUH2), crop fractions (from MIRCA2000, incl. irrigated vs. rainfed), fertilizer intensity (from AgGRID_nutrient_input_v1.1), and optionally manure application (from Zhang et al., 2017) into a form LPJ-GUESS can use. 

I refer to this as "remapping" because of the step that "re-maps" the MIRCA crop types to LPJ-GUESS crop types. This remapping, or translation, is defined in the call to `get_remapv2_keys()`. Note that, in addition to `remapVer`, you must specify another version name called `thisVer`, which can be more verbose. `thisVer` gets used in `get_remapv2_keys()` here and in various other functions in this codebase.

v11 and beyond: You must specify at least one gridlist file in `gridlist_files` at the top of the remap script. The output land use files, as well as the output `gridlist.remapv${remapVer}.txt`, will contain any gridcell present in any of the input gridlist files, excluding those which ended up not having data in the input land use files. In addition, the remap script will save copies of the original gridlist files, again without those that got excluded.

### Generating crop calendar files

This is pretty easy. Specify `remapVer` and `thisVer` at the top of `g2p_growing_seasons.m` and run it to produce crop calendar files that can be used by LPJ-GUESS for `file_sdates`, `file_hdates`, `file_growseaslength_in`, and `file_Nfertdate2_in`.

This script follows the GGCMI Phase 2 protocol where possible:

- Planting and harvest dates are taken from the Phase 1 data (`AGMIP_GROWING_SEASON.HARM.version1.25`), except for wheat which uses Phase 2 (`AGMIP_GROWING_SEASON.HARM.version2.0`). Any cells in the gridlist but missing from either of those use the Phase 3 growing seasons, as Christoph used in his emulation runs.
- The first fertilization happens at sowing. Second fertilization happens at 40 days after sowing for spring-planted crops, or as specified in `wwh_rf_2nd_fertilizer_days_disseminate_v2.nc4` for winter wheat.
- Crop calendars are identical for rainfed and irrigated crops, in accordance with Phase 2 protocol. The files produced only have one column for use in both rainfed and irrigated crop types; this requires the use of an LPJ-GUESS version that supports the specification of PFT parameter `cropphen_col` (or the use of "simplified PFT lists", where multiple crop stands can use the same PFT but still output separate yields).

### Producing alternate land use area files

`make_justCPB.m` produces a one-year `file_lu` that takes all the non-BARREN area and splits it evenly between CROPLAND and pasture. This is used during "potential" runs (the ones with the factorial crop\*Nfert\*irrigation stands).

`make_someOfAllCrops.m` produces `file_lu` and `file_lucrop` for use in all "actual" runs contributing to potential yield generation for PLUM. It ensures that there is at least a tiny amount of every crop in every gridcell.

### Producing ins-files

The script `print_stand_pft_lists_g2p.m` will save ins-files for stand and crop PFT specifications as needed for a given remap version. 

`print_stand_simplePFT_lists_g2p.m` does the same, but for use in a "simplified PFT lists" setup.

