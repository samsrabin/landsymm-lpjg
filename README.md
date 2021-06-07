## Calibrate LPJ-GUESS crop yields against FAO data

These MATLAB scripts and functions read LPJ-GUESS crop yield outputs and calibrate them against FAO yield data at the country level, according to a mapping of crop types specified by the user. The calibration factors are the result of a slope-only regression (i.e., with Y-intercept set to 0).

## Calibration versions

1. As Peter Alexander described
2. Looking at each crop separately. NOT ACTUALLY USEFUL FOR CALIBRATION
3. As 1, but with oil+cakes equivalent
4. As 1, but using oil crops as in Commodity Balance "Oilcrops"
5. As 4, but using Commodity Balance for Production
6. As 1, but with pre-processing as 5
7. As 6, but with oil+cakes equivalent
8. As 6, but with updated FAOSTAT data (2017-02-17)
9. As 5, but with Peter's old FAOSTAT data (2015-12)
10. As 9, but without rye
11. As 9, but with Stijn's new crops
12. As 10, but with PLUM6
13. As 9, but with only exact crop matches (for pure6 runs)
14. As 13, but without Soybeans or Pulses (for pure4 runs)
15. ~9 for FAO data, plus Rye, and with PLUM6 types (PLUM6xtra)
16. As 15, but correctly with Rye instead of Oats in Products


## Master script

In addition to the included files, you also need a master script. It should look something like this:

```matlab
%% Information about this calibration run

verName_calib = 'remap2_PLUM6xtra_WWSW' ;   % calib_ver = 16
filename_guess_yield = '/Volumes/WDMPP_Storage/PLUM/trunk_runs_external/calib.remap2.PLUM6xtra_WWSW.1901-2005/output-2018-02-16-185453/yield.out.gz' ;
filename_guess_landuse = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
filename_guess_cropfrac = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.assignWWorSW_cruncep.txt' ;


%% Other options

% Calibration version
calib_ver = 16 ;

% Years for calibration
year1 = 1995 ;
yearN = 2005 ;

% Paths to calibration code and data
dir_code = '/Users/Shared/PLUM/crop_calib_code/' ;
dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;

% Path to figure output directory
dir_outfigs = '/Users/Shared/PLUM/crop_calibration_for_PLUM/figures/' ;


%% Do it

% Add code files to path (just for this session)
addpath(genpath(dir_code))

crop_calibration


%% Save figure

out_file = [dir_outfigs verName_calib '_v' num2str(calib_ver) '.pdf'] ;
export_fig(out_file,'-r300')
```






