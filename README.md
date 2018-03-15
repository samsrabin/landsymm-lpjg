## Calibrate LPJ-GUESS crop yields against FAO data

These MATLAB scripts and functions read LPJ-GUESS crop yield outputs and calibrate them against FAO yield data at the country level, according to a mapping of crop types specified by the user. The calibration factors are the result of a slope-only regression (i.e., with Y-intercept set to 0).


## Master script

In addition to the included files, you also need a master script. It should look something like this:

```matlab
%% Information about this calibration run

version_name = 'remap2_PLUM6xtra_WWSW' ;   % calib_ver = 16
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

out_file = [dir_outfigs version_name '_v' num2str(calib_ver) '.pdf'] ;
export_fig(out_file,'-r300')
```






