%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process crop growing seasons for GGCMI-to-PLUM sims %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses rainfed seasons only, in accordance with Phase 2 protocol.

% remapVer = '11_g2p' ; thisVer = 'WithFruitVeg_sepSugar_sepOil_sepC3' ;
% remapVer = '12_g2p' ; thisVer = 'ggcmi5' ;
remapVer = '13_g2p' ; thisVer = 'ggcmi5_preBNF' ;


%% Setup

addpath(genpath('/Users/sam/Documents/git_repos/g2p_emulation'))

% Get crop list
cropList_lpjg = get_remapv2_keys(thisVer) ;
Ncrops = length(cropList_lpjg) ;

inDir_phase1 = '/Volumes/Reacher/GGCMI/AgMIP.input/other.inputs/AGMIP_GROWING_SEASON.HARM.version1.25' ;
inDir_phase2 = '/Volumes/Reacher/GGCMI/AgMIP.input/other.inputs/AGMIP_GROWING_SEASON.HARM.version2.0' ;
inDir_phase3 = '/Volumes/Reacher/GGCMI/AgMIP.input/phase3/ISIMIP3/crop_calendar' ;

thisDir = sprintf('/Volumes/Reacher/G2P/inputs/LU/remaps_v%s', remapVer) ;
if ~exist(thisDir, 'dir')
    error('Working directory not found: %s', thisDir)
end

gridlist = lpjgu_matlab_read2geoArray( ...
    sprintf('%s/gridlist.remapv%s.txt', thisDir, remapVer)) ;
Ncells = length(gridlist.list2map) ;


%% Import

disp('Importing...')

cropList_as_gsType = cell(Ncrops,1) ;
sdate_xc = nan(Ncells, Ncrops) ;
hdate_xc = nan(Ncells, Ncrops) ;
gslen_xc = nan(Ncells, Ncrops) ;
fertDate2_xc = nan(Ncells, Ncrops) ;
source_sdate_xc = ones(Ncells, Ncrops) ;
source_hdate_xc = ones(Ncells, Ncrops) ;
source_gslen_xc = ones(Ncells, Ncrops) ;
for c = 1:Ncrops
    thisCrop_lpjg = cropList_lpjg{c} ;
    
    % Get corresponding growing-season type
    thisCrop_whe = '' ;
    switch thisCrop_lpjg
        case {'CerealsC3', 'CerealsC3s', ...
                'StarchyRoots', 'FruitAndVeg', ...
                'Sugarbeet', 'OilOther', 'ExtraCrop'}
            thisCrop_ph1 = 'spring_wheat' ;
        case {'CerealsC3w'}
            thisCrop_ph1 = 'winter_wheat' ;
        case {'CerealsC4', 'Miscanthus', 'Sugarcane'}
            thisCrop_ph1 = 'Maize' ;
        case {'Oilcrops', 'Pulses', 'OilNfix'}
            thisCrop_ph1 = 'Soybeans' ;
        case {'Rice'}
            thisCrop_ph1 = 'Rice' ;
        otherwise
            error('What crop from inputs should I use for %s?', thisCrop_lpjg)
    end
    cropList_as_gsType{c} = thisCrop_ph1 ;
    
    % If previously imported, just use that
    I_prev = find(strcmp(cropList_as_gsType, thisCrop_ph1)) ;
    if length(I_prev) > 1
        sdate_xc(:,c) = sdate_xc(:,I_prev(1)) ;
        hdate_xc(:,c) = hdate_xc(:,I_prev(1)) ;
        gslen_xc(:,c) = gslen_xc(:,I_prev(1)) ;
        fertDate2_xc(:,c) = fertDate2_xc(:,I_prev(1)) ;
        source_sdate_xc(:,c) = source_sdate_xc(:,I_prev(1)) ;
        source_hdate_xc(:,c) = source_hdate_xc(:,I_prev(1)) ;
        source_gslen_xc(:,c) = source_gslen_xc(:,I_prev(1)) ;
        continue
    end
    
    % Import phase 1 (or phase 2, for wheat)
    sdate_x = import_phase12(thisCrop_ph1, 'sdate', ...
        inDir_phase1, inDir_phase2, gridlist) ;
    hdate_x = import_phase12(thisCrop_ph1, 'hdate', ...
        inDir_phase1, inDir_phase2, gridlist) ;
    gslen_x = import_phase12(thisCrop_ph1, 'gslen', ...
        inDir_phase1, inDir_phase2, gridlist) ;
    
    gslen_sanity_check(sdate_x, hdate_x, gslen_x)
    
    % Import 2nd fertilization dates (winter wheat only)
    if strcmp(thisCrop_ph1, 'winter_wheat')
        thisFile = sprintf('%s/wwh_rf_2nd_fertilizer_days_disseminate_v2.nc4', ...
            inDir_phase2) ;
        tmp_YX = transpose(ncread(thisFile, '2nd fert appl')) ;
        fertDate2_x = tmp_YX(gridlist.list2map) ;
        clear tmp_YX
    end
        
    % Source tracking
    if contains(thisCrop_ph1, 'wheat')
        source_sdate_xc(:,c) = 2 ;
        source_hdate_xc(:,c) = 2 ;
        source_gslen_xc(:,c) = 2 ;
    end
    source_x = source_sdate_xc(:,c) ;
    source_x(isnan(sdate_x)) = 3 ;
    source_sdate_xc(:,c) = source_x ;
    source_x = source_hdate_xc(:,c) ;
    source_x(isnan(hdate_x)) = 3 ;
    source_hdate_xc(:,c) = source_x ;
    source_x = source_gslen_xc(:,c) ;
    source_x(isnan(gslen_x)) = 3 ;
    source_gslen_xc(:,c) = source_x ;
    clear source_x
    
    % Gap-fill with phase 3, if needed
    gapfill_these_x = isnan(sdate_x + hdate_x + gslen_x) ;
    if any(gapfill_these_x)
        % Unify NaN mask
        sdate_x(gapfill_these_x) = NaN ;
        hdate_x(gapfill_these_x) = NaN ;
        gslen_x(gapfill_these_x) = NaN ;
        % Gapfill
        sdate_x = gapfill_with_phase3( ...
            thisCrop_ph1, 'sdate', sdate_x, ...
            inDir_phase3, gridlist) ;
        hdate_x = gapfill_with_phase3( ...
            thisCrop_ph1, 'hdate', hdate_x, ...
            inDir_phase3, gridlist) ;
        gslen_x = gapfill_with_phase3( ...
            thisCrop_ph1, 'gslen', gslen_x, ...
            inDir_phase3, gridlist) ;
    end
    
    gslen_sanity_check(sdate_x, hdate_x, gslen_x)
    
    % 2nd fert, all spring-planted crops: 40 days after sowing
    if ~strcmp(thisCrop_ph1, 'winter_wheat')
        fertDate2_x = rem(sdate_x + 40, 365) ;
    end
    
    % Make sure no NaNs remain
    if any(isnan(sdate_x))
        error('NaN remain in sdate_x')
    end
    if any(isnan(hdate_x))
        error('NaN remain in hdate_x')
    end
    if any(isnan(gslen_x))
        error('NaN remain in gslen_x')
    end
    if any(isnan(fertDate2_x))
        error('NaN remain in fertDate2_x')
    end
    
    % Handle non-integers
    sdate_x = handle_non_integers(sdate_x) ;
    hdate_x = handle_non_integers(hdate_x) ;
    gslen_x = handle_non_integers(gslen_x) ;
    fertDate2_x = handle_non_integers(fertDate2_x) ;
        
    % Convert dates to zero-based, for LPJ-GUESS
    sdate_x = sdate_x - 1 ;
    hdate_x = hdate_x - 1 ;
    fertDate2_x = fertDate2_x - 1 ;
    
    % Check for negatives
    if any(sdate_x < 0)
        error('Negative in sdate_x')
    end
    if any(hdate_x < 0)
        error('Negative in hdate_x')
    end
    if any(gslen_x < 0)
        error('Negative in gslen_x')
    end
    if strcmp(thisCrop_ph1, 'winter_wheat') && any(fertDate2_x(fertDate2_x~=-100) < 0)
        error('Negative in fertDate2_x')
    end
    
    % Save to main arrays
    sdate_xc(:,c) = sdate_x ;
    hdate_xc(:,c) = hdate_x ;
    gslen_xc(:,c) = gslen_x ;
    fertDate2_xc(:,c) = fertDate2_x ;
    
end

% Make sure no NaNs remain (might have been introduced by forgetting to
% include a new array in "If previously imported, just use that" step)
if any(any(isnan(sdate_xc)))
    error('NaN in sdate_xc')
end
if any(any(isnan(hdate_xc)))
    error('NaN in hdate_xc')
end
if any(any(isnan(gslen_xc)))
    error('NaN in gslen_xc')
end
if any(any(isnan(fertDate2_xc)))
    error('NaN in fertDate2_xc')
end

% Make sure no non-integer(s) remain
if any(any(sdate_xc ~= round(sdate_xc)))
    error('Non-integer(s) remain in sdate_xc')
end
if any(any(hdate_xc ~= round(hdate_xc)))
    error('Non-integer(s) remain in hdate_xc')
end
if any(any(gslen_xc ~= round(gslen_xc)))
    error('Non-integer(s) remain in gslen_xc')
end
if any(any(fertDate2_xc ~= round(fertDate2_xc)))
    error('Non-integer(s) remain in fertDate2_xc')
end

% Make sure no negatives exist
if any(any(sdate_xc < 0))
    error('Negative(s) in sdate_xc')
end
if any(any(hdate_xc < 0))
    error('Negative(s) in hdate_xc')
end
if any(any(gslen_xc < 0))
    error('Negative(s) in gslen_xc')
end
if any(fertDate2_xc(fertDate2_xc~=-100) < 0)
    error('Negative(s) in fertDate2_xc (other than no-data -100)')
end

disp('Done.')


%% Save

%%% Options %%%%%%%%%%%%%
outPrec = 0 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
do_gzip = false ;
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Saving...')

out_header_cell = [{'Lon', 'Lat'} cropList_lpjg] ;

S_out.list2map = gridlist.list2map ;
S_out.lonlats = gridlist.lonlats ;
S_out.varNames = cropList_lpjg ;

S_out.garr_xv = sdate_xc ;
out_file = sprintf('%s/cropcals_sdate.%s.txt', thisDir, remapVer) ;
lpjgu_matlab_saveTable(out_header_cell, S_out, out_file,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;

S_out.garr_xv = hdate_xc ;
out_file = sprintf('%s/cropcals_hdate.%s.txt', thisDir, remapVer) ;
lpjgu_matlab_saveTable(out_header_cell, S_out, out_file,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;

S_out.garr_xv = gslen_xc ;
out_file = sprintf('%s/cropcals_gslen.%s.txt', thisDir, remapVer) ;
lpjgu_matlab_saveTable(out_header_cell, S_out, out_file,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;

S_out.garr_xv = fertDate2_xc ;
out_file = sprintf('%s/cropcals_fertDate2.%s.txt', thisDir, remapVer) ;
lpjgu_matlab_saveTable(out_header_cell, S_out, out_file,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'gzip', do_gzip) ;

disp('All done!')


%% Save source maps

thisPos = figurePos ;
spacing = [0.05 0.05] ; % v h
ylim = 65:360 ;
fontSize = 14 ;

cropList_in = unique(cropList_as_gsType) ;
Ncrops_in = length(cropList_in) ;

varList = {'sdate', 'hdate', 'gslen'} ;
for v = 1:length(varList)
    thisVar = varList{v} ;
    eval(['source_xc = source_' thisVar '_xc ;']) ;
    
    if Ncrops_in > 6
        error('Set ny and nx for %d crops', Ncrops)
    else
        ny = 2 ;
        nx = 3 ;
    end
    
    figure('Position', figurePos, 'Color', 'w') ;
    for c = 1:Ncrops_in
        thisCrop = cropList_in{c} ;
        subplot_tight(ny, nx, c, spacing) ;
        cc = find(strcmp(cropList_as_gsType, thisCrop), 1) ;
        map_YX = lpjgu_vector2map(source_xc(:,cc), ...
            size(gridlist.mask_YX), gridlist.list2map) ;
        pcolor(map_YX(ylim,:)) ;
        shading flat; axis equal tight off
        caxis([1 3])
        hcb = colorbar('Location', 'SouthOutside') ;
        hcb.Ticks = 1:3 ;
        title(hcb, 'Source phase')
        title(strrep( ...
            sprintf('%s: %s', thisVar, thisCrop), ...
            '_', '\_')) ;
        set(gca, ...
            'FontSize', fontSize) ;
    end
    
    outFile = sprintf('%s/sourcePhase_%s.png', ...
        thisDir, thisVar) ;
    export_fig(outFile, '-r150')
    close
    
    
    
end


%% FUNCTIONS

function gslen_sanity_check(sdate_x, hdate_x, gslen_x)

tmp = hdate_x ;
tmp(hdate_x < sdate_x) = tmp(hdate_x < sdate_x) + 365 ;
gslen_from_shdates_x = tmp - sdate_x ;

diff_x = abs(gslen_from_shdates_x - gslen_x) ;
isbad = abs(diff_x) > 0 ;
if any(isbad)
    error('Inconsistency in growing season length read from file vs. that calculated from sowing and harvest dates')
end

end

function data_x = handle_non_integers(data_x)
%%% From https://ebi-forecast.igb.illinois.edu/ggcmi/projects/ggcmi/wiki/Phase_2_known_issues
% there are some non-integer values in growing seasons and 2nd fertilizer
% date. If your processing does not handle this automatically, just use the
% floor() function that truncates any decimal digits. #334 
% As there are cases where the date is 0.5, floor() does not work, so you 
% either need to check for that in the processing or you need to use ceil()
% but also check for values >365. (#336)
%
% I've decided to floor(), then add 1 to any zero values.

data_x = floor(data_x) ;
data_x(data_x ==0) = 1 ;

end

function out_x = gapfill_with_phase3( ...
    thisCrop, whichVar, in_x, ...
    inDir_phase3, gridlist)

% is_wheat = contains(thisCrop, 'wheat') ;
% if is_wheat
%     thisCrop = sprintf('%swh', thisCrop(1)) ;
% end

switch whichVar
    case 'sdate'
        fileVar = 'planting_day' ;
    case 'hdate'
        fileVar = 'maturity_day' ;
    case 'gslen'
        fileVar = 'growing_season_length' ;
    otherwise
        error('whichVar %s not recognized', whichVar)
end

switch thisCrop
    case 'spring_wheat'
        fileCrop = 'swh' ;
    case 'winter_wheat'
        fileCrop = 'wwh' ;
    case 'Rice'
        fileCrop = 'ri1' ;
    case 'Soybeans'
        fileCrop = 'soy' ;
    case 'Maize'
        fileCrop = 'mai' ;
    otherwise
        error('thisCrop %s not recognized', thisCrop)
end

thisIrr = 'rf' ;
filename = sprintf('%s/%s_%s_ggcmi_crop_calendar_phase3_v1.01.nc4', ...
    inDir_phase3, fileCrop, thisIrr) ;
if ~exist(filename, 'file')
    error('File not found: %s', filename) ;
end
map_YX = flipud(transpose(ncread(filename, fileVar))) ;
map_YX(map_YX<0) = NaN ;

if any(~isint(map_YX) & ~isnan(map_YX))
    error('Non-integers in %s', whichVar)
end

out_x = in_x ;
out_x(isnan(out_x)) = map_YX(gridlist.list2map(isnan(out_x))) ;

end



function out_x = import_phase12( ...
    thisCrop, whichVar, ...
    inDir_phase1, inDir_phase2, gridlist)

is_wheat = contains(thisCrop, 'wheat') ;
if is_wheat
    thisCrop = sprintf('%swh', thisCrop(1)) ;
end

switch whichVar
    case 'sdate'
        fileVar = 'planting day' ;
    case 'hdate'
        fileVar = 'harvest day' ;
    case 'gslen'
        fileVar = 'growing season length' ;
    otherwise
        error('whichVar %s not recognized', whichVar)
end

thisIrr = 'rf' ;
if is_wheat
    filename = sprintf('%s/%s_%s_growing_season_dates_v2.nc4', ...
        inDir_phase2, thisCrop, thisIrr) ;
else
    filename = sprintf('%s/%s_%s_growing_season_dates_v1.25.nc4', ...
        inDir_phase1, thisCrop, thisIrr) ;
end
if ~exist(filename, 'file')
    error('File not found: %s', filename) ;
end
map_YX = transpose(ncread(filename, fileVar)) ;
map_YX(map_YX<0) = NaN ;

if any(~isint(map_YX) & ~isnan(map_YX))
    error('Non-integers in %s', whichVar)
end

out_x = map_YX(gridlist.list2map) ;




end
