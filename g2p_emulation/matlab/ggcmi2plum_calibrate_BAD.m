%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Import and process GGCMI output for crop calibration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definitely bad
%%% thisModel = 'CARAIB' ; % Did not do N-factorial runs
%%% thisModel = 'EPIC-IIASA' ; % Did not do N-factorial runs
%%% thisModel = 'GEPIC' ; % Did not complete N-factorial runs with Winf.
%%% thisModel = 'JULES' ; % Did not do N-factorial runs
%%% thisModel = 'ORCHIDEE-crop' ; % For T0, only did W10 runs. No spring wheat.
%%% thisModel = 'PEPIC' ; % Does not appear to have finished WxN factorial experiments for T0.
%%% thisModel = 'PROMET' ; % Does not appear to have finished N factorial experiments.
%%% thisModel = 'PRYSBI2' ; % Did not do N-factorial runs
%%% thisModel = 'SIMPLACE' ; % Does not appear to have done any runs.

% Maybe bad
% thisModel = 'APSIM-UGOE' ; % Rice: Rainfed: It looks like maybe there is a missing value (-0.1) in 3/1011197 cells.
%                            % It's not actually a missing value, but
%                            % recommendation from Marian Koch is to assume
%                            % zero, because this behavior only happens in
%                            % inhospitable cells.
%                            % Maize (at least): Inconsistent NaN mask among N levels

% Apparently good
% thisModel = 'EPIC-TAMU' ; % Maize (at least): Inconsistent NaN mask among N levels
% thisModel = 'LPJ-GUESS' ; % Spring wheat (at least): Inconsistent NaN mask among N levels
thisModel = 'LPJmL' ;
% thisModel = 'pDSSAT' ;

warn_missingVal = true ;
warn_yieldLowNgtHigh = false ;


%% Get crops

in_dir_top = addslashifneeded(['/Volumes/WDMPP_Storage/GGCMI/AgMIP.output/' thisModel '/phase2/']) ;

tmp = dir(in_dir_top) ;
listCrops_ggcmi_long = {tmp.name} ;
listCrops_ggcmi_long(contains(listCrops_ggcmi_long,'.')) = [] ;

Ncrops = length(listCrops_ggcmi_long) ;

% For compatibility with crop_calibration.m etc.
listCrops_lpj_comb = cell(Ncrops,1) ;
for c = 1:Ncrops
    thisCrop_ggcmi_long = listCrops_ggcmi_long{c} ;
    switch thisCrop_ggcmi_long
        case 'winter_wheat' ; listCrops_lpj_comb{c} = 'CerealsC3w' ;
        case 'spring_wheat' ; listCrops_lpj_comb{c} = 'CerealsC3s' ;
        case 'maize' ; listCrops_lpj_comb{c} = 'CerealsC4' ;
        case 'rice' ; listCrops_lpj_comb{c} = 'Rice' ;
        case 'soy' ; listCrops_lpj_comb{c} = 'Oilcrops' ;
        otherwise ; error([thisCrop_ggcmi_long ' not recognized!'])
    end
end


%% Setup

% filename_landuse = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
% filename_cropfrac = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.txt' ;
% filename_nfert = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/nfert.remapv2.20180214.m0.txt' ;

filename_landuse = '/Users/Shared/PLUM/input/remaps_v5e/LU.remapv5e.txt' ;
filename_cropfrac = '/Users/Shared/PLUM/input/remaps_v5e/cropfracs.remapv5e.txt' ;
filename_nfert = '/Users/Shared/PLUM/input/remaps_v5e/nfert.remapv5e.txt' ;

cd '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/matlab'
addpath(genpath('/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/matlab'))

Nlevels = [10 60 200] ;
Nnlevels = length(Nlevels) ;

yearList_sim = 1980:2010 ;
Nyears_sim = length(yearList_sim) ;

yearList_cal = 1995:2005 ;
Nyears_cal = length(yearList_cal) ;
[~,calib_year_indices] = intersect(yearList_sim,yearList_cal) ;

% For (compatibility with) crop calibration code
% version_name = 'isimip2b_lpjg' ;
version_name = 'ggcmi2_lpjg' ;  
calib_ver = 16 ;
year1 = min(yearList_cal) ;
yearN = max(yearList_cal) ;
dir_code = '/Users/Shared/PLUM/crop_calib_code/' ;
dir_data = '/Users/Shared/PLUM/crop_calib_data/' ;
dir_outfigs = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/calibration_figs/' ;
dir_outmaps = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/Ninterpd_yield_maps/' ;
dir_outtables = '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/calibration_tables/' ;
addpath(genpath(dir_code))
filename_guess_landuse = filename_landuse ;
filename_guess_cropfrac = filename_cropfrac ;
filename_countriesMap = 'country_boundaries62892.noNeg99.extrapd.asc' ;


%% Import ancillary files, historical LU and Nfert

% Import land area (km2 to ha)
if ~exist('land_area_YX','var')
    landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
    gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
    land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
    land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
    %%%% Convert to half-degree
    tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
    land_area_km2_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
    % Convert to ha
    land_area_YX = land_area_km2_YX*1e2 ;
    clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd land_area_km2_YX
end
% Land uses
if ~exist('croparea_ha_YX1y','var')
    landuse_lpj = lpjgu_matlab_readTable_then2map(filename_landuse,'force_mat_save',true) ;
    [~,IA] = intersect(landuse_lpj.yearList,yearList_cal) ;
    landuse_lpj.maps_YXvy = landuse_lpj.maps_YXvy(:,:,:,IA) ;
    landuse_lpj.yearList = landuse_lpj.yearList(IA) ;
    croparea_ha_YX1y = landuse_lpj.maps_YXvy(:,:,strcmp(landuse_lpj.varNames,'CROPLAND'),:) ...
        .* repmat(land_area_YX,[1 1 1 Nyears_cal]) ;
    clear IA
end

% Crop fractions
if ~exist('cropfrac_lpj','var')
    cropfrac_lpj = lpjgu_matlab_readTable_then2map(filename_cropfrac,'force_mat_save',true) ;
end

% Fertilizer
if ~exist('nfert','var')
    nfert = lpjgu_matlab_readTable_then2map(filename_nfert,'force_mat_save',true) ;
    % Convert from kgN/m2 to kgN/ha
    nfert.maps_YXv = nfert.maps_YXv*1e4 ;
end

% Crop mask
min_crop_mask_YX = transpose(ncread('/Users/Shared/GGCMI/inputs/other.inputs/phase2.masks/boolean_cropmask_ggcmi_phase2.nc4',...
    'minimum cropland mask')) ;

disp('Done.')


%% Import GGCMI run yields

% Get this model's list of files
thisModel_files = dir([in_dir_top '*/A0/yield/*.nc4']);
thisModel_dates = sort([thisModel_files.datenum]) ;
thisModel_lastModDate_str = datestr(thisModel_dates(end),'yyyymmddHHMMSS') ;

% Set up yield_lpj structure
tmp = [listCrops_lpj_comb;strcat(listCrops_lpj_comb,'i')] ;
yield_lpj.varNames = tmp(:) ;
clear tmp
yield_lpj.yearList = yearList_cal ;
yield_lpj.maps_YXvy = nan(360,720,Ncrops*2,Nyears_cal) ;

for c = 1:Ncrops
    
    thisCrop_ggcmi_long = listCrops_ggcmi_long{c} ;
    disp(thisCrop_ggcmi_long)
    
    thisCrop_calib = listCrops_lpj_comb{c} ;
    
    [yieldAdj_rf_YXy, yieldAdj_ir_YXy] = import_and_interp(...
        thisCrop_ggcmi_long, thisCrop_calib, in_dir_top, ...
        nfert, Nlevels, calib_year_indices, croparea_ha_YX1y, thisModel, ...
        Nyears_cal, warn_missingVal, warn_yieldLowNgtHigh) ;

    % Save
%     iR = (c-1)*2 + 1 ;
%     iI = iR + 1 ;
    iR = find(strcmp(yield_lpj.varNames,thisCrop_calib)) ;
    iI = find(strcmp(yield_lpj.varNames,[thisCrop_calib 'i'])) ;
    disp(yield_lpj.varNames{iR})
    disp(yield_lpj.varNames{iI})
    yield_lpj.maps_YXvy(:,:,iR,:) = permute(yieldAdj_rf_YXy,[1 2 4 3]) ;
    yield_lpj.maps_YXvy(:,:,iI,:) = permute(yieldAdj_ir_YXy,[1 2 4 3]) ;
end
% clear croparea_ha_YX1y

if isempty(find(~isnan(yield_lpj.maps_YXvy),1))
    error('isempty(find(~isnan(yield_lpj.maps_YXvy),1))')
end

% shademap(mean(yield_lpj.maps_YXvy(:,:,iR,:),4));
disp('Done.')


%% Account for winter vs spring wheat

if any(contains(listCrops_lpj_comb,{'CerealsC3w','CerealsC3s'}))
    if any(contains(listCrops_lpj_comb,'CerealsC3w')) && any(contains(listCrops_lpj_comb,'CerealsC3s'))
        % Add new
        listCrops_lpj_comb = [listCrops_lpj_comb;'CerealsC3'] ;
        yield_lpj.varNames = [yield_lpj.varNames;'CerealsC3';'CerealsC3i'] ;
        yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames,'CerealsC3'),:) = ...
            max(cat(3,yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames,'CerealsC3w'),:),...
                      yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames,'CerealsC3s'),:)...
                   ),[],3) ;
        yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames,'CerealsC3i'),:) = ...
            max(cat(3,yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames,'CerealsC3wi'),:),...
                      yield_lpj.maps_YXvy(:,:,strcmp(yield_lpj.varNames,'CerealsC3si'),:)...
                   ),[],3) ;
        % Remove old
        listCrops_lpj_comb(contains(listCrops_lpj_comb,{'CerealsC3w','CerealsC3s'})) = [] ;
        to_remove = contains(yield_lpj.varNames,{'CerealsC3w','CerealsC3s'}) ;
        yield_lpj.maps_YXvy(:,:,to_remove,:) = [] ;
        yield_lpj.varNames(to_remove) = [] ;
        clear to_remove
    end
end


%% Add crops that were not simulated

listCrops_lpj_comb_all = {'CerealsC3','CerealsC4','Rice','Oilcrops','StarchyRoots','Pulses'} ;
listCrops_corresp_ggcmi = cell(size(listCrops_lpj_comb_all)) ;
for c = 1:length(listCrops_lpj_comb_all)
    thisCrop_calib = listCrops_lpj_comb_all{c} ;
    if any(contains(listCrops_lpj_comb,thisCrop_calib))
        switch thisCrop_calib
            case 'CerealsC3' ; listCrops_corresp_ggcmi{c} = 'max_wheat' ;
            case 'Oilcrops'  ; listCrops_corresp_ggcmi{c} = 'soy' ;
            case 'Pulses'  ; listCrops_corresp_ggcmi{c} = 'soy' ;
            case 'CerealsC4' ; listCrops_corresp_ggcmi{c} = 'maize' ;
            case 'Rice'      ; listCrops_corresp_ggcmi{c} = 'rice' ;
            otherwise ; error([thisCrop_ggcmi_long ' not recognized!'])
        end
        continue
    end
    
    % Match this crop with one that WAS simulated
    if contains(thisCrop_calib,{'Rice','Oilcrops','StarchyRoots','Pulses'})
        if contains(thisCrop_calib,{'Oilcrops','Pulses'}) ...
        && contains('soy',listCrops_ggcmi_long)
            thisCrop_ggcmi_long = 'soy' ;
        else
            thisCrop_ggcmi_long = 'spring_wheat' ;
        end
    else
        error([thisCrop_calib ' not matched with anything in listCrops_ggcmi_long!']) ;
    end
    listCrops_corresp_ggcmi{c} = thisCrop_ggcmi_long ;
    disp([thisCrop_calib ' (as ' thisCrop_ggcmi_long ')'])
    
    % Get N-interpolated yield values
    [yieldAdj_rf_YXy, yieldAdj_ir_YXy] = import_and_interp(...
        thisCrop_ggcmi_long, thisCrop_calib, in_dir_top, ...
        nfert, Nlevels, calib_year_indices, croparea_ha_YX1y, thisModel, ...
        Nyears_cal, warn_missingVal, warn_yieldLowNgtHigh) ;
    
    % Add new crops
    listCrops_lpj_comb = [listCrops_lpj_comb ; thisCrop_calib] ;
    yield_lpj.varNames = [yield_lpj.varNames ; thisCrop_calib ; [thisCrop_calib 'i']] ;
    yield_lpj.maps_YXvy = cat(3,yield_lpj.maps_YXvy,...
        permute(yieldAdj_rf_YXy,[1 2 4 3]), ...
        permute(yieldAdj_ir_YXy,[1 2 4 3])) ;
end


%% Do calibration

crop_calibration


%% Save results figure

disp('Exporting results figure...')

% Set size and position
set(gcf,'Position',figurePos)

% Get filename, appending last modified date/time of most recently modified
% file in this model's directory 
filename_outfig = [dir_outfigs thisModel '_' thisModel_lastModDate_str '_cv' num2str(calib_ver) '.v2.pdf'] ;

% Export and close
export_fig(filename_outfig) ;
close
disp('Done.')


%% Save results table

disp('Exporting results table...')

% Get filename, appending last modified date/time of most recently modified
% file in this model's directory 
filename_outtable = [dir_outtables thisModel '_' thisModel_lastModDate_str '_cv' num2str(calib_ver) '.csv'] ;

% Create table
table_out = table(shiftdim(listCrops_fa2o), ...
                  shiftdim(listCrops_corresp_ggcmi), ...
                  shiftdim(round(calib_factors_u,3))) ;
table_out.Properties.VariableNames = {'PLUM_crop','GGCMI_crop','calib_factor'} ;

% Export
writetable(table_out, filename_outtable) ;

disp('Done.')


%% Save interpolated yield maps
%%%%%%%%%%%%%%%
%%% Options
spacing = 0.01 ;
southern_limit = 35 ;
pngres = 300 ;
fontSize = 14 ;
position = [1 347 1440 458] ;
%%%%%%%%%%%%%%%

if length(listCrops_lpj_comb_all) ~= 6
    error(['Set up layout for ' num2str(length(listCrops_lpj_comb_all)) ' maps.'])
end

% Rainfed
figure('Color','w','Position',position) ;
for c = 1:length(listCrops_lpj_comb_all)
    thisCrop = listCrops_lpj_comb_all{c} ;
    thisInd = find(strcmp(yield_lpj.varNames,thisCrop)) ;
    if isempty(thisInd)
        error([thisCrop ' not found in yield_lpj.varNames.'])
    end
    subplot_tight(2,3,c,spacing) ;
    pcolor(mean(yield_lpj.maps_YXvy(southern_limit:end,:,thisInd,:),4)) ;
    shading flat ; axis equal tight off
    colorbar ;
    title([thisModel ': ' thisCrop ' (rf)'])
    set(gca,'FontSize',fontSize)
end
filename_outfig = [dir_outmaps thisModel '_rf_' thisModel_lastModDate_str '.png'] ;
export_fig(filename_outfig,['-r' num2str(pngres)]) ;
close


% Rainfed
figure('Color','w','Position',position) ;
for c = 1:length(listCrops_lpj_comb_all)
    thisCrop = listCrops_lpj_comb_all{c} ;
    thisInd = find(strcmp(yield_lpj.varNames,[thisCrop 'i'])) ;
    if isempty(thisInd)
        error([thisCrop ' not found in yield_lpj.varNames.'])
    end
    subplot_tight(2,3,c,spacing) ;
    pcolor(mean(yield_lpj.maps_YXvy(southern_limit:end,:,thisInd,:),4)) ;
    shading flat ; axis equal tight off
    colorbar ;
    title([thisModel ': ' thisCrop ' (ir)'])
    set(gca,'FontSize',fontSize)
end
filename_outfig = [dir_outmaps thisModel '_ir_' thisModel_lastModDate_str '.png'] ;
export_fig(filename_outfig,['-r' num2str(pngres)]) ;
close
    


%% TS: Original irrigated yield at each Nlevel

% figure('Color','w')
% clim2 = 0 ;
% for n = 1:Nnlevels
%     
%     hp{n} = subplot_tight(3,1,n,0.025) ;
%     
%     pcolor(mean(yield_ir_YXyn(:,:,:,n),3));
%     shading flat; axis equal tight off
%     colorbar
%     title(num2str(Nlevels(n)))
%     set(hp{n},'FontSize',14)
%     
%     clim2 = max(clim2,max(caxis)) ;
%     
% end
% for n = 1:Nnlevels
%     caxis(hp{n},[0 clim2]) ;
% end



