function [yieldAdj_rf_YXy, yieldAdj_ir_YXy] = import_and_interp(...
    thisCrop_ggcmi_long, thisCrop_calib, in_dir_top, ...
    nfert, Nlevels, calib_year_indices, croparea_ha_YX1y, thisModel, ...
    Nyears_cal, warn_missingVal, warn_yieldLowNgtHigh)

% Get short name of crop (for filename)
switch thisCrop_ggcmi_long
    case 'winter_wheat' ; thisCrop_ggcmi = 'wwh' ;
    case 'spring_wheat' ; thisCrop_ggcmi = 'swh' ;
    case 'maize' ; thisCrop_ggcmi = 'mai' ;
    case 'rice' ; thisCrop_ggcmi = 'ric' ;
    case 'soy' ; thisCrop_ggcmi = 'soy' ;
    otherwise ; error([thisCrop_ggcmi_long ' not recognized!'])
end

% Get input directory
in_dir =  addslashifneeded([in_dir_top thisCrop_ggcmi_long '/A0/yield']) ;

Nnlevels = length(Nlevels) ;
yield_rf_YXyn = nan(360,720,Nyears_cal,Nnlevels) ;
yield_ir_YXyn = nan(360,720,Nyears_cal,Nnlevels) ;
for n = 1:Nnlevels
    thisN = num2str(Nlevels(n)) ;
    
    %%%%%%%%%%%%%%%
    %%% Rainfed %%%
    %%%%%%%%%%%%%%%
    
    % Import
    in_file = [in_dir lower(thisModel) '_agmerra_fullharm_yield_' thisCrop_ggcmi '_global_annual_1980_2010_C360_T0_W0_N' thisN '_A0.nc4'] ;
    tmp = ncread(in_file,['yield_' thisCrop_ggcmi]) ;
    
    % Check NaN mask
    nan_mask = isnan(flipud(permute(tmp(:,:,calib_year_indices),[2 1 3]))) ;
    nan_mask = sum(nan_mask,3) ;
    if ~isequal(unique(nan_mask),[0 ; Nyears_cal])
        shademap(nan_mask);
        title(['# NaN over ' num2str(Nyears_cal) ' years'])
        warning(['Inconsistent NaN mask among years within ' thisCrop_ggcmi_long '_rf for N level ' num2str(Nlevels(n)) '.'])
    end
    if n==1
        nan_mask_n1 = nan_mask ;
    elseif ~isequal(nan_mask,nan_mask_n1)
        shademap(nan_mask - nan_mask_n1);
        title(['nanmask - nanmask1, ' thisCrop_ggcmi_long 'rf'])
        warning(['Inconsistent NaN mask among N levels within ' thisCrop_ggcmi_long '_rf (' num2str(Nlevels(1)) ' vs. ' num2str(Nlevels(n)) ').'])
    end
    
    % Save
    yield_rf_YXyn(:,:,:,n) = flipud(permute(tmp(:,:,calib_year_indices),[2 1 3])) ;
    clear tmp
    
    %%%%%%%%%%%%%%%%%
    %%% Irrigated %%%
    %%%%%%%%%%%%%%%%%
    
    % Import
    in_file = [in_dir lower(thisModel) '_agmerra_fullharm_yield_' thisCrop_ggcmi '_global_annual_1980_2010_C360_T0_Winf_N' thisN '_A0.nc4'] ;
    tmp = ncread(in_file,['yield_' thisCrop_ggcmi]) ;
    
    % Check NaN mask
    nan_mask = isnan(flipud(permute(tmp(:,:,calib_year_indices),[2 1 3]))) ;
    nan_mask = sum(nan_mask,3) ;
    if ~isequal(unique(nan_mask),[0 ; Nyears_cal])
%         shademap(nan_mask);
%         title(['# NaN over ' num2str(Nyears_cal) ' years'])
        warning(['Inconsistent NaN mask among years within ' thisCrop_ggcmi_long '_ir for N level ' num2str(Nlevels(n)) '.'])
    end
    if n==1
        nan_mask_n1 = nan_mask ;
    elseif ~isequal(nan_mask,nan_mask_n1)
%         shademap(nan_mask - nan_mask_n1);
%         title(['nanmask - nanmask1, ' thisCrop_ggcmi_long 'ir'])
        warning(['Inconsistent NaN mask among N levels within ' thisCrop_ggcmi_long '_ir (' num2str(Nlevels(1)) ' vs. ' num2str(Nlevels(n)) ').'])
    end
    
    % Save
    yield_ir_YXyn(:,:,:,n) = flipud(permute(tmp(:,:,calib_year_indices),[2 1 3])) ;
    clear tmp
end

% Interpolate to get yields under historical N fertilization
if strcmp(thisCrop_calib,'CerealsC3w') || strcmp(thisCrop_calib,'CerealsC3s')
    nfert_toUse_YX = nfert.maps_YXv(:,:,strcmp(nfert.varNames,'CerealsC3')) ;
else
    nfert_toUse_YX = nfert.maps_YXv(:,:,strcmp(nfert.varNames,thisCrop_calib)) ;
end
[yieldAdj_rf_YXy, yieldAdj_ir_YXy] = interp_yield_f(...
    nfert_toUse_YX, true, ...
    squeeze(croparea_ha_YX1y), Nlevels, ...
    yield_rf_YXyn, yield_ir_YXyn, ...
    warn_missingVal, warn_yieldLowNgtHigh) ;


end