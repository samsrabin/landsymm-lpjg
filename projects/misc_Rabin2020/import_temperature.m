%%%%%%%%%%%%%%%%%%%%%%%
%%% Get temperature %%%
%%%%%%%%%%%%%%%%%%%%%%%

disp('Importing temperature...')

years_bc = 1961:1990 ;
Nyears_bc = length(years_bc) ;
Nmonths_bc = Nyears_bc * 12 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% "For bias correction" part in CRU-NCEP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correctionmethod_bl = get_correctionmethod(baselineDir) ;
if ~strcmp(correctionmethod_bl, 'c3')
    error('correctionmethod (%s) not recognized', correctionmethod_bl)
end
[cruncep_temp_bc_XYm, cruncep_temp_bc_yearList] = get_cmip5_temp( ...
    baselineDir, 'file_cru', '*tair*monmean*.nc4', 'Temperature', years_bc) ;
if ~isequal(shiftdim(cruncep_temp_bc_yearList), shiftdim(years_bc))
    error('Years missing from CRU-NCEP temp for bias correction!')
elseif ~isequal([size(land_area_YX') Nmonths_bc], size(cruncep_temp_bc_XYm))
    error('Problem with shape of cruncep_temp_bc_XYm!')
end
clear cruncep_temp_bc_yearList


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% "For bias correction" part in "historical CMIP5" %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cmip5_temp_bc_XYm, cmip5_temp_bc_yearList] = get_cmip5_temp( ...
    baselineDir, 'path_cmip5hist', 'tas*halfDeg.nc4', 'tas', years_bc) ;
if ~isequal(shiftdim(cmip5_temp_bc_yearList), shiftdim(years_bc))
    error('Years missing from CMIP5 temp for bias correction!')
elseif ~isequal([size(land_area_YX') Nmonths_bc], size(cmip5_temp_bc_XYm))
    error('Problem with shape of cmip5_temp_bc_XYm!')
end
clear cmip5_temp_bc_yearList


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% "End-historical" part in "historical CMIP5" %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cmip5_temp_blRAW_XYm, cmip5_temp_bl_yearList] = ...
    get_cmip5_temp(baselineDir, 'path_cmip5hist', 'tas*halfDeg.nc4', 'tas', years_endh) ;
if ~isequal(shiftdim(cmip5_temp_bl_yearList), shiftdim(years_endh))
    [cmip5_temp2_blRAW_XYm, cmip5_temp2_bl_yearList] = ...
        get_cmip5_temp(baselineDir, 'path_cmip5scen', 'tas*halfDeg.nc4', 'tas', years_endh) ;
    cmip5_temp_blRAW_XYm = cat(3, cmip5_temp_blRAW_XYm, cmip5_temp2_blRAW_XYm) ;
    cmip5_temp_bl_yearList = cat(1, cmip5_temp_bl_yearList, cmip5_temp2_bl_yearList) ;
    if ~isequal(shiftdim(cmip5_temp_bl_yearList), shiftdim(years_endh))
        error('Can''t find enough years!')
    elseif ~isequal([size(land_area_YX') 12*length(years_endh)], size(cmip5_temp_blRAW_XYm))
        error('Problem with shape of cmip5_temp_bl_XYm!')
    end
    clear cmip5_temp2_blRAW_XYm cmip5_temp2_bl_yearList
end
clear cmip5_temp_bl_yearList


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% "End-future" part in "RCP CMIP5" %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmip5_temp_fuRAW_XYmr = nan([size(cmip5_temp_blRAW_XYm) Nruns]) ;
for r = 1:Nruns
    correctionmethod_thisRun = get_correctionmethod(runDirs{r}) ;
    if ~strcmp(correctionmethod_bl, correctionmethod_thisRun)
        error('Correction method differs between baseline and thisRun!')
    end
    [tmp_XYm, tmp_yearList] = ...
        get_cmip5_temp(runDirs{r}, 'path_cmip5scen', 'tas*halfDeg.nc4', 'tas', years_endf) ;
    if ~isequal(shiftdim(tmp_yearList), shiftdim(years_endf))
        error('Can''t find enough years!')
    elseif ~isequal([size(land_area_YX') 12*length(years_endf)], size(tmp_XYm))
        error('Problem with shape of tmp_XYm!')
    end
    cmip5_temp_fuRAW_XYmr(:,:,:,r) = tmp_XYm ;
    clear tmp* correctionmethod_thisRun
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform bias correction %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(correctionmethod_bl, 'c3')
    % Get mean climatologies
    cruncep_temp_bc_XYmy = reshape(cruncep_temp_bc_XYm, [size(land_area_YX') 12 Nyears_bc]) ;
    cruncep_temp_bc_XYmMean = mean(cruncep_temp_bc_XYmy,4) ;
    clear cruncep_temp_bc_XYmy
    cmip5_temp_bc_XYmy = reshape(cmip5_temp_bc_XYm, [size(land_area_YX') 12 Nyears_bc]) ;
    cmip5_temp_bc_XYmMean = mean(cmip5_temp_bc_XYmy,4) ;
    clear cmip5_temp_bc_XYmy

    % Correct
    correction_XYmMean = cmip5_temp_bc_XYmMean - cruncep_temp_bc_XYmMean ;
    clear cmip5_temp_bc_XYmMean cruncep_temp_bc_XYmMean
    cmip5_temp_blRAW_XYmy = reshape(cmip5_temp_blRAW_XYm, [size(land_area_YX') 12 length(years_endh)]) ;
    temp_bl_XYmy = cmip5_temp_blRAW_XYmy - repmat(correction_XYmMean,[1 1 1 length(years_endh)]) ;
    clear cmip5_temp_blRAW_XYmy
    cmip5_temp_fuRAW_XYmyr = reshape(cmip5_temp_fuRAW_XYmr, [size(land_area_YX') 12 length(years_endf) Nruns]) ;
    clear cmip5_temp_fuRAW_XYmr
    temp_fu_XYmyr = cmip5_temp_fuRAW_XYmyr - repmat(correction_XYmMean,[1 1 1 length(years_endf) Nruns]) ;
    clear cmip5_temp_fuRAW_XYmyr correction_XYmMean
else
    error('correctionmethod (%s) not recognized', correctionmethod_bl)
end
clear correctionmethod_bl


%%%%%%%%%%%%%%%%%
%%% Finish up %%%
%%%%%%%%%%%%%%%%%

% Reshape
temp_bl_YXmy = flip(permute(temp_bl_XYmy, [2 1 3 4]),1) ;
clear temp_bl_XYmy
temp_fu_YXmyr = flip(permute(temp_fu_XYmyr, [2 1 3 4 5]),1) ;
clear temp_fu_XYmyr

% Annual means
monLengths = [31 28 31 30 31 30 31 31 30 31 30 31] ; % Assumes noleap
monWeights = monLengths / sum(monLengths) ;
if abs(sum(monWeights) - 1) > 1e-9
    error('Problem with month weights!')
end
temp_bl_YXy = squeeze(sum(temp_bl_YXmy.*repmat(permute(monWeights,[1 3 2]), [size(land_area_YX) 1 length(years_endh)]), 3)) ;
temp_fu_YXyr = squeeze(sum(temp_fu_YXmyr.*repmat(permute(monWeights,[1 3 2]), [size(land_area_YX) 1 length(years_endf)]), 3)) ;
clear monLengths monWeights

% Area-weighted time series
land_area_weights_YX = land_area_YX ./ nansum(nansum(land_area_YX)) ;
if abs(nansum(nansum(land_area_weights_YX)) - 1) > 1e-9
    error('Problem with area weights!')
end
land_area_unmasked_weights_YX = land_area_unmasked_YX ./ nansum(nansum(land_area_unmasked_YX)) ;
if abs(nansum(nansum(land_area_unmasked_weights_YX)) - 1) > 1e-9
    error('Problem with area weights!')
end
ts_temp_bl = nan(Nyears_bl,1) ;
[~,IA] = intersect(yearList_baseline,years_endh) ;
ts_temp_bl(IA) = squeeze(nansum(nansum(temp_bl_YXy.*repmat(land_area_weights_YX,[1 1 length(years_endh)]),1),2)) ;
ts_temp_yr = nan(Nyears_fu,Nruns) ;
[~,IA] = intersect(yearList_future,years_endf) ;
ts_temp_yr(IA,:) = squeeze(nansum(nansum(temp_fu_YXyr.*repmat(land_area_weights_YX,[1 1 length(years_endf) Nruns]),1),2)) ;
