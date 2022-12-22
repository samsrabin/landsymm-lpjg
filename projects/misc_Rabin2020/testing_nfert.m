is_baseline = true ;

%% Import

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Import others
if is_baseline
    LU = lpjgu_matlab_readTable_then2map('/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt') ;
    cropfrac = lpjgu_matlab_readTable_then2map('/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.Misc0s.txt') ;
    nfert = lpjgu_matlab_readTable_then2map('/project/fh1-project-lpjgpi/lr8247/PLUM/input/Nfert/nfert_1700_2015_luh2_aggregate_sum2x2_midpoint_rescaled_v20.txt') ;
    if isfield(nfert,'maps_YXvy')
        whichIs = size(nfert.maps_YXvy,4) ;
    else
        whichIs = 1 ;
    end
else
    LU = lpjgu_matlab_readTable_then2map('/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v3.s1.forLPJG.MATLAB.20180426/landcover.txt') ;
    cropfrac = lpjgu_matlab_readTable_then2map('/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v3.s1.forLPJG.MATLAB.20180426/cropfractions.txt') ;
    nfert = lpjgu_matlab_readTable_then2map('/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v3.s1.forLPJG.MATLAB.20180426/nfert.txt') ;
    whichIs = 1 ;
end

if exist('nfert_orig','var')
    nfert = nfert_orig ;
end
if is_baseline
    nfert_orig = nfert ;
    nfert.varNames = cropfrac.varNames ;
    nfertIndices = nan(length(cropfrac.varNames),1) ;
    for c = 1:length(cropfrac.varNames)
        thisCrop = cropfrac.varNames{c} ;
        switch thisCrop
            case {'CerealsC3','Oilcrops','StarchyRoots','Pulses','Rice'}; thisnfertType = 'CC3ann' ;
            case {'CerealsC4','Miscanthus'}                             ; thisnfertType = 'CC4ann' ;
            case {'CerealsC3i','Oilcropsi','StarchyRootsi','Pulsesi','Ricei'}; thisnfertType = 'CC3ann' ;
            case {'CerealsC4i','Miscanthusi'}                             ; thisnfertType = 'CC4ann' ;
            otherwise; error(['Which nfertType should I use for ' thisCrop '?'])
        end
        thisnfertIndex = find(strcmp(nfert_orig.varNames,thisnfertType)) ;
        nfertIndices(c) = thisnfertIndex ;
    end
    nfert_maps_YXv = mean(nfert_orig.maps_YXvy(:,:,nfertIndices,whichIs),4) ;
    cropfrac_maps_YXv = cropfrac.maps_YXv ;
    croparea_YX = LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),end) ...
        .* land_area_YX_m2 ;
else
    nfert_maps_YXv = mean(nfert.maps_YXvy(:,:,:,whichIs),4) ;
    cropfrac_maps_YXv = mean(cropfrac.maps_YXvy(:,:,:,whichIs),4) ;
    croparea_YX = mean(LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),whichIs) ...
        .* repmat(land_area_YX_m2,[1 1 1 length(whichIs)]),4) ;
end


%%

cfts = sort(cropfrac.varNames) ;
Ncrops = length(cfts) ;

for c = 1:Ncrops
    thisCrop = cfts{c} ;
    thisSum = nansum(nansum(...
        croparea_YX ...
        .* cropfrac_maps_YXv(:,:,c) ...
        .* nfert_maps_YXv(:,:,c)...
        )) ;
    disp([thisCrop ': ' num2str(thisSum*1e-9) ' Mt'])
end