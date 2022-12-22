%% Import data

yield = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/trunk_runs/PLUM2LPJG_SSP1_RCP45_v3s1/output-2018-04-23-145614/yield.out') ;
yield.maps_YXvy(:,:,contains(yield.varNames,'_ic'),:) = [] ;
yield.varNames(contains(yield.varNames,'_ic')) = [] ;
lu = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v3.s1.forLPJG.MATLAB.20180416/landcover.txt') ;
cropfracs = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v3.s1.forLPJG.MATLAB.20180416/cropfractions.txt') ;
cropfracs.maps_YXvy(:,:,contains(cropfracs.varNames,'Miscanthus'),:) = [] ;
cropfracs.varNames(contains(cropfracs.varNames,'Miscanthus')) = [] ;

lu_old = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt') ;
cropfracs_old = lpjgu_matlab_readTable_then2map('/Users/Shared/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.NfertEmpties0-0200-1000.txt') ;
cropfracs_old.maps_YXvy(:,:,contains(cropfracs_old.varNames,'Miscanthus'),:) = [] ;
cropfracs_old.varNames(contains(cropfracs_old.varNames,'Miscanthus')) = [] ;
cropfracs_old.maps_YXv(:,:,contains(cropfracs_old.varNames,'0')) = [] ;
cropfracs_old.varNames(contains(cropfracs_old.varNames,'0')) = [] ;

% Get irrigation only
isirrfun = @(x) strcmp(x(end),'i') ;
yield.maps_YXvy(:,:,~cellfun(isirrfun,yield.varNames),:) = [] ;
yield.varNames(~cellfun(isirrfun,yield.varNames)) = [] ;
cropfracs.maps_YXvy(:,:,~cellfun(isirrfun,cropfracs.varNames),:) = [] ;
cropfracs.varNames(~cellfun(isirrfun,cropfracs.varNames)) = [] ;

% Get same order
[yield.varNames,I] = sort(yield.varNames) ;
yield.maps_YXvy = yield.maps_YXvy(:,:,I,:) ;
[cropfracs.varNames,I] = sort(cropfracs.varNames) ;
cropfracs.maps_YXvy = cropfracs.maps_YXvy(:,:,I,:) ;
if ~isequal(yield.varNames,cropfracs.varNames)
    error('~isequal(yield.varNames,cropfracs.varNames)')
end

% Import land area (km2)
nanmask = isnan(lu.maps_YXvy(:,:,1,1)) ;
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX(nanmask) = NaN ;
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd


%% Calculate production

area = cropfracs.maps_YXvy .* ...
    repmat(lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:),[1 1 length(yield.varNames) 1]) .* ...
    repmat(land_area_YX_m2,[1 1 size(yield.maps_YXvy,3) size(yield.maps_YXvy,4)]) ;

prod_YXvy = yield.maps_YXvy .* area ;
    
    



%%

yield1_YXc = yield.maps_YXvy(:,:,:,1) ;
yield2_YXc = yield.maps_YXvy(:,:,:,3) ;

area1_YXc = area(:,:,:,1) ;
area2_YXc = area(:,:,:,3) ;

%%

figure('Color','w','Position',figurePos) ;

for c = 1:6
    thisCrop = yield.varNames{c} ;
    subplot(3,2,c)
    pcolor(yield2_YXc(:,:,c)-yield1_YXc(:,:,c)) ; shading flat ; axis equal tight
    title(thisCrop)
    caxis([-max(abs(caxis)) max(abs(caxis))])
    colorbar; 
end

%%

figure('Color','w','Position',figurePos) ;

for c = 1:6
    thisCrop = yield.varNames{c} ;
    subplot(3,2,c)
    tmp = double(yield2_YXc(:,:,c)>0 & yield1_YXc(:,:,c)==0 & area2_YXc(:,:,c)>0 & area1_YXc(:,:,c)>0) ;
    tmp(nanmask) = NaN ;
    pcolor(tmp) ; shading flat ; axis equal tight
    title(thisCrop)
    caxis([-max(abs(caxis)) max(abs(caxis))])
    colorbar; 
end


%%

figure('Color','w','Position',figurePos) ;

plot(yield.yearList(1:5:end),squeeze(nansum(nansum(nansum(prod_YXvy(:,:,:,1:5:end),3),2),1)))
hold on
plot(yield.yearList(3:5:end),squeeze(nansum(nansum(nansum(prod_YXvy(:,:,:,3:5:end),3),2),1)))
hold off

