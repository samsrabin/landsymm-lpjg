LUfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
NfertFile = '/Users/Shared/PLUM/input/remaps_v6/nfert.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;

%%

NfluxFile = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6_testNoFertOnFake/output-2018-11-28-182230/nflux.out.gz' ;
% NfluxFile = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-11-05-085615/nflux.out.gz' ;


%% Import land area (km2 to m2)

disp('Importing land area...')

% Read PLUMout_gridlist
gridlist_file = '/Users/Shared/lpj-guess/gridlists/PLUMout_gridlist.txt' ;
PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;

landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
gcel_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
% Convert to m2
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

disp('Done.')

isNH = false(360,720) ;
isNH(181:360,:) = true ;


%% Import others

disp('Importing other files...')

lu = lpjgu_matlab_readTable_then2map(LUfile) ;
lu.maps_YXvy = lu.maps_YXvy(:,:,:,end) ;

cf = lpjgu_matlab_readTable_then2map(cropfile) ;
cf.maps_YXvy = cf.maps_YXvy(:,:,:,end) ;

nfert = lpjgu_matlab_readTable_then2map(NfertFile) ;
nfert.maps_YXvy = nfert.maps_YXvy(:,:,:,end) ;

nflux = lpjgu_matlab_readTable_then2map(NfluxFile) ;
nflux.maps_YXvy = nflux.maps_YXvy(:,:,:,end) ;


disp('Done.')


%%

if ~isequal(cf.varNames,nfert.varNames)
    error('~isequal(cf.varNames,nfert.varNames)')
end

tmpG = land_area_YX ...
    .* lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:) ...
    .* sum(cf.maps_YXvy .* nfert.maps_YXvy,3) ;
fprintf('Good: %0.3e\n',nansum(nansum(tmpG)))

tmpB = land_area_YX .* (-1e-4*nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'fert'),:)) ;
fprintf(' Bad: %0.3e\n',nansum(nansum(tmpB)))

tmpB = land_area_YX .* (-1e-4*nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'fert'),:)) ...
- land_area_YX ...
    .* lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:) ...
    .* cf.maps_YXvy(:,:,strcmp(cf.varNames,'ExtraCrop'),:) .* nfert.maps_YXvy(:,:,strcmp(cf.varNames,'CerealsC3'),:) ;
fprintf('Bad2: %0.3e\n',nansum(nansum(tmpB)))


%%

mask_YX = isNH & sum(cf.maps_YXvy(:,:,contains(cf.varNames,'CerealsC3'),:),3)>1e-5 & lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:)>5e-6 ;
shademap(~mask_YX & lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:)>5e-6) ;

tmpG = land_area_YX ...
    .* lu.maps_YXvy(:,:,strcmp(lu.varNames,'CROPLAND'),:) ...
    .* sum(cf.maps_YXvy .* nfert.maps_YXvy,3) ;
tmpG(mask_YX) = NaN ;
fprintf('Good: %0.3e\n',nansum(nansum(tmpG)))

tmpB = land_area_YX .* (-1e-4*nflux.maps_YXvy(:,:,strcmp(nflux.varNames,'fert'),:)) ;
tmpB(mask_YX) = NaN ;
fprintf(' Bad: %0.3e\n',nansum(nansum(tmpB)))


%%

for c = 1:size(nflux.maps_YXvy,3)
    thisField = nflux.varNames{c} ;
    tmp = land_area_YX .* (-1e-4*nflux.maps_YXvy(:,:,c,1)) ;
    fprintf(' %s:\t %0.3e\n',thisField,nansum(nansum(tmp)))
end


