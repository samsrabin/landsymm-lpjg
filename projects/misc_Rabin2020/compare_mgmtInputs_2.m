%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare 20171228 vs. 20180109 PLUMouts %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisRun = 'SSP1' ;

CFTnames = {'TeWWi';
            'TeSWi'
            'TeCoi';
            'TrRii';
            } ;
Ncfts = length(CFTnames) ;


%% Setup

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/')

% Get input directories
inDir_v1 = addslashifneeded(['/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG/' thisRun '.v2.forLPJG.MATLAB.20171201']) ;
% inDir_v2 = addslashifneeded(['/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG/' thisRun '.v2.forLPJG.MATLAB.20180109.BAD']) ;
inDir_v2 = addslashifneeded(['/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG/' thisRun '.v2.forLPJG.MATLAB.20180109']) ;

% Get year info
yearList = 2011:2100 ;
Nyears = length(yearList) ;

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calibration_for_PLUM/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd


%% Import

landcover_v1 = lpjgu_matlab_readTable_then2map([inDir_v1 'landcover.txt'],'verbose',false,'force_mat_save',true) ;
cropfracs_v1 = lpjgu_matlab_readTable_then2map([inDir_v1 'cropfractions.txt'],'verbose',false,'force_mat_save',true) ;
irrig_v1 = lpjgu_matlab_readTable_then2map([inDir_v1 'irrig.txt'],'verbose',false,'force_mat_save',true) ;
nfert_v1 = lpjgu_matlab_readTable_then2map([inDir_v1 'nfert.txt'],'verbose',false,'force_mat_save',true) ;

landcover_v2 = lpjgu_matlab_readTable_then2map([inDir_v2 'landcover.txt'],'verbose',false,'force_mat_save',true) ;
cropfracs_v2 = lpjgu_matlab_readTable_then2map([inDir_v2 'cropfractions.txt'],'verbose',false,'force_mat_save',true) ;
irrig_v2 = lpjgu_matlab_readTable_then2map([inDir_v2 'irrig.txt'],'verbose',false,'force_mat_save',true) ;
nfert_v2 = lpjgu_matlab_readTable_then2map([inDir_v2 'nfert.txt'],'verbose',false,'force_mat_save',true) ;


%% Compare total cropland area

% Options %%%%%%%%
fontSize = 14 ;
%%%%%%%%%%%%%%%%%%

tmp_v1 = squeeze(nansum(nansum(landcover_v1.maps_YXvy(:,:,strcmp(landcover_v1.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
tmp_v2 = squeeze(nansum(nansum(landcover_v2.maps_YXvy(:,:,strcmp(landcover_v1.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;

figure('Color','w') ;
plot(yearList,tmp_v1,'--r','LineWidth',3)
hold on
plot(yearList,tmp_v2,'-b','LineWidth',3)
plot(yearList,tmp_v1,'--r','LineWidth',3)
hold off
legend('20171201','20180109','Location','NorthWest')
title(['Cropland area (' thisRun ')'])
ylabel('km^2')
set(gca,'FontSize',fontSize)


%% Compare area of each crop

% Options %%%%%%%%
fontSize = 14 ;
spacing = [0.07 0.05] ;
%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',figurePos) ;

for c = 1:Ncfts
    thisCrop = CFTnames{c} ;
    tmp_v1 = squeeze(nansum(nansum(cropfracs_v1.maps_YXvy(:,:,strcmp(cropfracs_v1.varNames,thisCrop),:) .* landcover_v1.maps_YXvy(:,:,strcmp(landcover_v1.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
    tmp_v2 = squeeze(nansum(nansum(cropfracs_v2.maps_YXvy(:,:,strcmp(cropfracs_v2.varNames,thisCrop),:) .* landcover_v2.maps_YXvy(:,:,strcmp(landcover_v2.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
    subplot_tight(2,2,c,spacing) ;
    plot(yearList,tmp_v1,'--r','LineWidth',3) ;
    hold on
    plot(yearList,tmp_v2,'-b','LineWidth',3) ;
    plot(yearList,tmp_v1,'--r','LineWidth',3) ;
    hold off
    legend('20171201','20180109','Location','NorthWest')
    title([thisCrop ' area (' thisRun ')'])
    ylabel('km^2')
    set(gca,'FontSize',fontSize)
end


%% Compare fertilization of each crop

% Options %%%%%%%%
fontSize = 14 ;
spacing = [0.07 0.05] ;
%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',figurePos) ;

for c = 1:Ncfts
    thisCrop = CFTnames{c} ;
    tmp_v1 = squeeze(nansum(nansum(1e6*nfert_v1.maps_YXvy(:,:,strcmp(nfert_v1.varNames,thisCrop),:) .* cropfracs_v1.maps_YXvy(:,:,strcmp(cropfracs_v1.varNames,thisCrop),:) .* landcover_v1.maps_YXvy(:,:,strcmp(landcover_v1.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
    tmp_v2 = squeeze(nansum(nansum(1e6*nfert_v2.maps_YXvy(:,:,strcmp(nfert_v2.varNames,thisCrop),:) .* cropfracs_v2.maps_YXvy(:,:,strcmp(cropfracs_v2.varNames,thisCrop),:) .* landcover_v2.maps_YXvy(:,:,strcmp(landcover_v2.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
    subplot_tight(2,2,c,spacing) ;
    plot(yearList,tmp_v1,'--r','LineWidth',3) ;
    hold on
    plot(yearList,tmp_v2,'-b','LineWidth',3) ;
    plot(yearList,tmp_v1,'--r','LineWidth',3) ;
    hold off
    legend('20171201','20180109','Location','NorthWest')
    title([thisCrop ' N fertilizer (' thisRun ')'])
    ylabel('kg')
    set(gca,'FontSize',fontSize)
end


%% Compare irrigation of each crop

% Options %%%%%%%%
fontSize = 14 ;
spacing = [0.07 0.05] ;
%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',figurePos) ;

for c = 1:Ncfts
    thisCrop = CFTnames{c} ;
    tmp_v1 = squeeze(nansum(nansum(irrig_v1.maps_YXvy(:,:,strcmp(nfert_v1.varNames,thisCrop),:) .* cropfracs_v1.maps_YXvy(:,:,strcmp(cropfracs_v1.varNames,thisCrop),:) .* landcover_v1.maps_YXvy(:,:,strcmp(landcover_v1.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
    tmp_v2 = squeeze(nansum(nansum(irrig_v2.maps_YXvy(:,:,strcmp(nfert_v2.varNames,thisCrop),:) .* cropfracs_v2.maps_YXvy(:,:,strcmp(cropfracs_v2.varNames,thisCrop),:) .* landcover_v2.maps_YXvy(:,:,strcmp(landcover_v2.varNames,'CROPLAND'),:) .* repmat(land_area_YX,[1 1 1 Nyears]),1),2)) ;
    subplot_tight(2,2,c,spacing) ;
    plot(yearList,tmp_v1,'--r','LineWidth',3) ;
    hold on
    plot(yearList,tmp_v2,'-b','LineWidth',3) ;
    plot(yearList,tmp_v1,'--r','LineWidth',3) ;
    hold off
    legend('20171201','20180109','Location','NorthWest')
    title([thisCrop ' irrigation (' thisRun ')'])
    ylabel('Ignore units')
    set(gca,'FontSize',fontSize)
end


