%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make figures of LU data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lu_file = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.txt' ;
% cf_file = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.txt' ;
% nf_file = '/Users/Shared/PLUM/input/remaps_v6p3/nfert.remapv6p3.txt' ;
% outDir = '/Users/Shared/PLUM/input/remaps_v6p3/' ;
% yearList = 1850:2010 ;

lu_file = '/Users/Shared/PLUM/input/remaps_v6p6/LU.remapv6p6.txt' ;
cf_file = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.txt' ;
nf_file = '/Users/Shared/PLUM/input/remaps_v6p6/nfert.remapv6p6.txt' ;
outDir = '/Users/Shared/PLUM/input/remaps_v6p6/' ;
yearList = 1850:2010 ;

% lu_file = '/Users/Shared/PLUM/input/remaps_v7a/LU.remapv7a.txt' ;
% cf_file = '/Users/Shared/PLUM/input/remaps_v7a/cropfracs.remapv7a.txt' ;
% nf_file = '/Users/Shared/PLUM/input/remaps_v7a/nfert.remapv7a.txt' ;
% outDir = '/Users/Shared/PLUM/input/remaps_v7a/' ;
% yearList = 1850:2010 ;


%% Setup

Nyears = length(yearList) ;
LUlist = {'CROPLAND','NATURAL','PASTURE','BARREN'} ;
Nlu = length(LUlist) ;

addpath(genpath(landsymm_lpjg_path()))

landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;

% Get LUH2 land area (km2)
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd tmp

% Get mask and apply to land area
file_in = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v10.s1/2010/LandCoverFract.txt' ;
S = lpjgu_matlab_readTable_then2map(file_in,'verboseIfNoMat',false,'force_mat_nosave',true) ;
mask_YX = isnan(S.maps_YXv(:,:,1)) ...
    | sum(S.maps_YXv(:,:,contains(S.varNames,{'CROPLAND','PASTURE','NATURAL'})),3)==0 ...
    | landArea_YX==0 ;
clear S
landArea_YX(mask_YX) = 0 ;


%% Import data 

% Import land uses
lu = lpjgu_matlab_readTable_then2map(lu_file) ;
lu = trim_years(lu, yearList) ;

% Trim land uses
if ~isequal(lu.varNames,LUlist)
    [~,IA] = setdiff(lu.varNames,LUlist) ;
    if any(lu.maps_YXvy(:,:,IA,:)>0)
        error('LUs excluded from lu_out have area >0!')
    end
    clear IA
    if ~isempty(setdiff(LUlist,lu.varNames))
        error('lu_out missing something from LUlist!')
    end
    [~,~,IB] = intersect(LUlist,lu.varNames,'stable') ;
    lu.maps_YXvy = lu.maps_YXvy(:,:,IB,:) ;
    clear IB
    lu.varNames = LUlist ;
end

% Convert to land use areas
lu.maps_YXvy = lu.maps_YXvy .* repmat(landArea_YX,[1 1 Nlu Nyears]) ;

% Import crop fractions
cf = import_croppy(cf_file, yearList) ;
% cropList = cf.varNames ;
cropList = {} ;
for c = 1:length(cf.varNames)
    thisCrop = cf.varNames{c} ;
    if ~any(strcmp(cf.varNames,thisCrop(1:end-1)))
        cropList = [cropList thisCrop] ;
    end
end
Ncrops = length(cropList) ;

% Combine crop fractions, if needed
if sum(contains(cf.varNames,{'CC3G','CC4G'}))==2
    tmpFrac_3 = cf.maps_YXvy(:,:,strcmp(cf.varNames,'CC3G'),:) ;
    tmpFrac_4 = cf.maps_YXvy(:,:,strcmp(cf.varNames,'CC4G'),:) ;
    tmpWeight_3 = tmpFrac_3 ./ (tmpFrac_3 + tmpFrac_4) ;
    tmpWeight_4 = tmpFrac_4 ./ (tmpFrac_3 + tmpFrac_4) ;
    tmpWeight_3(tmpFrac_3+tmpFrac_4 == 0) = 0 ;
    tmpWeight_4(tmpFrac_3+tmpFrac_4 == 0) = 0 ;
    cf.maps_YXvy(:,:,strcmp(cf.varNames,'CC3G'),:) = sum(cf.maps_YXvy(:,:,contains(cf.varNames,{'CC3G','CC4G'}),:),3) ;
    cf.varNames{strcmp(cf.varNames,'CC3G')} = 'ExtraCrop' ;
    cf.maps_YXvy(:,:,strcmp(cf.varNames,'CC4G'),:) = [] ;
    cf.varNames(strcmp(cf.varNames,'CC4G')) = [] ;
elseif any(contains(cf.varNames,{'CC3G','CC4G'}))
	error('???')
end

% Convert to crop areas
cf.maps_YXvy = cf.maps_YXvy .* repmat(landArea_YX,[1 1 length(cf.varNames) Nyears]) ;

% Import nfert
nf = import_croppy(nf_file, yearList) ;
if sum(contains(nf.varNames,{'CC3G','CC4G'}))==2
    if any(any(any(nf.maps_YXvy(:,:,contains(nf.varNames,{'CC3G','CC4G'}),:)>0)))
        warning('Fertilizer applied to CC3G and/or CC4G.')
        tmpNfert_3 = nf.maps_YXvy(:,:,contains(nf.varNames,{'CC3G'}),:) ;
        tmpNfert_4 = nf.maps_YXvy(:,:,contains(nf.varNames,{'CC4G'}),:) ;
        nf.maps_YXvy(:,:,contains(nf.varNames,{'CC3G'}),:) = ...
            tmpWeight_3.*tmpNfert_3 + tmpWeight_4.*tmpNfert_4 ;
        clear tmp*
    end
    nf.varNames{strcmp(nf.varNames,'CC3G')} = 'ExtraCrop' ;
    nf.maps_YXvy(:,:,strcmp(nf.varNames,'CC4G'),:) = [] ;
    nf.varNames(strcmp(nf.varNames,'CC4G')) = [] ;
elseif any(contains(nf.varNames,{'CC3G','CC4G'}))
	error('???')
end

% Get Nfert totals (kg)
nf.maps_YXvy = nf.maps_YXvy .* (1e6*cf.maps_YXvy) ;


%% Plot change in land use area

spacing = [0.1 0.1] ;

figure('Color','w','Position',figurePos) ;
for v = 1:Nlu
    subplot_tight(2,2,v,spacing) ;
    tmp1 = squeeze(nansum(nansum(lu.maps_YXvy(:,:,v,:),2),1)) ;
    plot(yearList,1e-6*tmp1,'LineWidth',3) ;
    title(LUlist{v}) ;
    xlabel('Year')
    ylabel('Million km^2')
    set(gca,'FontSize',14)
end
sgt = sgtitle('Land use areas') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_lu.pdf'],'-r300')


%% Plot change in total crop areas

spacing = [0.1 0.1] ;
thisPos = [1    33   720   772] ;

figure('Color','w','Position',thisPos) ;
for c = 1:(Ncrops+1)
    subplot_tight(4,2,c,spacing) ;
    if c==8
        thisCrop = 'Total' ;
        tmp1 = squeeze(nansum(nansum(nansum(cf_areaT_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(cf2_areaT_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(cf_areaT_YXvy(:,:,c,:),2),1)) ;
        tmp2 = squeeze(nansum(nansum(cf2_areaT_YXvy(:,:,c,:),2),1)) ;
    end
    plot(yearList,1e-6*tmp1,'LineWidth',3) ;
    hold on
    plot(yearList,1e-6*tmp2,'--','LineWidth',3) ;
    hold off
    title(thisCrop) ;
    xlabel('Year')
    ylabel('Million km^2')
    set(gca,'FontSize',14)
    legend({'Orig','New'},'Location','NorthWest')
end
sgt = sgtitle('Crop area') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_cropareas.pdf'],'-r300')


%% Plot change in irrigated crop areas

spacing = [0.1 0.1] ;
thisPos = [1    33   720   772] ;

figure('Color','w','Position',thisPos) ;
for c = 1:(Ncrops+1)
    subplot_tight(4,2,c,spacing) ;
    if c==8
        thisCrop = 'Total' ;
        tmp1 = squeeze(nansum(nansum(nansum(cf_areaI_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(cf2_areaI_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(cf_areaI_YXvy(:,:,c,:),2),1)) ;
        tmp2 = squeeze(nansum(nansum(cf2_areaI_YXvy(:,:,c,:),2),1)) ;
    end
    plot(yearList,1e-6*tmp1,'LineWidth',3) ;
    hold on
    plot(yearList,1e-6*tmp2,'--','LineWidth',3) ;
    hold off
    title(thisCrop) ;
    xlabel('Year')
    ylabel('Million km^2')
    set(gca,'FontSize',14)
    legend({'Orig','New'},'Location','NorthWest')
end

sgt = sgtitle('Irrigated crop area') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_cropareasIrr.pdf'],'-r300')


%% Plot change in fertilization RATE

spacing = [0.1 0.1] ;
thisPos = [1    33   720   772] ;

figure('Color','w','Position',thisPos) ;
for c = 1:(Ncrops+1)
    subplot_tight(4,2,c,spacing) ;
    if c==8
        thisCrop = 'Total' ;
        tmp1 = squeeze(nansum(nansum(nansum(nfTot_YXvy,2),1),3)) ./ squeeze(nansum(nansum(nansum(cf_areaT_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(nfTot2_YXvy,2),1),3)) ./ squeeze(nansum(nansum(nansum(cf2_areaT_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(nfTot_YXvy(:,:,c,:),2),1)) ./ squeeze(nansum(nansum(cf_areaT_YXvy(:,:,c,:),2),1)) ;
        tmp2 = squeeze(nansum(nansum(nfTot2_YXvy(:,:,c,:),2),1)) ./ squeeze(nansum(nansum(cf2_areaT_YXvy(:,:,c,:),2),1));
    end
    plot(yearList,1e-5*tmp1,'LineWidth',3) ;
    hold on
    plot(yearList,1e-5*tmp2,'--','LineWidth',3) ;
    hold off
    title(thisCrop) ;
    xlabel('Year')
    ylabel('tons ha^{-1}')
    set(gca,'FontSize',14)
    legend({'Orig','New'},'Location','NorthWest')
end
sgt = sgtitle('N fertilization rate') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_nfert.pdf'],'-r300')


%% Plot change in fertilization TOTAL

spacing = [0.1 0.1] ;
thisPos = [1    33   720   772] ;

figure('Color','w','Position',thisPos) ;
for c = 1:(Ncrops+1)
    subplot_tight(4,2,c,spacing) ;
    if c==8
        thisCrop = 'Total' ;
        tmp1 = squeeze(nansum(nansum(nansum(nf.maps_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(nf.maps_YXvy(:,:,c,:),2),1)) ;
    end
    plot(yearList,1e-9*tmp1,'LineWidth',3) ;
    title(thisCrop) ;
    xlabel('Year')
    ylabel('TgN')
    set(gca,'FontSize',14)
end
sgt = sgtitle('N fertilization') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_nfertTot.pdf'],'-r300')
close
