%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare two versions of historical LU data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % lu_file_1 = '/Users/Shared/PLUM/input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
% % cf_file_1 = '/Users/Shared/PLUM/input/remaps_v4/cropfracs.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt' ;
% % nf_file_1 = '/Users/Shared/PLUM/input/remaps_v4/nfert.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt' ;
% 
% % lu_file_1 = '/Users/Shared/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
% % cf_file_1 = '/Users/Shared/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.txt' ;
% % nf_file_1 = '/Users/Shared/PLUM/input/remaps_v2/nfert.remapv2.20180214.m0.txt' ;
% 
% % lu_file_1 = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % cf_file_1 = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % nf_file_1 = '/Users/Shared/PLUM/input/remaps_v5/nfert.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% 
% % lu_file_1 = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v11.s1.forLPJG.MATLAB.20181031/landcover.txt' ;
% % cf_file_1 = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v11.s1.forLPJG.MATLAB.20181031/cropfractions.txt' ;
% % nf_file_1 = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v11.s1.forLPJG.MATLAB.20181031/nfert.txt' ;
% 
% 
% % lu_file_2 = '/Users/Shared/PLUM/input/remaps_v3/LU_xtraCROPtoPAST.remapv3.20180214.cgFertIrr0.m0.txt' ;
% % cf_file_2 = '/Users/Shared/PLUM/input/remaps_v3/cropfracs.remapv3.20180214.cgFertIrr0.m0.txt' ;
% % nf_file_2 = '/Users/Shared/PLUM/input/remaps_v3/nfert.remapv3.20180214.cgFertIrr0.m0.txt' ;
% 
% % lu_file_2 = '/Users/Shared/PLUM/input/remaps_v4/LU.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt' ;
% % cf_file_2 = '/Users/Shared/PLUM/input/remaps_v4/cropfracs.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt' ;
% % nf_file_2 = '/Users/Shared/PLUM/input/remaps_v4/nfert.remapv4.20180214.cgFertIrr0.setaside0103.m0.txt' ;
% 
% % lu_file_2 = '/Users/Shared/PLUM/input/remaps_v5/LU.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % cf_file_2 = '/Users/Shared/PLUM/input/remaps_v5/cropfracs.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % nf_file_2 = '/Users/Shared/PLUM/input/remaps_v5/nfert.remapv5.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % 
% % lu_file_2 = '/Users/Shared/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % cf_file_2 = '/Users/Shared/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% % nf_file_2 = '/Users/Shared/PLUM/input/remaps_v6/nfert.remapv6.20180214.ecFertIrr0.setaside0103.m4.txt' ;
% 
% % lu_file_2 = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v11.s1.harm.forLPJG/landcover.txt' ;
% % cf_file_2 = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v11.s1.harm.forLPJG/cropfractions.txt' ;
% % nf_file_2 = '/Users/Shared/PLUM/PLUM_outputs_for_LPJG/SSP1.v11.s1.harm.forLPJG/nfert.txt' ;
% 
% % outDir = '/Users/Shared/PLUM/input/remaps_v2versus3/' ;
% % outDir = '/Users/Shared/PLUM/input/remaps_v2versus4/' ;
% % outDir = '/Users/Shared/PLUM/input/remaps_v2versus5/' ;
% % outDir = '/Users/Shared/PLUM/input/remaps_v2versus6/' ;
% % outDir = '/Users/Shared/PLUM/input/remaps_v5versus6/' ;
% % outDir = '/Users/Shared/PLUM/input/remaps_v5/' ;
% % outDir = '/Users/Shared/PLUM/input/remaps_v6/' ;
% % outDir = '/Users/Shared/PLUM/input/v11_ssp1_origVersusHarm/' ;

lu_file_1 = '/Users/Shared/PLUM/input/remaps_v6p3/LU.remapv6p3.txt' ;
cf_file_1 = '/Users/Shared/PLUM/input/remaps_v6p3/cropfracs.remapv6p3.txt' ;
nf_file_1 = '/Users/Shared/PLUM/input/remaps_v6p3/nfert.remapv6p3.txt' ;
lu_file_2 = '/Users/Shared/PLUM/input/remaps_v6p6/LU.remapv6p6.txt' ;
cf_file_2 = '/Users/Shared/PLUM/input/remaps_v6p6/cropfracs.remapv6p6.txt' ;
nf_file_2 = '/Users/Shared/PLUM/input/remaps_v6p6/nfert.remapv6p6.txt' ;
outDir = '/Users/Shared/PLUM/input/remaps_v6p3_versus_v6p6/' ;
yearList = 1850:2015 ;


%% Setup

if ~exist(outDir,'dir')
    mkdir(outDir) ;
end

Nyears = length(yearList) ;
LUlist = {'CROPLAND','NATURAL','PASTURE','BARREN'} ;
Nlu = length(LUlist) ;
cropList = {'CerealsC3','CerealsC4','Rice','Oilcrops','Pulses','StarchyRoots','ExtraCrop'} ;
Ncrops = length(cropList) ;

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

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

% v1
disp('Importing v1...')
[lu1, lu1_area_YXvy, ...
    cf1, cf1_areaT_YXvy, cf1_areaI_YXvy, ...
    nf1, nfTot1_YXvy] = ...
    compare_hist_LUs_import( ...
    lu_file_1, cf_file_1, nf_file_1, ...
    LUlist, yearList, landArea_YX, cropList) ;
disp('Done')

% v2
disp('Importing v2...')
[lu2, lu2_area_YXvy, ...
    cf2, cf2_areaT_YXvy, cf2_areaI_YXvy, ...
    nf2, nfTot2_YXvy] = ...
    compare_hist_LUs_import( ...
    lu_file_2, cf_file_2, nf_file_2, ...
    LUlist, yearList, landArea_YX, cropList) ;
disp('Done.')



%% Plot change in land use area

spacing = [0.1 0.1] ;

figure('Color','w','Position',figurePos) ;
for v = 1:Nlu
    subplot_tight(2,2,v,spacing) ;
    tmp1 = squeeze(nansum(nansum(lu1_area_YXvy(:,:,v,:),2),1)) ;
    tmp2 = squeeze(nansum(nansum(lu2_area_YXvy(:,:,v,:),2),1)) ;
    plot(yearList,1e-6*tmp1,'LineWidth',3) ;
    hold on
    plot(yearList,1e-6*tmp2,'--','LineWidth',3) ;
    hold off
    title(LUlist{v}) ;
    xlabel('Year')
    ylabel('Million km^2')
    set(gca,'FontSize',14)
    legend({'Orig','New'},'Location','NorthWest')
end
sgt = sgtitle('Land use areas') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_lu.pdf'],'-r300')
close


%% Plot change in total crop areas

spacing = [0.1 0.1] ;
thisPos = [1    33   720   772] ;

figure('Color','w','Position',thisPos) ;
for c = 1:(Ncrops+1)
    subplot_tight(4,2,c,spacing) ;
    if c==8
        thisCrop = 'Total' ;
        tmp1 = squeeze(nansum(nansum(nansum(cf1_areaT_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(cf2_areaT_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(cf1_areaT_YXvy(:,:,c,:),2),1)) ;
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
close


%% Plot change in irrigated crop areas

spacing = [0.1 0.1] ;
thisPos = [1    33   720   772] ;

figure('Color','w','Position',thisPos) ;
for c = 1:(Ncrops+1)
    subplot_tight(4,2,c,spacing) ;
    if c==8
        thisCrop = 'Total' ;
        tmp1 = squeeze(nansum(nansum(nansum(cf1_areaI_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(cf2_areaI_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(cf1_areaI_YXvy(:,:,c,:),2),1)) ;
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
        tmp1 = squeeze(nansum(nansum(nansum(nfTot1_YXvy,2),1),3)) ./ squeeze(nansum(nansum(nansum(cf1_areaT_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(nfTot2_YXvy,2),1),3)) ./ squeeze(nansum(nansum(nansum(cf2_areaT_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(nfTot1_YXvy(:,:,c,:),2),1)) ./ squeeze(nansum(nansum(cf1_areaT_YXvy(:,:,c,:),2),1)) ;
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
        tmp1 = squeeze(nansum(nansum(nansum(nfTot1_YXvy,2),1),3)) ;
        tmp2 = squeeze(nansum(nansum(nansum(nfTot2_YXvy,2),1),3)) ;
    else
        thisCrop = cropList{c} ;
        tmp1 = squeeze(nansum(nansum(nfTot1_YXvy(:,:,c,:),2),1)) ;
        tmp2 = squeeze(nansum(nansum(nfTot2_YXvy(:,:,c,:),2),1)) ;
    end
    plot(yearList,1e-9*tmp1,'LineWidth',3) ;
    hold on
    plot(yearList,1e-9*tmp2,'--','LineWidth',3) ;
    hold off
    title(thisCrop) ;
    xlabel('Year')
    ylabel('TgN')
    set(gca,'FontSize',14)
    legend({'Orig','New'},'Location','NorthWest')
end
sgt = sgtitle('N fertilization') ;
sgt.FontSize = 20 ; sgt.FontWeight = 'bold' ;

export_fig([outDir 'ts_nfertTot.pdf'],'-r300')
close
