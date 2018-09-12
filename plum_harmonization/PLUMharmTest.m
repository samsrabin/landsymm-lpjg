%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing harmonized PLUM land use trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisVer = '' ;
% thisVer = 'orig.' ;
% thisVer = '2deg.' ;

runList = {'SSP1.v10.s1' ;
%             'SSP3.v10.s1' ;
%             'SSP4.v10.s1' ;
%             'SSP5.v10.s1';
            } ;

base_year = 2010 ;

yearList = 2011:2020 ;

norm2extra = 0.177 ;

% Method for inpaint_nans()
inpaint_method = 0 ;

out_dir = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/harmonization_figs_20180831/' ;


%% Setup

yearList_luh2 = 1971:2010 ;

if length(runList) == 1
    legend_ts = {'LUH2','Orig','Harm'} ;
else
    legend_ts = {'LUH2'} ;
    for s = 1:length(runList)
        legend_ts = [legend_ts {runList{s}(1:4)}] ;
    end
end

topDir = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;
PLUM_in_toptop = strcat(topDir,runList) ;
PLUM_base_in = [addslashifneeded(PLUM_in_toptop{1}) '2010/'] ;

Nyears = length(yearList) ;
Nruns = length(runList) ;

addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;

% % Make lower-left lat/lon map (for compat. with PLUM style)
% lons_map_2deg = repmat(-180:2:178,[90 1]) ;
% lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
% lons_map = repmat(-180:0.5:179.5,[360 1]) ;
% lats_map = repmat((-90:0.5:89.5)',[1 720]) ;

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;


%% Import reference data

doHarm = false ;
run('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/plum_harmonization/PLUMharm_importRefData.m') ;


%% Import PLUM (original + harmonized)

is2deg = strcmp(thisVer,'2deg.') ;

if is2deg
    ny = 90 ;
    nx = 180 ;
else
    ny = 360 ;
    nx = 720 ;
end

disp('Setting up PLUM*_YXvyr arrays...')
PLUMorig_YXvyr = nan(ny,nx,Nlu,Nyears,Nruns,'single') ;
PLUMorig_nfert_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears,Nruns,'single') ;
PLUMorig_irrig_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears,Nruns,'single') ;
PLUMharm_YXvyr = nan(ny,nx,Nlu,Nyears,Nruns,'single') ;
PLUMharm_nfert_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears,Nruns,'single') ;
PLUMharm_irrig_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears,Nruns,'single') ;

for r = 1:Nruns
    thisRun = removeslashifneeded(runList{r}) ;

    % Original
    fprintf('Importing %s...\n', thisRun) ;
    [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
        [topDir thisRun], base_year, yearList, ...
        landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, base_bareFrac_YX, norm2extra, inpaint_method, '') ;
    PLUMorig_YXvyr(:,:,:,:,r) = S_out.maps_YXvy ;
    PLUMorig_nfert_YXvyr(:,:,:,:,r) = S_nfert_out.maps_YXvy ;
    PLUMorig_irrig_YXvyr(:,:,:,:,r) = S_irrig_out.maps_YXvy ;
%     PLUMorig_YXvyr(:,:,:,:,r) = S_out.maps_YXvy(:,:,:,2011:2100<=max(yearList)) ;
%     PLUMorig_nfert_YXvyr(:,:,:,:,r) = S_nfert_out.maps_YXvy(:,:,:,2011:2100<=max(yearList)) ;
%     PLUMorig_irrig_YXvyr(:,:,:,:,r) = S_irrig_out.maps_YXvy(:,:,:,2011:2100<=max(yearList)) ;
    clear S_*out

    % Harmonized
    fprintf('Importing %s.harm...\n', thisRun) ;
    [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
        [topDir thisRun '.harm'],base_year,yearList, ...
        landArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, [], 0, [], thisVer) ;
    PLUMharm_YXvyr(:,:,:,:,r) = S_out.maps_YXvy ;
    PLUMharm_nfert_YXvyr(:,:,:,:,r) = S_nfert_out.maps_YXvy ;
    PLUMharm_irrig_YXvyr(:,:,:,:,r) = S_irrig_out.maps_YXvy ;
    clear S_*out

end

disp('Done reading PLUM.')


%% Time series of LUs

ts_base_cy = squeeze(nansum(nansum(base.maps_YXvy,1),2)) ;
ts_base_cy = cat(1,sum(ts_base_cy(isCrop,:),1),ts_base_cy(~isCrop,:)) ;
ts_orig_cyr = squeeze(nansum(nansum(PLUMorig_YXvyr,1),2)) ;
ts_orig_cyr = cat(1,sum(ts_orig_cyr(isCrop,:,:),1),ts_orig_cyr(~isCrop,:,:)) ;
ts_harm_cyr = squeeze(nansum(nansum(PLUMharm_YXvyr,1),2)) ;
ts_harm_cyr = cat(1,sum(ts_harm_cyr(isCrop,:,:),1),ts_harm_cyr(~isCrop,:,:)) ;
combinedLUs = [{'CROPLAND'} LUnames(~isCrop)] ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:length(combinedLUs)
    subplot_tight(2,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:)*1e-6*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList,squeeze(ts_orig_cyr(v,:,:))*1e-6*1e-6,'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList,squeeze(ts_harm_cyr(v,:,:))*1e-6*1e-6,'-','LineWidth',1)
    hold off
    title(combinedLUs{v})
    set(gca,'FontSize',14)
    ylabel('Million km2')
    legend(legend_ts,'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_landUse.pdf']) ;
close


%% Time series of crops

ts_base_cy = squeeze(nansum(nansum(base.maps_YXvy,1),2)) ;
ts_orig_cyr = squeeze(nansum(nansum(PLUMorig_YXvyr,1),2)) ;
ts_harm_cyr = squeeze(nansum(nansum(PLUMharm_YXvyr,1),2)) ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:Ncrops_lpjg
    subplot_tight(4,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:)*1e-6*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList,squeeze(ts_orig_cyr(v,:,:))*1e-6*1e-6,'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList,squeeze(ts_harm_cyr(v,:,:))*1e-6*1e-6,'-','LineWidth',1)
    hold off
    title(LPJGcrops{v})
    set(gca,'FontSize',14)
    ylabel('Million km2')
    legend(legend_ts,'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_crops.pdf']) ;
close


%% Time series of Nfert

if is2deg
    ts_base_cy = cf_kg2Mt .* squeeze(nansum(nansum(base_nfertTot_2deg.maps_YXvy(:,:,isCrop,:)))) ;
else
    ts_base_cy = cf_kg2Mt .* squeeze(nansum(nansum(base_nfertTot.maps_YXvy(:,:,isCrop,:)))) ;
end
ts_orig_cyr = cf_kg2Mt .* squeeze(nansum(nansum(PLUMorig_YXvyr(:,:,isCrop,:,:) .* PLUMorig_nfert_YXvyr,1),2)) ;
ts_harm_cyr = cf_kg2Mt .* squeeze(nansum(nansum(PLUMharm_YXvyr(:,:,isCrop,:,:) .* PLUMharm_nfert_YXvyr,1),2)) ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:Ncrops_lpjg
    subplot_tight(4,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:),'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList,squeeze(ts_orig_cyr(v,:,:)),'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList,squeeze(ts_harm_cyr(v,:,:)),'-','LineWidth',1)
    hold off
    title(LPJGcrops{v})
    set(gca,'FontSize',14)
    ylabel('Mt')
    legend(legend_ts,'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_nfert.pdf']) ;
close


%% Maps: At three years

theseYears = [2011 2050 2100] ;
spacing = [0.01 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 60 ;
fontSize = 14 ;
png_res = 150 ;

for r = 1:Nruns
    thisRun = runList{r} ;
    for v = 1:Nlu
        figure('Color','w','Position',figurePos) ;
        thisLU = LUnames{v} ;
        for y = 1:3
            thisYear = theseYears(y) ;
            h1 = subplot_tight(2,3,y,spacing) ;
            tmp = PLUMorig_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,3,y+3,spacing) ;
            tmp = PLUMharm_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s harm: %s, %d (km^2)',thisRun,thisLU,thisYear)) ;
            new_caxis = [0 max([caxis(h1) caxis(h2)])] ;
            caxis(h1,new_caxis) ;
            caxis(h2,new_caxis) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'maps_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end


%% Maps: Diffs between orig and harm at 3 years

theseYears = [2011 2050 2100] ;
spacing = [0.01 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 60 ;
fontSize = 14 ;
png_res = 150 ;
thisPos = [1         500        1440         305] ;

for r = 1:Nruns
    thisRun = runList{r} ;
    for v = 1:Nlu
        figure('Color','w','Position',thisPos) ;
        thisLU = LUnames{v} ;
        for y = 1:3
            thisYear = theseYears(y) ;
            h1 = subplot_tight(1,3,y,spacing) ;
            tmp1 = PLUMorig_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp2 = PLUMharm_YXvyr(:,:,v,yearList==thisYear,r) ;
            tmp = tmp2 - tmp1 ;
%             tmp = tmp2/sum(tmp2(:)) - tmp1/sum(tmp1(:)) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(flipud(brewermap(64,'rdbu_ssr'))) ;
            caxis([-1 1]*max(abs(caxis))) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'mapsOHdiffs_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end




%% Maps: Differences between two pairs of years

theseYears = [2011 2050 2100] ;
spacing = [0.01 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 60 ;
fontSize = 14 ;
png_res = 150 ;

for r = 1:Nruns
    thisRun = runList{r} ;
    for v = 1:Nlu
        figure('Color','w','Position',figurePos) ;
        thisLU = LUnames{v} ;
        for y = 1:2
            thisYear1 = theseYears(y) ;
            thisYear2 = theseYears(y+1) ;
            h1 = subplot_tight(2,2,y,spacing) ;
            tmp1 = PLUMorig_YXvyr(:,:,v,yearList==thisYear1,r) ;
            tmp2 = PLUMorig_YXvyr(:,:,v,yearList==thisYear2,r) ;
            tmp = tmp2 - tmp1 ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(brewermap(64,'rdbu_ssr')) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d-%d',thisRun,thisLU,thisYear1,thisYear2)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,2,y+2,spacing) ;
            tmp1 = PLUMharm_YXvyr(:,:,v,yearList==thisYear1,r) ;
            tmp2 = PLUMharm_YXvyr(:,:,v,yearList==thisYear2,r) ;
            tmp = tmp2 - tmp1 ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(brewermap(64,'rdbu_ssr')) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s harm: %s, %d-%d',thisRun,thisLU,thisYear1,thisYear2)) ;
            new_caxis = [-1 1] * max([caxis(h1) caxis(h2)]) ;
            caxis(h1,new_caxis) ;
            caxis(h2,new_caxis) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'mapsChgs_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end



%%

v = 3 ;
y1 = 1 ;
yN = 3 ;

orig1_YX = PLUMorig_YXvyr(:,:,v,y1,1) .* PLUMorig_nfert_YXvyr(:,:,v,y1,1) ;
origN_YX = PLUMorig_YXvyr(:,:,v,yN,1) .* PLUMorig_nfert_YXvyr(:,:,v,yN,1) ;
orig_YX = origN_YX - orig1_YX ;
harm1_YX = PLUMharm_YXvyr(:,:,v,y1,1) .* PLUMharm_nfert_YXvyr(:,:,v,y1,1) ;
harmN_YX = PLUMharm_YXvyr(:,:,v,yN,1) .* PLUMharm_nfert_YXvyr(:,:,v,yN,1) ;
harm_YX = harmN_YX - harm1_YX ;

figure ;
spacing = 0.025 ;
h1 = subplot_tight(2,1,1,spacing) ;
pcolor(orig_YX); shading flat; axis equal tight
title('orig')
colorbar
h2 = subplot_tight(2,1,2,spacing) ;
pcolor(harm_YX); shading flat; axis equal tight
title('harm')
colorbar
new_caxis = [-1 1] * max([abs(caxis(h1)) abs(caxis(h2))]) ;
caxis(h1,new_caxis) ;
caxis(h2,new_caxis) ;
