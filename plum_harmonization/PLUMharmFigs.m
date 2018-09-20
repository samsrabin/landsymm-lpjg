%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing harmonized PLUM land use trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

read_MATs = true ;

thisVer = '' ;
% thisVer = 'orig.' ;
% thisVer = '2deg.' ;

runList = {...
           'SSP1.v10.s1' ;
           'SSP3.v10.s1' ;
           'SSP4.v10.s1' ;
           'SSP5.v10.s1';
            } ;

base_year = 2010 ;

yearList_harm = 2011:2100 ;

norm2extra = 0.177 ;

% Method for inpaint_nans()
inpaint_method = 0 ;

out_dir = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/harmonization_figs_20180831/' ;


%% Setup

yearList_luh2 = 1971:2010 ;
yearList_orig = [yearList_harm(1)-1 yearList_harm] ;

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

Nyears_orig = length(yearList_orig) ;
Nyears_harm = length(yearList_harm) ;
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
    thisLandArea_YX = landArea_2deg_YX ;
else
    ny = 360 ;
    nx = 720 ;
    thisLandArea_YX = landArea_YX ;
end

disp('Setting up PLUM*_YXvyr arrays...')
PLUMorig_YXvyr = nan(ny,nx,Nlu,Nyears_orig,Nruns,'single') ;
PLUMorig_nfert_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears_orig,Nruns,'single') ;
PLUMorig_irrig_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears_orig,Nruns,'single') ;
PLUMharm_YXvyr = nan(ny,nx,Nlu,Nyears_harm,Nruns,'single') ;
PLUMharm_nfert_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears_harm,Nruns,'single') ;
PLUMharm_irrig_YXvyr = nan(ny,nx,Ncrops_lpjg,Nyears_harm,Nruns,'single') ;

for r = 1:Nruns
    thisRun = removeslashifneeded(runList{r}) ;

    % Original
    fprintf('Importing %s...\n', thisRun) ;
    [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
        [topDir thisRun], base_year, yearList_orig, ...
        thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, [], norm2extra, inpaint_method, '', false, true) ;
    [~,~,year_indices] = intersect(S_out.yearList,yearList_orig,'stable') ;
    if length(year_indices)~=length(yearList_orig)
        error('length(year_indices)~=length(yearList_orig)')
    end
    PLUMorig_YXvyr(:,:,:,:,r) = S_out.maps_YXvy(:,:,:,year_indices) ;
    PLUMorig_nfert_YXvyr(:,:,:,:,r) = S_nfert_out.maps_YXvy(:,:,:,year_indices) ;
    PLUMorig_irrig_YXvyr(:,:,:,:,r) = S_irrig_out.maps_YXvy(:,:,:,year_indices) ;
    clear S_*out

    % Harmonized
    fprintf('Importing %s.harm...\n', thisRun) ;
    if exist([topDir thisRun '.harm.forLPJG'],'dir')
        thisDir = [topDir thisRun '.harm.forLPJG/'] ;
        
        % Land use fractions
        S_lu = lpjgu_matlab_readTable_then2map([thisDir 'landcover.txt'],'force_mat_save',true) ;
        [~,year_indices,~] = intersect(S_lu.yearList,yearList_harm,'stable') ;
        if length(year_indices)~=length(yearList_harm)
            error('length(year_indices)~=length(yearList_harm)')
        end
        S_cropf = lpjgu_matlab_readTable_then2map([thisDir 'cropfractions.txt'],'force_mat_save',true) ;
        crops_tmp = strcat(LUnames(isCrop),'i') ;
        crops_tmp(strcmp(crops_tmp,'ExtraCropi')) = {'ExtraCrop'} ;
        [C_lu,~,indices_lu] = intersect(LUnames(~isCrop),S_lu.varNames,'stable') ;
        [   ~,~,indices_cf] = intersect(crops_tmp,S_cropf.varNames,'stable') ;
        [C_cf,~,         ~] = intersect(LUnames(isCrop),S_cropf.varNames,'stable') ;
        
        if ~isequal([C_cf C_lu],LUnames)
            error('~isequal([C_cf C_lu],LUnames)')
        end
        cropland_frac_YXvy = repmat(S_lu.maps_YXvy(:,:,strcmp(S_lu.varNames,'CROPLAND'),year_indices),[1 1 length(crops_tmp) 1]) ;
        PLUMharm_YXvyr(:,:,:,:,r) = repmat(thisLandArea_YX,[1 1 Nlu Nyears_harm]) ...
                                    .* cat(3, ...
                                           S_cropf.maps_YXvy(:,:,indices_cf,year_indices) .* cropland_frac_YXvy, ...
                                           S_lu.maps_YXvy(:,:,indices_lu,year_indices)) ;
        clear S_lu S_cropf cropland_frac_YXvy
        
        % Fertilization
        S = lpjgu_matlab_readTable_then2map([thisDir 'nfert.txt'],'force_mat_save',true) ;
        [~,IA] = intersect(S.varNames,crops_tmp,'stable') ;
        if length(IA) ~= Ncrops_lpjg
            error('length(IA)~=Ncrops_lpjg')
        end
        PLUMharm_nfert_YXvyr(:,:,:,:,r) = S.maps_YXvy(:,:,IA,year_indices) ;
        clear S
        
        % Irrigation
        S = lpjgu_matlab_readTable_then2map([thisDir 'irrig.txt'],'force_mat_save',true) ;
        [~,IA] = intersect(S.varNames,crops_tmp,'stable') ;
        if length(IA) ~= Ncrops_lpjg
            error('length(IA)~=Ncrops_lpjg')
        end
        PLUMharm_irrig_YXvyr(:,:,:,:,r) = S.maps_YXvy(:,:,IA,year_indices) ;
        clear S
    else
        [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
            [topDir thisRun '.harm'],base_year,yearList_harm, ...
            thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, [], 0, [], thisVer, read_MATs, false) ;
        PLUMharm_YXvyr(:,:,:,:,r) = S_out.maps_YXvy ;
        PLUMharm_nfert_YXvyr(:,:,:,:,r) = S_nfert_out.maps_YXvy ;
        PLUMharm_irrig_YXvyr(:,:,:,:,r) = S_irrig_out.maps_YXvy ;
        clear S_*out
    end

end

disp('Done reading PLUM.')


%% Time series of LUs

ts_base_cy = squeeze(nansum(nansum(base.maps_YXvy,1),2)) ;
ts_base_cy = cat(1,sum(ts_base_cy(isCrop,:),1),ts_base_cy(~isCrop,:)) ;
ts_orig_cyr = squeeze(nansum(nansum(PLUMorig_YXvyr,1),2)) ;
ts_orig_cyr = cat(1,sum(ts_orig_cyr(isCrop,:,:),1),ts_orig_cyr(~isCrop,:,:)) ;
ts_harm_cyr = squeeze(nansum(nansum(PLUMharm_YXvyr,1),2)) ;
ts_harm_cyr = cat(1,sum(ts_harm_cyr(isCrop,:,:),1),ts_harm_cyr(~isCrop,:,:)) ;
ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;

combinedLUs = [{'CROPLAND'} LUnames(~isCrop)] ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:length(combinedLUs)
    subplot_tight(2,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:)*1e-6*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList_orig,squeeze(ts_orig_cyr(v,:,:))*1e-6*1e-6,'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList_orig,squeeze(ts_harm_cyr(v,:,:))*1e-6*1e-6,'-','LineWidth',1)
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
ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:Ncrops_lpjg
    subplot_tight(4,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:)*1e-6*1e-6,'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList_orig,squeeze(ts_orig_cyr(v,:,:))*1e-6*1e-6,'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList_orig,squeeze(ts_harm_cyr(v,:,:))*1e-6*1e-6,'-','LineWidth',1)
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
    ts_base_cy = cf_kg2Mt .* squeeze(nansum(nansum(base_nfertTot_2deg.maps_YXvy))) ;
else
    ts_base_cy = cf_kg2Mt .* squeeze(nansum(nansum(base_nfertTot.maps_YXvy))) ;
end
ts_orig_cyr = cf_kg2Mt .* squeeze(nansum(nansum(PLUMorig_YXvyr(:,:,isCrop,:,:) .* PLUMorig_nfert_YXvyr,1),2)) ;
ts_harm_cyr = cf_kg2Mt .* squeeze(nansum(nansum(PLUMharm_YXvyr(:,:,isCrop,:,:) .* PLUMharm_nfert_YXvyr,1),2)) ;
ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:Ncrops_lpjg
    subplot_tight(4,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:),'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList_orig,squeeze(ts_orig_cyr(v,:,:)),'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList_orig,squeeze(ts_harm_cyr(v,:,:)),'-','LineWidth',1)
    hold off
    title(LPJGcrops{v})
    set(gca,'FontSize',14)
    ylabel('Mt')
    legend(legend_ts,'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_nfert.pdf']) ;
close


%% Time series of irrig

if is2deg
    ts_base_cy = squeeze(nansum(nansum(base_irrigTot_2deg.maps_YXvy))) ;
else
    ts_base_cy = squeeze(nansum(nansum(base_irrigTot.maps_YXvy))) ;
end
ts_orig_cyr = squeeze(nansum(nansum(PLUMorig_YXvyr(:,:,isCrop,:,:) .* PLUMorig_irrig_YXvyr,1),2)) ;
ts_harm_cyr = squeeze(nansum(nansum(PLUMharm_YXvyr(:,:,isCrop,:,:) .* PLUMharm_irrig_YXvyr,1),2)) ;
ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;

spacing = [0.05 0.1] ;

figure('Color','w','Position',figurePos)

for v = 1:Ncrops_lpjg
    subplot_tight(4,2,v,spacing) ;
    plot(yearList_luh2,ts_base_cy(v,:),'-k','LineWidth',2) ;
    set(gca,'ColorOrderIndex',1) ;
    hold on
    plot(yearList_orig,squeeze(ts_orig_cyr(v,:,:)),'--','LineWidth',1)
    set(gca,'ColorOrderIndex',1) ;
    plot(yearList_orig,squeeze(ts_harm_cyr(v,:,:)),'-','LineWidth',1)
    hold off
    title(LPJGcrops{v})
    set(gca,'FontSize',14)
    ylabel('Arbitrary units')
    legend(legend_ts,'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_irrig.pdf']) ;
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
            tmp = PLUMorig_YXvyr(:,:,v,yearList_harm==thisYear,r) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,3,y+3,spacing) ;
            tmp = PLUMharm_YXvyr(:,:,v,yearList_harm==thisYear,r) ;
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
            tmp1 = PLUMorig_YXvyr(:,:,v,yearList_harm==thisYear,r) ;
            tmp2 = PLUMharm_YXvyr(:,:,v,yearList_harm==thisYear,r) ;
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
            tmp1 = PLUMorig_YXvyr(:,:,v,yearList_harm==thisYear1,r) ;
            tmp2 = PLUMorig_YXvyr(:,:,v,yearList_harm==thisYear2,r) ;
            tmp = tmp2 - tmp1 ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(brewermap(64,'rdbu_ssr')) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d-%d',thisRun,thisLU,thisYear1,thisYear2)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,2,y+2,spacing) ;
            tmp1 = PLUMharm_YXvyr(:,:,v,yearList_harm==thisYear1,r) ;
            tmp2 = PLUMharm_YXvyr(:,:,v,yearList_harm==thisYear2,r) ;
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
