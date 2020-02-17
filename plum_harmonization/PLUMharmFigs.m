%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing harmonized PLUM land use trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Baseline version
% baseline_ver = 1 ;
% baseline_ver = 2 ;   % Based on remap_v6
baseline_ver = 3 ;   % Based on remap_v6p7

thisVer = '' ;
% thisVer = 'orig.' ;
% thisVer = '2deg.' ;

runList = {...
           'SSP1.v12.s1' ;
           'SSP3.v12.s1' ;
           'SSP4.v12.s1' ;
           'SSP5.v12.s1';
            } ;

base_year = 2010 ;

yearList_harm = 2011:2100 ;

norm2extra = 0.177 ;

% Method for inpaint_nans()
% (moved to PLUMharm_importRefData)


%% Setup

% Determine which system you're on
tmp = pwd ;
if strcmp(tmp(1:5),'/User')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp

if onMac
    out_dir = sprintf('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/harmonization_figs_v%d/',baseline_ver) ;
    topDir = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;
    addpath(genpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work')) ;
    PLUMharm_top = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/plum_harmonization/' ;
    inDir_protectedAreas = '/Users/Shared/PLUM/input/protected_areas/' ;
else
    out_dir = sprintf('/home/fh1-project-lpjgpi/lr8247/PLUM_harmonization_figs/harmonization_figs_v%d/',baseline_ver) ;
    topDir = addslashifneeded('/home/fh1-project-lpjgpi/lr8247/PLUM/input/PLUMouts_2011-2100') ;
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/paper02-matlab-work/')) ;
    PLUMharm_top = '/pfs/data1/home/kit/imk-ifu/lr8247/plum_harmonization/' ;
    inDir_protectedAreas = '/home/fh1-project-lpjgpi/lr8247/PLUM/input/protected_areas/' ;
end

if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

yearList_luh2 = 1971:2010 ;
% yearList_luh2 = 2001:2010 ;
yearList_orig = [yearList_harm(1)-1 yearList_harm] ;

if length(runList) == 1
    legend_ts = {'LUH2','Orig','Harm'} ;
else
    legend_ts = {'LUH2'} ;
    for s = 1:length(runList)
        legend_ts = [legend_ts {runList{s}(1:4)}] ;
    end
end

PLUM_in_toptop = strcat(topDir,runList) ;
PLUM_base_in = [addslashifneeded(PLUM_in_toptop{1}) '2010/'] ;

Nyears_orig = length(yearList_orig) ;
Nyears_harm = length(yearList_harm) ;
Nruns = length(runList) ;

% % Make lower-left lat/lon map (for compat. with PLUM style)
% lons_map_2deg = repmat(-180:2:178,[90 1]) ;
% lats_map_2deg = repmat((-90:2:88)',[1 180]) ;
% lons_map = repmat(-180:0.5:179.5,[360 1]) ;
% lats_map = repmat((-90:0.5:89.5)',[1 720]) ;

% Conversion factors
cf_kg2Mt = 1e-3*1e-6 ;


%% Import reference data

doHarm = false ;
run([PLUMharm_top 'PLUMharm_importRefData.m']) ;


%% Import PLUM (original + harmonized)

is2deg = strcmp(thisVer,'2deg.') ;

if is2deg
    ny = 90 ;
    nx = 180 ;
    thisLandArea_YX = landArea_2deg_YX ;
    thisLandArea_x = landArea_2deg_YX(list2map_2deg) ;
else
    ny = 360 ;
    nx = 720 ;
    thisLandArea_YX = landArea_YX ;
    thisLandArea_x = landArea_YX(list2map) ;
end

disp('Setting up PLUM*_xvyr arrays...')
Ncells = length(find(~mask_YX)) ;
PLUMorig_xvyr = nan(Ncells,Nlu,Nyears_orig,Nruns,'single') ;
PLUMorig_nfert_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_orig,Nruns,'single') ;
PLUMorig_irrig_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_orig,Nruns,'single') ;
PLUMharm_xvyr = nan(Ncells,Nlu,Nyears_harm,Nruns,'single') ;
PLUMharm_nfert_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_harm,Nruns,'single') ;
PLUMharm_irrig_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_harm,Nruns,'single') ;

for r = 1:Nruns
    thisRun = removeslashifneeded(runList{r}) ;

    % Original
    fprintf('Importing %s...\n', thisRun) ;
    tic
    [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
        [topDir thisRun], base_year, yearList_orig, ...
        thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
        is2deg, [], norm2extra, inpaint_method, '', true) ;
    [~,~,year_indices] = intersect(S_out.yearList,yearList_orig,'stable') ;
    if length(year_indices)~=length(yearList_orig)
        error('length(year_indices)~=length(yearList_orig)')
    end
    
    if length(year_indices) ~= size(S_out.maps_YXvy,4)
        S_out.maps_YXvy = S_out.maps_YXvy(:,:,:,year_indices) ;
    end
    incl_YXvy = repmat(~mask_YX, [1 1 size(S_out.maps_YXvy, 3:4)]) ;
    PLUMorig_xvyr(:,:,:,r) = reshape(S_out.maps_YXvy(incl_YXvy), size(PLUMorig_xvyr, 1:3)) ;
    clear S_out incl_YXvy
    
    if length(year_indices) ~= size(S_nfert_out.maps_YXvy,4)
        S_nfert_out.maps_YXvy = S_nfert_out.maps_YXvy(:,:,:,year_indices) ;
    end
    incl_YXvy = repmat(~mask_YX, [1 1 size(S_nfert_out.maps_YXvy, 3:4)]) ;
    PLUMorig_nfert_xvyr(:,:,:,r) = reshape(S_nfert_out.maps_YXvy(incl_YXvy), size(PLUMorig_nfert_xvyr, 1:3)) ;
    clear S_nfert_out
    
    if length(year_indices) ~= size(S_irrig_out.maps_YXvy,4)
        S_irrig_out.maps_YXvy = S_irrig_out.maps_YXvy(:,:,:,year_indices) ;
    end
    incl_YXvy = repmat(~mask_YX, [1 1 size(S_irrig_out.maps_YXvy, 3:4)]) ;
    PLUMorig_irrig_xvyr(:,:,:,r) = reshape(S_irrig_out.maps_YXvy(incl_YXvy), size(PLUMorig_irrig_xvyr, 1:3)) ;
    clear S_nfert_out incl_YXvy
    disp(toc_hms(toc))

    % Harmonized
    fprintf('Importing %s.harm...\n', thisRun) ;
    tic
    if exist([topDir thisRun '.harm.forLPJG'],'dir')
        thisDir = [topDir thisRun '.harm.forLPJG/'] ;
        
        % Land use fractions
        S_lu = lpjgu_matlab_readTable_then2map([thisDir 'landcover.txt'],'force_mat_save',true) ;
        [~,year_indices,~] = intersect(S_lu.yearList,yearList_harm,'stable') ;
        if length(year_indices)~=length(yearList_harm)
            error('length(year_indices)~=length(yearList_harm)')
        end
        if length(year_indices)~=size(S_lu.maps_YXvy, 4)
            S_lu.maps_YXvy = S_lu.maps_YXvy(:,:,:,year_indices) ;
        end
        S_cropf = lpjgu_matlab_readTable_then2map([thisDir 'cropfractions.txt'],'force_mat_save',true) ;
        if length(year_indices)~=size(S_cropf.maps_YXvy, 4)
            S_cropf.maps_YXvy = S_cropf.maps_YXvy(:,:,:,year_indices) ;
        end
        crops_tmp = strcat(LUnames(isCrop),'i') ;
        crops_tmp(strcmp(crops_tmp,'ExtraCropi')) = {'ExtraCrop'} ;
        [C_lu,~,indices_lu] = intersect(LUnames(~isCrop),S_lu.varNames,'stable') ;
        [   ~,~,indices_cf] = intersect(crops_tmp,S_cropf.varNames,'stable') ;
        [C_cf,~,         ~] = intersect(LUnames(isCrop),S_cropf.varNames,'stable') ;
        if ~isequal([C_cf C_lu],LUnames)
            error('~isequal([C_cf C_lu],LUnames)')
        end
        
        incl_YXvy = repmat(~mask_YX, [1 1 size(S_lu.maps_YXvy, 3:4)]) ;
        lu_xvy = reshape(S_lu.maps_YXvy(incl_YXvy), [Ncells size(S_lu.maps_YXvy, 3:4)]) ;
        lu_xvy = lu_xvy(:,indices_lu,:) ;
        cropland_frac_YXvy = repmat(S_lu.maps_YXvy(:,:,strcmp(S_lu.varNames,'CROPLAND'),:),[1 1 length(crops_tmp) 1]) ;
        clear S_lu incl_YXvy
        incl_YXvy = repmat(~mask_YX, [1 1 size(cropland_frac_YXvy, 3:4)]) ;
        cropland_frac_xvy = reshape(cropland_frac_YXvy(incl_YXvy), [Ncells size(cropland_frac_YXvy, 3:4)]) ;
        clear cropland_frac_YXvy incl_YXvy
        cropf_YXvy = S_cropf.maps_YXvy(:,:,indices_cf,:) ;
        clear S_cropf
        incl_YXvy = repmat(~mask_YX, [1 1 size(cropf_YXvy, 3:4)]) ;
        cropf_xvy = reshape(cropf_YXvy(incl_YXvy), [Ncells size(cropf_YXvy, 3:4)]) ;
        clear cropf_YXvy incl_YXvy
        PLUMharm_xvyr(:,:,:,r) = repmat(thisLandArea_x,[1 Nlu Nyears_harm]) ...
            .* cat(2, ...
                   cropf_xvy .* cropland_frac_xvy, ...
                   lu_xvy) ;
        clear cropf_xvy cropland_frac_xvy lu_xvy
        
        % Fertilization
        S = lpjgu_matlab_readTable_then2map([thisDir 'nfert.txt'],'force_mat_save',true) ;
        [~,~,IB] = intersect(crops_tmp,S.varNames,'stable') ;
        if length(IA) ~= Ncrops_lpjg
            error('length(IA)~=Ncrops_lpjg')
        end
        if ~isequal(IB, shiftdim(1:length(S.varNames))) ...
                || length(year_indices)~=size(S.maps_YXvy,4)
            S.maps_YXvy = S.maps_YXvy(:,:,IB,year_indices) ;
        end
        incl_YXvy = repmat(~mask_YX, [1 1 size(S.maps_YXvy, 3:4)]) ;
        PLUMharm_nfert_xvyr(:,:,:,r) = reshape( ...
            S.maps_YXvy(incl_YXvy), [Ncells size(S.maps_YXvy, 3:4)]) ;
        clear S
        
        % Irrigation
        S = lpjgu_matlab_readTable_then2map([thisDir 'irrig.txt'],'force_mat_save',true) ;
        [~,~,IB] = intersect(crops_tmp,S.varNames,'stable') ;
        if length(IA) ~= Ncrops_lpjg
            error('length(IA)~=Ncrops_lpjg')
        end
        if ~isequal(IB, shiftdim(1:length(S.varNames))) ...
                || length(year_indices)~=size(S.maps_YXvy,4)
            S.maps_YXvy = S.maps_YXvy(:,:,IB,year_indices) ;
        end
        incl_YXvy = repmat(~mask_YX, [1 1 size(S.maps_YXvy, 3:4)]) ;
        PLUMharm_irrig_xvyr(:,:,:,r) = reshape( ...
            S.maps_YXvy(incl_YXvy), [Ncells size(S.maps_YXvy, 3:4)]) ;
        clear S
    else
        [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
            [topDir thisRun '.harm'],base_year,yearList_harm, ...
            thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, [], 0, [], thisVer, false) ;
        
        if length(year_indices) ~= size(S_out.maps_YXvy,4)
            S_out.maps_YXvy = S_out.maps_YXvy(:,:,:,year_indices) ;
        end
        incl_YXvy = repmat(~mask_YX, [1 1 size(S_out.maps_YXvy, 3:4)]) ;
        PLUMharm_xvyr(:,:,:,r) = reshape(S_out.maps_YXvy(incl_YXvy), size(PLUMharm_xvyr, 1:3)) ;
        clear S_out incl_YXvy
        
        if length(year_indices) ~= size(S_nfert_out.maps_YXvy,4)
            S_nfert_out.maps_YXvy = S_nfert_out.maps_YXvy(:,:,:,year_indices) ;
        end
        incl_YXvy = repmat(~mask_YX, [1 1 size(S_nfert_out.maps_YXvy, 3:4)]) ;
        PLUMharm_nfert_xvyr(:,:,:,r) = reshape(S_nfert_out.maps_YXvy(incl_YXvy), size(PLUMharm_nfert_xvyr, 1:3)) ;
        clear S_nfert_out
        
        if length(year_indices) ~= size(S_irrig_out.maps_YXvy,4)
            S_irrig_out.maps_YXvy = S_irrig_out.maps_YXvy(:,:,:,year_indices) ;
        end
        PLUMharm_irrig_xvyr(:,:,:,r) = reshape(S_irrig_out.maps_YXvy(incl_YXvy), size(PLUMharm_irrig_xvyr, 1:3)) ;
        clear S_irrig_out incl_YXvy
        
    end
    disp(toc_hms(toc))

end

landArea_x = landArea_YX(list2map) ;
landArea_xv = repmat(landArea_x, [1 Nlu]) ;
landArea_xvr = repmat(landArea_xv, [1 1 Nruns]) ;

disp('Done reading PLUM.')


%% Time series of LUs

ts_base_cy = squeeze(nansum(nansum(base.maps_YXvy,1),2)) ;
ts_base_cy = cat(1,sum(ts_base_cy(isCrop,:),1),ts_base_cy(~isCrop,:)) ;
ts_orig_cyr = squeeze(nansum(PLUMorig_xvyr,1)) ;
ts_orig_cyr = cat(1,sum(ts_orig_cyr(isCrop,:,:),1),ts_orig_cyr(~isCrop,:,:)) ;
ts_harm_cyr = squeeze(nansum(PLUMharm_xvyr,1)) ;
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
    title(['Area: ' combinedLUs{v}])
    set(gca,'FontSize',14)
    ylabel('Million km2')
    legend(legend_ts,'Location','NorthWest')
end

% Save
export_fig([out_dir 'timeSeries_landUse.pdf']) ;
close


%% Time series of crops

ts_base_cy = squeeze(nansum(nansum(base.maps_YXvy,1),2)) ;
ts_orig_cyr = squeeze(nansum(PLUMorig_xvyr,1)) ;
ts_harm_cyr = squeeze(nansum(PLUMharm_xvyr,1)) ;
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
    title(['Area: ' LPJGcrops{v}])
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
ts_orig_cyr = cf_kg2Mt .* squeeze(nansum(PLUMorig_xvyr(:,isCrop,:,:) .* PLUMorig_nfert_xvyr,1)) ;
ts_harm_cyr = cf_kg2Mt .* squeeze(nansum(PLUMharm_xvyr(:,isCrop,:,:) .* PLUMharm_nfert_xvyr,1)) ;
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
    title(['Fert.: ' LPJGcrops{v}])
    set(gca,'FontSize',14)
    ylabel('Mt N')
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
ts_orig_cyr = squeeze(nansum(PLUMorig_xvyr(:,isCrop,:,:) .* PLUMorig_irrig_xvyr,1)) ;
ts_harm_cyr = squeeze(nansum(PLUMharm_xvyr(:,isCrop,:,:) .* PLUMharm_irrig_xvyr,1)) ;
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
    title(['Irrigation: ' LPJGcrops{v}])
    set(gca,'FontSize',14)
    ylabel('intensity \times area')
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
            tmp = 1e-6*lpjgu_vector2map(PLUMorig_xvyr(:,v,yearList_orig==thisYear,r), [ny nx], list2map) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,3,y+3,spacing) ;
            tmp = 1e-6*lpjgu_vector2map(PLUMharm_xvyr(:,v,yearList_harm==thisYear,r), [ny nx], list2map) ;
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
            tmp1 = 1e-6*lpjgu_vector2map(PLUMorig_xvyr(:,v,yearList_orig==thisYear,r), [ny nx], list2map) ;
            tmp2 = 1e-6*lpjgu_vector2map(PLUMharm_xvyr(:,v,yearList_harm==thisYear,r), [ny nx], list2map) ;
            tmp = tmp2 - tmp1 ;
%             tmp = tmp2/sum(tmp2(:)) - tmp1/sum(tmp1(:)) ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(flipud(brewermap(64,'rdbu_ssr'))) ;
            caxis([-1 1]*max(abs(caxis))) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('Harm-Orig, %s: %s, %d',thisRun,thisLU,thisYear)) ;
            set(gca,'FontSize',fontSize)
        end
        export_fig([out_dir 'mapsOHdiffs_' thisLU '_' strrep(num2str(theseYears),'  ','-') '_' thisRun '.png'],['-r' num2str(png_res)]) ;
        close
    end
end


%% Maps: Diffs between orig and harm at one year for each run

thisYear = 2011 ;
spacing = [0.05 0.025] ;
cbar_loc = 'SouthOutside' ;
y1 = 66 ;
fontSize = 14 ;
png_res = 150 ;
thisPos = figurePos ;
as_frac_land = true ;

tmpO_xvr = squeeze(PLUMorig_xvyr(:,:,yearList_orig==thisYear,:)) ;
tmpH_xvr = squeeze(PLUMharm_xvyr(:,:,yearList_harm==thisYear,:)) ;
if as_frac_land
    tmp_xvr = 100*(tmpH_xvr - tmpO_xvr) ./ landArea_xvr ;
else
    tmp_xvr = 1e-6*(tmpH_xvr - tmpO_xvr) ;
end
tmp_xvr(landArea_xvr==0) = NaN ;

for v = 1:Nlu
    thisLU = LUnames{v} ;
    figure('Color','w','Position',thisPos) ;
    if as_frac_land
%         new_caxis = [-100 100] ;
        new_caxis = [-1 1]*max(max(abs(tmp_xvr(:,v,:)))) ;
    else
        new_caxis = [-1 1]*max(max(abs(tmp_xvr(:,v,:)))) ;
    end
    for r = 1:Nruns
        thisRun = runList{r} ;
        h1 = subplot_tight(2,2,r,spacing) ;
        tmp = lpjgu_vector2map(tmp_xvr(:,v,r), [ny nx], list2map) ;
        pcolor(tmp(y1:end,:)) ;
        shading flat ; axis equal tight off
        colormap(flipud(brewermap(64,'rdbu_ssr'))) ;
        caxis(new_caxis) ;
        colorbar('Location',cbar_loc) ;
        title(thisRun) ;
        set(gca,'FontSize',fontSize)
    end
    if as_frac_land
        hsgt = sgtitle(sprintf('Harm-Orig (%%): %s, %d', thisLU, thisYear)) ;
    else
        hsgt = sgtitle(sprintf('Harm-Orig (km^2): %s, %d', thisLU, thisYear)) ;
    end
    set(hsgt, 'FontSize', fontSize+2, 'FontWeight', 'bold')
    filename = sprintf('%s/mapsOHdiffs_%s_%d.png', out_dir, thisLU, thisYear) ;
    if as_frac_land
        filename = strrep(filename, 'diffs', 'diffsFrac') ;
    end
    export_fig(filename,['-r' num2str(png_res)]) ;
    close
end
disp('Done')


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
            tmp1 = 1e-6*lpjgu_vector2map(PLUMorig_xvyr(:,v,yearList_orig==thisYear1,r), [ny nx], list2map) ;
            tmp2 = 1e-6*lpjgu_vector2map(PLUMorig_xvyr(:,v,yearList_orig==thisYear2,r), [ny nx], list2map) ;
            tmp = tmp2 - tmp1 ;
            tmp(landArea_YX==0) = NaN ;
            pcolor(tmp(y1:end,:)) ;
            shading flat ; axis equal tight off
            colormap(brewermap(64,'rdbu_ssr')) ;
            colorbar('Location',cbar_loc) ;
            title(sprintf('%s orig: %s, %d-%d',thisRun,thisLU,thisYear1,thisYear2)) ;
            set(gca,'FontSize',fontSize)
            h2 = subplot_tight(2,2,y+2,spacing) ;
            tmp1 = 1e-6*lpjgu_vector2map(PLUMharm_xvyr(:,v,yearList_harm==thisYear1,r), [ny nx], list2map) ;
            tmp2 = 1e-6*lpjgu_vector2map(PLUMharm_xvyr(:,v,yearList_harm==thisYear2,r), [ny nx], list2map) ;
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