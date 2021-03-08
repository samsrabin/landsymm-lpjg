%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing harmonized PLUM land use trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLUM_version = 'v12.s1' ;
% runList = strcat({...
%            'SSP1' ;
%            'SSP3' ;
%            'SSP4' ;
%            'SSP5';
%             }, ['.' PLUM_version]) ;
% yearList_harm = 2011:2100 ;
% fruitveg_sugar_2oil = false ;
% runList_legend = {'SSP1-45', 'SSP3-60', 'SSP4-60', 'SSP5-85'} ;
% % baseline_ver = 1 ;
% % baseline_ver = 2 ;   % Based on remap_v6
% baseline_ver = 3 ;   % Based on remap_v6p7

% runList = {...
%     'halfearth/HEoct/baseline/s1';
%     'halfearth/HEoct/halfearth/s1';
%     } ;
% yearList_harm = 2011:2060 ;
% fruitveg_sugar_2oil = true ;
% runList_legend = {'baseline', 'halfearth'} ;
% legend_ts = [{'LUH2'} runList_legend] ;
% baseline_ver = 4 ;   % Based on remap_v8b2oil

PLUM_version = 'v13.s1' ;
runList = {...
%     'ssp13/SSP1/s1' ;
    'ssp13/SSP2/s1' ;
    'ssp13/SSP3/s1' ;
    'ssp13/SSP4/s1' ;
    'ssp13/SSP5/s1' ;
    } ;
yearList_harm = 2011:2100 ;
fruitveg_sugar_2oil = false ;
runList_legend = strrep(strrep(runList, 'ssp13/',''),'/s1','')' ;
legend_ts = [{'LUH2'} runList_legend] ;
baseline_ver = 4 ;   % Based on remap_v8b2oil

% Use combined-crop version?
combineCrops = false ;

thisVer = '' ;
% thisVer = 'orig.' ;
% thisVer = '2deg.' ;

base_year = 2010 ;

norm2extra = 0.177 ;

% Method for inpaint_nans()
% (moved to PLUMharm_importRefData)


%% Setup

% Determine which system you're on
tmp = pwd ;
if strcmp(tmp(1:5),'/User') || strcmp(tmp(1:8),'/Volumes')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp

if onMac
    topDir = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;
    if any(contains(runList, 'halfearth/HEoct'))
        out_dir = sprintf('%s/halfearth/HEoct/harm_figs/', topDir) ;
    else
        out_dir = sprintf('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/harm%dFigs_%s/',baseline_ver,PLUM_version) ;
    end
    addpath(genpath('/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work')) ;
    PLUMharm_top = '/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/plum_harmonization/' ;
    addpath(genpath(PLUMharm_top))
    inDir_protectedAreas = '/Users/Shared/PLUM/input/protected_areas/' ;
else
    out_dir = sprintf('/home/fh1-project-lpjgpi/lr8247/PLUM_harmonization_figs/harm%dFigs_%s/',baseline_ver,PLUM_version) ;
    topDir = addslashifneeded('/home/fh1-project-lpjgpi/lr8247/PLUM/input/PLUMouts_2011-2100') ;
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/paper02-matlab-work/')) ;
    PLUMharm_top = '/pfs/data1/home/kit/imk-ifu/lr8247/plum_harmonization/' ;
    inDir_protectedAreas = '/home/fh1-project-lpjgpi/lr8247/PLUM/input/protected_areas/' ;
end
out_dir = get_harm_dir(out_dir, fruitveg_sugar_2oil, combineCrops) ;

if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

yearList_luh2 = 1971:2010 ;
% yearList_luh2 = 2001:2010 ;
yearList_orig = [yearList_harm(1)-1 yearList_harm] ;
add_baseline_to_harm = ...
    yearList_harm(1)==yearList_orig(1)+1 ...
    && yearList_orig(1)==base_year ...
    && any(yearList_luh2==base_year) ;
if ~exist('legend_ts', 'var')
    if length(runList) == 1
        legend_ts = {'LUH2','Orig','Harm'} ;
    else
        legend_ts = {'LUH2'} ;
        for s = 1:length(runList)
            legend_ts = [legend_ts {runList{s}(1:4)}] ;
        end
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

R = georasterref('RasterSize', [360 720], ...
    'RasterInterpretation', 'cells', ...
    'ColumnsStartFrom', 'north', ...
    'LatitudeLimits', [-90 90], ...
    'LongitudeLimits', [-180 180]) ;

% Get shortest distinct runList
shortened_runList = runList ;
tmp = shortened_runList ;
while length(unique(tmp)) == Nruns
    shortened_runList = tmp ;
    for r = 1:Nruns
        tmp{r} = tmp{r}(1:end-1) ;
    end
end


%% Import reference data

doHarm = false ;

run([PLUMharm_top 'PLUMharm_importRefData.m']) ;

biomeID_YX = flipud(imread( ...
    '/Users/sam/Geodata/General/WWF terrestrial ecosystems/wwf_terr_ecos_dissolveBiome_halfDeg_id.tif')) ;
biomeID_YX(biomeID_YX<0) = NaN ;
biomeID_x = biomeID_YX(list2map) ;
countries_YX = flipud(dlmread( ...
    '/Users/Shared/PLUM/crop_calib_data/countries/PLUM-specific/country_boundaries62892.noNeg99.extrapd.asc', ...
    '',6,0)) ;
countries_YX(countries_YX<=0) = NaN ;
countries_x = countries_YX(list2map) ;
countries_key = readtable('/Users/Shared/lpj-guess-plum/input/unifying_gridlist/country_boundaries_codes4.csv') ;



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
PLUMharm_xvyr = nan(Ncells,Nlu,Nyears_harm,Nruns,'single') ;
if ~combineCrops
    PLUMorig_nfert_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_orig,Nruns,'single') ;
    PLUMorig_irrig_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_orig,Nruns,'single') ;
    PLUMharm_nfert_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_harm,Nruns,'single') ;
    PLUMharm_irrig_xvyr = nan(Ncells,Ncrops_lpjg,Nyears_harm,Nruns,'single') ;
end

for r = 1:Nruns
    thisRun = removeslashifneeded(runList{r}) ;

    % Original
    fprintf('Importing %s...\n', thisRun) ;
    tic
    if combineCrops
        [S_out, ~, ~] = PLUMharm_pp_readPLUM(...
            [topDir thisRun], base_year, yearList_orig, ...
            thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, [], norm2extra, inpaint_method, '', true, ...
            fruitveg_sugar_2oil) ;
        S_out.maps_YXvy = cat(3, ...
            sum(S_out.maps_YXvy(:,:,~contains(S_out.varNames,{'PASTURE','NATURAL','BARREN'}),:), 3), ...
            S_out.maps_YXvy(:,:,contains(S_out.varNames,{'PASTURE','NATURAL','BARREN'}),:)) ;
        S_out.varNames = [LPJGcrops, S_out.varNames(contains(S_out.varNames,{'PASTURE','NATURAL','BARREN'}))] ;
    else
        [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
            [topDir thisRun], base_year, yearList_orig, ...
            thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
            is2deg, [], norm2extra, inpaint_method, '', true, ...
            fruitveg_sugar_2oil) ;
    end
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
    
    if ~combineCrops
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
    end
    disp(toc_hms(toc))

    % Harmonized
    fprintf('Importing %s.harm...\n', thisRun) ;
    tic
    thisDir = [topDir thisRun '.harm'] ;
    thisDir = get_harm_dir(thisDir, fruitveg_sugar_2oil, combineCrops) ;
    thisDir = [thisDir '.forLPJG/'] ;
    if exist(thisDir,'dir')
        if combineCrops
            error('Need to rework this code to work with combineCrops')
        end
        
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
        tmpDir = [topDir thisRun '.harm'] ;
        tmpDir = get_harm_dir(tmpDir, fruitveg_sugar_2oil, combineCrops) ;
        if combineCrops
            [S_out, ~, ~] = PLUMharm_pp_readPLUM(...
                tmpDir,base_year,yearList_harm, ...
                thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
                is2deg, [], 0, [], thisVer, false, fruitveg_sugar_2oil) ;
            S_out.maps_YXvy = cat(3, ...
                sum(S_out.maps_YXvy(:,:,~contains(S_out.varNames,{'PASTURE','NATURAL','BARREN'}),:), 3), ...
                S_out.maps_YXvy(:,:,contains(S_out.varNames,{'PASTURE','NATURAL','BARREN'}),:)) ;
            S_out.varNames = [LPJGcrops, S_out.varNames(contains(S_out.varNames,{'PASTURE','NATURAL','BARREN'}))] ;
        else
            [S_out, S_nfert_out, S_irrig_out] = PLUMharm_pp_readPLUM(...
                tmpDir,base_year,yearList_harm, ...
                thisLandArea_YX, LUnames, PLUMtoLPJG, LPJGcrops, ...
                is2deg, [], 0, [], thisVer, false, fruitveg_sugar_2oil) ;
        end
        clear tmpDir
        
        if length(year_indices) ~= size(S_out.maps_YXvy,4) && ~add_baseline_to_harm
            S_out.maps_YXvy = S_out.maps_YXvy(:,:,:,year_indices) ;
        end
        incl_YXvy = repmat(~mask_YX, [1 1 size(S_out.maps_YXvy, 3:4)]) ;
        PLUMharm_xvyr(:,:,:,r) = reshape(S_out.maps_YXvy(incl_YXvy), size(PLUMharm_xvyr, 1:3)) ;
        clear S_out incl_YXvy
        
        if ~combineCrops
            if length(year_indices) ~= size(S_nfert_out.maps_YXvy,4) && ~add_baseline_to_harm
                S_nfert_out.maps_YXvy = S_nfert_out.maps_YXvy(:,:,:,year_indices) ;
            end
            incl_YXvy = repmat(~mask_YX, [1 1 size(S_nfert_out.maps_YXvy, 3:4)]) ;
            PLUMharm_nfert_xvyr(:,:,:,r) = reshape(S_nfert_out.maps_YXvy(incl_YXvy), size(PLUMharm_nfert_xvyr, 1:3)) ;
            clear S_nfert_out
            
            if length(year_indices) ~= size(S_irrig_out.maps_YXvy,4) && ~add_baseline_to_harm
                S_irrig_out.maps_YXvy = S_irrig_out.maps_YXvy(:,:,:,year_indices) ;
            end
            PLUMharm_irrig_xvyr(:,:,:,r) = reshape(S_irrig_out.maps_YXvy(incl_YXvy), size(PLUMharm_irrig_xvyr, 1:3)) ;
            clear S_irrig_out incl_YXvy
        end
        
    end
    disp(toc_hms(toc))
    
end

gcelArea_x = gcelArea_YX(list2map) ;
landArea_x = landArea_YX(list2map) ;
landArea_xv = repmat(landArea_x, [1 Nlu]) ;
landArea_xvr = repmat(landArea_xv, [1 1 Nruns]) ;

% Import food production units
fpu_YX = flipud(dlmread('/Users/Shared/PLUM/food_production_units/FPU.asc',' ',6,0)) ;
fpu_YX(fpu_YX==-9999) = NaN ;
fpu_x = fpu_YX(list2map) ;
fpu_list = unique(fpu_x(~isnan(fpu_x))) ;
Nfpu = length(fpu_list) ;

map_size = size(landArea_YX) ;
if add_baseline_to_harm
    % Area
    tmp_base_YXv = base.maps_YXvy(:,:,:,yearList_luh2==base_year) ;
    tmp_base_xv = nan(Ncells, Nlu) ;
    tmp_list2map_all = lpjgu_get_list2map_all(list2map, map_size, Nlu) ;
    tmp_base_xv(:) = tmp_base_YXv(tmp_list2map_all) ;
    PLUMharm_xvyr = cat(3, ...
        repmat(tmp_base_xv, [1 1 1 Nruns]), ...
        PLUMharm_xvyr) ;
    clear tmp*
    
    if ~combineCrops
        % Fert
        tmp_base_YXv = base_nfert.maps_YXvy(:,:,:,yearList_luh2==base_year) ;
        tmp_base_xv = nan(Ncells, Ncrops_lpjg) ;
        tmp_list2map_all = lpjgu_get_list2map_all(list2map, map_size, Ncrops_lpjg) ;
        tmp_base_xv(:) = tmp_base_YXv(tmp_list2map_all) ;
        PLUMharm_nfert_xvyr = cat(3, ...
            repmat(tmp_base_xv, [1 1 1 Nruns]), ...
            PLUMharm_nfert_xvyr) ;
        clear tmp*
        
        % Irrig
        tmp_base_YXv = base_irrig.maps_YXvy(:,:,:,yearList_luh2==base_year) ;
        tmp_base_xv = nan(Ncells, Ncrops_lpjg) ;
        tmp_list2map_all = lpjgu_get_list2map_all(list2map, map_size, Ncrops_lpjg) ;
        tmp_base_xv(:) = tmp_base_YXv(tmp_list2map_all) ;
        PLUMharm_irrig_xvyr = cat(3, ...
            repmat(tmp_base_xv, [1 1 1 Nruns]), ...
            PLUMharm_irrig_xvyr) ;
        clear tmp*
    end
    
    yearList_harm = [base_year yearList_harm] ;
end

disp('Done reading PLUM.')


%% Scatter plots after Hurtt et al. (2011) Fig. 4 (crop, pre-harm)
do_save = true ;

% Options %%%
ny = 1 ;
nx = Nruns ;
spacing = [0.05 0.05] ; % v h
fontSize = 14 ;
thisPos = [0         324        1440         376] ;
%%%%%%%%%%%%%

y1 = 2010 ;

figure('Color','w','Position',thisPos) ;

for r = 1:Nruns
    
    % Establish axis
    subplot_tight(ny, nx, r, spacing) ;
    
    % Plot crop scatter
    plot( ...
        sum(PLUMharm_xvyr(:,isCrop,yearList_harm==y1,r),2) ./ gcelArea_x, ...
        sum(PLUMorig_xvyr(:,isCrop,yearList_orig==y1,r),2) ./ gcelArea_x, ...
        '.k')
    
    % Finish up
    axis equal tight ;
    set(gca,'XLim',[0 1],'YLim',[0 1], 'FontSize', fontSize)
    title(runList{r})
    xlabel(sprintf('Fraction of gridcell %d (LUH2)', y1))
    ylabel(sprintf('Fraction of gridcell %d (PLUM output)', y1))
    
end

if do_save
    export_fig([out_dir 'scatter_hurtt2011_fig4.png'], '-r300') ;
    close
end



%% Scatter plots after Hurtt et al. (2011) Fig. 5 (crop and past, post-harm)
do_save = true ;

% Options %%%
ny = 2 ;
nx = Nruns ;
spacing = [0.1 0.05] ; % v h
fontSize = 14 ;
%%%%%%%%%%%%%

y1 = 2010 ;
yN = yearList_harm(end) ;

thisGray = 0.65*ones(3,1) ;

diff_crop_orig_xr = ...
    squeeze(sum(PLUMorig_xvyr(:,isCrop,yearList_orig==yN,:) ...
    - PLUMorig_xvyr(:,isCrop,yearList_orig==y1,:),2)) ...
    ./ repmat(gcelArea_x, [1 Nruns]) ;
diff_past_orig_xr = ...
    squeeze(PLUMorig_xvyr(:,strcmp(LUnames,'PASTURE'),yearList_orig==yN,:) ...
    - PLUMorig_xvyr(:,strcmp(LUnames,'PASTURE'),yearList_orig==y1,:)) ...
    ./ repmat(gcelArea_x, [1 Nruns]) ;
diff_orig_xrL = cat(3, diff_crop_orig_xr, diff_past_orig_xr) ;
diff_crop_harm_xr = ...
    squeeze(sum(PLUMharm_xvyr(:,isCrop,yearList_harm==yN,:) ...
    - PLUMharm_xvyr(:,isCrop,yearList_harm==y1,:),2)) ...
    ./ repmat(gcelArea_x, [1 Nruns]) ;
diff_past_harm_xr = ...
    squeeze(PLUMharm_xvyr(:,strcmp(LUnames,'PASTURE'),yearList_harm==yN,:) ...
    - PLUMharm_xvyr(:,strcmp(LUnames,'PASTURE'),yearList_harm==y1,:)) ...
    ./ repmat(gcelArea_x, [1 Nruns]) ;
diff_harm_xrL = cat(3, diff_crop_harm_xr, diff_past_harm_xr) ;

diff2_orig_xrL = nan(length(list2map_2deg), Nruns, 2) ;
diff2_harm_xrL = nan(length(list2map_2deg), Nruns, 2) ;
map_size = size(landArea_YX) ;
tmp = ...
    gcelArea_YX(:,1:4:720) + gcelArea_YX(:,2:4:720) + ...
    gcelArea_YX(:,3:4:720) + gcelArea_YX(:,4:4:720) ;
gcelArea_2deg_YX = ...
    tmp(1:4:360,:) + tmp(2:4:360,:) + ...
    tmp(3:4:360,:) + tmp(4:4:360,:) ;
clear tmp
for r = 1:Nruns
    for L = 1:2
        % Orig
        tmp_YX = lpjgu_vector2map(diff_orig_xrL(:,r,L).*gcelArea_x, ...
            map_size, list2map) ;
        tmp = ...
            tmp_YX(:,1:4:720) + tmp_YX(:,2:4:720) + ...
            tmp_YX(:,3:4:720) + tmp_YX(:,4:4:720) ;
        tmp_2deg_YX = ...
            tmp(1:4:360,:) + tmp(2:4:360,:) + ...
            tmp(3:4:360,:) + tmp(4:4:360,:) ;
        clear tmp_YX tmp
        tmp_2deg_YX = tmp_2deg_YX ./ gcelArea_2deg_YX ;
        diff2_orig_xrL(:,r,L) = tmp_2deg_YX(list2map_2deg) ;
        clear tmp_2deg_YX
        
        % Harm
        tmp_YX = lpjgu_vector2map(diff_harm_xrL(:,r,L).*gcelArea_x, ...
            map_size, list2map) ;
        tmp = ...
            tmp_YX(:,1:4:720) + tmp_YX(:,2:4:720) + ...
            tmp_YX(:,3:4:720) + tmp_YX(:,4:4:720) ;
        tmp_2deg_YX = ...
            tmp(1:4:360,:) + tmp(2:4:360,:) + ...
            tmp(3:4:360,:) + tmp(4:4:360,:) ;
        clear tmp_YX tmp
        tmp_2deg_YX = tmp_2deg_YX ./ gcelArea_2deg_YX ;
        diff2_harm_xrL(:,r,L) = tmp_2deg_YX(list2map_2deg) ;
        clear tmp_2deg_YX
    end
end

diffFPU_orig_xrL = nan(Nfpu, Nruns, 2) ;
diffFPU_harm_xrL = nan(Nfpu, Nruns, 2) ;
for f = 1:Nfpu
    thisFPU = fpu_list(f) ;
    isThisFPU = fpu_x==thisFPU ;
    if ~any(isThisFPU)
        error('No cells in this FPU?')
    end
    gcelArea_thisFPU_x = gcelArea_x(isThisFPU) ;
    for r = 1:Nruns
        for L = 1:2
            diffFPU_orig_xrL(f,r,L) = sum(diff_orig_xrL(isThisFPU,r,L) ...
                .* gcelArea_thisFPU_x) ./ sum(gcelArea_thisFPU_x) ;
            diffFPU_harm_xrL(f,r,L) = sum(diff_harm_xrL(isThisFPU,r,L) ...
                .* gcelArea_thisFPU_x) ./ sum(gcelArea_thisFPU_x) ;
            clear tmp_2deg_YX
        end
    end
end

figure('Color','w','Position',figurePos) ;

for r = 1:Nruns
    for L = 1:2
        % Establish axis
        if L==1
            thisPlot = r ;
        else
            thisPlot = Nruns + r ;
        end
        subplot_tight(ny, nx, thisPlot, spacing) ;
        
        % Plot half-degree points
        plot(diff_orig_xrL(:,r,L), diff_harm_xrL(:,r,L), '.', ...
            'MarkerFaceColor', thisGray, 'MarkerEdgeColor', thisGray);
        
        % Plot two-degree points
        hold on
        plot(diff2_orig_xrL(:,r,L), diff2_harm_xrL(:,r,L), '.k') ;
        hold off
        
%         % Plot FPU points
%         hold on
%         plot(diffFPU_orig_xrL(:,r,L), diffFPU_harm_xrL(:,r,L), 'or') ;
%         hold off
        
        % Plot 1:1 line
        hold on
        plot([-1 1], [-1 1], '--k')
        hold off
        
        % Finish up
        axis equal tight ;
        set(gca,'XLim',[-1 1],'YLim',[-1 1], 'FontSize', fontSize)
        title(runList{r})
        xlabel('\Delta gridcell fraction (original)')
        ylabel('\Delta gridcell fraction (harmonized)')
    end
end

if do_save
    export_fig([out_dir 'scatter_hurtt2011_fig5.png'], '-r300') ;
    close
end


%% Map deltas for orig and harm

tmp_lu_list = {'CROPLAND','PASTURE','NATURAL'} ;
% y1_list = 2011 ;
% yN_list= 2012 ;
% y1_list = 2010 ;
% yN_list= yearList_harm(end) ;
% y1_list = 2011:1:2099 ;
% yN_list = 2012:1:2100 ;
% y1_list = 2011:5:2099 ;
% yN_list = 2015:5:2100 ;

yN_list = 2030:10:2060 ;
y1_list = 2010*ones(size(yN_list)) ;

% Options %%%%%%%%%
fontSize = 14 ;
% spacing = [0.02 0.02] - 0.0025*8 ;   % [vert, horz]
spacing = 0 ;
textX = 0.115 ;
textY_1 = 50/360 ;
textY_2 = 20/360 ;
shiftup = 15/360 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ; 
nx = 2 ;
ny = length(runList) ;
as_frac_land = false ;
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
units_total = 'Mkm^2' ;
do_caps = false ;
ntrl_colormap_name = 'PiYG_ssr' ;
% lines_overlay = 'landareas.shp' ;
lines_overlay = '/Users/sam/Geodata/General/continents_from_countries/continents_from_countries.shp' ;
%%%%%%%%%%%%%%%%%%%

if ny == 4
    thisPos = [1    33   770   772] ;
elseif ny==2
    thisPos = [1   308   770   497] ;
else
    error('Set thisPos for ny = %d', ny)
end
textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ; 

if length(y1_list) > 1
    this_outdir = sprintf('%s/maps_manyDeltas_beforeAfter_%d-%d_by%d', ...
        out_dir, min(y1_list), max(yN_list), yN_list(1)-y1_list(1)+1) ;
    pngres = '-r150' ;
else
    this_outdir = out_dir ;
    pngres = '-r300' ;
end
this_outdir_geo = [removeslashifneeded(out_dir) 'geo/'] ;

if ~exist(this_outdir, 'dir')
    mkdir(this_outdir) ;
end
if ~exist(this_outdir_geo, 'dir')
    mkdir(this_outdir_geo) ;
end

map_size = size(landArea_YX) ;
orig_diff_YXrH = nan([map_size Nruns]) ;
harm_diff_YXrH = nan([map_size Nruns]) ;

if as_frac_land
    bins_lowBnds = [-100:20:-20 -5 [-1 1]*0.1 5 20:20:80] ;
    conv_fact_map = 1 ;
    units_map = '%' ;
else
    bins_lowBnds = [] ;
%     conv_fact_map = 1e-6 ;   % m2 to km2
%     units_map = 'km^2' ;
    conv_fact_map = 1e-4*1e-6 ;   % m2 to Mha
    units_map = 'Mha' ;
end

for l = 1:length(tmp_lu_list)
    
    % Get long and short LU name
    thisLU = tmp_lu_list{l} ;
    thisLU_short = lower(thisLU(1:4)) ;
    if strcmp(thisLU_short, 'natu')
        thisLU_short = 'ntrl' ;
    end
    
    if strcmp(thisLU, 'CROPLAND')
        v = isCrop ;
    else
        v = contains(LUnames, thisLU) ;
    end
    if ~strcmp(thisLU, 'NATURAL')
        if strcmp(ntrl_colormap_name(1), '-')
            this_colormap_name = ntrl_colormap_name(2:end) ;
        else
            this_colormap_name = ['-' ntrl_colormap_name] ;
        end
    end
    
    % Get color map
    if strcmp(thisLU, 'CROPLAND')
        v = isCrop ;
    else
        v = contains(LUnames, thisLU) ;
    end
    this_colormap_name = ntrl_colormap_name ;
    if ~strcmp(thisLU, 'NATURAL')
        if strcmp(this_colormap_name(1), '-')
            this_colormap_name = this_colormap_name(2:end) ;
        else
            this_colormap_name = ['-' this_colormap_name] ;
        end
    end
    for y = 1:length(y1_list)
        
        y1 = y1_list(y) ;
        yN = yN_list(y) ;
        
        area_orig_bl_r = squeeze(sum(sum( ...
            PLUMorig_xvyr(:,v,yearList_orig==y1,:),1),2))*conv_fact_total ;
        area_harm_bl_r = squeeze(sum(sum( ...
            PLUMharm_xvyr(:,v,yearList_harm==y1,:),1),2))*conv_fact_total ;
        
        total_origDiff_r = squeeze(sum(sum( ...
            PLUMorig_xvyr(:,v,yearList_orig==yN,:),1),2))*conv_fact_total ...
            - area_orig_bl_r ;
        total_harmDiff_r = squeeze(sum(sum( ...
            PLUMharm_xvyr(:,v,yearList_harm==yN,:),1),2))*conv_fact_total ...
            - area_harm_bl_r ;
        
        % Get difference (%)
        diff_orig_xr = squeeze(nansum( ...
            PLUMorig_xvyr(:,v,yearList_orig==yN,:) ...
            - PLUMorig_xvyr(:,v,yearList_orig==y1,:), ...
            2)) ;
        diff_harm_xr = squeeze(nansum( ...
            PLUMharm_xvyr(:,v,yearList_harm==yN,:) ...
            - PLUMharm_xvyr(:,v,yearList_harm==y1,:), ...
            2)) ;
        for r = 1:Nruns
            orig_diff_YXrH(:,:,r) = lpjgu_vector2map(100*diff_orig_xr(:,r), map_size, list2map) ;
            harm_diff_YXrH(:,:,r) = lpjgu_vector2map(100*diff_harm_xr(:,r), map_size, list2map) ;
        end
        if as_frac_land
            orig_diff_YXrH = orig_diff_YXrH ./ repmat(gcelArea_YX, [1 1 Nruns]) ;
            harm_diff_YXrH = harm_diff_YXrH ./ repmat(gcelArea_YX, [1 1 Nruns]) ;
        end
        orig_diff_YXrH = orig_diff_YXrH * conv_fact_map ;
        harm_diff_YXrH = harm_diff_YXrH * conv_fact_map ;
        
        
%         if ~as_frac_land
%             error('This only works with as_frac_land TRUE')
%         end
        
        col_titles = {sprintf('Original %s %s, %d%s%d', '\Delta', thisLU, y1, char(8211), yN), ...
            sprintf('Harmonized %s %s, %d%s%d', '\Delta', thisLU, y1, char(8211), yN)} ;
        [diff_crop_YXr, diff_past_YXr] = make_LUdiff_fig_v5(...
            area_orig_bl_r, area_harm_bl_r, total_origDiff_r, total_harmDiff_r, ...
            orig_diff_YXrH, harm_diff_YXrH, ...
            y1, yN, runList_legend, ...
            spacing, fontSize, textX, textY_1, textY_2, ...
            nx, ny, ...
            Nruns, thisPos, units_map, units_total, do_caps, ...
            bins_lowBnds, this_colormap_name, col_titles, ...
            lines_overlay) ;

        if do_save
            filename = sprintf('%s/maps_deltas_%s_%d-%d_beforeAfter.png', ...
                this_outdir, thisLU, y1, yN) ;
            if ~as_frac_land
                filename = strrep(filename, '.png', sprintf('.%s.png', units_map)) ;
            end
            export_fig(filename, pngres) ;
            close
            if as_frac_land && length(y1_list) == 1
                disp('Saving GeoTIFFs...')
                for r = 1:Nruns
                    filename = sprintf('%s/D%s_%d-%d_%s_orig.tif', ...
                        removeslashifneeded(this_outdir_geo), ...
                        thisLU_short, y1_list, yN_list, runList_legend{r}) ;
                    geotiffwrite_ssr(filename,orig_diff_YXrH(:,:,r),R,-999)
                    filename = strrep(filename, 'orig', 'harm') ;
                    geotiffwrite_ssr(filename,harm_diff_YXrH(:,:,r),R,-999)
                end
            end
        end
    end
    
end
disp('Done')


%% Harmonization effects by the numbers

incl_years = 2010:2100 ;

harm_focus_regions = { ...
... Number  Super-region    General biome           Geog. restriction?      Description
    unique(biomeID_x), ...
            'World',        'World',                [],                     'World' ;
    1095,   'Amazon',       'Trop. rainforest',     [],                     'Amazon' ;
    429,    'N. America',   'Temp. grassland',      [],                     'Great Plains' ;
    437,    'N. America',   'Temp. forest',         [],                     'E US mixed for' ;
    306,    'N. America',   'Temp. forest',         [],                     'U. Midw US br/mix for' ;
    477,    'N. America',   'Temp. forest',         [],                     'E US conif for' ;
    436,    'N. America',   'Temp. forest',         [],                     'Texarkana conif for' ;
    229,    'Alaska',       'Tundra',               [],                     'Alaskan forest' ;
    [62 63 68 72 78 79 80 87 95:97 101 102 109:113 124:127 ...
     131 133 134:136 139 143 144 150 163 175 202 213 220 229]', ...
            'Alaska',       'Bor. forest',               {'United States of America'}, ...
                                                                            'Alaskan tundra' ;
    950,    'Sub-Sah. Afr.','Trop. rainforest',     [],                     'C Afr rainfor.' ;
    1251,   'Sub-Sah. Afr.','Savanna',              '0&N',                  'N Afr savanna' ;
    1251,   'Sub-Sah. Afr.','Savanna',              '0&S',                  'S Afr savanna' ;
    811,    'South Asia',   'Desert/xeric',         {'India', 'Pakistan', 'Afghanistan'}, ...
                                                                            'S Asia xeric/desert' ;
    661,    'South Asia',   '(Sub)trop. dry for.',  [],                     'C Ind subt dry for' ;
    712,    'South Asia',   '(Sub)trop. dry for.',  [],                     'S Ind subt dry for' ;
    766,    'South Asia',   '(Sub)trop. dry for.',  [],                     'S Ind scrub for' ;
    808,    'South Asia',   '(Sub)trop. dry for.',  [],                     'SriL subt dry for' ;
    708,    'South Asia',   '(Sub)trop. wet for.',  [],                     'W Ind subt wet for' ;
    830,    'South Asia',   '(Sub)trop. wet for.',  [],                     'SriL subt wet for' ;
    597,    'South Asia',   '(Sub)trop. wet for.',  [],                     'C Ind subt wet for' ;
    [551;567;632;638], ...
            'South Asia',   '(Sub)trop. wet for.',  [],                     'E Ind subt wet for' ;
    662,    'South Asia',   '(Sub)trop. wet for.',  {'India', 'Bangladesh'},'NWInd+Bangl subt wet for' ;
    480,    'East Asia',    'Temp. forest',         [],                     'E Asia temp for' ;
    487,    'East Asia',    'Montane gr/shr',       [],                     'Tibetan Plat. steppe' ;
    342,    'East Asia',    'Desert/xeric',         [],                     'E Asia xeric/desert' ;
    348,    'East Asia',    'Temp. grass/sav/shr',  [],                     'E Asia temp grass' ;
    480,    'China',        'Temp. forest',         {'China'},              'E China temp for' ;
    487,    'China',        'Montane gr/shr',       {'China'},              'China Tib. Plat. steppe' ;
    342,    'China',        'Desert/xeric',         {'China'},              'China xeric/desert' ;
    348,    'China',        'Temp. grass/sav/shr',  {'China'},              'China temp grass' ;
    662,    'China',        '(Sub)trop. wet for.',  {'China'},              'China subt wet for' ;
    setdiff(unique(biomeID_x(countries_x==countries_key.numCode(strcmp(countries_key.Country,'China')))), ...
            [480 487 342 348 662]), ...
            'China',        'China other',          {'China'},              'China other' ;
    [132;153;169;171;252;313], ...
            'Europe+Nafr',  'Temp. forest',         [],                     'Eur temp br/mix for' ;
    [168;287;288;347], ...
            'Europe+Nafr',  'Temp. forest',         [],                     'Eur temp conif for' ;
    349,    'Europe+Nafr',  'Mediterranean',        [],                     'Mediterr. mediterr.' ;
    } ;

[~, yi_orig] = intersect(yearList_orig, incl_years) ;
[~, yi_harm] = intersect(yearList_harm, incl_years) ;
if length(yi_orig) ~= length(incl_years)
    error('length(yi_orig) ~= length(incl_years)')
elseif length(yi_harm) ~= length(incl_years)
    error('length(yi_harm) ~= length(incl_years)')
end
y1 = min(incl_years) ;
yN = max(incl_years) ;

tmp_lu_list = {'NATURAL','CROPLAND','PASTURE'} ;

% Setup
list_regions = harm_focus_regions(:,5) ;
list_superRegs = unique(harm_focus_regions(:,2)) ;
Nregions = length(list_regions) ;
NsuperRegs = length(list_superRegs) ;
harmByNums_outdir = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/harmonization_qgis' ;
templatefile = sprintf('%s/harm_by_numbers.template.xlsx', harmByNums_outdir) ;
templatefile_relLandArea = sprintf('%s/harm_by_numbers.template.relLandArea.xlsx', harmByNums_outdir) ;

% Do it
for l = 1:length(tmp_lu_list)
    thisLU = tmp_lu_list{l} ;
    outfile = sprintf('%s/harm_by_numbers.%d-%d.%s.xlsx', ...
        harmByNums_outdir, ...
        y1, yN, thisLU) ;
    if combineCrops
        outfile = strrep(outfile, thisLU, ['combCrops.' thisLU]) ;
    end
    outfile_relLandArea = strrep(outfile, '.xlsx', '.relLandArea.xlsx') ;
    
    if strcmp(thisLU, 'CROPLAND')
        v = isCrop ;
    else
        v = strcmp(LUnames, thisLU) ;
    end
    
    % Set up Excel files based on template
    copyfile(templatefile, outfile) ;
    copyfile(templatefile_relLandArea, outfile_relLandArea) ;
    
    for r = 1:Nruns
        thisRun = runList{r} ;
        fprintf('%d-%d, %s, %s\n', y1, yN, thisRun, thisLU)
        landArea_byReg = [] ;
        for s = 1:NsuperRegs
            [table_orig, table_orig_relY1, landArea_byReg] = PLUMharmFigs_iterate_superReg( ...
                sum(PLUMorig_xvyr(:,v,yi_orig,r),2), s, 'orig', ...
                list_superRegs, list_regions, harm_focus_regions, ...
                biomeID_x, countries_x, countries_key, lats, lons, ...
                list2map, landArea_x) ;
            if s==1
                table_orig_out = table_orig ;
                table_orig_relY1_out = table_orig_relY1 ;
                landArea_byReg_out = landArea_byReg ;
            else
                table_orig_out = cat(1, table_orig_out, table_orig) ;
                table_orig_relY1_out = cat(1, table_orig_relY1_out, table_orig_relY1) ;
                landArea_byReg_out = cat(1, landArea_byReg_out, landArea_byReg) ;
            end
            [table_harm, table_harm_relY1] = PLUMharmFigs_iterate_superReg( ...
                sum(PLUMharm_xvyr(:,v,yi_harm,r),2), s, 'harm', ...
                list_superRegs, list_regions, harm_focus_regions, ...
                biomeID_x, countries_x, countries_key, lats, lons, ...
                list2map, landArea_x) ;
            if s==1
                table_harm_out = table_harm ;
                table_harm_relY1_out = table_harm_relY1 ;
            else
                table_harm_out = cat(1, table_harm_out, table_harm) ;
                table_harm_relY1_out = cat(1, table_harm_relY1_out, table_harm_relY1) ;
            end
        end
        
        % Calculate tables relative to land area in region
        table_orig_relLandArea_out = table_orig_out ;
        table_orig_relLandArea_out.Gain = ...
            table_orig_out.Gain ./ landArea_byReg_out ;
        table_orig_relLandArea_out.Loss = ...
            table_orig_out.Loss ./ landArea_byReg_out ;
        table_orig_relLandArea_out.Net = ...
            table_orig_out.Net ./ landArea_byReg_out ;
        table_harm_relLandArea_out = table_harm_out ;
        table_harm_relLandArea_out.Gain = ...
            table_harm_out.Gain ./ landArea_byReg_out ;
        table_harm_relLandArea_out.Loss = ...
            table_harm_out.Loss ./ landArea_byReg_out ;
        table_harm_relLandArea_out.Net = ...
            table_harm_out.Net ./ landArea_byReg_out ;
        
        % Save to Excel files
        thisRun_short = shortened_runList{r} ;
        thisRange = sprintf('A1:F%d', size(table_orig_out, 1) + 1) ;
        sheetName = [strrep(thisRun_short,'.','_') '_o2'] ;
        writetable(table_orig_out, outfile, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        writetable(table_orig_relLandArea_out, outfile_relLandArea, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        sheetName = [strrep(thisRun_short,'.','_') '_h2'] ;
        writetable(table_harm_out, outfile, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        writetable(table_harm_relLandArea_out, outfile_relLandArea, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        sheetName = [strrep(thisRun_short,'.','_') '_o2Ry1'] ;
        writetable(table_orig_relY1_out, outfile, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        writetable(table_orig_relY1_out, outfile_relLandArea, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        sheetName = [strrep(thisRun_short,'.','_') '_h2Ry1'] ;
        writetable(table_harm_relY1_out, outfile, ...
            'Sheet', sheetName, 'Range', thisRange) ;
        writetable(table_harm_relY1_out, outfile_relLandArea, ...
            'Sheet', sheetName, 'Range', thisRange) ;
    end
end
disp('Done.')


%% Time series of harmonization effect on change in non-ag area

% Options %%%%%%%%%
fontSize = 14 ;
lineWidth = 2 ;
thisPos = [1         455        1440         350] ;
%%%%%%%%%%%%%%%%%%%

if ~isequal(yearList_orig, yearList_harm)
    error('This code assumes original and harmonized have same yearList.')
end

ii = strcmp(LUnames, 'NATURAL') ;
PLUMorig_incr_xyr = squeeze(PLUMorig_xvyr(:,ii,2:end,:) - PLUMorig_xvyr(:,ii,1:end-1,:)) ;
PLUMharm_incr_xyr = squeeze(PLUMharm_xvyr(:,ii,2:end,:) - PLUMharm_xvyr(:,ii,1:end-1,:)) ;
PLUMorig_incr_xyr(PLUMorig_incr_xyr<0) = 0 ;
PLUMharm_incr_xyr(PLUMharm_incr_xyr<0) = 0 ;
PLUMorig_incr_yr = squeeze(sum(PLUMorig_incr_xyr,1)) ;
PLUMharm_incr_yr = squeeze(sum(PLUMharm_incr_xyr,1)) ;
harmEffect_yr = PLUMharm_incr_yr - PLUMorig_incr_yr ;

x = yearList_orig(2:end) ;
lms = cell(Nruns,1) ;
for r = 1:Nruns
    lms{r} = fitlm(x, harmEffect_yr(:,r)) ;
end

figure('Color', 'w', 'Position', thisPos) ;
plot(x, harmEffect_yr, ...
    'LineWidth', 1)
set(gca, 'FontSize', fontSize) ;
hold on
set(gca,'ColorOrderIndex',1) ;
for r = 1:Nruns
    plot(x, lms{r}.Fitted, '--', ...
        'LineWidth', 2)
end
hold off
legend(shortened_runList, 'Location', 'best')


%% Map one year, one LC for orig and harm
do_save = true ;

thisLU = 'NATURAL' ;
thisYear = 2010 ;

% Options %%%%%%%%%
fontSize = 14 ;
% spacing = [0.02 0.02] - 0.0025*8 ;   % [vert, horz]
spacing = 0 ;
textX = 0.115 ;
textY_1 = 50/360 ;
textY_2 = 20/360 ;
% shiftup = 0 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ;
shiftup = 15/360 ; textY_1 = textY_1 + shiftup ; textY_2 = textY_2 + shiftup - shiftup/3 ;
thisPos = [1    33   770   772] ;
nx = 2 ;
ny = 4 ;
as_frac_land = true ;
bins_lowBnds = [0:99] ;
conv_fact_total = 1e-6*1e-6 ;   % m2 to Mkm2
do_caps = false ;
this_colormap_name = 'parula' ;
% lines_overlay = 'landareas.shp' ;
lines_overlay = '/Users/sam/Geodata/General/continents_from_countries/continents_from_countries.shp' ;
%%%%%%%%%%%%%%%%%%%

v = contains(LUnames, thisLU) ;

if length(y1_list) > 1
    this_outdir = sprintf('%s/maps_manyDeltas_beforeAfter_%d-%d_by%d', ...
        out_dir, min(y1_list), max(yN_list), yN_list(1)-y1_list(1)+1) ;
    pngres = '-r150' ;
else
    this_outdir = out_dir ;
    pngres = '-r300' ;
end

if ~exist(this_outdir, 'dir')
    mkdir(this_outdir) ;
end

units_map = '%' ;
units_total = 'Mkm^2' ;

map_size = size(landArea_YX) ;
orig_frac_YXrH = nan([map_size Nruns]) ;
harm_frac_YXrH = nan([map_size Nruns]) ;

y1 = y1_list(y) ;
yN = yN_list(y) ;

% Get gridcell fraction that is this LU type (%)
area_orig_xr = squeeze(nansum( ...
    PLUMorig_xvyr(:,v,yearList_orig==thisYear,:), ...
    2)) ;
area_harm_xr = squeeze(nansum( ...
    PLUMharm_xvyr(:,v,yearList_harm==thisYear,:), ...
    2)) ;
for r = 1:Nruns
    orig_frac_YXrH(:,:,r) = lpjgu_vector2map(100*area_orig_xr(:,r)./gcelArea_x, map_size, list2map) ;
    harm_frac_YXrH(:,:,r) = lpjgu_vector2map(100*area_harm_xr(:,r)./gcelArea_x, map_size, list2map) ;
end
area_orig_r = sum(area_orig_xr,1)*conv_fact_total ;
area_harm_r = sum(area_harm_xr,1)*conv_fact_total ;

if ~as_frac_land
    error('This only works with as_frac_land TRUE')
end

col_titles = {sprintf('Original %s, %d', thisLU, thisYear), ...
    sprintf('Harmonized %s, d', thisLU, thisYear)} ;
make_LUfrac_fig_v5(...
    area_orig_r, area_harm_r, ...
    orig_frac_YXrH, harm_frac_YXrH, ...
    runList_legend, ...
    spacing, fontSize, textX, textY_1, textY_2, ...
    nx, ny, ...
    Nruns, thisPos, units_map, units_total, do_caps, ...
    bins_lowBnds, this_colormap_name, col_titles, ...
    lines_overlay) ;

if do_save
    filename = sprintf('%s/maps_%d_beforeAfter.png', ...
        this_outdir, thisYear) ;
    export_fig(filename, pngres) ;
    close
end
    

%% Time series of LUs

ts_base_cy = squeeze(nansum(nansum(base.maps_YXvy,1),2)) ;
ts_base_cy = cat(1,sum(ts_base_cy(isCrop,:),1),ts_base_cy(~isCrop,:)) ;
ts_orig_cyr = squeeze(nansum(PLUMorig_xvyr,1)) ;
ts_orig_cyr = cat(1,sum(ts_orig_cyr(isCrop,:,:),1),ts_orig_cyr(~isCrop,:,:)) ;
ts_harm_cyr = squeeze(nansum(PLUMharm_xvyr,1)) ;
ts_harm_cyr = cat(1,sum(ts_harm_cyr(isCrop,:,:),1),ts_harm_cyr(~isCrop,:,:)) ;
if ~add_baseline_to_harm
    ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;
end

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
if ~add_baseline_to_harm
    ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;
end

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
if ~add_baseline_to_harm
    ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;
end

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
if ~isequal(yearList_orig, yearList_harm)
    ts_harm_cyr = cat(2, ts_harm_cyr(:,1,:)-(ts_orig_cyr(:,2,:)-ts_orig_cyr(:,1,:)), ts_harm_cyr) ;
end

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

thisYear = 2010 ;
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
    tmp_xvr = 100*(tmpH_xvr - tmpO_xvr) ./ repmat(gcelArea_x, [1 Nlu Nruns]) ;
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
        tmp = lpjgu_vector2map(tmp_xvr(:,v,r), map_size, list2map) ;
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

