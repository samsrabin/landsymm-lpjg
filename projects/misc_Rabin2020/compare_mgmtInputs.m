%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compare management inputs from LPJG4 and PLUM6 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisDir = 'SSP1.v2' ;
y1 = 2011 ;
yN = 2015 ;
yStep = 1 ;


%% Setup

cf_ktNha_kgNm2 = 1e-4 ;
cf_tPha_kgPm2 = 0.1 ;

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/')

inDir = find_PLUM2LPJG_inputs(thisDir) ;
disp(inDir)

% Get year info
yearList = y1:yStep:yN ;
Nyears = length(yearList) ;

% % Get cells missing from climate
% lpjg_in = lpjgu_matlab_readTable('/Users/Shared/PLUM/input/gridlists/missing_climate_searchradius0.5.txt','dont_save_MAT',true,'verbose',false) ;

% Get cells present in previous LPJ-GUESS output
lpjg_in = lpjgu_matlab_readTable('/Users/Shared/PLUM/trunk_runs/runs_20170925/LPJG_PLUM_expt1.2_rcp26_ipsl_irrV7104_b20170728_bothDS_irrWaits_forED_20170925154752/2011-2015/yield.out.gz','dont_save_MAT',true,'verbose',false) ;
lons_lpjg = lpjg_in.Lon ;
lats_lpjg = lpjg_in.Lat ;
clear lpjg_in

% Get land uses from previous LPJ-GUESS run
LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
lpj_lu = lpjgu_matlab_readTable(LUfile) ;
lpj_lumap = lpjgu_matlab_readTable_then2map(LUfile) ;

outPrec_LC = 6 ;
outPrec_mgmtInputs = outPrec_LC ;

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


%% Import and save table from normal

for y = 1:Nyears
    thisYear = yearList(y) ;
    disp(['Reading ' num2str(thisYear) '...'])
    
    landcover_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
        '/LandCoverFract.txt'],'dont_save_MAT',true,'verbose',false) ;
    cropland_tmp = landcover_in_tmp.CROPLAND ;
    cropfracs_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
        '/CropFract.txt'],'dont_save_MAT',true,'verbose',false) ;
    nfert_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
        '/Fert.txt'],'dont_save_MAT',true,'verbose',false) ;
    irrig_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
        '/Irrig.txt'],'dont_save_MAT',true,'verbose',false) ;
    
    % Trim unneeded columns
    nfert_in_tmp.Pasture = [] ;
    irrig_in_tmp.Pasture = [] ;
    
    if y==1
        % Get variable names
        landcover_cols_v1 = landcover_in_tmp.Properties.VariableNames ;
        cropfracs_cols_v1 = cropfracs_in_tmp.Properties.VariableNames ;
        nfert_cols_v1 = nfert_in_tmp.Properties.VariableNames ;
        irrig_cols_v1 = irrig_in_tmp.Properties.VariableNames ;
        
        % Convert crop names
        % % % % %         cropfracs_cols_v1 = strrep(cropfracs_cols_v1,'irr','') ;
        cftList_v1 = strrep(cropfracs_cols_v1(4:end),'irr','') ;
        Ncfts_v1 = length(cftList_v1) ;
%         cropfracs_cols_v1 = [strrep(cropfracs_cols_v1,'irr','i') strrep(cropfracs_cols_v1(4:end),'irr','')] ;
%         nfert_cols_v1 = [strrep(nfert_cols_v1,'irr','i') strrep(nfert_cols_v1(4:end),'irr','')] ;
%         irrig_cols_v1 = [strrep(irrig_cols_v1,'irr','i') strrep(irrig_cols_v1(4:end),'irr','')] ;
        cropfracs_cols_v1 = strrep(cropfracs_cols_v1,'irr','i') ;
        nfert_cols_v1 = strrep(nfert_cols_v1,'irr','i') ;
        irrig_cols_v1 = strrep(irrig_cols_v1,'irr','i') ;
        
        % Get # variables
        Nvar_landcover_v1 = length(landcover_cols_v1) - 3 ;
        %             Nvar_landcover_v1 = Nvar_landcover_v1 - 1 ; % Will get rid of URBAN later
        Nvar_cropfracs_v1 = length(cropfracs_cols_v1) - 3 ;
        Nvar_nfert_v1 = length(nfert_cols_v1) - 3 ;
        Nvar_irrig_v1 = length(irrig_cols_v1) - 3 ;
    end
    
    % Get lat/lons_v1
    if y==1
        lons_v1 = cropfracs_in_tmp.Lon ;
        lats_v1 = cropfracs_in_tmp.Lat ;
        Ncells = length(lons_v1) ;
        lonlats_randomized = false ;
    end
    
%     % Fill CrOp values (i.e., rainfed versions) with 0
%     for c = 1:Ncfts_v1
%         eval(['cropfracs_in_tmp.' cftList_v1{c} ' = zeros(length(cropfracs_in_tmp.Lon),1) ;']);
%         eval(['nfert_in_tmp.' cftList_v1{c} ' = zeros(length(nfert_in_tmp.Lon),1) ;']);
%         eval(['irrig_in_tmp.' cftList_v1{c} ' = zeros(length(irrig_in_tmp.Lon),1) ;']);
%     end
    
    % Remove rows that sum to zero
    if y==1
        remove_these_cells = find(sum(table2array(landcover_in_tmp(:,4:end)),2)==0) ;
        lons_v1(remove_these_cells) = [] ;
        lats_v1(remove_these_cells) = [] ;
        Ncells = length(lons_v1) ;
    end
    landcover_in_tmp(remove_these_cells,:) = [] ;
    cropland_tmp(remove_these_cells) = [] ;
    cropfracs_in_tmp(remove_these_cells,:) = [] ;
    nfert_in_tmp(remove_these_cells,:) = [] ;
    irrig_in_tmp(remove_these_cells,:) = [] ;
    if any(sum(table2array(landcover_in_tmp(:,4:end)),2)==0)
        error('At least one row in landcover_in_tmp sums to zero!')
    end
    cropfracs_in_sum = sum(table2array(cropfracs_in_tmp(:,4:end)),2) ;
    if any(cropfracs_in_sum==0 & cropland_tmp>0)
        arebad = find(cropfracs_in_sum==0 & cropland_tmp>0) ;
        nbad = length(arebad) ;
        max_cropFrac_inBad = max(cropland_tmp(arebad)) ;
        warning(['cropfracs_in_tmp sums to zero in ' num2str(nbad) ' crop-containing row(s)! Setting to zero. (Max crop frac = ' num2str(max_cropFrac_inBad) ')'])
        if max_cropFrac_inBad>1e-4
            error('max_cropFrac_inBad is too large for me to be comfortable with this!')
        end
        landcover_in_tmp.CROPLAND(arebad) = 0 ;
        cropland_tmp = landcover_in_tmp.CROPLAND ;
        clear arebad nbad max_cropFrac_inBad
    end
    
    % Make sure lat/lon cols_v1 are equal
    if ~isequal(lons_v1,landcover_in_tmp.Lon) ...
            || ~isequal(lons_v1,cropfracs_in_tmp.Lon) ...
            || ~isequal(lons_v1,nfert_in_tmp.Lon) ...
            || ~isequal(lons_v1,irrig_in_tmp.Lon)
        error('Lons do not match!')
    end
    if ~isequal(lats_v1,landcover_in_tmp.Lat) ...
            || ~isequal(lats_v1,cropfracs_in_tmp.Lat) ...
            || ~isequal(lats_v1,nfert_in_tmp.Lat) ...
            || ~isequal(lats_v1,irrig_in_tmp.Lat)
        error('Lats do not match!')
    end
    
    if y==1
        % Set up empty arrays
        landcover_in_Xvy_v1 = nan(Ncells,Nvar_landcover_v1,Nyears) ;
        cropfracs_in_Xvy_v1 = nan(Ncells,Nvar_cropfracs_v1,Nyears) ;
        nfert_in_Xvy_v1 = nan(Ncells,Nvar_nfert_v1,Nyears) ;
        irrig_in_Xvy_v1 = nan(Ncells,Nvar_irrig_v1,Nyears) ;
    end
    
    % Check variable match
    if ~isequal(landcover_cols_v1,landcover_in_tmp.Properties.VariableNames) ...
            || ~isequal(cropfracs_cols_v1,strrep(cropfracs_in_tmp.Properties.VariableNames,'irr','i')) ...
            || ~isequal(nfert_cols_v1,strrep(nfert_in_tmp.Properties.VariableNames,'irr','i')) ...
            || ~isequal(irrig_cols_v1,strrep(irrig_in_tmp.Properties.VariableNames,'irr','i'))
        error('Columns do not match!')
    end
    
    % Convert to arrays
    landcover_in_tmp = table2array(landcover_in_tmp(:,4:end)) ;
    cropfracs_in_tmp = table2array(cropfracs_in_tmp(:,4:end)) ;
    nfert_in_tmp = table2array(nfert_in_tmp(:,4:end)) ;
    irrig_in_tmp = table2array(irrig_in_tmp(:,4:end)) ;
    
    % Check for NaNs
    if any(isnan(landcover_in_tmp))
        error('Some value of landcover_in_tmp is NaN (before ensuring sum to 1).')
    end
    if any(isnan(cropfracs_in_tmp))
        error('Some value of cropfracs_in_tmp is NaN (before ensuring sum to 1).')
    end
    if any(isnan(nfert_in_tmp))
        error('Some value of nfert_in_tmp is NaN.')
    end
    if any(isnan(irrig_in_tmp))
        error('Some value of irrig_in_tmp is NaN.')
    end
    
    % Make sure landcover and cropfracs sum to 1
    % (Previously checked to avoid rows that sum to 0)
    while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
        if y==1
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover_v1]) ;
        else
            % Because got rid of URBAN after doing this on y==1
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover_v1+1]) ;
        end
        if any(isnan(landcover_in_tmp))
            error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
        end
    end
    cropfracs_in_sum_whereCrop = sum(cropfracs_in_tmp(cropland_tmp>0,:),2) ;
    while any(abs(cropfracs_in_sum_whereCrop - 1) > 10^-(outPrec_LC+2))
        cropfracs_in_tmp(cropland_tmp>0,:) = cropfracs_in_tmp(cropland_tmp>0,:) ./ repmat(cropfracs_in_sum_whereCrop,[1 Nvar_cropfracs_v1]) ;
        if any(abs(cropfracs_in_sum_whereCrop - 1) > 0.01)
            error('Cropfrac difference too big to be comfortable!')
        end
        if any(isnan(cropfracs_in_tmp))
            error('While ensuring sum to 1 of cropfracs_in_tmp, NaN created!')
        end
        cropfracs_in_sum_whereCrop = sum(cropfracs_in_tmp(cropland_tmp>0,:),2) ;
    end
    
    % Deal with NATURAL/URBAN/BARREN
    if y==1
        barren_lpj_map = lpj_lumap.maps_YXvy(:,:,strcmp(lpj_lumap.varNames,'BARREN'),end) ;
        vegd_lpj_map = sum(lpj_lumap.maps_YXvy(:,:,strcmp(lpj_lumap.varNames,'NATURAL')|strcmp(lpj_lumap.varNames,'CROPLAND')|strcmp(lpj_lumap.varNames,'PASTURE'),end),3) ;
        lons4map = -0.25 + -179.75:0.5:179.75 ;
        lats4map = -0.25 + -89.75:0.5:89.75 ;
        lons_map = repmat(lons4map,[length(lats4map) 1]) ;
        lats_map = repmat(lats4map',[1 length(lons4map)]) ;
        plum_list2map = nan(Ncells,1) ;
        for c = 1:Ncells
            thisLon = lons_v1(c) ;
            thisLat = lats_v1(c) ;
            plum_list2map(c) = find(lons_map==thisLon & lats_map==thisLat) ;
        end
        barren_lpj = barren_lpj_map(plum_list2map) ;
        vegd_lpj = vegd_lpj_map(plum_list2map) ;
    end
    barren_old = landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'BARREN'))-3) ;
    natural_old = landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'NATURAL'))-3) ;
    urban_old = landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'URBAN'))-3) ;
    vegd_old = sum(landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'NATURAL')|strcmp(landcover_cols_v1,'CROPLAND')|strcmp(landcover_cols_v1,'PASTURE'))-3),2) ;
    mgmt_old = vegd_old - natural_old ;
    [barren_new, natural_new, urban_new] = ...
        BarNatUrb(barren_old, mgmt_old, natural_old, urban_old, ...
        barren_lpj, vegd_lpj) ;
    
    % Save
    landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'BARREN'))-3) = barren_new ;
    landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'NATURAL'))-3) = natural_new ;
    landcover_in_tmp(:,find(strcmp(landcover_cols_v1,'URBAN'))-3) = [] ;
    if y==1
        % Because we got rid of URBAN
        landcover_in_Xvy_v1(:,end,:) = [] ;
        Nvar_landcover_v1 = Nvar_landcover_v1 - 1 ;
    end
    
    % Again make sure landcover and cropfracs sum to 1
    % (Previously checked to avoid rows that sum to 0)
    while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
        landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover_v1]) ;
        if any(isnan(landcover_in_tmp))
            error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
        end
    end
    
    % Save to array
    landcover_in_Xvy_v1(:,:,y) = landcover_in_tmp ;
    cropfracs_in_Xvy_v1(:,:,y) = cropfracs_in_tmp ;
    nfert_in_Xvy_v1(:,:,y) = nfert_in_tmp ;
    irrig_in_Xvy_v1(:,:,y) = irrig_in_tmp ;
    
    clear *tmp
end

% Check for NaNs
if any(isnan(landcover_in_Xvy_v1))
    error('Some value of landcover_in_Xvy_v1 is NaN.')
end
if any(isnan(cropfracs_in_Xvy_v1))
    error('Some value of cropfracs_in_Xvy_v1 is NaN.')
end
if any(isnan(nfert_in_Xvy_v1))
    error('Some value of nfert_in_Xvy_v1 is NaN.')
end
if any(isnan(irrig_in_Xvy_v1))
    error('Some value of irrig_in_Xvy_v1 is NaN.')
end

% Check for nonsensical irr values
if min(irrig_in_Xvy_v1)<0
    warning('min(irrig_in_Xvy_v1)<0')
end
if max(irrig_in_Xvy_v1)>1
    warning('max(irrig_in_Xvy_v1)>1')
end

% Change lon/lat from lower-left to center
lons_v1 = lons_v1 + 0.25 ;
lats_v1 = lats_v1 + 0.25 ;

% Convert ktN/ha to kgN/m2
nfert_in_Xvy_v1 = nfert_in_Xvy_v1*cf_ktNha_kgNm2 ;

disp('Done.')


%% Import and save table from detailed

PLUMcrops = {'wheat' ...
    'maize', ...
    'rice' ...
    'oilcrops' ...
    'pulses' ...
    'starchyRoots' ...
    'energycrops' ...
    } ;
LPJGcrops = {'CerealsC3', ...
    'CerealsC4', ...
    'Rice', ...
    'Oilcrops', ...
    'Pulses', ...
    'StarchyRoots', ...
    'Miscanthus'} ;
Ncfts_v2 = length(PLUMcrops) ;
if length(LPJGcrops) ~= Ncfts_v2
    error('length(LPJGcrops) ~= Ncfts_v2') ;
end

for y = 1:Nyears
    thisYear = yearList(y) ;
    disp(['Reading ' num2str(thisYear) '...'])
    
    landcover_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
        '/LandCoverFract.txt'],'dont_save_MAT',true,'verbose',false) ;
    details_in_tmp = lpjgu_matlab_readTable([inDir num2str(thisYear) ...
        '/LandUse.txt'],'dont_save_MAT',true,'verbose',false) ;
    cropland_tmp = landcover_in_tmp.CROPLAND ;
    
    if y==1
        % Get variable names, part 1
        landcover_cols_v2 = landcover_in_tmp.Properties.VariableNames ;
        details_cols_v2 = details_in_tmp.Properties.VariableNames ;
        
        % Get columns indices from detailed table
        inds_PLUMcrops_cropfracs = find(...
            getCellsWithString(details_cols_v2,'_A') ...
            & getCellsWithoutString(details_cols_v2,'ruminants') ...
            & getCellsWithoutString(details_cols_v2,'monogastrics') ...
            & getCellsWithoutString(details_cols_v2,'pasture') ...
            ) ;
        inds_PLUMcrops_nfert = find(...
            getCellsWithString(details_cols_v2,'_FQ') ...
            & getCellsWithoutString(details_cols_v2,'ruminants') ...
            & getCellsWithoutString(details_cols_v2,'monogastrics') ...
            & getCellsWithoutString(details_cols_v2,'pasture') ...
            ) ;
        inds_PLUMcrops_irrig = find(...
            getCellsWithString(details_cols_v2,'_II') ...
            & getCellsWithoutString(details_cols_v2,'ruminants') ...
            & getCellsWithoutString(details_cols_v2,'monogastrics') ...
            & getCellsWithoutString(details_cols_v2,'pasture') ...
            ) ;
        
        % Get temporary tables from detailed table
        cropfracs_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_cropfracs]) ;
        nfert_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_nfert]) ;
        irrig_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_irrig]) ;
        
        % Get variable names, part 2
        cropfracs_cols_v2 = cropfracs_in_tmp.Properties.VariableNames ;
        nfert_cols_v2 = nfert_in_tmp.Properties.VariableNames ;
        irrig_cols_v2 = irrig_in_tmp.Properties.VariableNames ;
        
        % Convert crop names
        %             % % % % %         cropfracs_cols_v2 = strrep(cropfracs_cols_v2,'irr','') ;
        %             cftList_v2 = strrep(cropfracs_cols_v2(3:end),'irr','') ;
        %             Ncfts_v2 = length(cftList_v2) ;
        %             cropfracs_cols_v2 = [strrep(cropfracs_cols_v2,'irr','i') strrep(cropfracs_cols_v2(3:end),'irr','')] ;
        %             nfert_cols_v2 = [strrep(nfert_cols_v2,'irr','i') strrep(nfert_cols_v2(3:end),'irr','')] ;
        %             irrig_cols_v2 = [strrep(irrig_cols_v2,'irr','i') strrep(irrig_cols_v2(3:end),'irr','')] ;
        cftList_v2 = LPJGcrops ;
        for c = 1:Ncfts_v2
            thisCrop_plum = PLUMcrops{c} ;
            thisCrop_lpjg = LPJGcrops{c} ;
            cropfracs_cols_v2(strcmp(cropfracs_cols_v2,[thisCrop_plum '_A'])) = {[thisCrop_lpjg 'i']} ;
            nfert_cols_v2(strcmp(nfert_cols_v2,[thisCrop_plum '_FQ'])) = {[thisCrop_lpjg 'i']} ;
            irrig_cols_v2(strcmp(irrig_cols_v2,[thisCrop_plum '_II'])) = {[thisCrop_lpjg 'i']} ;
        end
%         % Add columns for rainfed
%         for c = 1:Ncfts_v2
%             thisCrop_lpjg = LPJGcrops{c} ;
%             cropfracs_cols_v2(end+1) = {thisCrop_lpjg} ;
%             nfert_cols_v2(end+1) = {thisCrop_lpjg} ;
%             irrig_cols_v2(end+1) = {thisCrop_lpjg} ;
%         end
        
        % Get # variables
        Nvar_landcover_v2 = length(landcover_cols_v2) - 3 ;
        %             Nvar_landcover_v2 = Nvar_landcover_v2 - 1 ; % Will get rid of URBAN later
        Nvar_cropfracs_v2 = length(cropfracs_cols_v2) - 2 ;
        Nvar_nfert_v2 = length(nfert_cols_v2) - 2 ;
        Nvar_irrig_v2 = length(irrig_cols_v2) - 2 ;
        
        % Get lat/lons_v2
        lons_v2 = cropfracs_in_tmp.Lon ;
        lats_v2 = cropfracs_in_tmp.Lat ;
        Ncells = length(lons_v2) ;
        lonlats_randomized = false ;
    else
        % Get temporary tables from detailed table
        cropfracs_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_cropfracs]) ;
        nfert_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_nfert]) ;
        irrig_in_tmp = details_in_tmp(:,[1:2 inds_PLUMcrops_irrig]) ;
    end
    
%     % Fill CrOp values (i.e., rainfed versions) with 0
%     for c = 1:Ncfts_v2
%         eval(['cropfracs_in_tmp.' cftList_v2{c} ' = zeros(length(cropfracs_in_tmp.Lon),1) ;']);
%         eval(['nfert_in_tmp.' cftList_v2{c} ' = zeros(length(nfert_in_tmp.Lon),1) ;']);
%         eval(['irrig_in_tmp.' cftList_v2{c} ' = zeros(length(irrig_in_tmp.Lon),1) ;']);
%     end
    
    % Remove rows that sum to zero
    if y==1
        remove_these_cells = find(sum(table2array(landcover_in_tmp(:,4:end)),2)==0) ;
        lons_v2(remove_these_cells) = [] ;
        lats_v2(remove_these_cells) = [] ;
        Ncells = length(lons_v2) ;
    end
    landcover_in_tmp(remove_these_cells,:) = [] ;
    cropland_tmp(remove_these_cells) = [] ;
    cropfracs_in_tmp(remove_these_cells,:) = [] ;
    nfert_in_tmp(remove_these_cells,:) = [] ;
    irrig_in_tmp(remove_these_cells,:) = [] ;
    if any(sum(table2array(landcover_in_tmp(:,4:end)),2)==0)
        error('At least one row in landcover_in_tmp sums to zero!')
    end
    cropfracs_in_sum = sum(table2array(cropfracs_in_tmp(:,3:end)),2) ;
    if any(cropfracs_in_sum==0 & cropland_tmp>0)
        arebad = find(cropfracs_in_sum==0 & cropland_tmp>0) ;
        nbad = length(arebad) ;
        max_cropFrac_inBad = max(cropland_tmp(arebad)) ;
        warning(['cropfracs_in_tmp sums to zero in ' num2str(nbad) ' crop-containing row(s)! Setting to zero. (Max crop frac = ' num2str(max_cropFrac_inBad) ')'])
        if max_cropFrac_inBad>1e-4
            error('max_cropFrac_inBad is too large for me to be comfortable with this!')
        end
        landcover_in_tmp.CROPLAND(arebad) = 0 ;
        cropland_tmp = landcover_in_tmp.CROPLAND ;
        clear arebad nbad max_cropFrac_inBad
    end
    
    % Make sure lat/lon cols_v2 are equal
    if ~isequal(lons_v2,landcover_in_tmp.Lon) ...
            || ~isequal(lons_v2,cropfracs_in_tmp.Lon) ...
            || ~isequal(lons_v2,nfert_in_tmp.Lon) ...
            || ~isequal(lons_v2,irrig_in_tmp.Lon)
        error('Lons do not match!')
    end
    if ~isequal(lats_v2,landcover_in_tmp.Lat) ...
            || ~isequal(lats_v2,cropfracs_in_tmp.Lat) ...
            || ~isequal(lats_v2,nfert_in_tmp.Lat) ...
            || ~isequal(lats_v2,irrig_in_tmp.Lat)
        error('Lats do not match!')
    end
    
    if y==1
        % Set up empty arrays
        landcover_in_Xvy_v2 = nan(Ncells,Nvar_landcover_v2,Nyears) ;
        cropfracs_in_Xvy_v2 = nan(Ncells,Nvar_cropfracs_v2,Nyears) ;
        nfert_in_Xvy_v2 = nan(Ncells,Nvar_nfert_v2,Nyears) ;
        irrig_in_Xvy_v2 = nan(Ncells,Nvar_irrig_v2,Nyears) ;
    end
    
    %         % Check variable match
    %         if ~isequal(landcover_cols_v2,landcover_in_tmp.Properties.VariableNames) ...
    %                 || ~isequal(cropfracs_cols_v2,strrep(cropfracs_in_tmp.Properties.VariableNames,'irr','i')) ...
    %                 || ~isequal(nfert_cols_v2,strrep(nfert_in_tmp.Properties.VariableNames,'irr','i')) ...
    %                 || ~isequal(irrig_cols_v2,strrep(irrig_in_tmp.Properties.VariableNames,'irr','i'))
    %             error('Columns do not match!')
    %         end
    
    % Convert to arrays
    landcover_in_tmp = table2array(landcover_in_tmp(:,4:end)) ;
    cropfracs_in_tmp = table2array(cropfracs_in_tmp(:,3:end)) ;
    nfert_in_tmp = table2array(nfert_in_tmp(:,3:end)) ;
    irrig_in_tmp = table2array(irrig_in_tmp(:,3:end)) ;
    
    % Check for NaNs
    if any(isnan(landcover_in_tmp))
        error('Some value of landcover_in_tmp is NaN (before ensuring sum to 1).')
    end
    if any(isnan(cropfracs_in_tmp))
        error('Some value of cropfracs_in_tmp is NaN (before ensuring sum to 1).')
    end
    if any(isnan(nfert_in_tmp))
        error('Some value of nfert_in_tmp is NaN.')
    end
    if any(isnan(irrig_in_tmp))
        error('Some value of irrig_in_tmp is NaN.')
    end
    
    % Make sure landcover and cropfracs sum to 1
    % (Previously checked to avoid rows that sum to 0)
    while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
        if y==1
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover_v2]) ;
        else
            % Because got rid of URBAN after doing this on y==1
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover_v2+1]) ;
        end
        if any(isnan(landcover_in_tmp))
            error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
        end
    end
    cropfracs_in_sum_whereCrop = sum(cropfracs_in_tmp(cropland_tmp>0,:),2) ;
    while any(abs(cropfracs_in_sum_whereCrop - 1) > 10^-(outPrec_LC+2))
        cropfracs_in_tmp(cropland_tmp>0,:) = cropfracs_in_tmp(cropland_tmp>0,:) ./ repmat(cropfracs_in_sum_whereCrop,[1 Nvar_cropfracs_v2]) ;
        if any(abs(cropfracs_in_sum_whereCrop - 1) > 0.01)
            error('Cropfrac difference too big to be comfortable!')
        end
        if any(isnan(cropfracs_in_tmp))
            error('While ensuring sum to 1 of cropfracs_in_tmp, NaN created!')
        end
        cropfracs_in_sum_whereCrop = sum(cropfracs_in_tmp(cropland_tmp>0,:),2) ;
    end
    
    % Deal with NATURAL/URBAN/BARREN
    if y==1
        barren_lpj_map = lpj_lumap.maps_YXvy(:,:,strcmp(lpj_lumap.varNames,'BARREN'),end) ;
        vegd_lpj_map = sum(lpj_lumap.maps_YXvy(:,:,strcmp(lpj_lumap.varNames,'NATURAL')|strcmp(lpj_lumap.varNames,'CROPLAND')|strcmp(lpj_lumap.varNames,'PASTURE'),end),3) ;
        lons4map = -0.25 + -179.75:0.5:179.75 ;
        lats4map = -0.25 + -89.75:0.5:89.75 ;
        lons_map = repmat(lons4map,[length(lats4map) 1]) ;
        lats_map = repmat(lats4map',[1 length(lons4map)]) ;
        plum_list2map = nan(Ncells,1) ;
        for c = 1:Ncells
            thisLon = lons_v2(c) ;
            thisLat = lats_v2(c) ;
            plum_list2map(c) = find(lons_map==thisLon & lats_map==thisLat) ;
        end
        barren_lpj = barren_lpj_map(plum_list2map) ;
        vegd_lpj = vegd_lpj_map(plum_list2map) ;
    end
    barren_old = landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'BARREN'))-3) ;
    natural_old = landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'NATURAL'))-3) ;
    urban_old = landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'URBAN'))-3) ;
    vegd_old = sum(landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'NATURAL')|strcmp(landcover_cols_v2,'CROPLAND')|strcmp(landcover_cols_v2,'PASTURE'))-3),2) ;
    mgmt_old = vegd_old - natural_old ;
    [barren_new, natural_new, urban_new] = ...
        BarNatUrb(barren_old, mgmt_old, natural_old, urban_old, ...
        barren_lpj, vegd_lpj) ;
    
    % Save
    landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'BARREN'))-3) = barren_new ;
    landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'NATURAL'))-3) = natural_new ;
    landcover_in_tmp(:,find(strcmp(landcover_cols_v2,'URBAN'))-3) = [] ;
    if y==1
        % Because we got rid of URBAN
        landcover_in_Xvy_v2(:,end,:) = [] ;
        Nvar_landcover_v2 = Nvar_landcover_v2 - 1 ;
    end
    
    % Again make sure landcover and cropfracs sum to 1
    % (Previously checked to avoid rows that sum to 0)
    while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
        landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover_v2]) ;
        if any(isnan(landcover_in_tmp))
            error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
        end
    end
    
    % Save to array
    landcover_in_Xvy_v2(:,:,y) = landcover_in_tmp ;
    cropfracs_in_Xvy_v2(:,:,y) = cropfracs_in_tmp ;
    nfert_in_Xvy_v2(:,:,y) = nfert_in_tmp ;
    irrig_in_Xvy_v2(:,:,y) = irrig_in_tmp ;
    
    clear *tmp
end

% Check for NaNs
if any(isnan(landcover_in_Xvy_v2))
    error('Some value of landcover_in_Xvy_v2 is NaN.')
end
if any(isnan(cropfracs_in_Xvy_v2))
    error('Some value of cropfracs_in_Xvy_v2 is NaN.')
end
if any(isnan(nfert_in_Xvy_v2))
    error('Some value of nfert_in_Xvy_v2 is NaN.')
end
if any(isnan(irrig_in_Xvy_v2))
    error('Some value of irrig_in_Xvy_v2 is NaN.')
end

% Check for nonsensical values
if min(nfert_in_Xvy_v2)<0
    warning('min(nfert_in_Xvy_v2)<0')
end
if min(irrig_in_Xvy_v2)<0
    warning('min(irrig_in_Xvy_v2)<0')
end
if max(irrig_in_Xvy_v2)>1
    warning('max(irrig_in_Xvy_v2)>1')
end

% Change lon/lat from lower-left to center
lons_v2 = lons_v2 + 0.25 ;
lats_v2 = lats_v2 + 0.25 ;

% Convert ktN/ha to kgN/m2
nfert_in_Xvy_v2 = nfert_in_Xvy_v2*cf_ktNha_kgNm2 ;

disp('Done.')


%% Get land area array

ind2map = nan(size(lons_v1)) ;
for i = 1:Ncells
    ind2map(i) = find(lons_map+0.25==lons_v1(i) & lats_map+0.25==lats_v1(i)) ;
end
landarea_m2_x_y = repmat(land_area_YX_m2(ind2map),[1 1 Nyears]) ;
% landarea_m2_x4y = repmat(landarea_m2_x_y,[1 4 1]) ;
% landarea_m2_x7y = repmat(landarea_m2_x_y,[1 7 1]) ;


%% Compare Nfert: TrRii and Ricei

% Options %%%%%%%%%%%%
thisVar_v1 = 'TrRii' ;
thisVar_v2 = 'Ricei' ;
% thisVar_v1 = 'TeWWi' ;
% thisVar_v2 = 'CerealsC3i' ;
fontSize = 14 ;
%%%%%%%%%%%%%%%%%%%%%%

v1 = squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v1(:,find(strcmp(landcover_cols_v1,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v1(:,find(strcmp(cropfracs_cols_v1,thisVar_v1))-3,:) ...
    .* nfert_in_Xvy_v1(:,find(strcmp(nfert_cols_v1,thisVar_v1))-3,:) ...
    , 1)) ;

v2 = squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v2(:,find(strcmp(landcover_cols_v2,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v2(:,find(strcmp(cropfracs_cols_v2,thisVar_v2))-2,:) ...
    .* nfert_in_Xvy_v2(:,find(strcmp(nfert_cols_v2,thisVar_v2))-2,:) ...
    , 1)) ;

figure('Color','w') ;
plot(yearList,v1,'-b','LineWidth',3)
hold on
plot(yearList,v2,'--r','LineWidth',3)
hold off
legend(thisVar_v1,thisVar_v2)
set(gca,'FontSize',fontSize,'XTick',yearList)
ylabel('kg N')


%% Compare crop area of CerealsC3/TeWW+TeSW

% Options %%%%%%%%%%%%
fontSize = 18 ;
%%%%%%%%%%%%%%%%%%%%%%

m2_to_Mha = 1e-4*1e-6 ;

area_TeWW = m2_to_Mha * squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v1(:,find(strcmp(landcover_cols_v1,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v1(:,find(strcmp(cropfracs_cols_v1,'TeWWi'))-3,:) ...
    , 1)) ;
area_TeSW = m2_to_Mha * squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v1(:,find(strcmp(landcover_cols_v1,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v1(:,find(strcmp(cropfracs_cols_v1,'TeSWi'))-3,:) ...
    , 1)) ;
area_CerealsC3 = m2_to_Mha * squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v2(:,find(strcmp(landcover_cols_v2,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v2(:,find(strcmp(cropfracs_cols_v2,'CerealsC3i'))-2,:) ...
    , 1)) ;
area_Oilcrops = m2_to_Mha * squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v2(:,find(strcmp(landcover_cols_v2,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v2(:,find(strcmp(cropfracs_cols_v2,'Oilcropsi'))-2,:) ...
    , 1)) ;
area_Pulses = m2_to_Mha * squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v2(:,find(strcmp(landcover_cols_v2,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v2(:,find(strcmp(cropfracs_cols_v2,'Pulsesi'))-2,:) ...
    , 1)) ;
area_StarchyRoots = m2_to_Mha * squeeze(sum( ...
landarea_m2_x_y .* ...
    landcover_in_Xvy_v2(:,find(strcmp(landcover_cols_v2,'CROPLAND'))-3,:) ...
    .* cropfracs_in_Xvy_v2(:,find(strcmp(cropfracs_cols_v2,'StarchyRootsi'))-2,:) ...
    , 1)) ;
area_TeSWcereal = area_TeSW - (area_Oilcrops+area_Pulses+area_StarchyRoots) ;

figure('Color','w') ;
plot(yearList,area_TeWW+area_TeSWcereal,'-k','LineWidth',3)
hold on
plot(yearList,[area_TeWW area_TeSW],'-','LineWidth',3)
plot(yearList,[area_CerealsC3 area_Oilcrops area_Pulses area_StarchyRoots],'--','LineWidth',3)
plot(yearList,area_TeSWcereal,':','LineWidth',3)
hold off
legend('TeWW + (TeSW - [Oilcrops+Pulses+StarchyRoots])','TeWW','TeSW','CerealsC3','Oilcrops','Pulses','StarchyRoots','(TeSW - [Oilcrops+Pulses+StarchyRoots])',...
       'Location','SouthOutside')
set(gca,'FontSize',fontSize,'XTick',yearList)
ylabel('Mha')
title('Crop area')

