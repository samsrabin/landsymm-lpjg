%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Separate agricultural input and output into different PLUM types %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inDir_tmp = 'PLUM2LPJGblIrr_SSP1v2_RCP45/output-2017-12-25-140114' ;
y1 = 2011 ;
yN = 2100 ;
yStep = 1 ;


%% Setup

cf_ktNha_kgNm2 = 1e-4 ;
outPrec_LC = 6 ;
% outPrec_mgmtInputs = 3 ;
% outPrec_mgmtInputs = 4 ;
outPrec_mgmtInputs = outPrec_LC ;
force_overwrite = true ;
fclose_every = 1000 ;
pct_progress = 25 ;

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/')

% Get year info
yearList = y1:yStep:yN ;
Nyears = length(yearList) ;

% Get cells present in previous LPJ-GUESS output
lpjg_in = lpjgu_matlab_readTable('/Users/Shared/PLUM/trunk_runs/runs_20170925/LPJG_PLUM_expt1.2_rcp26_ipsl_irrV7104_b20170728_bothDS_irrWaits_forED_20170925154752/2011-2015/yield.out.gz','dont_save_MAT',true,'verbose',false) ;
lons_lpjg = lpjg_in.Lon ;
lats_lpjg = lpjg_in.Lat ;
clear lpjg_in

% Get land uses from previous LPJ-GUESS run
LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
lpj_lu = lpjgu_matlab_readTable(LUfile) ;
lpj_lumap = lpjgu_matlab_readTable_then2map(LUfile) ;


%% Import PLUM outputs / LPJG inputs

inDir = find_PLUM2LPJG_run(inDir_tmp) ;
[~,LUcropDir_tmp] = unix(['/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work/get_lu_dir.sh '...
    ...'/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work']) ;
    inDir]) ;
LUcropDir_tmp = LUcropDir_tmp(1:end-1) ;
LUcropDir_matlab = find_PLUM2LPJG_inputs(LUcropDir_tmp) ;
LUcropDir_orig = find_PLUM2LPJG_inputs(strrep(LUcropDir_tmp,'.forLPJG.MATLAB.20171201','')) ;

%%%%%%%%%%%%%%
%%% Import %%%
%%%%%%%%%%%%%%

for y = 1:Nyears
    thisYear = yearList(y) ;
    disp(['Reading ' num2str(thisYear) '...'])
    
    landcover_in_tmp = lpjgu_matlab_readTable([LUcropDir_orig num2str(thisYear) ...
        '/LandCoverFract.txt'],'dont_save_MAT',true,'verbose',false) ;
    cropland_tmp = landcover_in_tmp.CROPLAND ;
    cropfracs_in_tmp = lpjgu_matlab_readTable([LUcropDir_orig num2str(thisYear) ...
        '/CropFract.txt'],'dont_save_MAT',true,'verbose',false) ;
    nfert_in_tmp = lpjgu_matlab_readTable([LUcropDir_orig num2str(thisYear) ...
        '/Fert.txt'],'dont_save_MAT',true,'verbose',false) ;
    irrig_in_tmp = lpjgu_matlab_readTable([LUcropDir_orig num2str(thisYear) ...
        '/Irrig.txt'],'dont_save_MAT',true,'verbose',false) ;
    detailedLU_in_tmp = lpjgu_matlab_readTable([LUcropDir_orig num2str(thisYear) ...
        '/LandUse.txt'],'dont_save_MAT',true,'verbose',false) ;
    
    % Trim unneeded columns
    nfert_in_tmp.Pasture = [] ;
    irrig_in_tmp.Pasture = [] ;
    
    if y==1
        % Get variable names
        landcover_cols = landcover_in_tmp.Properties.VariableNames ;
        cropfracs_cols = cropfracs_in_tmp.Properties.VariableNames ;
        nfert_cols = nfert_in_tmp.Properties.VariableNames ;
        irrig_cols = irrig_in_tmp.Properties.VariableNames ;
        dtlLU_cols = detailedLU_in_tmp.Properties.VariableNames ;
        
        % Convert crop names
        % % % % %         cropfracs_cols = strrep(cropfracs_cols,'irr','') ;
        cftList = strrep(cropfracs_cols(4:end),'irr','') ;
        Ncfts = length(cftList) ;
        cropfracs_cols = [strrep(cropfracs_cols,'irr','i') strrep(cropfracs_cols(4:end),'irr','')] ;
        nfert_cols = [strrep(nfert_cols,'irr','i') strrep(nfert_cols(4:end),'irr','')] ;
        irrig_cols = [strrep(irrig_cols,'irr','i') strrep(irrig_cols(4:end),'irr','')] ;
        
        % Get # variables
        Nvar_landcover = length(landcover_cols) - 3 ;
        %             Nvar_landcover = Nvar_landcover - 1 ; % Will get rid of URBAN later
        Nvar_cropfracs = length(cropfracs_cols) - 3 ;
        Nvar_nfert = length(nfert_cols) - 3 ;
        Nvar_irrig = length(irrig_cols) - 3 ;
        Nvar_dtlLU = 
    end
    
    % Get lat/lons
    if y==1
        lons = landcover_in_tmp.Lon ;
        lats = landcover_in_tmp.Lat ;
        Ncells = length(lons) ;
        lonlats_randomized = false ;
    end
    
    % Fill CrOp values (i.e., rainfed versions) with 0
    for c = 1:Ncfts
        eval(['cropfracs_in_tmp.' cftList{c} ' = zeros(length(cropfracs_in_tmp.Lon),1) ;']);
        eval(['nfert_in_tmp.' cftList{c} ' = zeros(length(nfert_in_tmp.Lon),1) ;']);
        eval(['irrig_in_tmp.' cftList{c} ' = zeros(length(irrig_in_tmp.Lon),1) ;']);
    end
    
    % Remove rows that sum to zero
    if y==1
        remove_these_cells = find(sum(table2array(landcover_in_tmp(:,4:end)),2)==0) ;
        lons(remove_these_cells) = [] ;
        lats(remove_these_cells) = [] ;
        Ncells = length(lons) ;
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
    
    % Make sure lat/lon cols are equal
    if ~isequal(lons,landcover_in_tmp.Lon) ...
            || ~isequal(lons,cropfracs_in_tmp.Lon) ...
            || ~isequal(lons,nfert_in_tmp.Lon) ...
            || ~isequal(lons,irrig_in_tmp.Lon)
        error('Lons do not match!')
    end
    if ~isequal(lats,landcover_in_tmp.Lat) ...
            || ~isequal(lats,cropfracs_in_tmp.Lat) ...
            || ~isequal(lats,nfert_in_tmp.Lat) ...
            || ~isequal(lats,irrig_in_tmp.Lat)
        error('Lats do not match!')
    end
    
    if y==1
        % Set up empty arrays
        landcover_in_Xvy = nan(Ncells,Nvar_landcover,Nyears) ;
        cropfracs_in_Xvy = nan(Ncells,Nvar_cropfracs,Nyears) ;
        nfert_in_Xvy = nan(Ncells,Nvar_nfert,Nyears) ;
        irrig_in_Xvy = nan(Ncells,Nvar_irrig,Nyears) ;
    end
    
    % Check variable match
    if ~isequal(landcover_cols,landcover_in_tmp.Properties.VariableNames) ...
            || ~isequal(cropfracs_cols,strrep(cropfracs_in_tmp.Properties.VariableNames,'irr','i')) ...
            || ~isequal(nfert_cols,strrep(nfert_in_tmp.Properties.VariableNames,'irr','i')) ...
            || ~isequal(irrig_cols,strrep(irrig_in_tmp.Properties.VariableNames,'irr','i'))
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
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover]) ;
        else
            % Because got rid of URBAN after doing this on y==1
            landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover+1]) ;
        end
        if any(isnan(landcover_in_tmp))
            error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
        end
    end
    cropfracs_in_sum_whereCrop = sum(cropfracs_in_tmp(cropland_tmp>0,:),2) ;
    while any(abs(cropfracs_in_sum_whereCrop - 1) > 10^-(outPrec_LC+2))
        cropfracs_in_tmp(cropland_tmp>0,:) = cropfracs_in_tmp(cropland_tmp>0,:) ./ repmat(cropfracs_in_sum_whereCrop,[1 Nvar_cropfracs]) ;
        if any(abs(cropfracs_in_sum_whereCrop - 1) > 0.01)
            error('Cropfrac difference too big to be comfortable!')
        end
        if any(isnan(cropfracs_in_tmp))
            error('While ensuring sum to 1 of cropfracs_in_tmp, NaN created!')
        end
        cropfracs_in_sum_whereCrop = sum(cropfracs_in_tmp(cropland_tmp>0),2) ;
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
            thisLon = lons(c) ;
            thisLat = lats(c) ;
            plum_list2map(c) = find(lons_map==thisLon & lats_map==thisLat) ;
        end
        barren_lpj = barren_lpj_map(plum_list2map) ;
        vegd_lpj = vegd_lpj_map(plum_list2map) ;
    end
    barren_old = landcover_in_tmp(:,find(strcmp(landcover_cols,'BARREN'))-3) ;
    natural_old = landcover_in_tmp(:,find(strcmp(landcover_cols,'NATURAL'))-3) ;
    urban_old = landcover_in_tmp(:,find(strcmp(landcover_cols,'URBAN'))-3) ;
    vegd_old = sum(landcover_in_tmp(:,find(strcmp(landcover_cols,'NATURAL')|strcmp(landcover_cols,'CROPLAND')|strcmp(landcover_cols,'PASTURE'))-3),2) ;
    mgmt_old = vegd_old - natural_old ;
    [barren_new, natural_new, urban_new] = ...
        BarNatUrb(barren_old, mgmt_old, natural_old, urban_old, ...
        barren_lpj, vegd_lpj) ;
    
    % Save
    landcover_in_tmp(:,find(strcmp(landcover_cols,'BARREN'))-3) = barren_new ;
    landcover_in_tmp(:,find(strcmp(landcover_cols,'NATURAL'))-3) = natural_new ;
    landcover_in_tmp(:,find(strcmp(landcover_cols,'URBAN'))-3) = [] ;
    if y==1
        % Because we got rid of URBAN
        landcover_in_Xvy(:,end,:) = [] ;
        Nvar_landcover = Nvar_landcover - 1 ;
    end
    
    % Again make sure landcover and cropfracs sum to 1
    % (Previously checked to avoid rows that sum to 0)
    while any(abs(sum(landcover_in_tmp,2) - 1) > 10^-(outPrec_LC+2))
        landcover_in_tmp = landcover_in_tmp ./ repmat(sum(landcover_in_tmp,2),[1 Nvar_landcover]) ;
        if any(isnan(landcover_in_tmp))
            error('While ensuring sum to 1 of landcover_in_tmp, NaN created!')
        end
    end
    
    % Save to array
    landcover_in_Xvy(:,:,y) = landcover_in_tmp ;
    cropfracs_in_Xvy(:,:,y) = cropfracs_in_tmp ;
    nfert_in_Xvy(:,:,y) = nfert_in_tmp ;
    irrig_in_Xvy(:,:,y) = irrig_in_tmp ;
    
    clear *tmp
end

% Check for NaNs
if any(isnan(landcover_in_Xvy))
    error('Some value of landcover_in_Xvy is NaN.')
end
if any(isnan(cropfracs_in_Xvy))
    error('Some value of cropfracs_in_Xvy is NaN.')
end
if any(isnan(nfert_in_Xvy))
    error('Some value of nfert_in_Xvy is NaN.')
end
if any(isnan(irrig_in_Xvy))
    error('Some value of irrig_in_Xvy is NaN.')
end

% Check for nonsensical irr values
if min(irrig_in_Xvy)<0
    warning('min(irrig_in_Xvy)<0')
end
if max(irrig_in_Xvy)>1
    warning('max(irrig_in_Xvy)>1')
end

% Change lon/lat from lower-left to center
lons = lons + 0.25 ;
lats = lats + 0.25 ;

% Convert ktN/ha to kgN/m2
nfert_in_Xvy = nfert_in_Xvy*cf_ktNha_kgNm2 ;

disp('Done.')



%%

