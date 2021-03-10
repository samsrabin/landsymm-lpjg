%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LUH1-style harmonization for PLUM outputs, at cropType level %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input directories and settings specific thereto %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Typical runs (2010-2100), v1-3
% dirList = {...
%                   'SSP1.v10.s1' ;
%                   'SSP3.v10.s1' ;
%                   'SSP4.v10.s1' ;
%                   'SSP5.v10.s1' ;
%                   'SSP1.v11.s1' ;
%                   'SSP3.v11.s1' ;
%                   'SSP4.v11.s1' ;
%                   'SSP5.v11.s1' ;
%                   'SSP1.v11.s1test' ;
%                   'ssp11/SSP5/s3' ;
%                   'SSP1.v12.s1' ;
%                   'SSP3.v12.s1' ;
%                   'SSP4.v12.s1' ;
%                   'SSP5.v12.s1' ;
%                   } ;
% base_year = 2010 ;
% year1 = 2011 ;
% yearN = 2100 ;
% fruitveg_sugar_2oil = false ;
% % baseline_ver = 1 ;
% % baseline_ver = 2 ;   % Based on remap_v6
% baseline_ver = 3 ;   % Based on remap_v6p7

%%% Typical runs (2010-2100), v4
dirList = {...
%     'SSP1/s1' ;
    'SSP2/s1' ;
    'SSP3/s1' ;
    'SSP4/s1' ;
    'SSP5/s1' ;
    } ;
base_year = 2010 ;
year1 = 2011 ;
yearN = 2100 ;
fruitveg_sugar_2oil = false ;
baseline_ver = 4 ;   % Based on remap_v8b(2oil, if fruitveg_sugar_2oil true)

%%% Half-Earth runs (2010-2060)
% dirList = {...
%                   'halfearth/HEoct/baseline/s1';
%                   'halfearth/HEoct/halfearth/s1';
%                   } ;
% base_year = 2010 ;
% year1 = 2011 ;
% yearN = 2060 ;
% fruitveg_sugar_2oil = false ;
% baseline_ver = 4 ;   % Based on remap_v8b(2oil, if fruitveg_sugar_2oil true)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General behavior options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save?
save_halfDeg_mat = true ;
save_2deg_mat = false ;
save_halfDeg_txt = false ;
save_2deg_txt = false ;

% Debugging outputs?
debug_areas = false ;
debug_nfert = false ;
debug_irrig = false ;

% Coordinates of 2-degree cell to debug (leave empty for no debug)
debugIJ_2deg = [] ;

% Land use of interest
dbCrop = '' ;

% Print verbose messages?
verbose = false ;

% Combine all crops? This is useful only as a test to see how different the
% area is this way as opposed to the right way (keeping every crop
% separate). This also may not actually work!
combineCrops = false ;


%% Setup

warning('on','all')

% Determine which system you're on and set up
thisSystem = get_system_name() ;
if strcmp(thisSystem, 'ssr_mac')
    addpath(genpath('/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/MATLAB_work')) ;
    plumharm_repo_path = '/Users/sam/Documents/Dropbox/2016_KIT/LandSyMM/plum_harmonization/' ;
elseif strcmp(thisSystem, 'ssr_keal')
    addpath(genpath('/pd/data/lpj/sam/paper02-matlab-work')) ;
    plumharm_repo_path = '/pd/data/lpj/sam/PLUM/plum_harmonization' ;
else
    error('thisSystem not recognized: %s', thisSystem)
end
addpath(genpath(plumharm_repo_path))

% Method for inpaint_nans()
% (moved to PLUMharm_importRefData)

% Use "latest PLUM management" in cells that don't have thisCrop thisYear
% but did in a previous year? FALSE = rely solely on interpolation of
% thisYear's thisMgmt.
useLatestPLUMmgmt = true ;

% Save details
save_halfDeg_any = save_halfDeg_mat || save_halfDeg_txt ;
save_2deg_any = save_2deg_mat || save_2deg_txt ;
save_any = save_halfDeg_any || save_2deg_any ;

% Output file details
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;

% How much divergence from PLUM-original transition is acceptable? (%)
% Also used for checking whether out_total > max_mgmt_thisRing.
conserv_tol_pct = 0.2 ;
% How much divergence from PLUM-original transition is acceptable? (m2)
% Also used for checking whether agricultural area exceeds vegetated area
% in each gridcell.
conserv_tol_area = 1e3 ;

% The fraction of real crops (by area) not included in PLUM's 6
% non-energyCrops commodities. Multiply PLUM area by (1-norm2extra) to
% get area of the crop itself without the norm2extra fraction included.
norm2extra = 0.177 ;

% Debug?
do_debug = ~isempty(debugIJ_2deg) && (debug_areas || debug_nfert || debug_irrig) ;
clear db2i db2j
if do_debug
    db2i = debugIJ_2deg(1) ;
    db2j = debugIJ_2deg(2) ;
end
debug_header = ':\t\tArea\t\tFrac\n' ;

if combineCrops
    warning('Combining all crops into one! Will skip harmonization of nfert and irrig.')
end


%% Import reference data

doHarm = true ;
PLUMharm_importRefData

% If specified crop name instead of index, find index.
if ~isempty(dbCrop) && ischar(dbCrop)
    dbCrop = find(strcmp(LPJGcrops,dbCrop)) ;
end


%% Do it

warning('on','all')

% The years we want to produce PLUM outputs for (will begin with
% transitions from years(1)-1 to years(1)
years = year1:yearN ;
Nyears = length(years) ;
if year1 <= base_year
    error('year1 (%d) must be > base_year (%d)!\n', year1, base_year)
elseif year1 > yearN
    error('year1 (%d) must be <= yearN (%d)!\n', year1, yearN)
end

for d = 1:length(dirList)

    %%%%%%%%%%%%%
    %%% Setup %%% 
    %%%%%%%%%%%%%

    % Get directories
    inDir = removeslashifneeded(dirList{d}) ;
    if ~exist(inDir, 'dir')
        error('inDir %s not found. Try changing MATLAB working directory to inDir''s parent.')
    end
    disp(inDir)
    outDir = [inDir '.harm'] ;
    outDir = get_harm_dir(outDir, fruitveg_sugar_2oil, combineCrops) ;
    outDir = addslashifneeded(outDir) ;
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end
    inDir = addslashifneeded(inDir) ;
    
    % Save diary file
    diary off
    diaryfile = sprintf('%s/matlab_log.txt', outDir) ;
    if exist(diaryfile, 'file')
        delete(diaryfile)
    end
    diary(diaryfile)
    diary on

    % Read fraction of VEGETATED protected by...
    %%% PLUM's minimum natural fraction ("rate")
    [~,r] = unix(['grep "MIN_NATURAL_RATE" ' inDir 'config.properties | sed "s/MIN_NATURAL_RATE=//"  | tr -d ''\n''']) ;
    min_natural_rate = str2double(r) ;
    resFrac_minN_YX = min_natural_rate * ones(size(resFrac_prot_YX)) ;
    %%% Total
    resFrac_YX = 1 - (1-resFrac_terr_YX).*(1-resFrac_prot_YX)*(1-min_natural_rate) ;
    % Convert to area
    resArea_YX = resFrac_YX .* base_vegd_YX ;
    resArea_2deg_YX = PLUMharm_aggregate(resArea_YX,0.5,2) ;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Process PLUM outputs %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Debugging output
    if do_debug
        fprintf('Center Lon %.2f, Lat %.2f, land area %0.4g\n\n',1+lons_map_2deg(db2i,db2j),1+lats_map_2deg(db2i,db2j),round(landArea_2deg_YX(db2i,db2j),4)) ;
        if year1==base_year+1
            if debug_areas && ~isempty(debugIJ_2deg)
                PLUMharm_debugOut('LUH2','areas',base_2deg.maps_YXv,landArea_2deg_YX,debugIJ_2deg,LUnames)
            end
            if debug_nfert && ~isempty(debugIJ_2deg)
                PLUMharm_debugOut('LUH2','Nfert',base_nfertTot_2deg.maps_YXv,base_2deg.maps_YXv(:,:,isCrop),debugIJ_2deg,LPJGcrops)
            end
            if debug_irrig && ~isempty(debugIJ_2deg)
                PLUMharm_debugOut('LUH2','irrig',base_irrigTot_2deg.maps_YXv,base_2deg.maps_YXv(:,:,isCrop),debugIJ_2deg,LPJGcrops)
            end
        end
    end

    bareFrac_y0_YX = base_bareFrac_YX ;
    if useLatestPLUMmgmt
        latestPLUMin_2deg_nfert_YXv = -1*ones([size(landArea_2deg_YX) Ncrops_lpjg]) ;
        latestPLUMin_2deg_irrig_YXv = -1*ones([size(landArea_2deg_YX) Ncrops_lpjg]) ;
    else
        latestPLUMin_2deg_nfert_YXv = [] ; %#ok<*UNRCH>
        latestPLUMin_2deg_irrig_YXv = [] ;
    end
    clear out_* in_*
    
    isFirstYearInLoop = true ;
    
    for y = 1:Nyears
        
        thisYear = years(y) ;
        disp(num2str(thisYear)) ;
        tic ;
                
        % Import previous harmonized year, if needed
        if thisYear-1 == base_year
            out_y0 = base ;
            if do_debug
                disp('out_y0 from luh2')
            end
            out_y0_2deg = base_2deg ;
            out_y0_vegd_YX = sum(base.maps_YXv(:,:,notBare),3) ;
            out_y0_2deg_vegd_YX = sum(base_2deg.maps_YXv(:,:,notBare),3) ;
            out_y0_agri_YXv = base.maps_YXv(:,:,isAgri) ;
            out_y0_2deg_agri_YXv = base_2deg.maps_YXv(:,:,isAgri) ;
            if ~combineCrops
                out_y0_nfert_YXv = base_nfert.maps_YXv ;
                out_y0_irrig_YXv = base_irrig.maps_YXv ;
                out_y0_2deg_nfert = base_2deg_nfert ;
                out_y0_2deg_irrig = base_2deg_irrig ;
                out_y0_2deg_nfert_YXv = base_2deg_nfert.maps_YXv ;
                out_y0_2deg_irrig_YXv = base_2deg_irrig.maps_YXv ;
            end
        elseif y==1
            file_in = [removeslashifneeded(inDir) '.harm/' num2str(thisYear-1) 'post.base' num2str(base_year) '.mat'] ;
            % Load previous MAT-file
            if do_debug
                disp(['*y0* from ' file_in])
            end
            load(file_in) ;
            if combineCrops
                clear *nfert* *irrig*
            else
                % Rename nfert vars, if necessary
                badVars = whos('*nfert_2deg*') ;
                if ~isempty(badVars)
                    for v = 1:length(badVars)
                        thisVar = badVars(v).name ;
                        newVar = strrep(thisVar,'nfert_2deg','2deg_nfert') ;
                        eval([newVar ' = ' thisVar ' ;']) ;
                        eval(['clear ' thisVar]) ;
                    end
                end
                % Rename irrig vars, if necessary
                badVars = whos('*irrig_2deg*') ;
                if ~isempty(badVars)
                    for v = 1:length(badVars)
                        thisVar = badVars(v).name ;
                        newVar = strrep(thisVar,'irrig_2deg','2deg_irrig') ;
                        eval([newVar ' = ' thisVar ' ;']) ;
                        eval(['clear ' thisVar]) ;
                    end
                end
            end
        end
        
        % Import and process previous year, if needed
        if ~exist('in_y0','var')
            file_in = [inDir num2str(thisYear-1) '/LandCoverFract.txt'] ;
            if do_debug
                disp(['in_y0 from ' file_in])
            end
            if combineCrops
                [in_y0, ~, ~, in_y0_2deg] = ...
                    PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                    LUnames, bareFrac_y0_YX, [], [], ...
                    PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, ...
                    fruitveg_sugar_2oil) ;
            else
                [in_y0, in_y0_nfert, in_y0_irrig, in_y0_2deg, in_y0_2deg_nfert, in_y0_2deg_irrig, ...
                    latestPLUMin_2deg_nfert_YXv, latestPLUMin_2deg_irrig_YXv, ...
                    max_orig_nfert_y0] = ...
                    PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                    LUnames, bareFrac_y0_YX, latestPLUMin_2deg_nfert_YXv, latestPLUMin_2deg_irrig_YXv, ...
                    PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, ...
                    fruitveg_sugar_2oil) ;
            end
            bareFrac_y0_YX = in_y0.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ./ landArea_YX ;
            in_y0_agri_YXv = in_y0.maps_YXv(:,:,isAgri) ;
            in_y0_2deg_agri_YXv = in_y0_2deg.maps_YXv(:,:,isAgri) ;
            in_y0_2deg_vegd_YX = sum(in_y0_2deg.maps_YXv(:,:,notBare),3) ;
            % Check for bad values
            PLUMharm_checkBadVals(in_y0.maps_YXv, [], [], landArea_YX, LUnames, 'in_y0') ;
        end
        
        % Check for extreme discrepancies and warn.
        % These are not necessarily problems! We expect SOME discrepancy.
        % Printing these warnings lets us do things like check whether the
        % baseline is way too high (PLUM discrepancy -100%). A very large
        % positive PLUM discrepancy can foreshadow an eventual drop to zero
        % after harmonization.
        if thisYear-1 == base_year
            for v = 1:length(in_y0.varNames)
                thisVar = in_y0.varNames{v} ;
                
                % Areas
                tmp_in_area_YX = in_y0.maps_YXv(:,:,v) ;
                tmp_in = squeeze(nansum(nansum(tmp_in_area_YX,2),1)) ;
                tmp_out_area_YX = out_y0.maps_YXv(:,:,v) ;
                tmp_out = squeeze(nansum(nansum(tmp_out_area_YX,2),1)) ;
                check_discrepancy(tmp_in, tmp_out, base_year, thisVar, 'area')
                
                % Managements
                if ~isCrop(v)
                    continue
                end
                tmp_in = squeeze(nansum(nansum( ...
                    tmp_in_area_YX .* ...
                    in_y0_nfert.maps_YXv(:,:,strcmp(LPJGcrops, thisVar)),1),2)) ;
                tmp_out = squeeze(nansum(nansum( ...
                    tmp_out_area_YX .* ...
                    out_y0_nfert_YXv(:,:,strcmp(LPJGcrops, thisVar)),1),2)) ;
                check_discrepancy(tmp_in, tmp_out, base_year, thisVar, 'fert')
                tmp_in = squeeze(nansum(nansum( ...
                    tmp_in_area_YX .* ...
                    in_y0_irrig.maps_YXv(:,:,strcmp(LPJGcrops, thisVar)),1),2)) ;
                tmp_out = squeeze(nansum(nansum( ...
                    tmp_out_area_YX .* ...
                    out_y0_irrig_YXv(:,:,strcmp(LPJGcrops, thisVar)),1),2)) ;
                check_discrepancy(tmp_in, tmp_out, base_year, thisVar, 'irrig')                
            end 
        end

        % Import this year and convert to area
        file_in = [inDir num2str(thisYear) '/LandCoverFract.txt'] ;
        if do_debug
            disp(['in_y1 from ' file_in])
        end
        if combineCrops
            [in_y1, ~, ~, in_y1_2deg] = ...
                PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                LUnames, bareFrac_y0_YX, [], [], ...
                PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, ...
                    fruitveg_sugar_2oil) ;
        else
            [in_y1, in_y1_nfert, in_y1_irrig, in_y1_2deg, in_y1_2deg_nfert, in_y1_2deg_irrig, ...
                latestPLUMin_2deg_nfert_YXv, latestPLUMin_2deg_irrig_YXv, ...
                max_orig_nfert_y1] = ...
                PLUMharm_processPLUMin_areaCrops(file_in, landArea_YX, landArea_2deg_YX, ...
                LUnames, bareFrac_y0_YX, latestPLUMin_2deg_nfert_YXv, latestPLUMin_2deg_irrig_YXv, ...
                PLUMtoLPJG, LPJGcrops, norm2extra, inpaint_method, ...
                    fruitveg_sugar_2oil) ;
        end
        in_y1_agri_YXv = in_y1.maps_YXv(:,:,isAgri) ;
        in_y1_2deg_agri_YXv = in_y1_2deg.maps_YXv(:,:,isAgri) ;
        in_y1_ntrl_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'NATURAL')) ;
        in_y1_bare_YX = in_y1.maps_YXv(:,:,strcmp(in_y1.varNames,'BARREN')) ;
        in_y1_vegd_YX = sum(in_y1.maps_YXv(:,:,notBare),3) ;
        in_y1_vegdFrac_YX = in_y1_vegd_YX ./ landArea_YX ;
        in_y1_bareFrac_YX = in_y1_bare_YX ./ landArea_YX ;
        in_y1_2deg_vegd_YX = sum(in_y1_2deg.maps_YXv(:,:,notBare),3) ;
        
        % Check for bad values
        PLUMharm_checkBadVals(in_y1.maps_YXv, [], [], landArea_YX, LUnames, 'in_y1') ;
        
        % Debugging
        if debug_areas
            debug_global_areas(in_y0_2deg.maps_YXv, in_y1_2deg.maps_YXv, ...
                'Initial import', 'in_2deg', 'in_2deg', ...
                LUnames, isCrop, isAgri, dbCrop, thisYear)
        end
        if ~combineCrops
            if debug_nfert
                debug_global_mgmts(...
                    in_y0_2deg_nfert.maps_YXv, in_y1_2deg_nfert.maps_YXv, 0, ...
                    in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                    'Mt', 1e-6*1e-3, ...
                    'Initial import', 'in_2deg', 'in_2deg', 'nfert', ...
                    LPJGcrops, dbCrop, thisYear)
            end
            if debug_irrig
                debug_global_mgmts(...
                    in_y0_2deg_irrig.maps_YXv, in_y1_2deg_irrig.maps_YXv, 0, ...
                    in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                    'arb', 1, ...
                    'Initial import', 'in_2deg', 'in_2deg', 'irrig', ...
                    LPJGcrops, dbCrop, thisYear)
            end
        end
        
        % Set up empty arrays for tracking global totals
        if isFirstYearInLoop
            out_ts.yearList = [(years(1)-1) years] ;
            out_ts.luNames = in_y0.varNames ;
            out_ts.cropNames = LPJGcrops ;
            out_ts.area_vy = nan(Nlu, Nyears+1) ;
            out_ts.nfert_vy = nan(length(find(isCrop)), Nyears+1) ;
            out_ts.irrig_vy = nan(length(find(isCrop)), Nyears+1) ;
            in_ts = out_ts ;
        end
        
        % Make sure that total global loss of a land use does not exceed
        % the total global area of that land use.
        % NOTE: v10.SSP4 has Miscanthus in 2010, and some gridcells lose
        % Miscanthus from 2010 to 2011. This cannot be satisfied using my
        % 2010 maps because I have no Miscanthus. However, it shouldn't be
        % a problem as long as the grid cells PLUM specifies as losing
        % Miscanthus lose ALL of their Miscanthus.
        setAreaToZero = false(size(LUnames)) ;
        for i = 1:length(LUnames)
            thisLU = LUnames{i} ;
            out0_thisLU_YX = out_y0.maps_YXv(:,:,i) ;
            in0_thisLU_YX = in_y0.maps_YXv(:,:,i) ;
            in1_thisLU_YX = in_y1.maps_YXv(:,:,i) ;
            inDiff_thisLU_YX = in1_thisLU_YX - in0_thisLU_YX ;
            globLoss_thisLU = sum(inDiff_thisLU_YX(inDiff_thisLU_YX<0)) ;
            globArea_thisLU = sum(out0_thisLU_YX(:)) ;
            NETglobChg_thisLU = sum(inDiff_thisLU_YX(:)) ;
            if -globLoss_thisLU > globArea_thisLU
                if do_debug
                    disp(LUnames{i})
                    fprintf('globArea_y0_thisLu = %0.4e\n',sum(in0_thisLU_YX(:))) ;
                    fprintf('globArea_y1_thisLu = %0.4e\n',sum(in1_thisLU_YX(:))) ;
                    fprintf('NETglobChg_thisLU  = %0.4e\n',NETglobChg_thisLU) ;
                    fprintf('globLoss_thisLU    = %0.4e\n',globLoss_thisLU) ;
                    fprintf('globArea_thisLU    = %0.4e\n',globArea_thisLU) ;
                end
                if strcmp(thisLU,'Miscanthus') ...
                && thisYear-1==base_year
                    if sum(in1_thisLU_YX(in0_thisLU_YX>0)) == 0
                        warning('PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. However, it shouldn''t be a problem because the grid cells PLUM specifies as losing Miscanthus lose ALL it.')
                    else
                        warning('PLUM has Miscanthus in 2010, and some gridcells lose Miscanthus from 2010 to 2011. This cannot be satisfied using my 2010 maps because I have no Miscanthus. But not all grid cells PLUM specifies as losing Miscanthus lose ALL of it! Difference of %0.1g m2. Ignoring.',sum(in1_thisLU_YX(in0_thisLU_YX>0)))
                    end
                    tmp = in_y0.maps_YXv(:,:,i) ;
                    tmp(inDiff_thisLU_YX<0) = in1_thisLU_YX(inDiff_thisLU_YX<0) ;
                    in_y0.maps_YXv(:,:,i) = tmp ;
                    in_y0_2deg.maps_YXv = PLUMharm_aggregate(in_y0.maps_YXv,0.5,2) ;
                    in_y0_agri_YXv = in_y0.maps_YXv(:,:,isAgri) ;
                    in_y0_2deg_agri_YXv = in_y0_2deg.maps_YXv(:,:,isAgri) ;
                elseif 0 <= NETglobChg_thisLU && NETglobChg_thisLU < conserv_tol_area
                    warning('Total global loss of %s (%0.4e m2) exceeds its y0 area (%0.4e m2). LIKELY WILL cause infinite loop in ringRedist.\n', ...
                        thisLU, -globLoss_thisLU, globArea_thisLU) ;
                elseif globArea_thisLU + sum(inDiff_thisLU_YX(inDiff_thisLU_YX>0)) > -globLoss_thisLU
                    warning('Total global loss of %s (%0.4e m2) exceeds its y0 area (%0.4e m2), but not y0 area + total global gain (%0.4e m2). (((May))) cause infinite loop in ringRedist.\n', ...
                        thisLU, -globLoss_thisLU, globArea_thisLU, globArea_thisLU + sum(inDiff_thisLU_YX(inDiff_thisLU_YX>0))) ;
                else
                    warning('Total global loss of %s (%0.4e m2) exceeds its y0 area (%0.4e m2) AND ALSO y0 area + total global gain (%0.4e m2). Setting area to zero.\n', ...
                        thisLU, -globLoss_thisLU, globArea_thisLU, globArea_thisLU + sum(inDiff_thisLU_YX(inDiff_thisLU_YX>0))) ;
                    setAreaToZero(i) = true ;
                end
                pause(0.1)
            end
        end; clear *_thisLU*
        
        % Check for bad values
        PLUMharm_checkBadVals(in_y0.maps_YXv, [], [], landArea_YX, LUnames, 'in_y0.2') ;
        PLUMharm_checkBadVals(in_y1.maps_YXv, [], [], landArea_YX, LUnames, 'in_y1.2') ;
        
        % Debugging output
        if debug_areas && ~isempty(debugIJ_2deg)
            PLUMharm_debugOut_deltas('iny0_to_iny1','areas',in_y0_2deg.maps_YXv,in_y1_2deg.maps_YXv,debugIJ_2deg,LUnames)
        end
        
        % Get maximum allowed N fertilization (maximum seen for this crop
        % in any PLUM output thus far, or in LUH2, or in any harmonized
        % output thus far [although no harmonized output should exceed max
        % anyway)
        if ~combineCrops
            max_harm_nfert_y0 = squeeze(max(max(out_y0_2deg_nfert.maps_YXv,[],1),[],2)) ;
            if ~exist('max_nfert_y0','var')
                max_nfert_y1 = max([max_orig_nfert_y0 max_orig_nfert_y1 max_harm_nfert_y0],[],2) ;
            else
                max_nfert_y1 = max([max_orig_nfert_y0 max_orig_nfert_y1 max_harm_nfert_y0 max_nfert_y0],[],2) ;
            end
            max_nfert_y1_YXv = repmat(permute(max_nfert_y1,[2 3 1]),[size(landArea_2deg_YX) 1]) ;
            tooMuchNfert_YXv = in_y0_2deg_nfert.maps_YXv > max_nfert_y1_YXv ;
            
            % Do not allow invalid management inputs. This should only affect
            % cells that did not have any of this crop in the original PLUM
            % outputs.
            in_y0_2deg_nfert.maps_YXv(in_y0_2deg_nfert.maps_YXv<0) = 0 ;
            in_y1_2deg_nfert.maps_YXv(in_y1_2deg_nfert.maps_YXv<0) = 0 ;
            in_y0_2deg_nfert.maps_YXv(tooMuchNfert_YXv) = max_nfert_y1_YXv(tooMuchNfert_YXv) ;
            in_y1_2deg_nfert.maps_YXv(tooMuchNfert_YXv) = max_nfert_y1_YXv(tooMuchNfert_YXv) ;
            in_y0_2deg_irrig.maps_YXv(in_y0_2deg_irrig.maps_YXv<0) = 0 ;
            in_y1_2deg_irrig.maps_YXv(in_y1_2deg_irrig.maps_YXv<0) = 0 ;
            in_y0_2deg_irrig.maps_YXv(in_y0_2deg_irrig.maps_YXv>1) = 1 ;
            in_y1_2deg_irrig.maps_YXv(in_y1_2deg_irrig.maps_YXv>1) = 1 ;
            
            % Check for bad values
            PLUMharm_checkBadVals([], in_y0_2deg_nfert.maps_YXv, in_y0_2deg_irrig.maps_YXv, [], LUnames, 'in_y0_2deg') ;
            PLUMharm_checkBadVals([], in_y1_2deg_nfert.maps_YXv, in_y1_2deg_irrig.maps_YXv, [], LUnames, 'in_y1_2deg') ;
            
            % Debugging
            if debug_nfert
                debug_global_mgmts(...
                    in_y0_2deg_nfert.maps_YXv, in_y1_2deg_nfert.maps_YXv, 0, ...
                    in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                    'Mt', 1e-6*1e-3, ...
                    'After limiting to [0, max]', 'in_2deg', 'in_2deg', 'nfert', ...
                    LPJGcrops, dbCrop, thisYear)
                if ~isempty(debugIJ_2deg)
                    PLUMharm_debugOut_deltas('iny0_to_iny1','Nfert', ...
                        in_y0_2deg_nfert.maps_YXv .* in_y0_2deg.maps_YXv(:,:,isCrop), ...
                        in_y1_2deg_nfert.maps_YXv .* in_y1_2deg.maps_YXv(:,:,isCrop), ...
                        debugIJ_2deg, LPJGcrops)
                end
            end
            if debug_irrig
                debug_global_mgmts(...
                    in_y0_2deg_irrig.maps_YXv, in_y1_2deg_irrig.maps_YXv, 0, ...
                    in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                    'arb', 1, ...
                    'After limiting to [0, max]', 'in_2deg', 'in_2deg', 'irrig', ...
                    LPJGcrops, dbCrop, thisYear)
                if ~isempty(debugIJ_2deg)
                    PLUMharm_debugOut_deltas('iny0_to_iny1','irrig', ...
                        in_y0_2deg_irrig.maps_YXv .* in_y0_2deg.maps_YXv(:,:,isCrop), ...
                        in_y1_2deg_irrig.maps_YXv .* in_y1_2deg.maps_YXv(:,:,isCrop), ...
                        debugIJ_2deg, LPJGcrops)
                end
            end
        end
        
        % Calculate changes in PLUM agri area grids at 2 degrees
        %%% Negative indicates LOSS of thisAgri area or mgmt
        agri_d_YXv = in_y1_2deg_agri_YXv - in_y0_2deg_agri_YXv ;
        for i = find(setAreaToZero) % Skips if all false
            agri_d_YXv(:,:,i) = -out_y0_2deg_agri_YXv(:,:,i) ;
        end
        
        % Apply area changes to previous grid (@2deg)
        mid1_y1_2deg_agri_YXv = out_y0_2deg_agri_YXv + agri_d_YXv ;
        mid1_y1_2deg_ntrl_YX = out_y0_2deg.maps_YXv(:,:,strcmp(out_y0_2deg.varNames,'NATURAL')) ...
                             - sum(agri_d_YXv,3) ;
        
        % Debugging output
        if debug_areas && ~isempty(debugIJ_2deg)
            PLUMharm_debugOut_deltas('outy0_to_mid1y1', 'areas',...
                out_y0_2deg.maps_YXv, ...
                cat(3,mid1_y1_2deg_agri_YXv,mid1_y1_2deg_ntrl_YX,base_2deg_bare_YX), ...
                debugIJ_2deg,LUnames) ;
        end

        % Get area available for conversion to agriculture,
        % as max(NATURAL-RESERVED,0)
        out_y0_2deg_ntrl_YX = out_y0_2deg.maps_YXv(:,:,strcmp(out_y0_2deg.varNames,'NATURAL')) ;
        nonResNtrl_YX = max(out_y0_2deg_ntrl_YX - resArea_2deg_YX, 0) ;

        % Check that neither crop nor pasture exceeds nonBare area;
        % check that neither crop nor pasture is negative;
        % check that nonBare area is conserved;
        % compute the total amount of crop or pasture increase / decrease that
        % is not able to be met within the 2 degree gridcells
        %%% Negative total_unmet indicates TOO MUCH LOSS of thisAgri area
        if debug_areas
            thisDBij = debugIJ_2deg ;
        else
            thisDBij = [] ;
        end
        [unm_y1_2deg_agri_YXv, ...
            mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX] = ...
            PLUMharm_getUnmet_cropAreaRes(...
            mid1_y1_2deg_agri_YXv, out_y0_2deg_vegd_YX, ...
            resArea_2deg_YX, sum(out_y0_2deg_agri_YXv,3), thisDBij, dbCrop) ;
        mid_y1_2deg_vegd_YX = sum(mid_y1_2deg_agri_YXv,3) + mid_y1_2deg_ntrl_YX ;
        mid_y1_2deg_bare_YX = landArea_2deg_YX - mid_y1_2deg_vegd_YX ;
                        
        % Rounding errors can result in small negative values. Fix.
        % Check for bad values
        tmp_YXv = cat(3, mid_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX, mid_y1_2deg_bare_YX) ;
        tmp_YXv = PLUMharm_fixTinyNegs(tmp_YXv, repmat(landArea_2deg_YX,[1 1 Nlu])) ;
        PLUMharm_checkBadVals(tmp_YXv, [], [], landArea_2deg_YX, LUnames, 'mid_y1_2deg') ;
        mid_y1_2deg_agri_YXv = tmp_YXv(:,:,1:end-2) ;
        mid_y1_2deg_ntrl_YX = tmp_YXv(:,:,end-1) ;
        mid_y1_2deg_bare_YX = tmp_YXv(:,:,end) ;
        clear tmp_YXv
        
        % Check 2: Check that global area changes are (mostly) conserved
        PLUMharm_checkCons_area(...
            out_y0_2deg_agri_YXv, mid_y1_2deg_agri_YXv, ...
            in_y0_2deg_agri_YXv, in_y1_2deg_agri_YXv, ...
            unm_y1_2deg_agri_YXv, LUnames(isAgri), ...
            conserv_tol_pct, conserv_tol_area, '2') ;
        
        % Check 2: Check that cells' vegetated area does not change much
        nonCons_pct_YX = 100*(mid_y1_2deg_vegd_YX - out_y0_2deg_vegd_YX) ./ out_y0_2deg_vegd_YX ;
        if max(abs(nonCons_pct_YX(:))) > conserv_tol_pct
            error('By-cell vegetated area not conserved! Worst errors %0.1e and %0.1e.\n', ...
                min(nonCons_pct_YX(:)), max(nonCons_pct_YX(:))) ;
        end
        
        % Debugging
        if debug_areas
            debug_global_areas(out_y0_2deg.maps_YXv, ...
                cat(3,mid_y1_2deg_agri_YXv+unm_y1_2deg_agri_YXv, mid_y1_2deg_ntrl_YX, mid_y1_2deg_bare_YX), ...
                'After addition of deltas', 'out_2deg', 'mid_2deg', ...
                LUnames, isCrop, isAgri, dbCrop, thisYear)
            if ~isempty(debugIJ_2deg)
                PLUMharm_debugOut_deltas('outy0_to_midy1', 'areas', ...
                    out_y0_2deg.maps_YXv, ...
                    cat(3,mid_y1_2deg_agri_YXv,mid_y1_2deg_ntrl_YX,mid_y1_2deg_bare_YX), ...
                    debugIJ_2deg,LUnames)
            end
        end
        
        % Loop through every 2-degree gridcell. If a gridcell has unmet crop
        % or pasture, look for place to put this unmet amount in neighboring
        % rings, starting with gridcells that are 1 unit away, then 2, etc.
        % until all unmet has been displaced to new 2 degree cells. Track the
        % displaced crop and pasture and the # of "rings" needed for each
        % 2-degree gridcell.
        [out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX] = ...
            PLUMharm_ringRedist_areaCropsRes(...
            mid_y1_2deg_agri_YXv, ...
            unm_y1_2deg_agri_YXv, ...
            debugIJ_2deg, conserv_tol_pct, '2b areas', ...
            in_y0_2deg_agri_YXv, in_y1_2deg_agri_YXv, ...
            out_y0_2deg_agri_YXv, ...
            mid_y1_2deg_ntrl_YX, resArea_2deg_YX, LUnames_agri) ;
        out_y1_2deg_vegd_YX = sum(out_y1_2deg_agri_YXv,3) + out_y1_2deg_ntrl_YX ;
        out_y1_2deg_bare_YX = landArea_2deg_YX - out_y1_2deg_vegd_YX ;
        
        % Rounding errors can result in small negative values. Fix.
        % Check for bad values
        tmp_YXv = cat(3, out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX, out_y1_2deg_bare_YX) ;
        tmp_YXv = PLUMharm_fixTinyNegs(tmp_YXv, repmat(landArea_2deg_YX,[1 1 Nlu])) ;
        PLUMharm_checkBadVals(tmp_YXv, [] ,[], landArea_2deg_YX, LUnames, 'out_y1_2deg') ;
        out_y1_2deg_agri_YXv = tmp_YXv(:,:,1:end-2) ;
        out_y1_2deg_ntrl_YX = tmp_YXv(:,:,end-1) ;
        out_y1_2deg_bare_YX = tmp_YXv(:,:,end) ;
        clear tmp_YXv
        
        % Debugging
        if debug_areas
            debug_global_areas(out_y0_2deg.maps_YXv, ...
                cat(3,out_y1_2deg_agri_YXv,out_y1_2deg_ntrl_YX,out_y1_2deg_bare_YX), ...
                'After ringRedist', 'out_2deg', 'out_2deg', ...
                LUnames, isCrop, isAgri, dbCrop, thisYear)
            if ~isempty(debugIJ_2deg)
                PLUMharm_debugOut_deltas('outy0_to_outy1', 'areas', ...
                    out_y0_2deg.maps_YXv, ...
                    cat(3,out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX, out_y1_2deg_bare_YX), ...
                    debugIJ_2deg,LUnames) ;
            end
        end
        
        % Check 3: Check that global area changes are (mostly) conserved
        PLUMharm_checkCons_area(...
            out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, ...
            in_y0_2deg_agri_YXv, in_y1_2deg_agri_YXv, ...
            zeros(size(out_y0_2deg_agri_YXv)), LUnames(isAgri), ...
            conserv_tol_pct, conserv_tol_area, ...
            '3') ;
        
        % Apply management deltas.
        if ~combineCrops
            % Make sure that managements are within acceptable bounds,
            % considering (a) how much each cell had available to lose and (b)
            % how much room each cell has to gain before exceeding max rate.
            if debug_nfert ; thisDBij = debugIJ_2deg ;
            else ; thisDBij = [] ;
            end
            [unm_y1_2deg_nfert_YXv, mid_y1_2deg_nfert_YXv] = ...
                PLUMharm_getUnmet_mgmt(...
                out_y0_2deg_nfert.maps_YXv, out_y0_2deg.maps_YXv(:,:,isCrop), ...
                in_y0_2deg_nfert.maps_YXv, in_y1_2deg_nfert.maps_YXv, ...
                in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), max_nfert_y1, thisDBij, 'Nfert', dbCrop) ;
            if debug_irrig ; thisDBij = debugIJ_2deg ;
            else ; thisDBij = [] ;
            end
            [unm_y1_2deg_irrig_YXv, mid_y1_2deg_irrig_YXv] = ...
                PLUMharm_getUnmet_mgmt(...
                out_y0_2deg_irrig.maps_YXv, out_y0_2deg.maps_YXv(:,:,isCrop), ...
                in_y0_2deg_irrig.maps_YXv, in_y1_2deg_irrig.maps_YXv, ...
                in_y0_2deg.maps_YXv(:,:,isCrop), in_y1_2deg.maps_YXv(:,:,isCrop), ...
                out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ones(size(max_nfert_y1)), thisDBij, 'irrig', dbCrop) ;
            
            % Do not allow invalid management inputs.
            PLUMharm_checkBadVals([], mid_y1_2deg_nfert_YXv, mid_y1_2deg_irrig_YXv, [], LUnames, 'mid_y1_2deg') ;
            
            % Debugging
            if debug_nfert
                debug_global_mgmts(...
                    out_y0_2deg_nfert.maps_YXv, mid_y1_2deg_nfert_YXv, unm_y1_2deg_nfert_YXv, ...
                    out_y0_2deg_agri_YXv(:,:,1:end-1), out_y1_2deg_agri_YXv(:,:,1:end-1), ...
                    'Mt', 1e-6*1e-3, ...
                    'After applying deltas', 'out_2deg', 'mid_2deg', 'nfert', ...
                    LPJGcrops, dbCrop, thisYear)
                if ~isempty(debugIJ_2deg)
                    PLUMharm_debugOut_deltas('outy0_to_midy1', 'Nfert', ...
                        out_y0_2deg_nfert_YXv.*out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                        mid_y1_2deg_nfert_YXv.*out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                        debugIJ_2deg,LPJGcrops) ;
                end
            end
            if debug_irrig
                debug_global_mgmts(...
                    out_y0_2deg_irrig.maps_YXv, mid_y1_2deg_irrig_YXv, unm_y1_2deg_irrig_YXv, ...
                    out_y0_2deg_agri_YXv(:,:,1:end-1), out_y1_2deg_agri_YXv(:,:,1:end-1), ...
                    'arb', 1, ...
                    'After applying deltas', 'out_2deg', 'mid_2deg', 'irrig', ...
                    LPJGcrops, dbCrop, thisYear)
                if ~isempty(debugIJ_2deg)
                    PLUMharm_debugOut_deltas('outy0_to_midy1', 'irrig', ...
                        out_y0_2deg_irrig_YXv.*out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                        mid_y1_2deg_irrig_YXv.*out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                        debugIJ_2deg,LPJGcrops) ;
                end
            end
            
            % Check 2b: Check that management changes are (mostly) conserved
            PLUMharm_checkCons_mgmt(...
                out_y0_2deg_nfert.maps_YXv, out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                mid_y1_2deg_nfert_YXv, out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_nfert.maps_YXv, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_nfert.maps_YXv, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                unm_y1_2deg_nfert_YXv, LPJGcrops, conserv_tol_pct, false(Ncrops_lpjg,1), ...
                '2b nfert', true) ;
            PLUMharm_checkCons_mgmt(...
                out_y0_2deg_irrig.maps_YXv, out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                mid_y1_2deg_irrig_YXv, out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_irrig.maps_YXv, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_irrig.maps_YXv, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                unm_y1_2deg_irrig_YXv, LPJGcrops, conserv_tol_pct, false(Ncrops_lpjg,1), ...
                '2b irrig', true) ;
            
            % Do ring redistribution for management inputs
            if debug_nfert ; thisDBij = debugIJ_2deg ;
            else ; thisDBij = [] ;
            end
            [out_y1_2deg_nfert_YXv, unm2_y1_2deg_nfert_YXv, notEnough_nfert] = ...
                PLUMharm_ringRedist_mgmt(...
                mid_y1_2deg_nfert_YXv, out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                unm_y1_2deg_nfert_YXv, max_nfert_y1, ...
                LPJGcrops, thisDBij, ...
                out_y0_2deg_nfert, out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_nfert, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_nfert, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                conserv_tol_pct, 'ringRedist nfert', dbCrop) ;
            if debug_irrig ; thisDBij = debugIJ_2deg ;
            else ; thisDBij = [] ;
            end
            
            [out_y1_2deg_irrig_YXv, unm2_y1_2deg_irrig_YXv, notEnough_irrig] = ...
                PLUMharm_ringRedist_mgmt(...
                mid_y1_2deg_irrig_YXv, out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                unm_y1_2deg_irrig_YXv, ones(size(max_nfert_y1)), ...
                LPJGcrops, thisDBij, ...
                out_y0_2deg_irrig, out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_irrig, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_irrig, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                conserv_tol_pct, 'ringRedist irrig', dbCrop) ;
            
            % Do not allow invalid management inputs.
            PLUMharm_checkBadVals([], out_y1_2deg_nfert_YXv, out_y1_2deg_irrig_YXv, [], LUnames, 'out_y1_2deg') ;
            
            % Debugging
            if debug_nfert
                debug_global_mgmts(...
                    out_y0_2deg_nfert.maps_YXv, out_y1_2deg_nfert_YXv, 0, ...
                    out_y0_2deg_agri_YXv(:,:,1:end-1), out_y1_2deg_agri_YXv(:,:,1:end-1), ...
                    'Mt', 1e-6*1e-3, ...
                    'After ringRedist', 'out_2deg', 'out_2deg', 'nfert', ...
                    LPJGcrops, dbCrop, thisYear)
            end
            if debug_irrig
                debug_global_mgmts(...
                    out_y0_2deg_irrig.maps_YXv, out_y1_2deg_irrig_YXv, 0, ...
                    out_y0_2deg_agri_YXv(:,:,1:end-1), out_y1_2deg_agri_YXv(:,:,1:end-1), ...
                    'arb', 1, ...
                    'After ringRedist', 'out_2deg', 'out_2deg', 'irrig', ...
                    LPJGcrops, dbCrop, thisYear)
            end
            
            % Check 3b: Check that management changes are (mostly) conserved
            PLUMharm_checkCons_mgmt(...
                out_y0_2deg_nfert.maps_YXv, out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                out_y1_2deg_nfert_YXv, out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_nfert.maps_YXv, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_nfert.maps_YXv, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                unm2_y1_2deg_nfert_YXv, LPJGcrops, conserv_tol_pct, notEnough_nfert, ...
                '3b nfert', true) ;
            
            PLUMharm_checkCons_mgmt(...
                out_y0_2deg_irrig.maps_YXv, out_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                out_y1_2deg_irrig_YXv, out_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_irrig.maps_YXv, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_irrig.maps_YXv, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                unm2_y1_2deg_irrig_YXv, LPJGcrops, conserv_tol_pct, notEnough_irrig, ...
                '3b irrig', true) ;
        end
        
        % Loop through all 2-degree gridcells and distribute the new crop and
        % pasture area changes to the half-degree gridcells within.
        %
        % If the crop or pasture change is a decrease, apply the same
        % PERCENTAGE decrease to all half-degree crop or pasture gridcells 
        % within the 2-degree cell.
        %
        % If the crop or pasture change is an increase, apply this increase to 
        % all half-degree gridcells within the 2-degree cell proportionally 
        % according to available land area. Make sure that total area within a 
        % half-degree gridcell does not exceed 1 or 0.
        if do_debug
            debug_tmp = [debugIJ_2deg dbCrop] ;
        else
            debug_tmp = [] ;
        end
        out_y1_agri_YXv = ...
            PLUMharm_distDeltas_areaCrops_recursive( ...
            landArea_YX, landArea_2deg_YX, ...
            out_y0_2deg_agri_YXv, out_y1_2deg_agri_YXv, out_y0_agri_YXv, ...
            out_y0_vegd_YX, conserv_tol_pct, conserv_tol_area, LUnames, ...
            debug_tmp) ; clear debug_tmp
        out_y1_ntrl_YX = out_y0_vegd_YX - sum(out_y1_agri_YXv,3) ;
        out_y1_vegd_YX = sum(out_y1_agri_YXv,3) + out_y1_ntrl_YX ;
        out_y1_bare_YX = landArea_YX - out_y1_vegd_YX ;
        
        % Debugging
        if debug_areas
            debug_global_areas(out_y0_2deg.maps_YXv, ...
                cat(3,out_y1_agri_YXv, out_y1_ntrl_YX, out_y1_bare_YX), ...
                'Now at half-degree but before negative value fix', 'out', 'out', ...
                LUnames, isCrop, isAgri, dbCrop, thisYear)
        end

        % Rounding errors can result in small negative values. Fix.
        % Check for bad values.
        tmp_YXv = cat(3, out_y1_agri_YXv, out_y1_ntrl_YX, out_y1_bare_YX) ;
        tmp_YXv = PLUMharm_fixTinyNegs(tmp_YXv, repmat(landArea_YX,[1 1 Nlu])) ;
        PLUMharm_checkBadVals(tmp_YXv, [] ,[], landArea_YX, LUnames, 'out_y1_2deg') ;
        out_y1_agri_YXv = tmp_YXv(:,:,1:end-2) ;
        out_y1_ntrl_YX = tmp_YXv(:,:,end-1) ;
        out_y1_bare_YX = tmp_YXv(:,:,end) ;        
        
        % Debugging
        if debug_areas
            debug_global_areas(out_y0_2deg.maps_YXv, ...
                cat(3,out_y1_agri_YXv, out_y1_ntrl_YX, out_y1_bare_YX), ...
                'At half-degree after negative value fix', 'out', 'out', ...
                LUnames, isCrop, isAgri, dbCrop, thisYear)
        end        
                
        % Check 5: Check that global area changes are (mostly) conserved
        PLUMharm_checkCons_area(...
            out_y0_agri_YXv, out_y1_agri_YXv, ...
            in_y0_agri_YXv, in_y1_agri_YXv, ...
            zeros(size(out_y0_agri_YXv)), LUnames(isAgri), ...
            conserv_tol_pct, conserv_tol_area, ...
            '5') ;
        
        if ~combineCrops
            % Distribute management inputs from 2-degree to half-degree cells.
            out_y1_nfert_YXv = PLUMharm_distMgmt(out_y1_2deg_nfert_YXv,2,0.5) ;
            out_y1_irrig_YXv = PLUMharm_distMgmt(out_y1_2deg_irrig_YXv,2,0.5) ;
            out_y1_nfert_YXv(out_y1_agri_YXv(:,:,~isAgri_isPast)==0) = 0 ;
            out_y1_irrig_YXv(out_y1_agri_YXv(:,:,~isAgri_isPast)==0) = 0 ;
            
            % If no thisCrop, no mgmt additions to thisCrop
            out_y1_nfert_YXv(out_y1_agri_YXv(:,:,~isAgri_isPast)==0) = 0 ;
            out_y1_irrig_YXv(out_y1_agri_YXv(:,:,~isAgri_isPast)==0) = 0 ;
            
            % Do not allow invalid management inputs.
            PLUMharm_checkBadVals([], out_y1_nfert_YXv, out_y1_irrig_YXv, [], LUnames, 'out_y1') ;
            
            % Debugging
            if debug_nfert
                debug_global_mgmts(...
                    out_y0_nfert_YXv, out_y1_nfert_YXv, 0, ...
                    out_y0_agri_YXv(:,:,1:end-1), out_y1_agri_YXv(:,:,1:end-1), ...
                    'Mt', 1e-6*1e-3, ...
                    'Now at half-degree', 'out', 'out', 'nfert', ...
                    LPJGcrops, dbCrop, thisYear)
            end
            if debug_irrig
                debug_global_mgmts(...
                    out_y0_irrig_YXv, out_y1_irrig_YXv, 0, ...
                    out_y0_agri_YXv(:,:,1:end-1), out_y1_agri_YXv(:,:,1:end-1), ...
                    'arb', 1, ...
                    'Now at half-degree', 'out', 'out', 'irrig', ...
                    LPJGcrops, dbCrop, thisYear)
            end
            
            % Check 5: Check that global mgmt changes are (mostly) conserved
            PLUMharm_checkCons_mgmt(...
                out_y0_nfert_YXv, out_y0_agri_YXv(:,:,~isAgri_isPast), ...
                out_y1_nfert_YXv, out_y1_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_nfert.maps_YXv, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_nfert.maps_YXv, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
...                 in_y0_nfert.maps_YXv, in_y0_agri_YXv(:,:,~isAgri_isPast), ...
...                 in_y1_nfert.maps_YXv, in_y1_agri_YXv(:,:,~isAgri_isPast), ...
                zeros(size(out_y0_nfert_YXv)), LPJGcrops, conserv_tol_pct, notEnough_nfert, ...
                '5 nfert', false) ;
            PLUMharm_checkCons_mgmt(...
                out_y0_irrig_YXv, out_y0_agri_YXv(:,:,~isAgri_isPast), ...
                out_y1_irrig_YXv, out_y1_agri_YXv(:,:,~isAgri_isPast), ...
                in_y0_2deg_irrig.maps_YXv, in_y0_2deg_agri_YXv(:,:,~isAgri_isPast), ...
                in_y1_2deg_irrig.maps_YXv, in_y1_2deg_agri_YXv(:,:,~isAgri_isPast), ...
...                 in_y0_irrig.maps_YXv, in_y0_agri_YXv(:,:,~isAgri_isPast), ...
...                 in_y1_irrig.maps_YXv, in_y1_agri_YXv(:,:,~isAgri_isPast), ...
                zeros(size(out_y0_irrig_YXv)), LPJGcrops, conserv_tol_pct, notEnough_irrig, ...
                '5 irrig', false) ;
            
            % Make sure that ExtraCrop receives no management
            if max(max(out_y1_nfert_YXv(:,:,strcmp(LPJGcrops,'ExtraCrop'))))>0
                error('Some ExtraCrop Nfert >0! (%0.3e)',max(max(out_y1_nfert_YXv(:,:,strcmp(LPJGcrops,'ExtraCrop')))))
            elseif max(max(out_y1_irrig_YXv(:,:,strcmp(LPJGcrops,'ExtraCrop'))))>0
                error('Some ExtraCrop irrig >0! (%0.3e)',max(max(out_y1_irrig_YXv(:,:,strcmp(LPJGcrops,'ExtraCrop')))))
            end
        end

        % Get land use areas
        out_y1_past_YX = out_y1_agri_YXv(:,:,isAgri_isPast) ;
        out_y1_crop_YX = sum(out_y1_agri_YXv(:,:,~isAgri_isPast),3) ;
        out_y1_2deg_past_YX = out_y1_2deg_agri_YXv(:,:,isAgri_isPast) ;
        out_y1_2deg_crop_YX = sum(out_y1_2deg_agri_YXv(:,:,~isAgri_isPast),3) ;
        if ~exist('out_y1_2deg','var')
            out_y1_2deg.varNames = out_y0_2deg.varNames ;
        end
        out_y1_2deg.maps_YXv = cat(3, out_y1_2deg_agri_YXv, out_y1_2deg_ntrl_YX, out_y1_2deg_bare_YX) ;

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Write output files %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if save_any
            disp(['  Done processing (' toc_hms(toc) '). Now writing.'])
        end
        
        outDir_thisYear = sprintf('%s%d', outDir, thisYear) ;
        if save_halfDeg_any
            if ~exist(outDir_thisYear, 'dir')
                mkdir(outDir_thisYear)
            end
            % Save new LandCoverFract.txt (0.5-degree)
            out_y1.varNames = {'PASTURE','CROPLAND','NATURAL','BARREN'} ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'CROPLAND')) = out_y1_crop_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'PASTURE')) = out_y1_past_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'NATURAL')) = out_y1_ntrl_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'BARREN')) = out_y1_bare_YX ;
            out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_YX,[1 1 length(out_y1.varNames)]) ;
            
            if save_halfDeg_mat
                file_out = [outDir_thisYear '/LandCoverFract.base' num2str(base_year) '.mat'] ;
                save(file_out,'out_y1') ;
            end
            if save_halfDeg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                file_out = [outDir_thisYear '/LandCoverFract.base' num2str(base_year) '.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1
            
            % Save new CropFract.txt (0.5-degree)
            if ~combineCrops
                out_y1.maps_YXv = out_y1_agri_YXv(:,:,~isAgri_isPast) ./ repmat(out_y1_crop_YX,[1 1 length(LPJGcrops)]) ;
                out_y1.maps_YXv(repmat(out_y1_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
                out_y1.varNames = LPJGcrops ;
                if save_halfDeg_mat
                    file_out = [outDir_thisYear '/CropFract.base' num2str(base_year) '.mat'] ;
                    save(file_out,'out_y1') ;
                end
                if save_halfDeg_txt
                    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                    file_out = [outDir_thisYear '/CropFract.base' num2str(base_year) '.txt'] ;
                    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                        'outPrec', outPrec, ...
                        'outWidth', outWidth, ...
                        'delimiter', delimiter, ...
                        'overwrite', overwrite, ...
                        'fancy', fancy, ...
                        'progress_step_pct', 20, ...
                        'verbose', false) ;
                end
                clear out_y1
            end
            
            % Save new Fert.txt (0.5-degree)
            if ~combineCrops
                % Convert from kg/m2 to kg/ha for compatibility with original
                % PLUM style.
                out_y1.maps_YXv = 1e4*out_y1_nfert_YXv ;
                out_y1.varNames = LPJGcrops ;
                if save_halfDeg_mat
                    file_out = [outDir_thisYear '/Fert.base' num2str(base_year) '.mat'] ;
                    save(file_out,'out_y1') ;
                end
                if save_halfDeg_txt
                    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                    file_out = [outDir_thisYear '/Fert.base' num2str(base_year) '.txt'] ;
                    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                        'outPrec', outPrec, ...
                        'outWidth', outWidth, ...
                        'delimiter', delimiter, ...
                        'overwrite', overwrite, ...
                        'fancy', fancy, ...
                        'progress_step_pct', 20, ...
                        'verbose', false) ;
                end
                clear out_y1
            end
            
            % Save new Irrig.txt (0.5-degree)
            if ~combineCrops
                out_y1.maps_YXv = out_y1_irrig_YXv ;
                out_y1.varNames = LPJGcrops ;
                if save_halfDeg_mat
                    file_out = [outDir_thisYear '/Irrig.base' num2str(base_year) '.mat'] ;
                    save(file_out,'out_y1') ;
                end
                if save_halfDeg_txt
                    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map) ;
                    file_out = [outDir_thisYear '/Irrig.base' num2str(base_year) '.txt'] ;
                    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                        'outPrec', outPrec, ...
                        'outWidth', outWidth, ...
                        'delimiter', delimiter, ...
                        'overwrite', overwrite, ...
                        'fancy', fancy, ...
                        'progress_step_pct', 20, ...
                        'verbose', false) ;
                end
                clear out_y1
            end
            
        end

        if save_2deg_any
            if ~exist(outDir_thisYear, 'dir')
                mkdir(outDir_thisYear)
            end
            % Save new LandCoverFract.txt (2-degree)
            out_y1.varNames = {'PASTURE','CROPLAND','NATURAL','BARREN'} ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'CROPLAND')) = out_y1_2deg_crop_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'PASTURE')) = out_y1_2deg_past_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'NATURAL')) = out_y1_2deg_ntrl_YX ;
            out_y1.maps_YXv(:,:,strcmp(out_y1.varNames,'BARREN')) = out_y1_2deg_bare_YX ;
            out_y1.maps_YXv = out_y1.maps_YXv ./ repmat(landArea_2deg_YX,[1 1 length(out_y1.varNames)]) ;
            if save_2deg_mat
                file_out = [outDir_thisYear '/LandCoverFract.base' num2str(base_year) '.2deg.mat'] ;
                save(file_out,'out_y1') ;
            end
            if save_2deg_txt
                [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                file_out = [outDir_thisYear '/LandCoverFract.base' num2str(base_year) '.2deg.txt'] ;
                lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                    'outPrec', outPrec, ...
                    'outWidth', outWidth, ...
                    'delimiter', delimiter, ...
                    'overwrite', overwrite, ...
                    'fancy', fancy, ...
                    'progress_step_pct', 20, ...
                    'verbose', false) ;
            end
            clear out_y1

            % Save new CropFract.txt (2-degree)
            if ~combineCrops
                out_y1.maps_YXv = out_y1_2deg_agri_YXv(:,:,~isAgri_isPast) ./ repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)]) ;
                out_y1.maps_YXv(repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
                out_y1.varNames = LPJGcrops ;
                if save_2deg_mat
                    file_out = [outDir_thisYear '/CropFract.base' num2str(base_year) '.2deg.mat'] ;
                    save(file_out,'out_y1') ;
                end
                if save_2deg_txt
                    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                    file_out = [outDir_thisYear '/CropFract.base' num2str(base_year) '.2deg.txt'] ;
                    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                        'outPrec', outPrec, ...
                        'outWidth', outWidth, ...
                        'delimiter', delimiter, ...
                        'overwrite', overwrite, ...
                        'fancy', fancy, ...
                        'progress_step_pct', 20, ...
                        'verbose', false) ;
                end
                clear out_y1
            end
            
            % Save new Fert.txt (2-degree)
            if ~combineCrops
                % Convert from kg/m2 to kg/ha for compatibility with original
                % PLUM style.
                out_y1.maps_YXv = 1e4*out_y1_2deg_nfert_YXv ;
                out_y1.varNames = LPJGcrops ;
                if save_2deg_mat
                    file_out = [outDir_thisYear '/Fert.base' num2str(base_year) '.2deg.mat'] ;
                    save(file_out,'out_y1') ;
                end
                if save_2deg_txt
                    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                    file_out = [outDir_thisYear '/Fert.base' num2str(base_year) '.2deg.txt'] ;
                    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                        'outPrec', outPrec, ...
                        'outWidth', outWidth, ...
                        'delimiter', delimiter, ...
                        'overwrite', overwrite, ...
                        'fancy', fancy, ...
                        'progress_step_pct', 20, ...
                        'verbose', false) ;
                    clear out_y1
                end
            end
            
            % Save new Irrig.txt (2-degree)
            if ~combineCrops
                out_y1.maps_YXv = out_y1_2deg_irrig_YXv ;
                out_y1.maps_YXv(repmat(out_y1_2deg_crop_YX,[1 1 length(LPJGcrops)])==0) = 0 ;
                out_y1.varNames = LPJGcrops ;
                if save_2deg_mat
                    file_out = [outDir_thisYear '/Irrig.base' num2str(base_year) '.2deg.mat'] ;
                    save(file_out,'out_y1') ;
                end
                if save_2deg_txt
                    [out_y1_array, out_header_cell] = lpjgu_matlab_maps2table(out_y1,list2map_2deg) ;
                    file_out = [outDir_thisYear '/Irrig.base' num2str(base_year) '.2deg.txt'] ;
                    lpjgu_matlab_saveTable(out_header_cell, out_y1_array, file_out,...
                        'outPrec', outPrec, ...
                        'outWidth', outWidth, ...
                        'delimiter', delimiter, ...
                        'overwrite', overwrite, ...
                        'fancy', fancy, ...
                        'progress_step_pct', 20, ...
                        'verbose', false) ;
                end
                clear out_y1
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Save last year's global total to time series %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isFirstYearInLoop
            in_ts = save_to_timeseries_struct(in_ts, ...
                in_y0, in_y0_nfert.maps_YXv, in_y0_irrig.maps_YXv, 1) ;
            out_ts = save_to_timeseries_struct(out_ts, ...
                out_y0, out_y0_nfert_YXv, out_y0_irrig_YXv, 1) ;
        end
        
        in_ts = save_to_timeseries_struct(in_ts, ...
            in_y1, in_y1_nfert.maps_YXv, in_y1_irrig.maps_YXv, y+1) ;
        
        out_y1 = out_y0 ;
        out_y1.maps_YXv = nan(size(out_y1.maps_YXv)) ;
        [~,IA,IB] = intersect(out_y1.varNames,LUnames_agri,'stable') ;
        out_y1.maps_YXv(:,:,IA) = out_y1_agri_YXv(:,:,IB) ;
        clear IA IB
        out_y1.maps_YXv(:,:,strcmp(LUnames,'NATURAL')) = out_y1_ntrl_YX ;
        out_y1.maps_YXv(:,:,strcmp(LUnames,'BARREN')) = out_y1_bare_YX ;
        out_ts = save_to_timeseries_struct(out_ts, ...
            out_y1, out_y1_nfert_YXv, out_y1_irrig_YXv, y+1) ;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Prepare for next iteration %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clear out_y1_past_YX
        if y < Nyears
            in_y0 = in_y1 ;
            in_y0_2deg = in_y1_2deg ;
            in_y0_agri_YXv = in_y1_agri_YXv ;
            in_y0_2deg_agri_YXv = in_y1_2deg_agri_YXv ;
            out_y0 = out_y1 ;
            out_y0_2deg = out_y1_2deg ;
            out_y0_2deg_agri_YXv = out_y1_2deg_agri_YXv ;
            out_y0_agri_YXv = out_y1_agri_YXv ;
            bareFrac_y0_YX = in_y1_bareFrac_YX ;
            in_y0_2deg_vegd_YX = in_y1_2deg_vegd_YX ;
            out_y0_vegd_YX = out_y1_vegd_YX ;
            out_y0_2deg_vegd_YX = out_y1_2deg_vegd_YX ;
            
            if ~combineCrops
                out_y0_2deg_nfert.maps_YXv = out_y1_2deg_nfert_YXv ;
                out_y0_2deg_irrig.maps_YXv = out_y1_2deg_irrig_YXv ;
                in_y0_2deg_nfert = in_y1_2deg_nfert ;
                in_y0_2deg_irrig = in_y1_2deg_irrig ;
                max_orig_nfert_y0 = max_orig_nfert_y1 ;
                max_nfert_y0 = max_nfert_y1 ;
                out_y0_nfert_YXv = out_y1_nfert_YXv ;
                out_y0_irrig_YXv = out_y1_irrig_YXv ;
                out_y0_2deg_nfert_YXv = out_y1_2deg_nfert_YXv ;
                out_y0_2deg_irrig_YXv = out_y1_2deg_irrig_YXv ;
                in_y0_nfert = in_y1_nfert ;
                in_y0_irrig = in_y1_irrig ;
            end
            clear *y1*
        end
        
        % Save full-precision outputs for use in restarting
        thisMATfile = [outDir_thisYear 'post.base' num2str(base_year) '.mat'] ;
        save(thisMATfile, ...
            '*y0*','latestPLUMin_*','-v7.3') ;
        if ~combineCrops
            save(thisMATfile, ...
                'latestPLUMin_2deg_nfert_YXv','latestPLUMin_2deg_irrig_YXv',...
                '-append') ;
        end

        disp(['  Done (' toc_hms(toc) ').'])
        isFirstYearInLoop = false ;


    end % years loop
    
    
    %%%%%%%%%%%%%%%%%%%%
    %%% Make figures %%%
    %%%%%%%%%%%%%%%%%%%%
    
    figDir = addslashifneeded(sprintf('%s_figs', removeslashifneeded(outDir))) ;
    if ~exist(figDir, 'dir')
        mkdir(figDir)
    end
    legend_ts = {'From PLUM', 'Harmonized'} ;
    
    % LU areas
    tmp_names = [{'CROPLAND'} in_ts.luNames(~isCrop)] ;
    ts_orig_cy = nan(length(tmp_names), length(in_ts.yearList)) ;
    ts_orig_cy(strcmp(tmp_names, 'CROPLAND'),:) = sum(in_ts.area_vy(isCrop,:),1) ;
    ts_orig_cy(~strcmp(tmp_names, 'CROPLAND'),:) = in_ts.area_vy(~isCrop,:) ;
    ts_harm_cy = ts_orig_cy ;
    ts_harm_cy(strcmp(tmp_names, 'CROPLAND'),:) = sum(out_ts.area_vy(isCrop,:),1) ;
    ts_harm_cy(~strcmp(tmp_names, 'CROPLAND'),:) = out_ts.area_vy(~isCrop,:) ;
    units = 'Million km2' ;
    ts_orig_cy = ts_orig_cy*1e-6*1e-6 ;
    ts_harm_cy = ts_harm_cy*1e-6*1e-6 ;
    make_crops_timeseries_fig([], ts_orig_cy, ts_harm_cy, ...
        tmp_names, legend_ts, [], in_ts.yearList, units, ...
        'Area', 'landUse', figDir)
    
    % Crops
    units = 'Million km2' ;
    ts_orig_cy = in_ts.area_vy(isCrop,:)*1e-6*1e-6 ;
    ts_harm_cy = out_ts.area_vy(isCrop,:)*1e-6*1e-6 ;
    make_crops_timeseries_fig([], ts_orig_cy, ts_harm_cy, ...
        in_ts.cropNames, legend_ts, [], in_ts.yearList, units, ...
        'Area', 'crops', figDir)
    
    % Fertilization
    units = 'Mt N' ;
    ts_orig_cy = in_ts.nfert_vy*1e-3*1e-6 ;
    ts_harm_cy = out_ts.nfert_vy*1e-3*1e-6 ;
    make_crops_timeseries_fig([], ts_orig_cy, ts_harm_cy, ...
        in_ts.cropNames, legend_ts, [], in_ts.yearList, units, ...
        'Fert.', 'nfert', figDir)
    
    % Irrigation
    units = 'intensity \times area' ;
    ts_orig_cy = in_ts.irrig_vy ;
    ts_harm_cy = out_ts.irrig_vy ;
    make_crops_timeseries_fig([], ts_orig_cy, ts_harm_cy, ...
        in_ts.cropNames, legend_ts, [], in_ts.yearList, units, ...
        'Irrigation', 'irrig', figDir)
    
    disp('Done')
    diary off

end


%% FUNCTIONS

function check_discrepancy(tmp_in, tmp_out, base_year, thisVar, thisName)

errPct = (tmp_in - tmp_out) / tmp_out * 100 ;
errSign = '' ;
if errPct > 0
    errSign = '+' ;
end
if tmp_out > 0 && abs(errPct) > 50
    warning('PLUM baseline (%d) discrepancy: %s: %s: %s%0.1f%%', ...
        base_year, thisVar, thisName, errSign, errPct)
end
if errPct < -99
    keyboard
end

end



function S_ts = save_to_timeseries_struct(S_ts, ...
    S_area, nfert_YXv, irrig_YXv, y)

S_ts.area_vy(:,y) = squeeze(nansum(nansum(S_area.maps_YXv, 1), 2)) ;

[C,IA,IB] = intersect(S_ts.cropNames, S_ts.luNames, 'stable') ;
if ~isequal(C, S_ts.cropNames)
    error('Crop list mismatch')
end

S_ts.nfert_vy(:,y) = squeeze(nansum(nansum( ...
    nfert_YXv(:,:,IA) .* S_area.maps_YXv(:,:,IB), 1), 2)) ;
S_ts.irrig_vy(:,y) = squeeze(nansum(nansum( ...
    irrig_YXv(:,:,IA) .* S_area.maps_YXv(:,:,IB), 1), 2)) ;


end


