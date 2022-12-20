%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for tradeoffs paper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_save_figs = true ;

do_save.LU              = true ;
do_save.crops           = false ;
do_save.yield           = false ;
do_save.yield_exp       = false ;
do_save.yield_map       = false ;
do_save.yield_exp_map   = false ;
do_save.irrig           = false ;
do_save.water           = false ;
do_save.carbon          = true ;
do_save.mrunoff         = false ;
do_save.albedo          = true ;
do_save.bvocs           = true ;
do_save.Nflux           = false ;
do_save.Nfert           = false ;
do_save.fpc             = false ;

inDir_list = {...
%     'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-11-03-233845' ;
%     'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-11-03-235127' ;
%     'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-11-03-222204' ;
%     'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-11-03-233706' ;
%     'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-11-05-085615' ;
%     'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-11-05-071233' ;
%     'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-11-05-071814' ;
%     'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-11-05-110842' ;
%     'LPJGPLUM_2011-2100_harm2_SSP1_constClimCO2/output-2018-11-03-233941' ;
%     'LPJGPLUM_2011-2100_harm2_SSP3_constClimCO2/output-2018-11-03-221433' ;
%     'LPJGPLUM_2011-2100_harm2_SSP4_constClimCO2/output-2018-11-03-221734' ;
%     'LPJGPLUM_2011-2100_harm2_SSP5_constClimCO2/output-2018-11-03-233913' ;
%     'LPJGPLUM_2011-2100_harm2_constLU_RCP45/output-2018-11-05-003042' ;
%     'LPJGPLUM_2011-2100_harm2_constLU_RCP60/output-2018-11-05-003344' ;
%     'LPJGPLUM_2011-2100_harm2_constLU_RCP85/output-2018-11-05-000411' ;
    } ;



%% Setup

test_cropfracs_20170108 = false ;

outDir = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper04_Sam_ESintxns/test_outs/' ;

thisDir = addslashifneeded('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work') ;
if ~exist(thisDir,'dir')
    error('thisDir does not exist')
end
cd(thisDir) ;
addpath(genpath(pwd))

is_baseline_list = false(length(inDir_list),1) ;
for d = 1:length(inDir_list)
    if contains(inDir_list{d},'1850-')
        is_baseline_list(d) = true ;
    end
end

pftList_noCrops = {'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE',...
            'TrBR','C3G','C4G','PC3G','PC4G',...
            'CC3G_ic','CC4G_ic','ExtraCrop','Total','Crop_sum','Pasture_sum',...
            'Natural_sum','Barren_sum'} ;
LUlist = {'CROPLAND','PASTURE','NATURAL','BARREN'} ;

% Anonymous function for checking whether any version of a file exists
anyfileexist = @(in_file) ...
    exist(in_file,'file') ...
    || exist([in_file '.maps.mat'],'file') ...
    || exist([in_file '.mat'],'file') ...
    || exist([in_file '.gz'],'file') ;

% Get calibration factors
LPJGcrops_2_PLUM = readtable([thisDir 'PLUM6xtra_calib_20180423.csv']) ;
cropTypes_conv = LPJGcrops_2_PLUM.Crop ;
calibFactors = LPJGcrops_2_PLUM.calibFactor ;
clear LPJGcrops_2_PLUM

% Import bare-soil albedo
if do_save.albedo
    baresoil_albedo_file = [thisDir 'soilmap.txt'] ;
    baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true,'verbose',false) ;
    baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
    clear baresoil_albedo
end

% Read PLUMout_gridlist
if any(is_baseline_list) || contains(inDir_list{d},'LPJGPLUM_2011-2100_harm2_constLU_')
    gridlist_file = '/Users/Shared/lpj-guess/gridlists/PLUMout_gridlist.txt' ;
    PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;
end

% Import land area (km2 to m2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
if any(is_baseline_list) || contains(inDir_list{d},'LPJGPLUM_2011-2100_harm2_constLU_')
    land_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
end
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
if any(is_baseline_list) || contains(inDir_list{d},'LPJGPLUM_2011-2100_harm2_constLU_')
    gcel_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
end
% Convert to m2
land_area_YXqd = land_area_YXqd*1e6 ;
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd

% Import Koeppen-Geiger map and set up interpretations
kg_file = '/Users/sam/Geodata/KoeppenGeiger/Koeppen-Geiger_v1.1/koeppen-geiger_GTiffs/koeppen-geiger_0.5.tif' ;
kg_YX = flipud(imread(kg_file)) ;
kg_YX(kg_YX==0) = NaN ;
kg_names = {'tropical','dry','temperate','continental','polar'} ;
kg_ranges = [10 19 ; 20 29 ; 30 39 ; 40 59 ; 60 69] ;
Nbiomes_kg = length(kg_names) ;
kg_maps_YXv = false([size(kg_YX) Nbiomes_kg]) ;
for b = 1:Nbiomes_kg
    kg_maps_YXv(:,:,b) = kg_YX>=kg_ranges(b,1) & kg_YX<=kg_ranges(b,2) ;
end

% Ice/water

file_luh2_etc = '/Users/sam/Geodata/LUH2/supporting/staticData_quarterdeg.nc' ;
icwtr_YX_qd = flipud(transpose(ncread(file_luh2_etc,'icwtr'))) .* land_area_YXqd ;
%%%% Convert to half-degree
tmp = icwtr_YX_qd(:,1:2:1440) + icwtr_YX_qd(:,2:2:1440) ;
icwtr_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
if any(is_baseline_list) || contains(inDir_list{d},'LPJGPLUM_2011-2100_harm2_constLU_')
    icwtr_YX(~PLUMout_gridlist.mask_YX) = NaN ;
end
icwtr_YX = icwtr_YX ./ land_area_YX ;

% Conversion factors
%%% All masses in kg
%%% All areas in m2
%%% All volumes in m3
cf_cpool = 1 ;   % LPJ-GUESS outputs kgC/m2
cf_water = 1e-3 ; % LPJ-GUESS outputs mm (will be multiplied by land_area_m2 --> m3)
cf_nflux = 1e-4 ;   % LPJ-GUESS outputs kgN/ha
cf_bvoc = 1e-6 ;   % LPJ-GUESS outputs mgC/m2


%% Loop through inDir_list

for d = 1:length(inDir_list)
    
    %%%%%%%%%%%%%
    %%% Setup %%%
    %%%%%%%%%%%%%
    
    % Get full inDir
    inDir = find_PLUM2LPJG_run(inDir_list{d}) ;
    is_baseline = is_baseline_list(d) ;
    if is_baseline
        disp('is_baseline')
    end
    
    timeseries_out = [inDir 'timeseries.mat'] ;
    timeseries_PLUMexp_out = [inDir 'timeseries_PLUMexp.mat'] ;
    firstdecade_out = [inDir 'first_decade.mat'] ;
    lastdecade_out = [inDir 'last_decade.mat'] ;
    
    disp(inDir)
    
    % Get baseline/not info
    if is_baseline
        yearList = 1850:2010 ;
        if ~exist('PLUMout_mask_YX1y','var')
            PLUMout_mask_YX1y = repmat(PLUMout_gridlist.mask_YX,[1 1 1 length(yearList)]) ;
        end
    else
        yearList = 2011:2100 ;
    end
    yearList = transpose(yearList) ;
    Nyears = length(yearList) ;
        
    % Get land use and crop fractions input files
    clear LUfile cropfile
    if strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_2006-2010_PLUM6xtra/output-2018-04-21-061836') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_PLUM6xtra/matlab_merge_20180424050342/')
        LUfile = '/Users/Shared/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
        cropfile = '/Users/Shared/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.NfertEmpties0-0200-1000.txt' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180502104840') ...
        || strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180605160616/')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
        cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.Misc0s.txt' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_PLUM6xtraMisc_cg/matlab_merge_20180801105737')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v3/LU_xtraCROPtoPAST.remapv3.20180214.cgFertIrr0.m0.txt' ;
        cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v3/cropfracs.remapv3.20180214.cgFertIrr0.m0.txt' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_remap6/output-2018-10-27-073916') ...
        || strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/LU.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
        cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v6/cropfracs.remapv6.20180214.ecFertIrr0.setaside0103.m4.someOfEachCrop.txt' ;
    else
        [x,LUfile_tmp] = unix([thisDir 'get_lu_file.sh "' inDir '"']) ;
        if x~=0
            error(['get_lu_file.sh failed with error ' num2str(x)])
        end
        LUfile_tmp = regexprep(LUfile_tmp,'[\n\r]+','') ; % Remove extraneous newline
        if contains(LUfile_tmp,'param "file_lu" (str "')
            warning('LUfile_tmp badly parsed. Trying crude substitution.')
            LUfile = crude_file_find(LUfile_tmp) ;
        else
            LUfile = find_PLUM2LPJG_input_file(LUfile_tmp) ;
        end
        [x,cropfile_tmp] = unix([thisDir 'get_cropfrac_file.sh "' inDir '"']) ;
        if x~=0
            error(['get_cropfrac_file.sh failed with error ' num2str(x)])
        end
        cropfile_tmp = regexprep(cropfile_tmp,'[\n\r]+','') ; % Remove extraneous newline
        if contains(cropfile_tmp,'param "file_lucrop" (str "')
            warning('cropfile_tmp badly parsed. Trying crude substitution.')
            cropfile = crude_file_find(cropfile_tmp) ;
        else
            cropfile = find_PLUM2LPJG_input_file(cropfile_tmp) ;
            if do_save.yield_exp || do_save.yield_exp_map
                expyieldfile_tmp = strrep(cropfile_tmp,'cropfractions','yield') ;
                expyieldfile = find_PLUM2LPJG_input_file(expyieldfile_tmp,false) ;
            end
        end
        
    end
    % Do expected yields exist? If not, make sure do_save.yield_exp* are
    % false. Reset back to true at end of this for loop (before going to
    % next directory)
    flip_do_save_yield_exp = false ;
    if do_save.yield_exp && (~exist('expyieldfile','var') || isempty(expyieldfile)) && ~exist('LUcropDir_plum7','var')
        warning('No expected yields exist; skipping do_save.yield_exp.')
        flip_do_save_yield_exp = true ;
        do_save.yield_exp = ~do_save.yield_exp ;
    end
    flip_do_save_yield_exp_map = false ;
    if do_save.yield_exp_map && (~exist('expyieldfile','var') || isempty(expyieldfile)) && ~exist('LUcropDir_plum7','var')
        warning('No expected yields exist; skipping do_save.yield_exp_map.')
        flip_do_save_yield_exp_map = true ;
        do_save.yield_exp_map = ~do_save.yield_exp_map ;
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get important info %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yield = lpjgu_matlab_readTable_then2map([inDir 'yield.out'],'force_mat_save',true) ;
    list_to_map = yield.list_to_map ;
    cropTypes = yield.varNames ;
    cropTypes(contains(cropTypes,{'CC3G_ic','CC4G_ic','ExtraCrop'})) = [] ;
    matches = [] ;
    for c = 1:length(cropTypes)
        thisCrop = cropTypes{c} ;
        thisMatch = find(strcmp(cropTypes,[thisCrop 'i'])) ;
        if length(thisMatch)==1
            matches = [matches thisMatch] ;
        elseif length(thisMatch)>1
            error('How is length(thisMatch)>1 ?')
        end
    end
    cropTypes(matches) = [] ;
    if any(strcmp(cropTypes,'CC3G'))
        cropTypes_conv = [cropTypes_conv ; 'CC3G'] ;
        calibFactors = [calibFactors ; 1] ;
    end
    if any(strcmp(cropTypes,'CC4G'))
        cropTypes_conv = [cropTypes_conv ; 'CC4G'] ;
        calibFactors = [calibFactors ; 1] ;
    end
    if any(strcmp(cropTypes,'ExtraCrop'))
        cropTypes_conv = [cropTypes_conv ; 'ExtraCrop'] ;
        calibFactors = [calibFactors ; 1] ;
    end
    Ncrops = length(cropTypes) ;
    clear matches c thisCrop thisMatch
    pftList = [pftList_noCrops cropTypes] ;
    Npad = 0 ;
    if min(yield.yearList) == min(yearList) - 5 && max(yield.yearList)==max(yearList)
        warning('Adjusting yearList to account for 5 years'' padding at beginning.')
        Npad = 5 ;
        yearList = yield.yearList ;
        Nyears = length(yearList) ;
    elseif ~isequal(yield.yearList,yearList)
        error('Rework so yearLists match!')
    end
    if ~(do_save.yield || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map)
        clear yield
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Import land use %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.LU || do_save.crops || do_save.irrig || do_save.fpc ...
    || do_save.albedo || do_save.Nfert || do_save.yield ...
    || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map
%%%%|| do_save.water || do_save.carbon || do_save.bvocs || do_save.Nflux || do_save.mrunoff

        if contains(LUfile,'someOfEachCrop')
            warning('RUN USED SOMEOFEACH LUFILE; IGNORING')
            LUfile = strrep(LUfile,'.someOfEachCrop','') ;
        end

        LU = lpjgu_matlab_readTable_then2map(LUfile,'force_mat_save',true) ;
        if contains(LUfile,'LU_xtraCROPtoPAST')
            ignoredCropFracs_file = strrep(strrep(LUfile,...
                'LU_xtraCROPtoPAST', 'cropfracsIGNORED'),...
                '.txt',              '.mat') ;
            if exist(ignoredCropFracs_file,'file')
                load(ignoredCropFracs_file) ;
            else
                warning('LUfile has extra crop in pasture, but ignoredCropFracs_file does not exist.')
            end
        end
        if ~isfield(LU,'yearList') && size(LU.maps_YXv,4)==1
            % One-year LU file
            LU.maps_YXvy = repmat(LU.maps_YXv,[1 1 1 Nyears]) ;
            LU.yearList = yearList ;
            if exist('ignored_LUCROParea_YX_y','var')
                ignored_LUCROParea_YX_y = repmat(ignored_LUCROParea_YX_y,[1 1 1 Nyears]) ;
            end
        end
        % Check LU types, removing extraneous if necessary
        if any(strcmp(LU.varNames,'URBAN'))
            thisMax = max(max(max(LU.maps_YXvy(:,:,strcmp(LU.varNames,'URBAN'),:)))) ;
            if thisMax > 0
                warning(['Removing LU: URBAN even though up to ' num2str(thisMax) ' is present.'])
            end
            LU.maps_YXvy(:,:,strcmp(LU.varNames,'URBAN'),:) = [] ;
            LU.varNames(strcmp(LU.varNames,'URBAN')) = [] ;
        end
        if any(strcmp(LU.varNames,'PEATLAND'))
            thisMax = max(max(max(LU.maps_YXvy(:,:,strcmp(LU.varNames,'PEATLAND'),:)))) ;
            if thisMax > 0
                warning(['Removing LU: PEATLAND even though up to ' num2str(thisMax) ' is present.'])
            end
            LU.maps_YXvy(:,:,strcmp(LU.varNames,'PEATLAND'),:) = [] ;
            LU.varNames(strcmp(LU.varNames,'PEATLAND')) = [] ;
        end
        if any(~(strcmp(LU.varNames,'PASTURE') | strcmp(LU.varNames,'CROPLAND') | strcmp(LU.varNames,'NATURAL') | strcmp(LU.varNames,'BARREN')))
            error('Something is wrong with LU.varList.')
        end

        % Make LU0, if accounting for CROPLAND moved to PASTURE
        if exist('ignored_LUCROParea_YX_y','var')
            if size(ignored_LUCROParea_YX_y,4)~=size(LU.maps_YXvy,4)
                error('Number of years in ignored_LUCROParea_YX_y does not match that of LU.maps_YXvy!')
            end
            LU0 = LU ;
            LU0.maps_YXvy(:,:,strcmp(LU0.varNames,'CROPLAND'),:) =  ignored_LUCROParea_YX_y + LU0.maps_YXvy(:,:,strcmp(LU0.varNames,'CROPLAND'),:) ;
            LU0.maps_YXvy(:,:,strcmp(LU0.varNames,'PASTURE'),:)  = -ignored_LUCROParea_YX_y + LU0.maps_YXvy(:,:,strcmp(LU0.varNames,'PASTURE'),:) ;
            clear ignored_LUCROParea_YX_y
        end

        % Align yearLists, if needed
        if ~isequal(LU.yearList,yearList)
            if min(LU.yearList) <= min(yearList) && max(LU.yearList) >= max(yearList)
                LU.maps_YXvy = LU.maps_YXvy(:,:,:,LU.yearList>=min(yearList) & LU.yearList<=max(yearList)) ;
                LU.yearList = transpose(LU.yearList(1):max(yearList)) ;
            end
            if min(LU.yearList) == min(yearList) + 5 && max(LU.yearList)==max(yearList)
                warning('Adjusting LU to account for 5 years'' padding at beginning.')
                LU.yearList = yearList ;
                LU.maps_YXvy = cat(4,repmat(LU.maps_YXvy(:,:,:,1),[1 1 1 5]),LU.maps_YXvy) ;
            elseif ~isequal(LU.yearList,yearList)
    %             error('Rework LU so yearLists match!')
                warning('Yearlists don''t match. Just doing intersect.')
                [~,IA] = intersect(LU.yearList,yearList) ;
                LU.maps_YXvy = LU.maps_YXvy(:,:,:,IA) ;
                LU.yearList = yearList ;
            end
        end
        if exist('LU0','var')
            if min(LU0.yearList) <= min(yearList) && max(LU0.yearList) >= max(yearList)
                LU0.maps_YXvy = LU0.maps_YXvy(:,:,:,LU0.yearList>=min(yearList) & LU0.yearList<=max(yearList)) ;
                LU0.yearList = transpose(LU0.yearList(1):max(yearList)) ;
            end
            if min(LU0.yearList) == min(yearList) + 5 && max(LU0.yearList)==max(yearList)
                warning('Adjusting LU0 to account for 5 years'' padding at beginning.')
                LU0.yearList = yearList ;
                LU0.maps_YXvy = cat(4,repmat(LU0.maps_YXvy(:,:,:,1),[1 1 1 5]),LU0.maps_YXvy) ;
            elseif ~isequal(LU0.yearList,yearList)
                %                 error('Rework LU0 so yearLists match!')
                warning('Yearlists don''t match. Just doing intersect.')
                [~,IA] = intersect(LU0.yearList,yearList) ;
                LU0.maps_YXvy = LU0.maps_YXvy(:,:,:,IA) ;
                LU0.yearList = yearList ;
            end
        end

        % Align NaN mask, if needed
        if is_baseline
            LU.list_to_map = PLUMout_gridlist.list_to_map ;
            LU.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(LU.varNames) 1])) = NaN ;
            if exist('LU0','var')
                LU0.list_to_map = PLUMout_gridlist.list_to_map ;
                LU0.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(LU0.varNames) 1])) = NaN ;
            end
        end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Import crop fractions %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if contains(cropfile,'someOfEachCrop')
            warning('RUN USED SOMEOFEACH CROPFILE; IGNORING')
            LUfile = strrep(cropfile,'.someOfEachCrop','') ;
        end

        cropfracs = lpjgu_matlab_readTable_then2map(cropfile,'force_mat_save',true) ;
        cropfracs = adjust_cropinput_yearLists(cropfracs, yearList) ;
        cropfracs_orig = cropfracs ;

        % Will we be merging rainfed and irrigated? Or just replacing rainfed
        % with irrigated?
        is_rainfed = false(length(cropfracs_orig.varNames),1) ;
        for c = 1:length(cropfracs_orig.varNames)
            thisCrop = cropfracs_orig.varNames{c} ;
            if any(strcmp(cropfracs_orig.varNames,[thisCrop 'i']))
                is_rainfed(c) = true ;
            end
        end
        has_some_rainfed = ~isempty(find(cropfracs_orig.maps_YXvy(:,:,is_rainfed,:)>0,1)) ;
        if has_some_rainfed
            merge_or_replace = 'merge' ;
        else
            merge_or_replace = 'replace' ;
        end

        [cropfracs, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig, inds_cropTypes_nonCGs_cropFracsOrig] = ...
            CrOp_and_CrOpi(cropfracs, 'cropfracs', cropTypes, merge_or_replace, cropfile) ;

        if is_baseline
            cropfracs.list_to_map = PLUMout_gridlist.list_to_map ;
            cropfracs.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs.varNames) 1])) = NaN ;
            if exist('cropfracs_orig','var')
                cropfracs_orig.list_to_map = PLUMout_gridlist.list_to_map ;
                cropfracs_orig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs_orig.varNames) 1])) = NaN ;
            end
        end
        % Get area of each crop (m2)
        cropareas = cropfracs ;
        cropareas.maps_YXvy = ...
            cropfracs.maps_YXvy ...
            .* repmat(LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:),[1 1 length(cropfracs.varNames) 1]) ...
            .* repmat(land_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ;
        if exist('LU0','var')
            cropareas0 = cropfracs ;
            cropareas0.maps_YXvy = ...
                cropfracs.maps_YXvy ...
                .* repmat(LU0.maps_YXvy(:,:,strcmp(LU0.varNames,'CROPLAND'),:),[1 1 length(cropfracs.varNames) 1]) ...
                .* repmat(land_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ;
        end
    
    end
    
    
    %%%%%%%%%%%%%%%%%%
    %%% Import FPC %%%
    %%%%%%%%%%%%%%%%%%
    
    if do_save.fpc || do_save.albedo
        fpc = lpjgu_matlab_readTable_then2map([inDir 'fpc.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('Processing FPC...')
        [fpc, ~, ~] = CrOp_and_CrOpi(fpc, 'fpc', cropTypes, merge_or_replace) ;
        % Normalize to 1
        fpc_total_YX1y = fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Total'),:) ;
        isbad_YX1y = fpc_total_YX1y > 1 ;
        i = 0 ;
        while any(isbad_YX1y(:))
            i = i + 1 ;
            if i > 5
                error('Too many iterations!')
            end
            isbad_YXvy = repmat(isbad_YX1y,[1 1 size(fpc.maps_YXvy,3) 1]) ;
            fpc_total_YXvy = repmat(fpc_total_YX1y,[1 1 size(fpc.maps_YXvy,3) 1]) ;
            fpc.maps_YXvy(isbad_YXvy) = fpc.maps_YXvy(isbad_YXvy) ./ fpc_total_YXvy(isbad_YXvy) ;
            clear fpc_total_YXvy
            fpc_total_YX1y = fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Total'),:) ;
            isbad_YX1y = fpc_total_YX1y > 1 ;
        end
        clear fpc_total_YX1y isbad* i
        if is_baseline
            list_to_map = PLUMout_gridlist.list_to_map ;
            fpc.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(fpc.varNames) 1])) = NaN ;
        end
    end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save yield (kgDM/m2 --> kgDM) (i.e., actually PRODUCTION) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_save.yield || do_save.yield_map
        
        [yield, ~, ~] = CrOp_and_CrOpi(yield, 'yield', cropTypes, merge_or_replace, ...
            '', cropfracs_orig, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig, inds_cropTypes_nonCGs_cropFracsOrig) ;
        yield.maps_YXvy(:,:,contains(yield.varNames,{'CC3G_ic','CC4G_ic','ExtraCrop'}),:) = [] ;
        tmp = yield.varNames ;
        tmp(contains(tmp,{'CC3G_ic','CC4G_ic','ExtraCrop'})) = [] ;
        if ~isequal(sort(tmp),sort(cropTypes))
            error('yield doesn''t have the same CFTs as cropTypes!')
        end
        yield.varNames = tmp ;
        
        
    end
    clear yield
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save EXPECTED yield (kgDM/m2 --> kgDM) (i.e., actually PRODUCTION) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_save.yield_exp || do_save.yield_exp_map
        expyield = lpjgu_matlab_readTable_then2map(expyieldfile,'force_mat_save',true) ;
        if ~isfield(expyield,'maps_YXvy')
            warning('This appears to be a constant-LU run. Skipping expected yield.')
        else
        end
    end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save irrigation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.irrig
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of irrigation?')
        end
        irrig = lpjgu_matlab_readTable_then2map([inDir 'irrigation.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        gsirrig = lpjgu_matlab_readTable_then2map([inDir 'gsirrigation.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving irrigation...')
        [gsirrig, ~, ~] = CrOp_and_CrOpi(gsirrig, 'gsirrig', cropTypes, merge_or_replace, ...
            '', cropfracs_orig, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig, inds_cropTypes_nonCGs_cropFracsOrig) ;
        gsirrig.maps_YXvy(:,:,contains(gsirrig.varNames,{'CC3G_ic','CC4G_ic','ExtraCrop'}),:) = [] ;
        gsirrig.varNames(contains(gsirrig.varNames,{'CC3G_ic','CC4G_ic','ExtraCrop'})) = [] ;
        if is_baseline
            irrig.list_to_map = PLUMout_gridlist.list_to_map ;
            irrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(irrig.varNames) 1])) = NaN ;
            gsirrig.list_to_map = PLUMout_gridlist.list_to_map ;
            gsirrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(gsirrig.varNames) 1])) = NaN ;
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save other annual water %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.water && anyfileexist([inDir 'awater.out'])
        awater = lpjgu_matlab_readTable_then2map([inDir 'awater.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
    elseif do_save.water
        % Evaporation
        if anyfileexist([inDir 'aevap.out']) || anyfileexist([inDir 'mevap.out'])
            if anyfileexist([inDir 'aevap.out'])
                aevap = lpjgu_matlab_readTable_then2map([inDir 'aevap.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
            elseif anyfileexist([inDir 'mevap.out'])
                mevap = lpjgu_matlab_readTable_then2map([inDir 'mevap.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
                monthlengths = [31 28 31 30 31 30 31 31 30 31 30 31] ;
                aevap.list_to_map = mevap.list_to_map ;
                aevap.varNames = {'Total'} ;
                aevap.yearList = mevap.yearList ;
                aevap.maps_YXvy = sum(mevap.maps_YXvy .* repmat(permute(monthlengths,[4 3 2 1]),size(PLUMout_mask_YX1y)) / 365, 3) ;
            else
                error('???') ;
            end
            if is_baseline
                aevap.list_to_map = PLUMout_gridlist.list_to_map ;
                aevap.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aevap.varNames) 1])) = NaN ;
            end
        else
            warning('Evaporation missing! Skipping.')
        end
        
        % Evapotranspiration (mm/yr)
        aaet = lpjgu_matlab_readTable_then2map([inDir 'aaet.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving aaet...')
        if is_baseline
            aaet.list_to_map = PLUMout_gridlist.list_to_map ;
            aaet.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aaet.varNames) 1])) = NaN ;
        end
        
        % Runoff (mm/yr)
        tot_runoff = lpjgu_matlab_readTable_then2map([inDir 'tot_runoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving runoff...')
        if is_baseline
            tot_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            tot_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(tot_runoff.varNames) 1])) = NaN ;
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save carbon (kgC/m2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.carbon
        cpool = lpjgu_matlab_readTable_then2map([inDir 'cpool.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving cpool...')
        if is_baseline
            cpool.list_to_map = PLUMout_gridlist.list_to_map ;
            cpool.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cpool.varNames) 1])) = NaN ;
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save monthly runoff (mm/month) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.mrunoff
        mon_runoff = lpjgu_matlab_readTable_then2map([inDir 'mrunoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        if is_baseline
            mon_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            mon_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(mon_runoff.varNames) 1])) = NaN ;
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save albedo %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.albedo
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of albedo?')
        end
        snowdepth = lpjgu_matlab_readTable_then2map([inDir 'msnowdepth.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Getting and saving albedo...') 
        if is_baseline
            snowdepth.list_to_map = PLUMout_gridlist.list_to_map ;
            snowdepth.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(snowdepth.varNames) 1])) = NaN ;
        end
        
        % Remove ice/water from barren
        warning('Removing ice/water from barren')
        LU_tmp = LU ;
        tmp = LU_tmp.maps_YXvy(:,:,strcmp(LU_tmp.varNames,'BARREN'),:) ;
        tmp = tmp - repmat(icwtr_YX,[1 1 1 Nyears]) ;
        if min(min(min(tmp,[],4),[],2),[],1)<-1e-6
            error('Removing ice/water from BARREN results in too-negative values!')
        end
        tmp(tmp<0) = 0 ;
        LU_tmp.maps_YXvy(:,:,strcmp(LU_tmp.varNames,'BARREN'),:) = tmp ;
        LU_tmp.maps_YXvy = LU_tmp.maps_YXvy ./ repmat(sum(LU_tmp.maps_YXvy,3),[1 1 length(LU_tmp.varNames) 1]) ;
        
        [albedo_jan_YXy,albedo_jul_YXy] = ...
            get_albedo(fpc, snowdepth, LU_tmp, baresoil_albedo_YX, land_area_YX, pftList) ;
        
        clear LU_tmp
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save BVOCs (isoprene, monoterpenes) (mgC/m2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.bvocs
        if anyfileexist([inDir 'aiso_smry.out'])
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        if is_baseline
            aiso.list_to_map = PLUMout_gridlist.list_to_map ;
            aiso.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aiso.varNames) 1])) = NaN ;
        end
        
        if anyfileexist([inDir 'amon_smry.out'])
            amon = lpjgu_matlab_readTable_then2map([inDir 'amon_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            amon = lpjgu_matlab_readTable_then2map([inDir 'amon.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        if is_baseline
            amon.list_to_map = PLUMout_gridlist.list_to_map ;
            amon.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(amon.varNames) 1])) = NaN ;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save N flux (kgN/ha) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.Nflux
        
        thisArea_YX = land_area_YX ;
%         thisArea_YX = gcel_area_YX ;
        
        % From nflux.out
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of Nflux?')
        end
        nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        if is_baseline
            nflux.list_to_map = PLUMout_gridlist.list_to_map ;
            nflux.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(nflux.varNames) 1])) = NaN ;
        end
    end
    
    if do_save.Nfert
        % Fertilizer for each crop
%         [x,NfertFile_tmp] = unix([thisDir 'get_nfert_file.sh "' inDir '"']) ;
        cmd = 'grep "file_nfert" %s/landcover.ins | sed ''s@param "file_nfert" (str "/project/fh1-project-lpjgpi/lr8247@/Users/Shared@'' | sed ''s@")@@''' ;
        cmd_str = sprintf(cmd,removeslashifneeded(inDir)) ;
        [x,NfertFile_tmp] = unix(cmd_str) ;
        if x~=0
            error(['get_nfert_file.sh failed with error ' num2str(x)])
        end
        NfertFile_tmp = regexprep(NfertFile_tmp,'[\n\r]+','') ; % Remove extraneous newline
        if is_baseline
            if isempty(dir([NfertFile_tmp '*']))
                NfertFile = ['/project/fh1-project-lpjgpi/lr8247/PLUM/input/Nfert/' NfertFile_tmp] ;
                if ~exist(NfertFile,'file')
                    error('NfertFile not found!')
                end
            else
                NfertFile = NfertFile_tmp ;
            end
        elseif strcmp(NfertFile_tmp,'nfert_2010_luh2_aggregate_sum2x2_midpoint_rescaled_v20.txt')
            NfertFile = '/Users/Shared/PLUM/input/Nfert/LUH2/nfert_2010_luh2_aggregate_sum2x2_midpoint_rescaled_v20.txt' ;
        else
%             NfertFile = find_PLUM2LPJG_input_file(NfertFile_tmp) ;
            NfertFile = crude_file_find(NfertFile_tmp) ;
        end
        
        Nfert = lpjgu_matlab_readTable_then2map(NfertFile,'force_mat_save',true) ;
        Nfert = adjust_cropinput_yearLists(Nfert, yearList) ;
        if isequal(sort(intersect(Nfert.varNames,{'CC3ann','CC3per','CC3nfx','CC4ann','CC4per'})),sort(Nfert.varNames))
            Nfert_orig = Nfert ;
            Nfert = rmfield(Nfert,'maps_YXvy') ;
            Nfert.varNames = cropTypes ;
            NfertIndices = nan(Ncrops,1) ;
            for c = 1:Ncrops
                thisCrop = cropTypes{c} ;
                switch thisCrop
                    case {'CerealsC3', 'Oilcrops', 'StarchyRoots', 'Pulses', 'Rice', ...
                          'CerealsC3i','Oilcropsi','StarchyRootsi','Pulsesi','Ricei'}
                      thisNfertType = 'CC3ann' ;
                    case {'CerealsC4', 'Miscanthus', ...
                          'CerealsC4i','Miscanthusi'}
                      thisNfertType = 'CC4ann' ;
                    case {'CC3G', 'CC4G', 'ExtraCrop'}
                        thisNfertType = 'none' ;
                    otherwise
                        error(['Which NfertType should I use for ' thisCrop '?'])
                end
                if strcmp(thisNfertType,'none') && ~any(strcmp(Nfert_orig.varNames,'none'))
                    Nfert_orig.varNames = [Nfert_orig.varNames 'none'] ;
                    zeros_YX1y = zeros(size(Nfert_orig.maps_YXvy,1), size(Nfert_orig.maps_YXvy,2), 1, size(Nfert_orig.maps_YXvy,4)) ;
                    Nfert_orig.maps_YXvy = cat(3,Nfert_orig.maps_YXvy,zeros_YX1y) ;
                end
                thisNfertIndex = find(strcmp(Nfert_orig.varNames,thisNfertType)) ;
                if isempty(thisNfertIndex)
                    error('isempty(thisNfertIndex)')
                end
                NfertIndices(c) = thisNfertIndex ;
            end
            Nfert.maps_YXvy = Nfert_orig.maps_YXvy(:,:,NfertIndices,:) ;
            clear Nfert_orig
        else
            Nfert = CrOp_and_CrOpi(Nfert, 'Nfert', cropTypes, 'replace') ;
        end
        
    end
        
    
    %%%%%%%%%%%%%%%%%%%%
    %%% Housekeeping %%%
    %%%%%%%%%%%%%%%%%%%%
    
    if flip_do_save_yield_exp
        do_save.yield_exp = ~do_save.yield_exp ;
    end
    if flip_do_save_yield_exp_map
        do_save.yield_exp_map = ~do_save.yield_exp_map ;
    end
    
    disp('Done.')
    
end


%% Testing figure: Veg C vs. Isoprene 

thisPos = figurePos ;

decade = 2090 ;
yr1 = decade ;
yrN = decade + 9 ;

isCrop = strcmp(cpool.varNames,'VegC') ;
tmp1_YX = mean(cpool.maps_YXvy(:,:,isCrop,cpool.yearList>=yr1 & cpool.yearList<=yrN),4) ;
tmp2_YX = max(max(aiso.maps_YXvy(:,:,:,aiso.yearList>=yr1 & aiso.yearList<=yrN),[],4),[],3) ;

figure('Color','w','Position',thisPos) ;
ESintxns_subplots_biomes( ...
    tmp1_YX, tmp2_YX, kg_names, kg_maps_YXv, ...
    'marker', '.', ...
    'markerSize', 25, ...
    'markerColor', 'b', ...
    'markerAlpha', 0.05, ...
    'xlab', 'Vegetation C pool', ...
    'ylab', 'Isoprene emissions', ...
    'titl_prefix', 'Veg C vs. Isoprene', ...
    'bestFitType', 'linear') ;

if do_save_figs
    outFile = sprintf('%sVegC_vs_Isoprene_%ds.png',outDir,decade) ;
    export_fig(outFile,'-r150') ;
    close
end


%% Testing figure: Veg C vs. Monoterpenes 

thisPos = figurePos ;

decade = 2090 ;
yr1 = decade ;
yrN = decade + 9 ;

isCrop = strcmp(cpool.varNames,'VegC') ;
tmp1_YX = mean(cpool.maps_YXvy(:,:,isCrop,cpool.yearList>=yr1 & cpool.yearList<=yrN),4) ;
tmp2_YX = max(max(amon.maps_YXvy(:,:,:,yearList>=yr1 & yearList<=yrN),[],4),[],3) ;

figure('Color','w','Position',thisPos) ;
ESintxns_subplots_biomes( ...
    tmp1_YX, tmp2_YX, kg_names, kg_maps_YXv, ...
    'marker', '.', ...
    'markerSize', 25, ...
    'markerColor', 'b', ...
    'markerAlpha', 0.05, ...
    'xlab', 'Vegetation C pool', ...
    'ylab', 'Monoterpene emissions', ...
    'titl_prefix', 'Veg C vs. Monoterpene emissions', ...
    'bestFitType', 'linear') ;

if do_save_figs
    outFile = sprintf('%sVegC_vs_Monoterpenes_%ds.png',outDir,decade) ;
    export_fig(outFile,'-r150') ;
    close
end


%% Testing figure: Veg C vs. Max(Jan., Jul.) Albedo

thisPos = figurePos ;

decade = 2090 ;
yr1 = decade ;
yrN = decade + 9 ;

isCrop = strcmp(cpool.varNames,'VegC') ;
tmp1_YX = mean(cpool.maps_YXvy(:,:,isCrop,yearList>=yr1 & yearList<=yrN),4) ;
tmp2_YX = max(cat(3, ...
    mean(albedo_jan_YXy(:,:,yearList>=yr1 & yearList<=yrN),3), ...
    mean(albedo_jul_YXy(:,:,yearList>=yr1 & yearList<=yrN),3)), [], 3) ;

figure('Color','w','Position',thisPos) ;
ESintxns_subplots_biomes( ...
    tmp1_YX, tmp2_YX, kg_names, kg_maps_YXv, ...
    'marker', '.', ...
    'markerSize', 25, ...
    'markerColor', 'b', ...
    'markerAlpha', 0.05, ...
    'xlab', 'Vegetation C pool', ...
    'ylab', 'Max albedo', ...
    'titl_prefix', 'Veg C vs. Max Albedo', ...
    'bestFitType', 'linear') ;

if do_save_figs
    outFile = sprintf('%sVegC_vs_Albedo_%ds.png',outDir,decade) ;
    export_fig(outFile,'-r150') ;
    close
end


%% Testing figure: Natural area vs. Max(Jan., Jul.) Albedo

thisPos = figurePos ;

decade = 2090 ;
yr1 = decade ;
yrN = decade + 9 ;

tmp1_YX = mean(sum(LU.maps_YXvy(:,:,strcmp(LU.varNames,'NATURAL'),yearList>=yr1 & yearList<=yrN),3),4) ;
tmp2_YX = max(cat(3, ...
    mean(albedo_jan_YXy(:,:,yearList>=yr1 & yearList<=yrN),3), ...
    mean(albedo_jul_YXy(:,:,yearList>=yr1 & yearList<=yrN),3)), [], 3) ;

figure('Color','w','Position',thisPos) ;
ESintxns_subplots_biomes( ...
    tmp1_YX, tmp2_YX, kg_names, kg_maps_YXv, ...
    'marker', '.', ...
    'markerSize', 25, ...
    'markerColor', 'b', ...
    'markerAlpha', 0.05, ...
    'xlab', 'Natural frac.', ...
    'ylab', 'Max albedo', ...
    'titl_prefix', 'Natural frac. vs. Max Albedo', ...
    'bestFitType', 'linear') ;

if do_save_figs
    outFile = sprintf('%sNaturalFrac_vs_Albedo_%ds.png',outDir,decade) ;
    export_fig(outFile,'-r150') ;
    close
end

