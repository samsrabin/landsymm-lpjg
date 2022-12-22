%% Setup

test_cropfracs_20170108 = false ;

onMac = strcmp(thisSystem, 'ssr_mac') ;
addpath(genpath(landsymm_lpjg_path()))
gridlist_file = sprintf('%s/%s', paper02_repo_path, gridlist_file) ;
biomes_map_file = sprintf('%s/input_data/wwf_terr_ecos_UnpackClip.halfDeg.tif', paper02_repo_path)  ;
biomes_key_file = sprintf('%s/input_data/wwf_terr_ecos.codes.csv', paper02_repo_path)  ;
landarea_file = sprintf('%s/input_data/staticData_quarterdeg.nc', plumharm_repo_path) ;
baresoil_albedo_file = sprintf('%s/input_data/soilmap.txt', paper02_repo_path) ;

is_baseline_list = false(length(inDir_list),1) ;
do_PLUMout_gridlist_adjust_list = false(length(inDir_list),1) ;
for d = 1:length(inDir_list)
    if contains(inDir_list{d},'1850-')
        is_baseline_list(d) = true ;
        do_PLUMout_gridlist_adjust_list(d) = true ;
    elseif contains(inDir_list{d},'constLU')
        do_PLUMout_gridlist_adjust_list(d) = true ;
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
LPJGcrops_2_PLUM = readtable(calib_file) ;
cropTypes_conv = LPJGcrops_2_PLUM.Crop ;
calibFactors = LPJGcrops_2_PLUM.calibFactor ;
clear LPJGcrops_2_PLUM

% Import bare-soil albedo
if do_save.albedo
    baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true,'verbose',false) ;
    baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
    clear baresoil_albedo
end

% Read PLUMout_gridlist, if needed
do_PLUMout_gridlist_adjust = any(do_PLUMout_gridlist_adjust_list) ;
if do_PLUMout_gridlist_adjust
    PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;
end

% Import land area (km2 to m2)
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX_orig = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX_orig = tmp(1:2:720,:) + tmp(2:2:720,:) ;
% Convert to m2
land_area_YX_orig = land_area_YX_orig*1e6 ;
gcel_area_YX_orig = gcel_area_YX_orig*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Import biomes
biomes_YX = flipud(imread(biomes_map_file)) ;
biomes_YX(biomes_YX<0) = NaN ;
biomes_key = readtable(biomes_key_file) ;

% Conversion factors
%%% All masses in kg
%%% All areas in m2
%%% All volumes in m3
cf_cpool = 1 ;   % LPJ-GUESS outputs kgC/m2
cf_water = 1e-3 ; % LPJ-GUESS outputs mm (will be multiplied by land_area_m2 --> m3)
cf_nflux = 1e-4 ;   % LPJ-GUESS outputs kgN/ha
cf_bvoc = 1e-6 ;   % LPJ-GUESS outputs mgC/m2

% Set up for "last 30 years of century" stats
last30_statList = { ...
    'last30_mean' ; ...
    'last30_median' ; ...
    'last30_std' ; ...
    'last30_prctile05' ; ...
%     'last30_prctile25' ; ...
%     'last30_prctile75' ; ...
    'last30_prctile95' ; ...
    'last30_min' ; ...
    'last30_max' ; ...
    } ;
Nstats = length(last30_statList) ;
last30_mean = @(array_YXvy) mean(array_YXvy,4) ;
last30_median = @(array_YXvy) median(array_YXvy,4) ;
last30_std = @(array_YXvy) std(array_YXvy,1,4) ;
last30_prctile05 = @(array_YXvy) prctile(array_YXvy,5,4) ;
last30_prctile25 = @(array_YXvy) prctile(array_YXvy,25,4) ;
last30_prctile75 = @(array_YXvy) prctile(array_YXvy,75,4) ;
last30_prctile95 = @(array_YXvy) prctile(array_YXvy,95,4) ;
last30_min = @(array_YXvy) min(array_YXvy,[],4) ;
last30_max = @(array_YXvy) max(array_YXvy,[],4) ;
last30_statHandles = cell(Nstats,1) ;
for s = 1:Nstats
    last30_statHandles{s} = eval(last30_statList{s}) ;
end; clear s


%% Do loop

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
    do_PLUMout_gridlist_adjust = do_PLUMout_gridlist_adjust_list(d) ;
    if do_PLUMout_gridlist_adjust
        disp('do_PLUMout_gridlist_adjust')
    end
    
    timeseries_out = [inDir 'timeseries.mat'] ;
    timeseries_PLUMexp_out = [inDir 'timeseries_PLUMexp.mat'] ;
    firstdecade_out = [inDir 'first_decade.mat'] ;
    lastdecade_out = [inDir 'last_decade.mat'] ;
    last30_out = [inDir 'last_30yrs.mat'] ;
    
    disp(inDir)
    
    % Get baseline/not info
    if is_baseline
        yearList = 1850:2010 ;
        last30_yearList = 1971:2000 ;
    else
        yearList = 2011:2100 ;
        last30_yearList = 2071:2100 ;
    end
    land_area_YX = land_area_YX_orig ;
    gcel_area_YX = gcel_area_YX_orig ;
    if do_PLUMout_gridlist_adjust
        land_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
        gcel_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
        if ~exist('PLUMout_mask_YX1y','var')
            PLUMout_mask_YX1y = repmat(PLUMout_gridlist.mask_YX,[1 1 1 length(yearList)]) ;
        end
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
        if onMac
            % Get LU file
            [x,LUfile_tmp] = unix([addslashifneeded(paper02_repo_path) 'get_lu_file.sh "' inDir '"']) ;
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
            % Get crop fractions file
            [x,cropfile_tmp] = unix([addslashifneeded(paper02_repo_path) 'get_cropfrac_file.sh "' inDir '"']) ;
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
        else
            % Get LU file
            cmd = sprintf('grep ''param "file_lu"'' %s/landcover.ins | grep -v -e "^[[:blank:]]!" | sed ''s@param "file_lu"@@'' | sed ''s@(str@@'' | sed ''s@)@@'' | sed ''s@"@@g''', ...
                inDir) ;
            [x,LUfile_tmp] = unix(cmd) ;
            if x~=0
                error(['Failed when trying to find LU file, with error ' num2str(x)])
            end
            LUfile_tmp = regexprep(LUfile_tmp,'[\n\r]+','') ; % Remove extraneous newline
            LUfile = strrep(LUfile_tmp,' ','') ; % Remove extraneous spaces
            % Get crop fractions file
            cmd = sprintf('grep ''param "file_lucrop"'' %s/landcover.ins | grep -v -e "^[[:blank:]]!" | sed ''s@param "file_lucrop"@@'' | sed ''s@(str@@'' | sed ''s@)@@'' | sed ''s@"@@g''', ...
                inDir) ;
            [x,cropfile_tmp] = unix(cmd) ;
            if x~=0
                error(['Failed when trying to find crop fractions file, with error ' num2str(x)])
            end
            cropfile_tmp = regexprep(cropfile_tmp,'[\n\r]+','') ; % Remove extraneous newline
            cropfile = strrep(cropfile_tmp,' ','') ; % Remove extraneous spaces
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
    || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map ...
    || do_save.water || do_save.carbon || do_save.bvocs || do_save.Nflux || do_save.mrunoff

        if contains(LUfile,'someOfEachCrop')
            warning('RUN USED SOMEOFEACH LUFILE; IGNORING')
            LUfile = strrep(LUfile,'.someOfEachCrop','') ;
        elseif contains(LUfile, 'harm.forLPJG') && ~contains(LUfile, 'noMinCropFrac')
            warning('RUN USED SOMEOFEACH LUFILE; IGNORING')
            LUfile = strrep(LUfile,'.txt','.noMinCropFrac.txt') ;
        end
        fprintf('LUfile = %s\n', LUfile)

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
                LU.yearList = yearList ;
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
        if do_PLUMout_gridlist_adjust
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
            cropfile = strrep(cropfile,'.someOfEachCrop','') ;
        elseif contains(cropfile, 'harm.forLPJG') && ~contains(cropfile, 'noMinCropFrac')
            warning('RUN USED SOMEOFEACH LUFILE; IGNORING')
            cropfile = strrep(cropfile,'.txt','.noMinCropFrac.txt') ;
        end
        fprintf('cropfile = %s\n', cropfile)
        
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

        if do_PLUMout_gridlist_adjust
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
            ....* repmat(land_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ;
            .* repmat(gcel_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ; %land_gcel_fix
        if exist('LU0','var')
            cropareas0 = cropfracs ;
            cropareas0.maps_YXvy = ...
                cropfracs.maps_YXvy ...
                .* repmat(LU0.maps_YXvy(:,:,strcmp(LU0.varNames,'CROPLAND'),:),[1 1 length(cropfracs.varNames) 1]) ...
                ....* repmat(land_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ;
                .* repmat(gcel_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ; %land_gcel_fix
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
        
        if do_save.yield || do_save.yield_map
            disp('   Saving yield/cropprod...')
            if do_PLUMout_gridlist_adjust
                yield.list_to_map = PLUMout_gridlist.list_to_map ;
                yield.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(yield.varNames) 1])) = NaN ;
            end
            if do_save.yield
                for c = 1:length(cropTypes)
                    thisCrop = cropTypes{c} ;
                    thisCalibFactor = unique(calibFactors(strcmp(cropTypes_conv,thisCrop))) ;
                    eval(['cropprod_ts_' thisCrop ' = thisCalibFactor*getTS(yield,thisCrop,cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
                    save(timeseries_out,['cropprod_ts_' thisCrop],v73_or_append(timeseries_out)) ;
                    if exist('LU0','var')
                        eval(['cropprod0_ts_' thisCrop ' = thisCalibFactor*getTS(yield,thisCrop,cropareas0.maps_YXvy(:,:,strcmp(cropareas0.varNames,thisCrop),:)) ;']) ;
                        save(timeseries_out,['cropprod0_ts_' thisCrop],v73_or_append(timeseries_out)) ;
                    end
                    clear thisCrop
                end; clear c
                clear cropprod_ts*
            end
            
        end
        
        if do_save.yield_map
            yield_d1 = yield ;
            [~,~,IB] = intersect(yield_d1.varNames,cropTypes_conv,'stable') ;
            if ~isequal(shiftdim(cropTypes_conv(IB)),shiftdim(yield_d1.varNames))
                error('You screwed up an intersect(): ~isequal(shiftdim(cropTypes_conv(IB)),shiftdim(yield_d1.varNames))')
            end
            calibFactors_YXvy = repmat(permute(calibFactors(IB),[3 2 1]),...
                [size(yield_d1.maps_YXvy,1) size(yield_d1.maps_YXvy,2) 1 10]) ;
            yield_d1.maps_YXvy = yield_d1.maps_YXvy(:,:,:,Npad+(1:10)).*calibFactors_YXvy ;
            yield_d1.yearList = yield_d1.yearList(1:10) ;
            save(firstdecade_out,'yield_d1',v73_or_append(firstdecade_out)) ;
            clear yield_d1
            yield_d9 = yield ;
            yield_d9.maps_YXvy = yield_d9.maps_YXvy(:,:,:,end-9:end).*calibFactors_YXvy ;
            yield_d9.yearList = yield_d9.yearList(end-9:end) ;
            save(lastdecade_out,'yield_d9',v73_or_append(lastdecade_out)) ;
            clear yield_d9
        end
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
            disp('   Saving EXPECTED yield/cropprod...')
            if do_PLUMout_gridlist_adjust
                expyield.list_to_map = PLUMout_gridlist.list_to_map ;
                expyield.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(expyield.varNames) 1])) = NaN ;
            end
            if do_save.yield_exp
                for c = 1:length(cropTypes)
                    thisCrop = cropTypes{c} ;
                    % Do not multiply by calibration factor, because the CF
                    % that PLUM expects is already included in the PLUM yields
                    eval(['cropprodExp_ts_' thisCrop ' = getTS(expyield,thisCrop,cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
                    save(timeseries_out,['cropprodExp_ts_' thisCrop],v73_or_append(timeseries_out)) ;
                    clear thisCrop
                end; clear c
                clear cropprod_ts*
            end
            if do_save.yield_exp_map
                expyield_d1 = expyield ;
                expyield_d1.maps_YXvy = expyield_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
                expyield_d1.yearList = expyield_d1.yearList(1:10) ;
                save(firstdecade_out,'expyield_d1',v73_or_append(firstdecade_out)) ;
                clear expyield_d1
                expyield_d9 = expyield ;
                expyield_d9.maps_YXvy = expyield_d9.maps_YXvy(:,:,:,end-9:end) ;
                expyield_d9.yearList = expyield_d9.yearList(end-9:end) ;
                save(lastdecade_out,'expyield_d9',v73_or_append(lastdecade_out)) ;
                clear expyield_d9
            end
        end
    end
    clear expyield
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save LU and crop time series %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.LU
        disp('   Saving land use time series...')
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
        LUarea_ts_ntrl = getTS(LU,'NATURAL',this_area_YX) ;
        LUarea_ts_bare = getTS(LU,'BARREN',this_area_YX) ;
        LUarea_ts_crop = getTS(LU,'CROPLAND',this_area_YX) ;
        LUarea_ts_past = getTS(LU,'PASTURE',this_area_YX) ;
        save(timeseries_out,'LUarea_ts_ntrl','LUarea_ts_bare','LUarea_ts_crop','LUarea_ts_past',v73_or_append(timeseries_out)) ;
        if exist('LU0','var')
            LUarea_ts_crop0 = getTS(LU0,'CROPLAND',this_area_YX) ;
            LUarea_ts_past0 = getTS(LU0,'PASTURE',this_area_YX) ;
            save(timeseries_out,'LUarea_ts_crop0','LUarea_ts_past0',v73_or_append(timeseries_out)) ;
        end
        clear *_ts_*
    end
    if do_save.crops
        disp('   Saving crop area time series...')
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['croparea_ts_' thisCrop ' = getTS(cropareas,thisCrop,ones(size(this_area_YX))) ;']) ;
            save(timeseries_out,['croparea_ts_' thisCrop],v73_or_append(timeseries_out)) ;
            clear thisCrop
        end; clear c
        if exist('LU0','var')
            for c = 1:length(cropTypes)
                thisCrop = cropTypes{c} ;
                eval(['croparea0_ts_' thisCrop ' = getTS(cropareas0,thisCrop,ones(size(this_area_YX))) ;']) ;
                save(timeseries_out,['croparea0_ts_' thisCrop],v73_or_append(timeseries_out)) ;
                clear thisCrop
            end; clear c
        end
        clear *_ts_*
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
        if do_PLUMout_gridlist_adjust
            irrig.list_to_map = PLUMout_gridlist.list_to_map ;
            irrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(irrig.varNames) 1])) = NaN ;
            gsirrig.list_to_map = PLUMout_gridlist.list_to_map ;
            gsirrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(gsirrig.varNames) 1])) = NaN ;
        end
%         irrig_ts = getTS(irrig,'Total',land_area_YX) * cf_water ;
        irrig_ts = getTS(irrig,'Total',gcel_area_YX) * cf_water ; %land_gcel_fix
        save(timeseries_out,'irrig_ts',v73_or_append(timeseries_out)) ;
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['gsirrig_ts_' thisCrop ' = getTS(gsirrig,thisCrop,cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) * cf_water ;']) ;
            save(timeseries_out,['gsirrig_ts_' thisCrop],v73_or_append(timeseries_out)) ;
            clear thisCrop gsirrig_ts_*
        end; clear c
        
        irrig_d1 = irrig ;
        irrig_d1.maps_YXvy = irrig_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        irrig_d1.yearList = irrig_d1.yearList(1:10) ;
        if ~exist(firstdecade_out,'file')
            save(firstdecade_out,'irrig_d1') ;
        else
            save(firstdecade_out,'irrig_d1',v73_or_append(firstdecade_out)) ;
        end
        clear irrig_d1
        irrig_d9 = irrig ;
        irrig_d9.maps_YXvy = irrig_d9.maps_YXvy(:,:,:,end-9:end) ;
        irrig_d9.yearList = irrig_d9.yearList(end-9:end) ;
        if ~exist(lastdecade_out,'file')
            save(lastdecade_out,'irrig_d9') ;
        else
            save(lastdecade_out,'irrig_d9',v73_or_append(lastdecade_out)) ;
        end
        clear irrig_d9
        clear irrig
        gsirrig_d1 = gsirrig ;
        gsirrig_d1.maps_YXvy = gsirrig_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        gsirrig_d1.yearList = gsirrig_d1.yearList(1:10) ;
        save(firstdecade_out,'gsirrig_d1',v73_or_append(firstdecade_out)) ;
        clear gsirrig_d1
        gsirrig_d9 = gsirrig ;
        gsirrig_d9.maps_YXvy = gsirrig_d9.maps_YXvy(:,:,:,end-9:end) ;
        gsirrig_d9.yearList = gsirrig_d9.yearList(end-9:end) ;
        save(lastdecade_out,'gsirrig_d9',v73_or_append(lastdecade_out)) ;
        clear gsirrig_d9
        clear gsirrig
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save other annual water %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.water && anyfileexist([inDir 'awater.out'])
        awater = lpjgu_matlab_readTable_then2map([inDir 'awater.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving evaporation, transpiration, and runoff...')
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
        aevap_ts = getTS(awater,'Evap',this_area_YX) * cf_water ;
        save(timeseries_out,'aevap_ts',v73_or_append(timeseries_out)) ;
        aaet_ts = getTS(awater,'Transp',this_area_YX) * cf_water ;
        save(timeseries_out,'aaet_ts',v73_or_append(timeseries_out)) ;
        tot_runoff_ts = getTS(awater,'Runoff',this_area_YX) * cf_water ;
        save(timeseries_out,'tot_runoff_ts',v73_or_append(timeseries_out)) ;
        clear *_ts
        awater_d1 = awater ;
        awater_d1.maps_YXvy = awater_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        awater_d1.yearList = awater_d1.yearList(1:10) ;
        save(firstdecade_out,'awater_d1',v73_or_append(firstdecade_out)) ;
        clear awater_d1
        awater_d9 = awater ;
        awater_d9.maps_YXvy = awater_d9.maps_YXvy(:,:,:,end-9:end) ;
        awater_d9.yearList = awater_d9.yearList(end-9:end) ;
        save(lastdecade_out,'awater_d9',v73_or_append(lastdecade_out)) ;
        clear awater_d9
        awater_last30.list_to_map = awater.list_to_map ;
        awater_last30.years_incl = last30_yearList ;
        [~,IA] = intersect(awater.yearList, last30_yearList) ;
        if length(IA) ~= length(last30_yearList)
            error('length(IA) ~= length(last30_yearList)')
        end
        awater_last30.varNames = awater.varNames ;
        awater_last30.statHandles = last30_statHandles ;
        awater_last30.statList = last30_statList ;
        awater_last30.maps_YXvs = nan([size(land_area_YX) length(awater.varNames) Nstats],'single') ;
        for s = 1:Nstats
            awater_last30.maps_YXvs(:,:,:,s) = eval(sprintf('%s(awater.maps_YXvy(:,:,:,IA)) ;', last30_statList{s})) ;
        end; clear s
        save(last30_out,'awater_last30',v73_or_append(last30_out)) ;
        clear awater_last30
        clear awater
    elseif do_save.water
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
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
                aevap.maps_YXvy = sum(mevap.maps_YXvy .* repmat(permute(monthlengths,[4 3 2 1]),[size(land_area_YX) 1 length(yearList)]) / 365, 3) ;
            else
                error('???') ;
            end
            disp('   Saving evaporation...')
            if do_PLUMout_gridlist_adjust
                aevap.list_to_map = PLUMout_gridlist.list_to_map ;
                aevap.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aevap.varNames) 1])) = NaN ;
            end
            aevap_ts = getTS(aevap,'Total',this_area_YX) * cf_water ;
            save(timeseries_out,'aevap_ts',v73_or_append(timeseries_out)) ;
            aevap_d1 = aevap ;
            aevap_d1.maps_YXvy = aevap_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
            aevap_d1.yearList = aevap_d1.yearList(1:10) ;
            save(firstdecade_out,'aevap_d1',v73_or_append(firstdecade_out)) ;
            clear aevap_d1
            aevap_d9 = aevap ;
            aevap_d9.maps_YXvy = aevap_d9.maps_YXvy(:,:,:,end-9:end) ;
            aevap_d9.yearList = aevap_d9.yearList(end-9:end) ;
            save(lastdecade_out,'aevap_d9',v73_or_append(lastdecade_out)) ;
            clear aevap_d9
            clear aevap            
        else
            warning('Evaporation missing! Skipping.')
        end
        
        % Evapotranspiration (mm/yr)
        aaet = lpjgu_matlab_readTable_then2map([inDir 'aaet.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving aaet...')
        if do_PLUMout_gridlist_adjust
            aaet.list_to_map = PLUMout_gridlist.list_to_map ;
            aaet.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aaet.varNames) 1])) = NaN ;
        end
        aaet_ts = getTS(aaet,'Total',this_area_YX) * cf_water ;
        save(timeseries_out,'aaet_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        aaet_d1 = aaet ;
        aaet_d1.maps_YXvy = aaet_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        aaet_d1.yearList = aaet_d1.yearList(1:10) ;
        save(firstdecade_out,'aaet_d1',v73_or_append(firstdecade_out)) ;
        clear aaet_d1
        aaet_d9 = aaet ;
        aaet_d9.maps_YXvy = aaet_d9.maps_YXvy(:,:,:,end-9:end) ;
        aaet_d9.yearList = aaet_d9.yearList(end-9:end) ;
        save(lastdecade_out,'aaet_d9',v73_or_append(lastdecade_out)) ;
        clear aaet_d9
        clear aaet
        
        % Runoff (mm/yr)
        tot_runoff = lpjgu_matlab_readTable_then2map([inDir 'tot_runoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving runoff...')
        if do_PLUMout_gridlist_adjust
            tot_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            tot_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(tot_runoff.varNames) 1])) = NaN ;
        end
        tot_runoff_ts = getTS(tot_runoff,'Total',this_area_YX) * cf_water ;
        save(timeseries_out,'tot_runoff_ts',v73_or_append(timeseries_out)) ;
        clear tot_runoff_ts
        tot_runoff_d1 = tot_runoff ;
        tot_runoff_d1.maps_YXvy = tot_runoff_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        tot_runoff_d1.yearList = tot_runoff_d1.yearList(1:10) ;
        save(firstdecade_out,'tot_runoff_d1',v73_or_append(firstdecade_out)) ;
        clear tot_runoff_d1
        tot_runoff_d9 = tot_runoff ;
        tot_runoff_d9.maps_YXvy = tot_runoff_d9.maps_YXvy(:,:,:,end-9:end) ;
        tot_runoff_d9.yearList = tot_runoff_d9.yearList(end-9:end) ;
        save(lastdecade_out,'tot_runoff_d9',v73_or_append(lastdecade_out)) ;
        clear tot_runoff_d9
        clear tot_runoff
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save carbon (kgC/m2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.carbon
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
        cpool = lpjgu_matlab_readTable_then2map([inDir 'cpool.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving cpool...')
        if do_PLUMout_gridlist_adjust
            cpool.list_to_map = PLUMout_gridlist.list_to_map ;
            cpool.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cpool.varNames) 1])) = NaN ;
        end
        cpool_ts_VegC = getTS(cpool,'VegC',this_area_YX) * cf_cpool ;
        save(timeseries_out,'cpool_ts_VegC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_LitterSoilC = getTS(cpool,{'LitterC','SoilC'},this_area_YX) * cf_cpool ;
        save(timeseries_out,'cpool_ts_LitterSoilC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_HarvSlowC = getTS(cpool,'HarvSlowC',this_area_YX) * cf_cpool ;
        save(timeseries_out,'cpool_ts_HarvSlowC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_Total = getTS(cpool,'Total',this_area_YX) * cf_cpool ;
        save(timeseries_out,'cpool_ts_Total',v73_or_append(timeseries_out)) ;
        clear *_ts_*        
        cpool_d1 = cpool ;
        cpool_d1.maps_YXvy = cpool_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        cpool_d1.yearList = cpool_d1.yearList(1:10) ;
        save(firstdecade_out,'cpool_d1',v73_or_append(firstdecade_out)) ;
        clear cpool_d1
        cpool_d9 = cpool ;
        cpool_d9.maps_YXvy = cpool_d9.maps_YXvy(:,:,:,end-9:end) ;
        cpool_d9.yearList = cpool_d9.yearList(end-9:end) ;
        save(lastdecade_out,'cpool_d9',v73_or_append(lastdecade_out)) ;
        clear cpool_d9
        clear cpool
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save monthly runoff (mm/month) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.mrunoff
        mon_runoff = lpjgu_matlab_readTable_then2map([inDir 'mrunoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        if do_PLUMout_gridlist_adjust
            mon_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            mon_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(mon_runoff.varNames) 1])) = NaN ;
        end
        
        % Save peak runoff timeseries
        pkrunoff_ts_ym = nan(Nyears,12) ;
        for m = 1:12
%             pkrunoff_ts_ym(:,m) = getTS(mon_runoff,mon_runoff.varNames{m},land_area_YX) ...
            pkrunoff_ts_ym(:,m) = getTS(mon_runoff,mon_runoff.varNames{m},gcel_area_YX) ... %land_gcel_fix
                                * cf_water ;
        end
        pkrunoff_ts = max(pkrunoff_ts_ym,[],2) ;
        save(timeseries_out,'pkrunoff_ts',v73_or_append(timeseries_out)) ;
        clear pkrunoff_ts*
        
        % Save peak runoff maps
        mon_runoff_d1 = mon_runoff ;
        mon_runoff_d1.maps_YXvy = cf_water * mon_runoff_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        mon_runoff_d1.yearList = mon_runoff_d1.yearList(1:10) ;
        save(firstdecade_out,'mon_runoff_d1',v73_or_append(firstdecade_out)) ;
        clear mon_runoff_d1
        mon_runoff_d9 = mon_runoff ;
        mon_runoff_d9.maps_YXvy = cf_water * mon_runoff_d9.maps_YXvy(:,:,:,end-9:end) ;
        mon_runoff_d9.yearList = mon_runoff_d9.yearList(end-9:end) ;
        save(lastdecade_out,'mon_runoff_d9',v73_or_append(lastdecade_out)) ;
        clear mon_runoff_d9
        mon_runoff_last30.list_to_map = mon_runoff.list_to_map ;
        mon_runoff_last30.years_incl = last30_yearList ;
        [~,IA] = intersect(mon_runoff.yearList, last30_yearList) ;
        if length(IA) ~= length(last30_yearList)
            error('length(IA) ~= length(last30_yearList)')
        end
        mon_runoff_last30.varNames = 'allmonths' ;
        mon_runoff_last30.statList = last30_statList ;
        mon_runoff_last30.statHandles = last30_statHandles ;
        mon_runoff_last30.maps_YXvs = nan([size(land_area_YX) 1 Nstats],'single') ;
        mon_runoff_maps_YX1y = reshape(mon_runoff.maps_YXvy(:,:,:,IA), [size(land_area_YX) 1 length(mon_runoff_last30.years_incl)*length(mon_runoff.varNames)]) ;
        for s = 1:Nstats
            mon_runoff_last30.maps_YXvs(:,:,:,s) = eval(sprintf('%s(mon_runoff_maps_YX1y) ;', last30_statList{s})) ;
        end; clear s
        save(last30_out,'mon_runoff_last30',v73_or_append(last30_out)) ;
        clear mon_runoff_last30
        clear mon_runoff
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import (and save, if doing so) FPC %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        if do_PLUMout_gridlist_adjust
            list_to_map = PLUMout_gridlist.list_to_map ;
            fpc.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(fpc.varNames) 1])) = NaN ;
        end
        
        if do_save.fpc
            disp('   Saving FPC...')
            fpc_d1 = fpc ;
            fpc_d1.maps_YXvy = fpc_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
            fpc_d1.yearList = fpc_d1.yearList(1:10) ;
            save(firstdecade_out,'fpc_d1',v73_or_append(firstdecade_out))
            clear fpc_d1
            fpc_d9 = fpc ;
            fpc_d9.maps_YXvy = fpc_d9.maps_YXvy(:,:,:,end-9:end) ;
            fpc_d9.yearList = fpc_d9.yearList(end-9:end) ;
            save(lastdecade_out,'fpc_d9',v73_or_append(lastdecade_out))
            clear fpc_d9
        end
        if ~do_save.albedo
            clear fpc
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
        if do_PLUMout_gridlist_adjust
            snowdepth.list_to_map = PLUMout_gridlist.list_to_map ;
            snowdepth.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(snowdepth.varNames) 1])) = NaN ;
        end
        
        [albedo_jan_YXy,albedo_jul_YXy] = ...
            get_albedo(fpc, snowdepth, LU, baresoil_albedo_YX, pftList) ;
        clear fpc
        
        % Global
%         weighting_YXy = repmat(land_area_YX,[1 1 Nyears]) ;
        weighting_YXy = repmat(gcel_area_YX,[1 1 Nyears]) .* (1-squeeze(LU.maps_YXvy(:,:,strcmp(LU.varNames,'BARREN'),:))) ; %land_gcel_fix
        albedo1_ts = squeeze(nansum(nansum(albedo_jan_YXy.*(weighting_YXy / nansum(nansum(weighting_YXy(:,:,1)))), 2), 1)) ;
        albedo7_ts = squeeze(nansum(nansum(albedo_jul_YXy.*(weighting_YXy / nansum(nansum(weighting_YXy(:,:,1)))), 2), 1)) ;
        save(timeseries_out,'albedo1_ts',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'albedo7_ts',v73_or_append(timeseries_out)) ;
        
        % Boreal
        biome_borfor_YX = biomes_YX==biomes_key.Code(strcmp(biomes_key.Biome,'Boreal Forests/Taiga')) ;
        biome_borfor_YXy = repmat(biome_borfor_YX, [1 1 size(weighting_YXy,3)]) ;
        weighting_YXy_borfor = weighting_YXy ;
        weighting_YXy_borfor(~biome_borfor_YXy) = NaN ;
        albedo_jan_YXy_borfor = albedo_jan_YXy ;
        albedo_jan_YXy_borfor(~biome_borfor_YXy) = NaN ;
        albedo1_borfor_ts = squeeze(nansum(nansum(albedo_jan_YXy_borfor.*(weighting_YXy_borfor / nansum(nansum(weighting_YXy_borfor(:,:,1)))), 2), 1)) ;
        albedo_jul_YXy_borfor = albedo_jul_YXy ;
        albedo_jul_YXy_borfor(~biome_borfor_YXy) = NaN ;
        albedo7_borfor_ts = squeeze(nansum(nansum(albedo_jul_YXy_borfor.*(weighting_YXy_borfor / nansum(nansum(weighting_YXy_borfor(:,:,1)))), 2), 1)) ;
        save(timeseries_out,'albedo1_borfor_ts',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'albedo7_borfor_ts',v73_or_append(timeseries_out)) ;
        
        % Tundra
        biome_tundra_YX = biomes_YX==biomes_key.Code(strcmp(biomes_key.Biome,'Tundra')) ;
        biome_tundra_YXy = repmat(biome_tundra_YX, [1 1 size(weighting_YXy,3)]) ;
        weighting_YXy_tundra = weighting_YXy ;
        weighting_YXy_tundra(~biome_tundra_YXy) = NaN ;
        albedo_jan_YXy_tundra = albedo_jan_YXy ;
        albedo_jan_YXy_tundra(~biome_tundra_YXy) = NaN ;
        albedo1_tundra_ts = squeeze(nansum(nansum(albedo_jan_YXy_tundra.*(weighting_YXy_tundra / nansum(nansum(weighting_YXy_tundra(:,:,1)))), 2), 1)) ;
        albedo_jul_YXy_tundra = albedo_jul_YXy ;
        albedo_jul_YXy_tundra(~biome_tundra_YXy) = NaN ;
        albedo7_tundra_ts = squeeze(nansum(nansum(albedo_jul_YXy_tundra.*(weighting_YXy_tundra / nansum(nansum(weighting_YXy_tundra(:,:,1)))), 2), 1)) ;
        save(timeseries_out,'albedo1_tundra_ts',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'albedo7_tundra_ts',v73_or_append(timeseries_out)) ;
        
        clear albedo*_ts albedo_*_YXy_* weighting_YXy* biome_*_YX
        
        albedo.list_to_map = snowdepth.list_to_map ;
        albedo.varNames = {'January','July'} ;
        tmp_YXyv = nan([size(baresoil_albedo_YX) Nyears 2],'single') ;
        tmp_YXyv(:,:,:,1) = albedo_jan_YXy ;
        tmp_YXyv(:,:,:,2) = albedo_jul_YXy ;
        
        albedo_d1 = albedo ;
        tmp_YXyv1 = tmp_YXyv(:,:,1:10,:) ;
        albedo_d1.maps_YXvy = permute(tmp_YXyv1,[1 2 4 3]) ;
        albedo_d1.yearList = snowdepth.yearList(1:10) ;
        save(firstdecade_out,'albedo_d1',v73_or_append(firstdecade_out))
        clear albedo_d1 tmp_YXyv1
        albedo_d9 = albedo ;
        tmp_YXyv9 = tmp_YXyv(:,:,end-9:end,:) ;
        albedo_d9.maps_YXvy = permute(tmp_YXyv9,[1 2 4 3]) ;
        clear tmp_YXyv
        albedo_d9.yearList = snowdepth.yearList(end-9:end) ;
        save(lastdecade_out,'albedo_d9',v73_or_append(lastdecade_out))
        clear albedo_d9 tmp_YXyv9
        
        clear snowdepth albedo tmp_YXyv
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save BVOCs (isoprene, monoterpenes) (mgC/m2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.bvocs
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
        if anyfileexist([inDir 'aiso_smry.out'])
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        disp('   Saving aiso...')
        if do_PLUMout_gridlist_adjust
            aiso.list_to_map = PLUMout_gridlist.list_to_map ;
            aiso.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aiso.varNames) 1])) = NaN ;
        end
        aiso_ts = getTS(aiso,'Total',this_area_YX) * cf_bvoc ;
        save(timeseries_out,'aiso_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        aiso_d1 = aiso ;
        aiso_d1.maps_YXvy = aiso_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        aiso_d1.yearList = aiso_d1.yearList(1:10) ;
        save(firstdecade_out,'aiso_d1',v73_or_append(firstdecade_out))
        clear aiso_d1
        aiso_d9 = aiso ;
        aiso_d9.maps_YXvy = aiso_d9.maps_YXvy(:,:,:,end-9:end) ;
        aiso_d9.yearList = aiso_d9.yearList(end-9:end) ;
        save(lastdecade_out,'aiso_d9',v73_or_append(lastdecade_out))
        clear aiso_d9
        clear aiso
        
        if anyfileexist([inDir 'amon_smry.out'])
            amon = lpjgu_matlab_readTable_then2map([inDir 'amon_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            amon = lpjgu_matlab_readTable_then2map([inDir 'amon.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        disp('   Saving amon...')
        if do_PLUMout_gridlist_adjust
            amon.list_to_map = PLUMout_gridlist.list_to_map ;
            amon.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(amon.varNames) 1])) = NaN ;
        end
        amon_ts = getTS(amon,'Total',this_area_YX) * cf_bvoc ;
        save(timeseries_out,'amon_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        amon_d1 = amon ;
        amon_d1.maps_YXvy = amon_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        amon_d1.yearList = amon_d1.yearList(1:10) ;
        save(firstdecade_out,'amon_d1',v73_or_append(firstdecade_out))
        clear amon_d1
        amon_d9 = amon ;
        amon_d9.maps_YXvy = amon_d9.maps_YXvy(:,:,:,end-9:end) ;
        amon_d9.yearList = amon_d9.yearList(end-9:end) ;
        save(lastdecade_out,'amon_d9',v73_or_append(lastdecade_out))
        clear amon_d9
        clear amon
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save N flux (kgN/ha) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.Nflux
        
%         this_area_YX = land_area_YX ;
        this_area_YX = gcel_area_YX ; %land_gcel_fix
        
        % From nflux.out
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of Nflux?')
        end
        nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving nflux...')
        if do_PLUMout_gridlist_adjust
            nflux.list_to_map = PLUMout_gridlist.list_to_map ;
            nflux.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(nflux.varNames) 1])) = NaN ;
        end
        nflux_ts_fert = getTS(nflux,'fert',this_area_YX) * cf_nflux ;
        nflux_ts_flux = getTS(nflux,'flux',this_area_YX) * cf_nflux ;
        nflux_ts_leach = getTS(nflux,'leach',this_area_YX) * cf_nflux ;
        nflux_ts_harvest = getTS(nflux,'harvest',this_area_YX) * cf_nflux ;
        nflux_ts_LU_ch = getTS(nflux,'LU_ch',this_area_YX) * cf_nflux ;
        save(timeseries_out,'nflux_ts_fert',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'nflux_ts_flux',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'nflux_ts_leach',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'nflux_ts_harvest',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'nflux_ts_LU_ch',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        nflux_d1 = nflux ;
        nflux_d1.maps_YXvy = nflux_d1.maps_YXvy(:,:,:,Npad+(1:10)) * cf_nflux ;
        nflux_d1.yearList = nflux_d1.yearList(1:10) ;
        save(firstdecade_out,'nflux_d1',v73_or_append(firstdecade_out))
        clear nflux_d1
        nflux_d9 = nflux ;
        nflux_d9.maps_YXvy = nflux_d9.maps_YXvy(:,:,:,end-9:end) * cf_nflux ;
        nflux_d9.yearList = nflux_d9.yearList(end-9:end) ;
        save(lastdecade_out,'nflux_d9',v73_or_append(lastdecade_out))
        clear nflux_d9
        nflux_last30.list_to_map = nflux.list_to_map ;
        nflux_last30.years_incl = last30_yearList ;
        [~,IA] = intersect(nflux.yearList, last30_yearList) ;
        if length(IA) ~= length(last30_yearList)
            error('length(IA) ~= length(last30_yearList)')
        end
        nflux_last30.varNames = nflux.varNames ;
        nflux_last30.statList = last30_statList ;
        nflux_last30.statHandles = last30_statHandles ;
        nflux_last30.maps_YXvs = nan([size(land_area_YX) length(nflux.varNames) Nstats],'single') ;
        for s = 1:Nstats
            nflux_last30.maps_YXvs(:,:,:,s) = eval(sprintf('%s(nflux.maps_YXvy(:,:,:,IA)) ;', last30_statList{s})) ;
        end; clear s
        save(last30_out,'nflux_last30',v73_or_append(last30_out)) ;
        clear nflux_last30
        clear nflux
    end
    
    if do_save.Nfert
        % Fertilizer for each crop
        cmd = sprintf('grep -i ''param "file_nfert"'' %s/landcover.ins | grep -v -E "^\\s*!" | tail -n 1 | grep -oE "str \\".+\\"" | sed ''s/str //'' | sed ''s/"//g''', ...
            inDir) ;
        [x, NfertFile] = unix(cmd) ;
        if x~=0
            warning('Error %d when trying to find Nfert file', x)
        end
        NfertFile = regexprep(NfertFile,'[\n\r]+','') ; % Remove extraneous newline
        if ~exist(NfertFile, 'file')
            if onMac
                cmd = 'grep -i "file_nfert" %s/landcover.ins | sed ''s@param "file_nfert" (str "/project/fh1-project-lpjgpi/lr8247@/Users/Shared@'' | sed ''s@")@@''' ;
                cmd_str = sprintf(cmd,removeslashifneeded(inDir)) ;
                [x,NfertFile_tmp] = unix(cmd_str) ;
                if x~=0
                    error(['get_nfert_file.sh failed with error ' num2str(x)])
                end
            else
                NfertFile_tmp= 'grep -i "file_nfert" %s/landcover.ins | sed ''s@param "file_nfert" (str "@@'' | sed ''s@")@@''' ;
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
            elseif onMac
                NfertFile = crude_file_find(NfertFile_tmp) ;
            else
                cmd = sprintf('grep -i ''param "file_nfert"'' %s/landcover.ins | grep -v -e "^[[:blank:]]!" | sed ''s@param "file_nfert"@@'' | sed ''s@(str@@'' | sed ''s@)@@'' | sed ''s@"@@g''', ...
                    inDir) ;
                [x,NfertFile_tmp] = unix(cmd) ;
                if x~=0
                    error(['Failed when trying to find Nfert file, with error ' num2str(x)])
                end
                NfertFile_tmp = regexprep(NfertFile_tmp,'[\n\r]+','') ; % Remove extraneous newline
                NfertFile = strrep(NfertFile_tmp,' ','') ; % Remove extraneous spaces
            end
        end
        fprintf('NfertFile = %s\n', NfertFile)
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
%             Nfert = CrOp_and_CrOpi(Nfert, 'Nfert', cropTypes, 'replace') ;
            Nfert = CrOp_and_CrOpi( ...
                Nfert, 'Nfert', cropTypes, merge_or_replace, ...
                '', cropfracs_orig, inds_cropTypes_cropFracsOrig, ...
                inds_cropTypesI_cropFracsOrig, ...
                inds_cropTypes_nonCGs_cropFracsOrig) ;
        end
        % Save time series (kg/m2 --> kg)
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['nflux_ts_fert_' thisCrop ' = getTS(Nfert,thisCrop,cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
            save(timeseries_out,['nflux_ts_fert_' thisCrop],v73_or_append(timeseries_out)) ;
            clear thisCrop
        end; clear c
        if exist('LU0','var')
            for c = 1:length(cropTypes)
                thisCrop = cropTypes{c} ;
                eval(['nflux0_ts_fert_' thisCrop ' = getTS(Nfert,thisCrop,cropareas0.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
                save(timeseries_out,['nflux0_ts_fert_' thisCrop],v73_or_append(timeseries_out)) ;
                clear thisCrop
            end; clear c
        end
        clear *_ts_*
        
        % Save maps (kg/m2)
        Nfert_d1 = Nfert ;
        Nfert_d1.maps_YXvy = Nfert_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        Nfert_d1.yearList = Nfert_d1.yearList(1:10) ;
        save(firstdecade_out,'Nfert_d1',v73_or_append(firstdecade_out))
        clear Nfert_d1
        Nfert_d9 = Nfert ;
        Nfert_d9.maps_YXvy = Nfert_d9.maps_YXvy(:,:,:,end-9:end) ;
        Nfert_d9.yearList = Nfert_d9.yearList(end-9:end) ;
        save(lastdecade_out,'Nfert_d9',v73_or_append(lastdecade_out))
        clear Nfert_d9
        clear Nfert
        
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%
%     %%% Save FPC stuff %%%
%     %%%%%%%%%%%%%%%%%%%%%%
% 
%     if do_save.fpc
%         disp('   Saving FPC...')
%         fpc_d1 = fpc ;
%         fpc_d1.maps_YXvy = fpc_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
%         fpc_d1.yearList = fpc_d1.yearList(1:10) ;
%         save(firstdecade_out,'fpc_d1',v73_or_append(firstdecade_out))
%         clear fpc_d1
%         fpc_d9 = fpc ;
%         fpc_d9.maps_YXvy = fpc_d9.maps_YXvy(:,:,:,end-9:end) ;
%         fpc_d9.yearList = fpc_d9.yearList(end-9:end) ;
%         save(lastdecade_out,'fpc_d9',v73_or_append(lastdecade_out))
%         clear fpc_d9
%         clear fpc
%     end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save land use and crop maps %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.crops
        disp('   Saving crop maps...')
        cropfracs_d1 = cropfracs ;
        cropfracs_d1.maps_YXvy = cropfracs_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        cropfracs_d1.yearList = cropfracs_d1.yearList(1:10) ;
        save(firstdecade_out,'cropfracs_d1',v73_or_append(firstdecade_out))
        clear cropfracs_d1
        cropfracs_d9 = cropfracs ;
        cropfracs_d9.maps_YXvy = cropfracs_d9.maps_YXvy(:,:,:,end-9:end) ;
        cropfracs_d9.yearList = cropfracs_d9.yearList(end-9:end) ;
        save(lastdecade_out,'cropfracs_d9',v73_or_append(lastdecade_out))
        clear cropfracs_d9
    end
    if do_save.LU
        disp('   Saving land use maps...')
        LUorder = nan(1,length(LUlist),'single') ;
        for L = 1:length(LU.varNames)
            thisLU = LUlist{L} ;
            LUorder(L) = find(strcmp(LU.varNames,thisLU)) ;
        end
        LU.maps_YXvy = LU.maps_YXvy(:,:,LUorder,:) ;
        LU.varNames = LUlist ;
        LU_d1 = LU ;
        LU_d1.maps_YXvy = LU_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        LU_d1.yearList = LU_d1.yearList(1:10) ;
        save(firstdecade_out,'LU_d1',v73_or_append(firstdecade_out))
        clear LU_d1
        LU_d9 = LU ;
        LU_d9.maps_YXvy = LU_d9.maps_YXvy(:,:,:,end-9:end) ;
        LU_d9.yearList = LU_d9.yearList(end-9:end) ;
        save(lastdecade_out,'LU_d9',v73_or_append(lastdecade_out))
        clear LU_d9
    end
    if do_save.LU && exist('LU0','var')
        disp('   Saving land use maps (from before moving extra CROP to PAST)...')
        LUorder = nan(1,length(LUlist),'single') ;
        for L = 1:length(LU0.varNames)
            thisLU0 = LUlist{L} ;
            LUorder(L) = find(strcmp(LU0.varNames,thisLU0)) ;
        end
        LU0.maps_YXvy = LU0.maps_YXvy(:,:,LUorder,:) ;
        LU0.varNames = LUlist ;
        LU0_d1 = LU0 ;
        LU0_d1.maps_YXvy = LU0_d1.maps_YXvy(:,:,:,Npad+(1:10)) ;
        LU0_d1.yearList = LU0_d1.yearList(1:10) ;
        save(firstdecade_out,'LU0_d1',v73_or_append(firstdecade_out))
        clear LU0_d1
        LU0_d9 = LU0 ;
        LU0_d9.maps_YXvy = LU0_d9.maps_YXvy(:,:,:,end-9:end) ;
        LU0_d9.yearList = LU0_d9.yearList(end-9:end) ;
        save(lastdecade_out,'LU0_d9',v73_or_append(lastdecade_out))
        clear LU0_d9
    end
    clear LU LU0 cropfracs* cropareas*
    
    
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

disp('All done!')
