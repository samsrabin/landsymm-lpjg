%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Overlay nflux.out vs. input data               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


inDir_list = {...
    %     'LPJGPLUM_1850-2010_remap6/output-2018-11-03-234931' ;
    %     'LPJGPLUM_1850-2010_remap6_testNoFertOnFake/output-2018-11-28-182230' ;
    %     'LPJGPLUM_2011-2100_harm2_SSP1_RCP45/output-2018-11-05-085615' ;
    %     'LPJGPLUM_2011-2100_harm2_SSP3_RCP60/output-2018-11-05-071233' ;
    %     'LPJGPLUM_2011-2100_harm2_SSP4_RCP60/output-2018-11-05-071814' ;
    %     'LPJGPLUM_2011-2100_harm2_SSP5_RCP85/output-2018-11-05-110842' ;
    'LPJGPLUM_1850-2010_remap6_testNoFertOnFake_ExtraCropFix/output-2018-11-29-195751' ;
    'LPJGPLUM_2011-2100_harm2_SSP1_RCP45_testExtraCropFix/output-2018-11-29-213540' ;
    'LPJGPLUM_2011-2100_harm2_SSP3_RCP60_testExtraCropFix/output-2018-11-29-213508' ;
    'LPJGPLUM_2011-2100_harm2_SSP4_RCP60_testExtraCropFix/output-2018-11-29-223009' ;
    'LPJGPLUM_2011-2100_harm2_SSP5_RCP85_testExtraCropFix/output-2018-11-29-223921' ;
    } ;



%% Setup

addpath(genpath(landsymm_lpjg_path()))

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

% Read PLUMout_gridlist
gridlist_file = '/Users/Shared/lpj-guess/gridlists/PLUMout_gridlist.txt' ;
PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;

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
land_area_YX = land_area_YX*1e6 ;
gcel_area_YX = gcel_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Conversion factors
%%% All masses in kg
%%% All areas in m2
%%% All volumes in m3
cf_cpool = 1 ;   % LPJ-GUESS outputs kgC/m2
cf_water = 1e-3 ; % LPJ-GUESS outputs mm (will be multiplied by land_area_m2 --> m3)
cf_nflux = 1e-4 ;   % LPJ-GUESS outputs kgN/ha
cf_bvoc = 1e-6 ;   % LPJ-GUESS outputs mgC/m2

ts_fert_bl_good = zeros(length(1850:2010),1) ;
ts_fert_yr_good = zeros(length(2011:2100),4) ;
ts_fert_yr_bad_land = nan(length(2011:2100),4) ;
ts_fert_yr_bad_crop = nan(length(2011:2100),4) ;


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
        end
        
    end
    % Do expected yields exist? If not, make sure do_save.yield_exp* are
    % false. Reset back to true at end of this for loop (before going to
    % next directory)
    flip_do_save_yield_exp = false ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get important info %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Reading yield...') ; tic
    yield = lpjgu_matlab_readTable_then2map([inDir 'yield.out'],'force_mat_save',true) ;
    disp(toc_hms(toc)) ;
    
    disp('Processing yield...') ; tic
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
    clear yield
    disp(toc_hms(toc)) ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Import land use %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    if contains(LUfile,'someOfEachCrop')
        warning('RUN USED SOMEOFEACH LUFILE; IGNORING')
        LUfile = strrep(LUfile,'.someOfEachCrop','') ;
    end
    
    disp('Reading LU...') ; tic
    LU = lpjgu_matlab_readTable_then2map(LUfile,'force_mat_save',true) ;
    disp(toc_hms(toc)) ;
    
    disp('Processing LU...') ; tic
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
    disp(toc_hms(toc)) ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import crop fractions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if contains(cropfile,'someOfEachCrop')
        warning('RUN USED SOMEOFEACH CROPFILE; IGNORING')
        % % % % % %         LUfile = strrep(cropfile,'.someOfEachCrop','') ;
        cropfile = strrep(cropfile,'.someOfEachCrop','') ;
    end
    
    disp('Reading cropfracs...') ; tic
    cropfracs = lpjgu_matlab_readTable_then2map(cropfile,'force_mat_save',true) ;
    disp(toc_hms(toc)) ;
    
    disp('Processing cropfracs...') ; tic
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
    disp(toc_hms(toc)) ;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save N flux (kgN/ha) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % From nflux.out
    % THIS IS WHERE THE PROBLEM IS
    if exist('LU0','var')
        warning('Should you do an "LU0" adjusted version of Nflux?')
    end
    
    disp('Reading nflux...') ; tic
    nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
    disp(toc_hms(toc)) ;
    
    disp('Processing nflux...') ; tic
    if is_baseline
        nflux.list_to_map = PLUMout_gridlist.list_to_map ;
        nflux.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(nflux.varNames) 1])) = NaN ;
    end
    
    thisVar = 'fert' ;
    %         thisVar = 'seed' ;
    if is_baseline
        ts_fert_bl_bad_land = getTS(nflux,thisVar,land_area_YX) * cf_nflux ;
        ts_fert_bl_bad_crop = getTS(nflux,thisVar,sum(cropareas.maps_YXvy(:,:,:,end),3)) * cf_nflux ;
    else
        ts_fert_yr_bad_land(:,d-1) = getTS(nflux,thisVar,land_area_YX) * cf_nflux ;
        ts_fert_yr_bad_crop(:,d-1) = getTS(nflux,thisVar,sum(cropareas.maps_YXvy(:,:,:,end),3)) * cf_nflux ;
    end
    clear nflux
    disp(toc_hms(toc)) ;
    
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
    
    disp('Reading Nfert...') ; tic
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
        if is_baseline
            eval(['ts_fert_bl_good = ts_fert_bl_good + nflux_ts_fert_' thisCrop ' ;']) ;
        else
            eval(['ts_fert_yr_good(:,d-1) = ts_fert_yr_good(:,d-1) + nflux_ts_fert_' thisCrop ' ;']) ;
        end
        clear thisCrop
    end; clear c
    clear Nfert
    disp(toc_hms(toc)) ;
    
    
    disp('Done.')
    
end


%% Make figure

lineWidth = 2 ;
fontSize = 14 ;

figure('Color','w','Position',figurePos) ;
plot(1850:2010,ts_fert_bl_good,'-k','LineWidth',1) ;
hold on
% plot(1850:2010,-ts_fert_bl_bad_crop,'-k','LineWidth',lineWidth) ;
plot(1850:2010,-ts_fert_bl_bad_land,'--k','LineWidth',lineWidth) ;
plot(2011:2100,ts_fert_yr_good,'-','LineWidth',1) ;
set(gca,'ColorOrderIndex',1) ;
plot(2011:2100,-ts_fert_yr_bad_land,'--','LineWidth',lineWidth) ;
% set(gca,'ColorOrderIndex',1) ;
% plot(2011:2100,-ts_fert_yr_bad_crop,'-','LineWidth',lineWidth) ;
hold off
xlabel('Year')
ylabel('kg N')
set(gca,'FontSize',fontSize) ;


