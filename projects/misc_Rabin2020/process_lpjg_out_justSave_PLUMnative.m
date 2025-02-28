%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_cropfracs_20170108 = false ;

do_save.LU              = false ;
do_save.crops           = false ;
do_save.irrig           = false ;
do_save.water           = false ;
do_save.carbon          = false ;
do_save.mrunoff         = false ;
do_save.albedo          = false ;
do_save.bvocs           = false ;
do_save.Nflux           = true ;
do_save.yield           = false ;
do_save.yield_exp       = false ;
do_save.yield_map       = false ;
do_save.yield_exp_map   = false ;
do_save.fpc             = false ;

% inDir_list = {...
% %     'PLUM2LPJG_SSP1_RCP45_v3s1/output-2018-04-23-145614' ;
% %     'PLUM2LPJG_SSP3_RCP60_v3s1/output-2018-04-23-145614' ;
% %     'PLUM2LPJG_SSP4_RCP60_v3s1/output-2018-04-23-145614' ;
% %     'PLUM2LPJG_SSP5_RCP85_v3s1/output-2018-04-23-145614' ;
% %     'LPJGPLUM_1850-2010_PLUM6xtra/matlab_merge_20180424050342' ;
% %     'PLUM2LPJG_SSP1_RCP45_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-134335' ;
% %     'PLUM2LPJG_SSP3_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-131116' ;
% %     'PLUM2LPJG_SSP4_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-134415' ;
% %     'PLUM2LPJG_SSP5_RCP85_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-24-132213' ;
%     } ;
% Ncrops = 6 ;

inDir_list = {...
    'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180502104840' ;
    'PLUM2LPJG_SSP1_RCP45_v3s1_v20180426/output-2018-05-01-164708' ;
    'PLUM2LPJG_SSP3_RCP60_v3s1_v20180426/output-2018-04-30-125214' ;
    'PLUM2LPJG_SSP4_RCP60_v3s1_v20180426/output-2018-04-30-125218' ;
    'PLUM2LPJG_SSP5_RCP85_v3s1_v20180426/output-2018-05-01-024615' ;
    } ;
Ncrops = 7 ;


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
            'CC3G_ic','CC4G_ic','Total','Crop_sum','Pasture_sum',...
            'Natural_sum','Barren_sum'} ;
LUlist = {'CROPLAND','PASTURE','NATURAL','BARREN'} ;

% Anonymous function for checking whether any version of a file exists
anyfileexist = @(in_file) ...
    exist(in_file,'file') ...
    || exist([in_file '.maps.mat'],'file') ...
    || exist([in_file '.mat'],'file') ...
    || exist([in_file '.gz'],'file') ;

% Get CFTs and PFTs
if Ncrops==6
    cropTypes = {'CerealsC3','Oilcrops',...
                 'StarchyRoots','Pulses','CerealsC4',...
                 'Rice'} ;
    LPJGcrops_2_PLUM = readtable([thisDir 'PLUM6xtra_calib_20180423.csv']) ;
    cropTypes_conv = LPJGcrops_2_PLUM.Crop ;
    calibFactors = LPJGcrops_2_PLUM.calibFactor ;
    clear LPJGcrops_2_PLUM
elseif Ncrops==7
    cropTypes = {'CerealsC3','Oilcrops',...
                 'StarchyRoots','Pulses','CerealsC4',...
                 'Rice','Miscanthus'} ;
    LPJGcrops_2_PLUM = readtable([thisDir 'PLUM6xtra_calib_20180423.csv']) ;
    cropTypes_conv = LPJGcrops_2_PLUM.Crop ;
    calibFactors = LPJGcrops_2_PLUM.calibFactor ;
    clear LPJGcrops_2_PLUM
else
    error(['Rework CFT/PFT code for Ncrops==' num2str(Ncrops)]) ;
end
pftList = [pftList_noCrops cropTypes] ;

% Import bare-soil albedo
if do_save.albedo
    baresoil_albedo_file = [thisDir 'soilmap.txt'] ;
    baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true,'verbose',false) ;
    baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
    clear baresoil_albedo
end

% Read PLUMout_gridlist
if any(is_baseline_list)
    gridlist_file = '/Users/Shared/lpj-guess/gridlists/PLUMout_gridlist.txt' ;
    PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;
end

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
if any(is_baseline_list)
    land_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
end
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Conversion factors
cf_cpool = 1e-3*1e-9 ;
cf_water = 1e-3*1e-3*1e-3 ;
cf_nflux = 1e-9 ;
cf_bvoc = 1e-15 ;
cf_lu = 1e-6 ;   % km2 to Mkm2



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
        replace_CrOp_with_CrOpi = false ;
        merge_CrOpi_into_CrOp = true ;
        if ~exist('PLUMout_mask_YX1y','var')
            PLUMout_mask_YX1y = repmat(PLUMout_gridlist.mask_YX,[1 1 1 length(yearList)]) ;
        end
    else
        yearList = 2011:2100 ;
        replace_CrOp_with_CrOpi = true ;
        merge_CrOpi_into_CrOp = false ;
    end
    yearList = transpose(yearList) ;
    Nyears = length(yearList) ;
        
    % Get land use and crop fractions input files
    clear LUfile cropfile
    if strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_2006-2010_PLUM6xtra/output-2018-04-21-061836') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_PLUM6xtra/matlab_merge_20180424050342/')
        LUfile = '/Users/Shared/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
        cropfile = '/Users/Shared/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.NfertEmpties0-0200-1000.txt' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'LPJGPLUM_1850-2010_PLUM6xtraMisc/matlab_merge_20180502104840')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/LU_xtraCROPtoPAST.remapv2.20180214.m0.txt' ;
        cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/remaps_v2/cropfracs.remapv2.20180214.m0.Misc0s.txt' ;
    else
        [x,LUfile_tmp] = unix([thisDir 'get_lu_file.sh ' inDir]) ;
        if x~=0
            error(['get_lu_file.sh failed with error ' num2str(x)])
        end
        LUfile_tmp = regexprep(LUfile_tmp,'[\n\r]+','') ; % Remove extraneous newline
        LUfile = find_PLUM2LPJG_input_file(LUfile_tmp) ;
        [x,cropfile_tmp] = unix([thisDir 'get_cropfrac_file.sh ' inDir]) ;
        if x~=0
            error(['get_cropfrac_file.sh failed with error ' num2str(x)])
        end
        cropfile_tmp = regexprep(cropfile_tmp,'[\n\r]+','') ; % Remove extraneous newline
        cropfile = find_PLUM2LPJG_input_file(cropfile_tmp) ;
    end
    
    % Do expected yields exist? If not, make sure do_save.yield_exp* are
    % false. Reset back to true at end of this for loop (before going to
    % next directory)
    flip_do_save_yield_exp = false ;
    if do_save.yield_exp && ~exist('LUcropDir_plum7','var')
        warning('No expected yields exist; skipping do_save.yield_exp.')
        flip_do_save_yield_exp = true ;
        do_save.yield_exp = ~do_save.yield_exp ;
    end
    flip_do_save_yield_exp_map = false ;
    if do_save.yield_exp_map && ~exist('LUcropDir_plum7','var')
        warning('No expected yields exist; skipping do_save.yield_exp_map.')
        flip_do_save_yield_exp_map = true ;
        do_save.yield_exp_map = ~do_save.yield_exp_map ;
    end
    
    % Setup if processing expected yields
    if do_save.yield_exp
        if ~exist('expyieldfile_plum7','var')
            expyieldfile_plum7 = [LUcropDir_plum7 'yield.txt'] ;
        end
        if ~exist(expyieldfile_plum7,'file')
            warning('Expected yield file not found! Setting do_save.yield_exp to false.')
            do_save.yield_exp = false ;
        end
    end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get important info %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yield = lpjgu_matlab_readTable_then2map([inDir 'yield.out'],'force_mat_save',true) ;
    list_to_map = yield.list_to_map ;
    if min(yield.yearList) == min(yearList) - 5 && max(yield.yearList)==max(yearList)
        warning('Adjusting yearList to account for 5 years'' padding at beginning.')
        yearList = yield.yearList ;
        Nyears = length(yearList) ;
    elseif ~isequal(yield.yearList,yearList)
        error('Rework so yearLists match!')
    end
    if ~(do_save.yield || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map)
        clear yield
    end
    
    
    %%%%%%%%%%%%%%%%%%
    %%% Import FPC %%%
    %%%%%%%%%%%%%%%%%%
    
    if do_save.fpc || do_save.albedo
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of FPC?')
        end
        fpc = lpjgu_matlab_readTable_then2map([inDir 'fpc.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('Processing FPC...')
        [fpc, ~, ~] = CrOp_and_CrOpi(fpc, 'fpc', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp) ;
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Import land use %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    cropfracs = lpjgu_matlab_readTable_then2map(cropfile,'force_mat_save',true) ;
    cropfracs = adjust_cropinput_yearLists(cropfracs, yearList) ;
    cropfracs_orig = cropfracs ;
    [cropfracs, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig] = ...
        CrOp_and_CrOpi(cropfracs, 'cropfracs', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, cropfile) ;
    if is_baseline
        cropfracs.list_to_map = PLUMout_gridlist.list_to_map ;
        cropfracs.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs.varNames) 1])) = NaN ;
        if exist('cropfracs_orig','var')
            cropfracs_orig.list_to_map = PLUMout_gridlist.list_to_map ;
            cropfracs_orig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs_orig.varNames) 1])) = NaN ;
        end
    end
    % Get area of each crop (km2)
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save yield (kgDM/m2 --> kgDM) (i.e., actually PRODUCTION) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_save.yield || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map
        [yield, ~, ~] = CrOp_and_CrOpi(yield, 'yield', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, ...
            '', cropfracs_orig, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig) ;
        yield.maps_YXvy(:,:,strcmp(yield.varNames,'CC3G_ic')|strcmp(yield.varNames,'CC4G_ic'),:) = [] ;
        tmp = yield.varNames ;
        tmp(strcmp(tmp,'CC3G_ic')) = [] ;
        tmp(strcmp(tmp,'CC4G_ic')) = [] ;
        yield.varNames = tmp ;
        
        if do_save.yield || do_save.yield_map
            disp('   Saving yield/cropprod...')
            if is_baseline
                yield.list_to_map = PLUMout_gridlist.list_to_map ;
                yield.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(yield.varNames) 1])) = NaN ;
            end
            if do_save.yield
%                 disp('Regular:')
                for c = 1:length(cropTypes)
                    thisCrop = cropTypes{c} ;
                    thisCalibFactor = unique(calibFactors(strcmp(cropTypes_conv,thisCrop))) ;
                    eval(['cropprod_ts_' thisCrop ' = thisCalibFactor*getTS(yield,thisCrop,1e6*cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
                    save(timeseries_out,['cropprod_ts_' thisCrop],v73_or_append(timeseries_out)) ;
                    if exist('LU0','var')
                        eval(['cropprod0_ts_' thisCrop ' = thisCalibFactor*getTS(yield,thisCrop,1e6*cropareas0.maps_YXvy(:,:,strcmp(cropareas0.varNames,thisCrop),:)) ;']) ;
                        save(timeseries_out,['cropprod0_ts_' thisCrop],v73_or_append(timeseries_out)) ;
                    end
                    clear thisCrop
                end; clear c
                clear cropprod_ts*
            end
            
        end

        if do_save.yield_map
            yield_d1 = yield ;
            yield_d1.maps_YXvy = yield_d1.maps_YXvy(:,:,:,1:10) ;
            yield_d1.yearList = yield_d1.yearList(1:10) ;
            save(firstdecade_out,'yield_d1',v73_or_append(firstdecade_out)) ;
            clear yield_d1
            yield_d9 = yield ;
            yield_d9.maps_YXvy = yield_d9.maps_YXvy(:,:,:,end-9:end) ;
            yield_d9.yearList = yield_d9.yearList(end-9:end) ;
            save(lastdecade_out,'yield_d9',v73_or_append(lastdecade_out)) ;
            clear yield_d9
        end
        
    end
    clear yield
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save LU and crop time series %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.LU
        LUarea_ts_ntrl = getTS(LU,'NATURAL',land_area_YX)*cf_lu ;
        LUarea_ts_bare = getTS(LU,'BARREN',land_area_YX)*cf_lu ;
        LUarea_ts_crop = getTS(LU,'CROPLAND',land_area_YX)*cf_lu ;
        LUarea_ts_past = getTS(LU,'PASTURE',land_area_YX)*cf_lu ;
        save(timeseries_out,'LUarea_ts_ntrl','LUarea_ts_bare','LUarea_ts_crop','LUarea_ts_past',v73_or_append(timeseries_out)) ;
        if exist('LU0','var')
            LUarea_ts_crop0 = getTS(LU0,'CROPLAND',land_area_YX)*cf_lu ;
            LUarea_ts_past0 = getTS(LU0,'PASTURE',land_area_YX)*cf_lu ;
            save(timeseries_out,'LUarea_ts_crop0','LUarea_ts_past0',v73_or_append(timeseries_out)) ;
        end
        clear *_ts_*
    end
    if do_save.crops
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['croparea_ts_' thisCrop ' = getTS(cropareas,thisCrop,ones(size(land_area_YX)))*cf_lu ;']) ;
            save(timeseries_out,['croparea_ts_' thisCrop],v73_or_append(timeseries_out)) ;
            clear thisCrop
        end; clear c
        if exist('LU0','var')
            for c = 1:length(cropTypes)
                thisCrop = cropTypes{c} ;
                eval(['croparea0_ts_' thisCrop ' = getTS(cropareas0,thisCrop,ones(size(land_area_YX)))*cf_lu ;']) ;
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
        [gsirrig, ~, ~] = CrOp_and_CrOpi(gsirrig, 'gsirrig', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, ...
            '', cropfracs_orig, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig) ;
        gsirrig.maps_YXvy(:,:,strcmp(gsirrig.varNames,'CC3G_ic'),:) = [] ;
        gsirrig.varNames(strcmp(gsirrig.varNames,'CC3G_ic')) = [] ;
        gsirrig.maps_YXvy(:,:,strcmp(gsirrig.varNames,'CC4G_ic'),:) = [] ;
        gsirrig.varNames(strcmp(gsirrig.varNames,'CC4G_ic')) = [] ;
        if is_baseline
            irrig.list_to_map = PLUMout_gridlist.list_to_map ;
            irrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(irrig.varNames) 1])) = NaN ;
            gsirrig.list_to_map = PLUMout_gridlist.list_to_map ;
            gsirrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(gsirrig.varNames) 1])) = NaN ;
        end
        irrig_ts = getTS(irrig,'Total',land_area_YX) * cf_water ;
        save(timeseries_out,'irrig_ts',v73_or_append(timeseries_out)) ;
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['gsirrig_ts_' thisCrop ' = getTS(gsirrig,thisCrop,cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) * cf_water ;']) ;
            save(timeseries_out,['gsirrig_ts_' thisCrop],v73_or_append(timeseries_out)) ;
            clear thisCrop gsirrig_ts_*
        end; clear c
        
        irrig_d1 = irrig ;
        irrig_d1.maps_YXvy = irrig_d1.maps_YXvy(:,:,:,1:10) ;
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
        gsirrig_d1.maps_YXvy = gsirrig_d1.maps_YXvy(:,:,:,1:10) ;
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
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of annual water?')
        end
        awater = lpjgu_matlab_readTable_then2map([inDir 'awater.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving evaporation, transpiration, and runoff...')
        aevap_ts = getTS(awater,'Evap',land_area_YX) * cf_water ;
        save(timeseries_out,'aevap_ts',v73_or_append(timeseries_out)) ;
        aaet_ts = getTS(awater,'Transp',land_area_YX) * cf_water ;
        save(timeseries_out,'aaet_ts',v73_or_append(timeseries_out)) ;
        tot_runoff_ts = getTS(awater,'Runoff',land_area_YX) * cf_water ;
        save(timeseries_out,'tot_runoff_ts',v73_or_append(timeseries_out)) ;
        clear *_ts
        awater_d1 = awater ;
        awater_d1.maps_YXvy = awater_d1.maps_YXvy(:,:,:,1:10) ;
        awater_d1.yearList = awater_d1.yearList(1:10) ;
        save(firstdecade_out,'awater_d1',v73_or_append(firstdecade_out)) ;
        clear awater_d1
        awater_d9 = awater ;
        awater_d9.maps_YXvy = awater_d9.maps_YXvy(:,:,:,end-9:end) ;
        awater_d9.yearList = awater_d9.yearList(end-9:end) ;
        save(lastdecade_out,'awater_d9',v73_or_append(lastdecade_out)) ;
        clear awater_d9
        clear awater
    elseif do_save.water
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of annual water?')
        end
        % Evaporation
        if anyfileexist([inDir 'aevap.out'])
            aevap = lpjgu_matlab_readTable_then2map([inDir 'aevap.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
            disp('   Saving evaporation...')
            if is_baseline
                aevap.list_to_map = PLUMout_gridlist.list_to_map ;
                aevap.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aevap.varNames) 1])) = NaN ;
            end
            aevap_ts = getTS(aevap,'Total',land_area_YX) * cf_water ;
            save(timeseries_out,'aevap_ts',v73_or_append(timeseries_out)) ;
            aevap_d1 = aevap ;
            aevap_d1.maps_YXvy = aevap_d1.maps_YXvy(:,:,:,1:10) ;
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
        if is_baseline
            aaet.list_to_map = PLUMout_gridlist.list_to_map ;
            aaet.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aaet.varNames) 1])) = NaN ;
        end
        aaet_ts = getTS(aaet,'Total',land_area_YX) * cf_water ;
        save(timeseries_out,'aaet_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        aaet_d1 = aaet ;
        aaet_d1.maps_YXvy = aaet_d1.maps_YXvy(:,:,:,1:10) ;
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
        if is_baseline
            tot_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            tot_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(tot_runoff.varNames) 1])) = NaN ;
        end
        tot_runoff_ts = getTS(tot_runoff,'Total',land_area_YX) * cf_water ;
        save(timeseries_out,'tot_runoff_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        tot_runoff_d1 = tot_runoff ;
        tot_runoff_d1.maps_YXvy = tot_runoff_d1.maps_YXvy(:,:,:,1:10) ;
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
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of carbon?')
        end
        cpool = lpjgu_matlab_readTable_then2map([inDir 'cpool.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving cpool...')
        if is_baseline
            cpool.list_to_map = PLUMout_gridlist.list_to_map ;
            cpool.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cpool.varNames) 1])) = NaN ;
        end
        % cpool_table = lpjgu_matlab_readTable([inDir 'cpool.out'],'do_save_mat',true) ;
        cpool_ts_VegC = getTS(cpool,'VegC',land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_VegC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_LitterSoilC = getTS(cpool,{'LitterC','SoilC'},land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_LitterSoilC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_HarvSlowC = getTS(cpool,'HarvSlowC',land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_HarvSlowC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_Total = getTS(cpool,'Total',land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_Total',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_d1 = cpool ;
        cpool_d1.maps_YXvy = cpool_d1.maps_YXvy(:,:,:,1:10) ;
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
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of monthly runoff?')
        end
        mon_runoff = lpjgu_matlab_readTable_then2map([inDir 'mrunoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        if is_baseline
            mon_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            mon_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(mon_runoff.varNames) 1])) = NaN ;
        end
        mon_runoff_d1 = mon_runoff ;
        mon_runoff_d1.maps_YXvy = mon_runoff_d1.maps_YXvy(:,:,:,1:10) ;
        mon_runoff_d1.yearList = mon_runoff_d1.yearList(1:10) ;
        save(firstdecade_out,'mon_runoff_d1',v73_or_append(firstdecade_out)) ;
        clear mon_runoff_d1
        mon_runoff_d9 = mon_runoff ;
        mon_runoff_d9.maps_YXvy = mon_runoff_d9.maps_YXvy(:,:,:,end-9:end) ;
        mon_runoff_d9.yearList = mon_runoff_d9.yearList(end-9:end) ;
        save(lastdecade_out,'mon_runoff_d9',v73_or_append(lastdecade_out)) ;
        clear mon_runoff_d9
        clear mon_runoff
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
        [albedo_jan_YXy,albedo_jul_YXy] = ...
            get_albedo(fpc, snowdepth, LU, cropfracs, baresoil_albedo_YX, land_area_YX, cropTypes, pftList) ;
        albedo1_ts = squeeze(nansum(nansum(albedo_jan_YXy.*(repmat(land_area_YX,[1 1 size(albedo_jan_YXy,3)]) / nansum(nansum(land_area_YX))), 2), 1)) ;
        albedo7_ts = squeeze(nansum(nansum(albedo_jul_YXy.*(repmat(land_area_YX,[1 1 size(albedo_jan_YXy,3)]) / nansum(nansum(land_area_YX))), 2), 1)) ;
        save(timeseries_out,'albedo1_ts',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'albedo7_ts',v73_or_append(timeseries_out)) ;
        clear albedo*_ts
        albedo.list_to_map = snowdepth.list_to_map ;
        albedo.varNames = {'January','July'} ;
        tmp_YXyv = nan([size(baresoil_albedo_YX) Nyears 2],'single') ;
        tmp_YXyv(:,:,:,1) = albedo_jan_YXy ;
        tmp_YXyv(:,:,:,2) = albedo_jul_YXy ;
        
        albedo_d1 = albedo ;
        tmp_YXyv1 = tmp_YXyv(:,:,1:10,:) ;
        albedo_d1.maps_YXvy = permute(tmp_YXyv1,[1 2 4 3]) ;
        clear tmp_YXyv
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
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of BVOCs?')
        end
        if anyfileexist([inDir 'aiso_smry.out'])
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        disp('   Saving aiso...')
        if is_baseline
            aiso.list_to_map = PLUMout_gridlist.list_to_map ;
            aiso.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aiso.varNames) 1])) = NaN ;
        end
        aiso_ts = getTS(aiso,'Total',land_area_YX_m2) * cf_bvoc ;
        save(timeseries_out,'aiso_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        aiso_d1 = aiso ;
        aiso_d1.maps_YXvy = aiso_d1.maps_YXvy(:,:,:,1:10) ;
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
        if is_baseline
            amon.list_to_map = PLUMout_gridlist.list_to_map ;
            amon.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(amon.varNames) 1])) = NaN ;
        end
        amon_ts = getTS(amon,'Total',land_area_YX_m2) * cf_bvoc ;
        save(timeseries_out,'amon_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        amon_d1 = amon ;
        amon_d1.maps_YXvy = amon_d1.maps_YXvy(:,:,:,1:10) ;
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
        
%         % From nflux.out
%         if exist('LU0','var')
%             warning('Should you do an "LU0" adjusted version of Nflux?')
%         end
%         nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
%         disp('   Saving nflux...')
%         if is_baseline
%             nflux.list_to_map = PLUMout_gridlist.list_to_map ;
%             nflux.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(nflux.varNames) 1])) = NaN ;
%         end
%         nflux_ts_fert = getTS(nflux,'fert',land_area_YX_m2*1e-4) * cf_nflux ;
%         save(timeseries_out,'nflux_ts_fert',v73_or_append(timeseries_out)) ;
%         clear *_ts_*
%         nflux_ts_flux = getTS(nflux,'flux',land_area_YX_m2*1e-4) * cf_nflux ;
%         save(timeseries_out,'nflux_ts_flux',v73_or_append(timeseries_out)) ;
%         clear *_ts_*
%         nflux_ts_leach = getTS(nflux,'leach',land_area_YX_m2*1e-4) * cf_nflux ;
%         save(timeseries_out,'nflux_ts_leach',v73_or_append(timeseries_out)) ;
%         clear *_ts_*
%         nflux_ts_harvest = getTS(nflux,'harvest',land_area_YX_m2*1e-4) * cf_nflux ;
%         save(timeseries_out,'nflux_ts_harvest',v73_or_append(timeseries_out)) ;
%         clear *_ts_*
%         nflux_ts_LU_ch = getTS(nflux,'LU_ch',land_area_YX_m2*1e-4) * cf_nflux ;
%         save(timeseries_out,'nflux_ts_LU_ch',v73_or_append(timeseries_out)) ;
%         clear *_ts_*
%         nflux_d1 = nflux ;
%         nflux_d1.maps_YXvy = nflux_d1.maps_YXvy(:,:,:,1:10) ;
%         nflux_d1.yearList = nflux_d1.yearList(1:10) ;
%         save(firstdecade_out,'nflux_d1',v73_or_append(firstdecade_out))
%         clear nflux_d1
%         nflux_d9 = nflux ;
%         nflux_d9.maps_YXvy = nflux_d9.maps_YXvy(:,:,:,end-9:end) ;
%         nflux_d9.yearList = nflux_d9.yearList(end-9:end) ;
%         save(lastdecade_out,'nflux_d9',v73_or_append(lastdecade_out))
%         clear nflux_d9
%         clear nflux
        
        % Fertilizer for each crop
        [x,NfertFile_tmp] = unix([thisDir 'get_nfert_file.sh ' inDir]) ;
        if x~=0
            error(['get_nfert_file.sh failed with error ' num2str(x)])
        end
        NfertFile_tmp = regexprep(NfertFile_tmp,'[\n\r]+','') ; % Remove extraneous newline
        if is_baseline
            NfertFile = ['/project/fh1-project-lpjgpi/lr8247/PLUM/input/Nfert/' NfertFile_tmp] ;
            if ~exist(NfertFile,'file')
                error('NfertFile not found!')
            end
        else
            NfertFile = find_PLUM2LPJG_input_file(NfertFile_tmp) ;
        end
        Nfert = lpjgu_matlab_readTable_then2map(NfertFile,'force_mat_save',true) ;
        Nfert = adjust_cropinput_yearLists(Nfert, yearList) ;
        Nfert_orig = Nfert ;
        if is_baseline
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
                    otherwise
                        error(['Which NfertType should I use for ' thisCrop '?'])
                end
                thisNfertIndex = find(strcmp(Nfert_orig.varNames,thisNfertType)) ;
                NfertIndices(c) = thisNfertIndex ;
            end
            Nfert.maps_YXvy = Nfert_orig.maps_YXvy(:,:,NfertIndices,:) ;
        else
            Nfert = CrOp_and_CrOpi(Nfert, 'Nfert', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp) ;
        end
        % Convert kg/m2 to tons/ha (kg/m2 * tons/kg * m2/ha)
        Nfert.maps_YXvy = Nfert.maps_YXvy * 1e-3 * 1e4 ;
        % Save time series (tons/ha --> tons)
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['nflux_ts_fert_' thisCrop ' = getTS(Nfert,thisCrop,1e2*cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
            save(timeseries_out,['nflux_ts_fert_' thisCrop],v73_or_append(timeseries_out)) ;
            clear thisCrop
        end; clear c
        if exist('LU0','var')
            for c = 1:length(cropTypes)
                thisCrop = cropTypes{c} ;
                eval(['nflux0_ts_fert_' thisCrop ' = getTS(Nfert,thisCrop,1e2*cropareas0.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
                save(timeseries_out,['nflux0_ts_fert_' thisCrop],v73_or_append(timeseries_out)) ;
                clear thisCrop
            end; clear c
        end
        clear *_ts_*
        % Save maps (tons/ha)
        Nfert_d1 = Nfert ;
        Nfert_d1.maps_YXvy = Nfert_d1.maps_YXvy(:,:,:,1:10) ;
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

    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Save FPC stuff %%%
    %%%%%%%%%%%%%%%%%%%%%%

    if do_save.fpc
        if exist('LU0','var')
            warning('Should you do an "LU0" adjusted version of FPC?')
        end
        disp('   Saving FPC...')
        fpc_d1 = fpc ;
        fpc_d1.maps_YXvy = fpc_d1.maps_YXvy(:,:,:,1:10) ;
        fpc_d1.yearList = fpc_d1.yearList(1:10) ;
        save(firstdecade_out,'fpc_d1',v73_or_append(firstdecade_out))
        clear fpc_d1
        fpc_d9 = fpc ;
        fpc_d9.maps_YXvy = fpc_d9.maps_YXvy(:,:,:,end-9:end) ;
        fpc_d9.yearList = fpc_d9.yearList(end-9:end) ;
        save(lastdecade_out,'fpc_d9',v73_or_append(lastdecade_out))
        clear fpc_d9
        clear fpc
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save land use and crop maps %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.crops
        disp('   Saving crop maps...')
        cropfracs_d1 = cropfracs ;
        cropfracs_d1.maps_YXvy = cropfracs_d1.maps_YXvy(:,:,:,1:10) ;
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
        LU_d1.maps_YXvy = LU_d1.maps_YXvy(:,:,:,1:10) ;
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
        LU0_d1.maps_YXvy = LU0_d1.maps_YXvy(:,:,:,1:10) ;
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
