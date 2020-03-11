%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate emulator-projected yields as deltas from LPJ-GUESS baseline %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

which_file = 'yield' ;
% which_file = 'gsirrigation' ;

excl_lowBL_agmerra = true ;
excl_lowBL_emu = true ;
interp_infs = true ;
remove_outliers = false ;

gcm = 'IPSL-CM5A-MR_r1i1p1' ;
ggcm = 'LPJ-GUESS' ;
rcp = 'rcp45' ;
thisVer = '20200310' ;

baseline_y1 = 2001 ;
baseline_yN = 2010 ;
future_ts = 10 ; % Number of years in future time step
future_yN_lpj = 2100 ;
future_yN_emu = 2099 ;

topDir_phase2 = sprintf('/Users/Shared/GGCMI/AgMIP.output/%s/phase2', ggcm) ;

topDir_lpj = '/Users/Shared/GGCMI2PLUM_sh/lpj-guess_runs/GGCMIPLUM_2001-2100_remap6p7_forPotYields_rcp45_forED_20191104133630' ;
topDir_emu_bl = sprintf('/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcropsBaseline_%s', thisVer) ;
topDir_emu_fu = sprintf('/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcrops_%s', thisVer) ;
outDir = sprintf('/Users/Shared/GGCMI2PLUM_sh/send_to_plum/emulator_%s', thisVer) ;


%% Setup 

cd '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/matlab/emu2plum'

getbasename = @(x) regexprep(x,'i?\d\d\d$','') ;
getbasenamei = @(x) regexprep(x,'\d\d\d$','') ;
getN = @(x) x(end-2:end) ;
get_unneeded = @(x) cellfun(@isempty,regexp(x, '.*\d\d+')) ;

outDir_lpj = sprintf('%s/LPJ-GUESS', outDir) ;
outDir_ggcm = sprintf('%s/emul_%s', outDir, ggcm) ;

if ~exist(outDir, 'dir')
    mkdir(outDir) ;
end
if ~exist(outDir_lpj, 'dir')
    mkdir(outDir_lpj) ;
end
if ~exist(outDir_ggcm, 'dir')
    mkdir(outDir_ggcm) ;
end


%% Import LPJ-GUESS

disp('Importing LPJ-GUESS...')

data_bl_lpj = e2p_import_bl_lpj(baseline_y1, baseline_yN, topDir_lpj, ...
    which_file, get_unneeded) ;
e2p_check_correct_zeros(data_bl_lpj.garr_xv, which_file, getbasenamei(data_bl_lpj.varNames))

data_fu_lpj = e2p_import_fu_lpj(baseline_yN, future_ts, future_yN_lpj, topDir_lpj, ...
    which_file, data_bl_lpj.varNames, get_unneeded) ;
e2p_check_correct_zeros(data_fu_lpj.garr_xvt, which_file, getbasenamei(data_fu_lpj.varNames))

[varNames_lpj, cropList_lpj, ...
    varNames_lpj_basei, cropList_lpj_basei, ...
    Nlist_lpj, ~] = ...
    e2p_get_names(data_bl_lpj.varNames, data_fu_lpj.varNames, ...
    getbasename, getbasenamei, getN) ;

disp('Done.')


%% Import emulator outputs

if contains(which_file, {'yield', 'gsirrigation'})

    disp('Importing emulator outputs...')
    
    [data_bl_emu, data_fu_emu, Ntpers] = e2p_import_emu( ...
        topDir_emu_bl, topDir_emu_fu, gcm, ggcm, rcp, which_file, ...
        baseline_yN, future_yN_emu) ;
    
    e2p_check_correct_zeros(data_bl_emu.garr_xv, which_file, getbasenamei(data_bl_emu.varNames))
    e2p_check_correct_zeros(data_fu_emu.garr_xvt, which_file, getbasenamei(data_fu_emu.varNames))
    
    [varNames_emu, cropList_emu, ...
        varNames_emu_basei, cropList_emu_basei, ...
        Nlist_emu, ~] = ...
        e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
        getbasename, getbasenamei, getN) ;
    
    disp('Done.')

elseif ~strcmp(which_file, 'anpp')
    error('which_file (%s) not recognized', which_file)
end


%% Get and apply exclusions, if doing so

if contains(which_file, {'yield','gsirrigation'}) && (excl_lowBL_agmerra || excl_lowBL_emu)
    if strcmp(which_file, 'yield')
        
        % Set up exclusion array
        exclude_xc = false(length(data_bl_emu.list2map),length(cropList_emu_basei)) ;
        
        % Where do we exclude based on low AgMERRA yield (or existing exclusions)?
        if excl_lowBL_agmerra
            disp('Excluding based on low AgMERRA yield...')
            % Get AgMERRA yield
            yield_agmerraBL_xv = e2p_get_agmerra_yield(...
                varNames_emu, topDir_phase2, ggcm, data_bl_emu.list2map, getbasenamei, getN) ;
            if ~any(any(~isnan(yield_agmerraBL_xv)))
                error('yield_agmerraBL_xv is all NaN')
            end
            % Exclude
            exclude_xc = e2p_exclude_lowBLyield_atMaxN( ...
                varNames_emu, cropList_emu_basei, Nlist_emu, ...
                yield_agmerraBL_xv, 0.01, exclude_xc) ;
        end
        
        if ~any(any(~exclude_xc))
            error('All cells excluded because of low AgMERRA yield')
        end
        
        % Where do we exclude based on low baseline-year emulated yield (or
        % existing exclusions)?
        if excl_lowBL_emu
            disp('Excluding based on low baseline-year emulated yield...')
            exclude_xc = e2p_exclude_lowBLyield_atMaxN( ...
                varNames_emu, cropList_emu_basei, Nlist_emu, ...
                data_bl_emu.garr_xv, 0.01, exclude_xc) ;
        end
        
        if ~any(any(~exclude_xc))
            error('All cells excluded because of low AgMERRA yield and/or low baseline-year emulated yield')
        end
        
    % Error checks
    elseif strcmp(which_file, 'gsirrigation')
        disp('Applying yield-based exclusions...')
        exclude_xc_file = sprintf('%s/exclude_xc.mat', outDir_ggcm) ;
        load(exclude_xc_file) ;
    else
        error('which_file (%s) not recognized', which_file)
    end
    
    for c = 1:length(cropList_emu_basei)
        thisCrop = cropList_emu_basei{c} ;
        thisCrop_i = find(strcmp(varNames_emu_basei,thisCrop)) ;
        if length(thisCrop_i)~=length(Nlist_emu)
            error('Error finding isThisCrop (%d found)', length(thisCrop_i))
        end
        
        % Apply to emulated baseline
        tmp = data_bl_emu.garr_xv(:,thisCrop_i) ;
        tmp(exclude_xc(:,c),:) = NaN ;
        data_bl_emu.garr_xv(:,thisCrop_i) = tmp ;
        clear tmp
        
        % Apply to emulated future
        tmp = data_fu_emu.garr_xvt(:,thisCrop_i,:) ;
        tmp(exclude_xc(:,c),:,:) = NaN ;
        data_fu_emu.garr_xvt(:,thisCrop_i,:) = tmp ;
        clear tmp
    end
    
    disp('Done.')

elseif ~strcmp(which_file,'anpp')
    error('which_file (%s) not recognized', which_file)
end

% Sanity checks
if ~any(any(~isnan(data_bl_emu.garr_xv)))
    error('data_bl_emu.garr_xv is all NaN!')
end
if ~any(any(any(~isnan(data_fu_emu.garr_xvt))))
    error('data_fu_emu.garr_xvt is all NaN!')
end


%% Get max wheats

disp('Getting max wheats...')

refresh_vars = true ;
if strcmp(which_file, 'yield')
    [data_bl_emu.garr_xv, data_bl_emu.varNames, is_ww_max_bl_gW, winter_wheats] = ...
        e2p_get_max_wheat(data_bl_emu.garr_xv, data_bl_emu.varNames) ;
    [data_fu_emu.garr_xvt, data_fu_emu.varNames, is_ww_max_fu_gWt, winter_wheats_test] = ...
        e2p_get_max_wheat(data_fu_emu.garr_xvt, data_fu_emu.varNames) ;
    if ~isequal(winter_wheats, winter_wheats_test)
        error('Mismatch between winter wheat lists')
    end
elseif strcmp(which_file, 'gsirrigation')    
    data_bl_emu = e2p_apply_max_wheat(data_bl_emu, outDir) ;
    data_fu_emu = e2p_apply_max_wheat(data_fu_emu, outDir) ;
elseif ~strcmp(which_file, 'anpp')
    error('which_file (%s) not recognized', which_file)
else
    refresh_vars = false ;
end

e2p_check_correct_zeros(data_bl_emu.garr_xv, which_file, getbasenamei(data_bl_emu.varNames))
e2p_check_correct_zeros(data_fu_emu.garr_xvt, which_file, getbasenamei(data_fu_emu.varNames))

% Refresh variable and crop lists
if refresh_vars
    [varNames_emu, cropList_emu, ...
        varNames_emu_basei, cropList_emu_basei, ...
        Nlist_emu, ~] = ...
        e2p_get_names(data_bl_emu.varNames, data_fu_emu.varNames, ...
        getbasename, getbasenamei, getN) ; %#ok<ASGLU>
end

disp('Done.')


%% Harmonize LPJ-GUESS and emulator data

disp('Harmonizing LPJ-GUESS and emulator data...')

% Align gridlists
[data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu, list2map] = ...
    e2p_align_gridlists(data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu) ;

e2p_check_correct_zeros(data_bl_lpj.garr_xv, which_file, getbasenamei(data_bl_lpj.varNames))
e2p_check_correct_zeros(data_fu_lpj.garr_xvt, which_file, getbasenamei(data_fu_lpj.varNames))
e2p_check_correct_zeros(data_bl_emu.garr_xv, which_file, getbasenamei(data_bl_emu.varNames))
e2p_check_correct_zeros(data_fu_emu.garr_xvt, which_file, getbasenamei(data_fu_emu.varNames))

% Make sure that N lists match
if ~isequal(Nlist_lpj,Nlist_emu)
    error('Mismatch in N levels between LPJ-GUESS and emulator')
end

% Translate crop names
[cropList_lpj_asEmu, used_emuCrops] = e2p_translate_crops( ...
    cropList_lpj, cropList_emu) ;

disp('Done.')


%% Get and apply deltas (and remove outliers, if doing so)

if contains(which_file, {'yield','gsirrigation'})
    disp('Getting deltas...')
    
    deltas_emu_xvt = e2p_get_deltas(...
        data_bl_emu, data_fu_emu, interp_infs, cropList_emu, ...
        getbasename, getbasenamei, which_file, ...
        used_emuCrops, list2map) ;
    e2p_check_correct_zeros(deltas_emu_xvt, which_file, getbasenamei(data_fu_emu.varNames))
    
    disp('Done.')
    
elseif ~strcmp(which_file, 'anpp')
    error('which_file (%s) not recognized', which_file)
end

disp('Applying deltas...')

[data_fu_lpj, data_fu_out] = e2p_apply_deltas( ...
    data_bl_lpj, data_fu_lpj, data_bl_emu, data_fu_emu, deltas_emu_xvt, ...
    cropList_lpj, cropList_lpj_asEmu, varNames_lpj, ...
    list2map, getbasename, getbasenamei, which_file) ;
e2p_check_correct_zeros(data_fu_out.garr_xvt, which_file, getbasenamei(data_fu_out.varNames))

if remove_outliers
    disp('Removing outliers...')
    
    if strcmp(which_file,'yield')
        smad_mult = 3 ;
    elseif strcmp(which_file,'gsirrigation')
        warning('Might have to change gsirrigation smad_mult once properly pre-thresholding')
        smad_mult = 3 ;
    else
        error('which_file (%s) not recognized', which_file)
    end

    [data_fu_lpj, outlier_info_lpj] = e2p_remove_outliers(data_fu_lpj, smad_mult) ;
    e2p_check_correct_zeros(data_fu_lpj.garr_xvt, which_file, getbasenamei(data_fu_lpj.varNames))
    
    [data_fu_out, outlier_info_out] = e2p_remove_outliers(data_fu_out, smad_mult) ;
    e2p_check_correct_zeros(data_fu_out.garr_xvt, which_file, getbasenamei(data_fu_out.varNames))
    
    % Save info
    outlier_info_cols = string([ ...
        repmat('y',[Ntpers 1]) ...
        num2str(data_fu_out.y1s) ...
        repmat('_',[Ntpers 1]) ...
        num2str(data_fu_out.yNs)]) ;
    e2p_save_outlier_info(outlier_info_lpj, outDir_lpj, which_file, outlier_info_cols)
    e2p_save_outlier_info(outlier_info_out, outDir_ggcm, which_file, outlier_info_cols)
    clear outlier_info_cols
end

disp('Done.')


%% Save outputs

disp('Saving...')

% Save exclusion info
out_file = sprintf('%s/exclude_xc.mat', outDir_ggcm) ;
save(out_file, 'exclude_xc')

% Save outputs
out_header_cell = [{'Lon', 'Lat'} data_fu_out.varNames] ;
for t = 1:Ntpers
    fprintf('    %d/%d...\n', t, Ntpers)
    y1 = data_fu_out.y1s(t) ;
    yN = data_fu_out.yNs(t) ;
        
    e2p_save(outDir_lpj, y1, yN, out_header_cell, ...
        data_fu_lpj.lonlats, data_fu_lpj.garr_xvt(:,:,t), which_file, ...
        interp_infs, remove_outliers, false)
    e2p_save(outDir_ggcm, y1, yN, out_header_cell, ...
        data_fu_out.lonlats, data_fu_out.garr_xvt(:,:,t), which_file, ...
        interp_infs, remove_outliers, true)
    
end

% For yield, save is_ww_max_bl_gW*
if strcmp(which_file, 'yield')
    out_file = sprintf('%s/is_ww_max.mat', outDir_ggcm) ;
    save(out_file, 'is_ww_max_bl_gW', 'is_ww_max_fu_gWt', 'winter_wheats') ;
elseif ~contains(which_file, {'anpp', 'gsirrigation'})
    error('which_file (%s) not recognized', which_file)
end

disp('All done!')


%% UNUSED: Applying deltas to AgMERRA baseline

% % Set up output structure
% data_fu_out.list2map = list2map ;
% data_fu_out.lonlats = data_bl_emu.lonlats ;
% data_fu_out.varNames = varNames_emu ;
% data_fu_out.y1s = data_fu_emu.y1s ;
% data_fu_out.yNs = data_fu_emu.yNs ;
% Actually apply the deltas
% data_fu_out.garr_xvt = deltas_emu_xvt .* ...
%     repmat(yield_agmerraBL_xv, [1 1 Ntpers]) ;


