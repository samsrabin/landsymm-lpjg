%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate emulator-projected yields as deltas from LPJ-GUESS baseline %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

whichfile_list = {'yield', 'gsirrigation'} ;

% Development vs. production
figure_visibility = 'off' ; % 'off' or 'on'. Determines whether figures are shown on screen
figure_extension = 'png' ; % fig or png
save_excl_figs = true ;
save_interp_figs = true ;
save_out_figs = true ;
which_out_figs = {'max'} ; % {'max', 'first', 'first0', '4th', '4th0'}
save_txt_files = true ;
load_existing_file = false ;

% Behaviors
excl_lowBL_agmerra = true ;
excl_lowBL_emu = true ;
interp_infs = true ;
when_remove_outliers = 'end' ; % end, before_interp, off
fake1k = true ;
overwrite_existing_txt = true ;
overwrite_existing_figs = true ;
use_lpjg_baseline = false ;
use_ph2_baseline = true ;

% Run info
gcm_list = {'UKESM1-0-LL'} ;
ggcm_list = {'GEPIC', 'EPIC-TAMU', 'pDSSAT'} ;
ssp_list = {'ssp126', 'ssp585'} ;
thisVer = '20210226' ;
emuVer = 'v2.5' ;
adaptation = 1 ;

baseline_y1 = 2001 ;
baseline_yN = 2010 ;
future_ts = 10 ; % Number of years in future time step
future_yN_emu = 2084 ;

tmp_rcp = 'rcp45' ;
warning('Arbitrarily using %s LPJ-GUESS run! Fix this once you''ve done the CMIP6 runs.', ...
    tmp_rcp)


%% Setup

current_dir = pwd ;
if strcmp(current_dir(1:6), '/Users') || strcmp(current_dir(1:6), '/Volum')
    topdir_db = '/Users/sam/Documents/Dropbox/2016_KIT/GGCMI/GGCMI2PLUM_DB' ;
    topDir_lpj = sprintf('%s/GGCMIPLUM_2001-2100_remap6p7_forPotYields_%s', ...
        '/Volumes/Reacher/GGCMI/g2p/lpj-guess_runs', ...
        tmp_rcp) ;
    topDir_agmipout = '/Volumes/Reacher/GGCMI/AgMIP.output' ;
    topDir_emu = sprintf('%s/CMIP_emulated', topDir_agmipout) ;
    topDir_lpj0 = sprintf( ...
        '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_%s_forED', ...
        tmp_rcp) ;
elseif strcmp(current_dir(1:3), '/pd')
    topdir_db = '/pd/data/lpj/sam/ggcmi2plum' ;
    topDir_lpj = sprintf('%s/GGCMIPLUM_2001-2100_remap6p7_forPotYields_%s', ...
        topdir_db, tmp_rcp) ;
    topDir_agmipout = sprintf('%s/AgMIP.output', topdir_db) ;
    topDir_emu = '/pd/data/lpj/sam/ggcmi2plum/CMIP_emulated' ;
    topDir_lpj0 = sprintf( ...
        '/pd/data/lpj/sam/ggcmi2plum/lpj-guess_runs/LPJGPLUM_2001-2100_remap6p7_forPotYields_%s_forED', ...
        tmp_rcp) ;
else
    error('Failed to interpret what system you''re on')
end

% If started with -nodisplay option, you have to specify '-painters' in the
% export_fig() call.
if ~usejava('desktop')
    renderer = '-painters' ;
else
    renderer = '-opengl' ;
end

cd(sprintf('%s/emulation/matlab/emu2plum', topdir_db))

save_out_figs_Nth0 = save_out_figs & any( ...
    strcmp(which_out_figs, 'first0') ...
    | strcmp(which_out_figs, '4th0')) ;

% Check outlier removal setting
if ~any(strcmp({'off', '', 'end', 'before_interp'}, when_remove_outliers))
    error('when_remove_outliers not recognized: %s', when_remove_outliers)
end
remove_outliers = ~strcmp(when_remove_outliers, 'off') ...
    & ~isempty(when_remove_outliers) ;

% Check that directories exist
if ~contains(topDir_lpj, tmp_rcp)
    error('~contains(topDir_lpj, rcp)')
elseif ~exist(topDir_lpj, 'dir')
    error('topDir_lpj does not exist:\n %s', topDir_lpj)
elseif ~exist(topDir_lpj0, 'dir')
    error('topDir_lpj0 does not exist:\n %s', topDir_lpj0)
elseif ~exist(topDir_emu, 'dir')
    error('topDir_emu does not exist:\n %s', topDir_emu)
elseif ~exist(topDir_agmipout, 'dir')
    error('topDir_agmipout does not exist:\n %s', topDir_agmipout)
end

% getbasename = @(x) regexprep(x,'i?\d\d\d$','') ;
% getbasenamei = @(x) regexprep(x,'\d\d\d$','') ;
% getbasename0 = @(x) regexprep(regexprep(regexprep(x,'i?\d\d\d\d$',''),'i0$',''),'0$','') ;
% getbasename0i = @(x) regexprep(regexprep(x,'\d\d\d\d$',''),'0$','') ;
getN = @(x) regexprep(regexprep(x, 'CerealsC[34]', ''), '^[a-zA-Z_]+', '') ;
get_unneeded = @(x)cellfun(@isempty, ...
    regexp(regexprep(x,'CerealsC[34]','CerealsC'),'.*\d+')) | contains(x,'G_ic') ;

Nlist = [10 60 200] ;
NN = length(Nlist) ;

irrList_in = {'rf', 'ir'} ;
irrList_out = {'', 'i'} ;

gridlist_file = sprintf('%s/../gridlist_62892.runAEclimOK.txt', ...
    topDir_lpj) ;
gridlist = lpjgu_matlab_readTable_then2map(gridlist_file, ...
    'force_mat_nosave', true, 'force_mat_save', false) ;
gridlist_target = {gridlist.lonlats gridlist.list_to_map} ;

% Figure out timesteps
Nyears_ts = future_ts ;
ts1_y1 = ceil(baseline_yN/Nyears_ts)*Nyears_ts+1 ;
tsN_y1 = floor(future_yN_emu/Nyears_ts)*Nyears_ts ;
ts1_list = ts1_y1:Nyears_ts:tsN_y1 ;
tsN_list = ts1_list + Nyears_ts - 1 ;
Ntpers = length(ts1_list) ;

% What is the last year of the last time period that should be read of the
% LPJ-GUESS run?
future_yN_lpj = max(tsN_list) ;

if save_out_figs && isempty(which_out_figs)
    warning('save_out_figs is true but which_out_figs is empty. Will not make figures.')
    save_out_figs = false ;
end

% Ensure you have these options set correctly
if use_lpjg_baseline == use_ph2_baseline
    error('Choose either use_lpjg_baseline or use_ph2_baseline')
end

% Exclude cells with <0.1 t/ha yield in Phase 2 baseline simulations,
% as Jim did
low_yield_threshold_kgm2 = 0.01 ;


%% Set up crop lists

cropList_in = {'spring_wheat', 'winter_wheat', 'maize', 'soy', 'rice'} ;
cropList_out = {'CerealsC3', 'CerealsC4', 'Rice', 'Oilcrops', 'Pulses', 'StarchyRoots'} ;
Ncrops_in = length(cropList_in) ;
Ncrops_out = length(cropList_out) ;

% Get output file header
cropIrrList_out = [cropList_out strcat(cropList_out, 'i')] ;
header_out = 'Lon\tLat' ;
format_out = '%0.2f\t%0.2f' ;
cropIrrNlist_out = {} ;
for c = 1:length(cropIrrList_out)
    thisCropIrr = cropIrrList_out{c} ;
    for n = 1:NN
        thisN = pad(num2str(Nlist(n)), 3, 'left', '0') ;
        thisCropIrrN = [thisCropIrr thisN] ;
        cropIrrNlist_out{end+1} = thisCropIrrN ; %#ok<SAGROW>
        header_out = [header_out '\t' thisCropIrrN] ; %#ok<AGROW>
        format_out = [format_out '\t%0.3f'] ; %#ok<AGROW>
    end
end
header_out = [header_out '\n'] ;
format_out = [format_out '\n'] ;


%% Import LPJ-GUESS yield and irrigation

disp('Importing LPJ-GUESS yield...')

which_file = 'yield' ;

data_bl_lpj_yield = e2p_import_bl_lpj(baseline_y1, baseline_yN, topDir_lpj, ...
    which_file, get_unneeded, gridlist_target) ;
e2p_check_correct_zeros(data_bl_lpj_yield.garr_xv, ...
    which_file, data_bl_lpj_yield.varNames, ...
    'Baseline', @getbasenamei)

data_fu_lpj_yield = e2p_import_fu_lpj(baseline_yN, future_ts, future_yN_lpj, topDir_lpj, ...
    which_file, data_bl_lpj_yield.varNames, get_unneeded, gridlist_target) ;
e2p_check_correct_zeros(data_fu_lpj_yield.garr_xvt, ...
    which_file, data_fu_lpj_yield.varNames, ...
    'Future', @getbasenamei)

[varNames_lpj, cropList_lpj, ...
    varNames_lpj_basei, cropList_lpj_basei, ...
    Nlist_lpj, ~] = ...
    e2p_get_names(data_bl_lpj_yield.varNames, data_fu_lpj_yield.varNames, ...
    getN, get_unneeded) ;
varNames_out = varNames_lpj ;
cropList_out = cropList_lpj ;
varNames_out_basei = varNames_lpj_basei ;
cropList_out_basei = cropList_lpj_basei ;

if ~isequal(sort(cropIrrNlist_out), sort(varNames_lpj))
    error('Mismatch between cropIrrNlist_out and varNames_lpj')
end


disp('Importing LPJ-GUESS irrigation...')

which_file = 'gsirrigation' ;

data_bl_lpj_irrig = e2p_import_bl_lpj(baseline_y1, baseline_yN, topDir_lpj, ...
    which_file, get_unneeded, gridlist_target) ;
e2p_check_correct_zeros(data_bl_lpj_irrig.garr_xv, ...
    which_file, data_bl_lpj_irrig.varNames, ...
    'Baseline', @getbasenamei)

data_fu_lpj_irrig = e2p_import_fu_lpj(baseline_yN, future_ts, future_yN_lpj, topDir_lpj, ...
    which_file, data_bl_lpj_irrig.varNames, get_unneeded, gridlist_target) ;
e2p_check_correct_zeros(data_fu_lpj_irrig.garr_xvt, ...
    which_file, data_fu_lpj_irrig.varNames, ...
    'Future', @getbasenamei)

[varNames_lpj, cropList_lpj2, ...
    varNames_lpj_basei2, cropList_lpj_basei2, ...
    Nlist_lpj2, ~] = ...
    e2p_get_names(data_bl_lpj_irrig.varNames, data_fu_lpj_irrig.varNames, ...
    getN, get_unneeded) ;

if ~isequal(sort(cropIrrNlist_out), sort(varNames_lpj))
    error('Mismatch between cropIrrNlist_out and varNames_lpj')
elseif ~isequal(cropList_lpj,cropList_lpj2)
    error('Mismatch between cropList_lpj for yield vs. gsirrigation')
elseif ~isequal(varNames_lpj_basei,varNames_lpj_basei2)
    error('Mismatch between varNames_lpj_basei for yield vs. gsirrigation')
elseif ~isequal(cropList_lpj_basei,cropList_lpj_basei2)
    error('Mismatch between cropList_lpj_basei for yield vs. gsirrigation')
elseif ~isequal(Nlist_lpj,Nlist_lpj2)
    error('Mismatch between Nlist_lpj2 for yield vs. gsirrigation')
end
clear cropList_lpj2 varNames_lpj_basei2 cropList_lpj_basei2 Nlist_lpj2

disp('Done.')


%% Import LPJ-GUESS-0 yield and irrigation

if save_out_figs_Nth0 || fake1k
    disp('Importing LPJ-GUESS-0 yield...')
    
    which_file = 'yield' ;
    
    data_bl_lpj0_yield = e2p_import_bl_lpj(baseline_y1, baseline_yN, topDir_lpj0, ...
        which_file, get_unneeded, gridlist_target) ;
    e2p_check_correct_zeros(data_bl_lpj0_yield.garr_xv, ...
        which_file, data_bl_lpj0_yield.varNames, ...
        'Baseline', @getbasenamei)
    
    data_fu_lpj0_yield = e2p_import_fu_lpj(baseline_yN, future_ts, future_yN_lpj, topDir_lpj0, ...
        which_file, data_bl_lpj0_yield.varNames, get_unneeded, gridlist_target) ;
    e2p_check_correct_zeros(data_fu_lpj0_yield.garr_xvt, ...
        which_file, data_fu_lpj0_yield.varNames, ...
        'Future', @getbasenamei)
    
    [varNames_lpj0, cropList_lpj0, ...
        varNames_lpj0_basei, cropList_lpj0_basei, ...
        Nlist_lpj0, ~] = ...
        e2p_get_names(data_bl_lpj0_yield.varNames, data_fu_lpj0_yield.varNames, ...
        getN, get_unneeded) ;
    
else
    data_bl_lpj0_yield = [] ;
    data_fu_lpj0_yield = [] ;
end


% Import LPJ-GUESS-0 irrigation

if save_out_figs_Nth0 || fake1k
    disp('Importing LPJ-GUESS-0 irrigation...')
    
    which_file = 'gsirrigation' ;
    
    data_bl_lpj0_irrig = e2p_import_bl_lpj(baseline_y1, baseline_yN, topDir_lpj0, ...
        which_file, get_unneeded, gridlist_target) ;
    e2p_check_correct_zeros(data_bl_lpj0_irrig.garr_xv, ...
        which_file, data_bl_lpj0_irrig.varNames, ...
        'Baseline', @getbasenamei)
    
    data_fu_lpj0_irrig = e2p_import_fu_lpj(baseline_yN, future_ts, future_yN_lpj, topDir_lpj0, ...
        which_file, data_bl_lpj0_irrig.varNames, get_unneeded, gridlist_target) ;
    e2p_check_correct_zeros(data_fu_lpj0_irrig.garr_xvt, ...
        which_file, data_fu_lpj0_irrig.varNames, ...
        'Future', @getbasenamei)
    
    [varNames_lpj0, cropList_lpj02, ...
        varNames_lpj0_basei2, cropList_lpj0_basei2, ...
        Nlist_lpj02, ~] = ...
        e2p_get_names(data_bl_lpj0_irrig.varNames, data_fu_lpj0_irrig.varNames, ...
        getN, get_unneeded) ;
    
    if ~isequal(sort(cropIrrNlist_out), sort(varNames_lpj0))
        warning('Mismatch between cropIrrNlist_out and varNames_lpj0')
    elseif ~isequal(cropList_lpj0,cropList_lpj02)
        error('Mismatch between cropList_lpj0 for yield vs. gsirrigation')
    elseif ~isequal(varNames_lpj0_basei,varNames_lpj0_basei2)
        error('Mismatch between varNames_lpj0_basei for yield vs. gsirrigation')
    elseif ~isequal(cropList_lpj0_basei,cropList_lpj0_basei2)
        error('Mismatch between cropList_lpj0_basei for yield vs. gsirrigation')
    elseif ~isequal(Nlist_lpj0,Nlist_lpj02)
        error('Mismatch between Nlist_lpj02 for yield vs. gsirrigation')
    end
    clear cropList_lpj02 varNames_lpj0_basei2 cropList_lpj0_basei2 Nlist_lpj02
    
    disp('Done.')
else
    data_bl_lpj0_irrig = [] ;
    data_fu_lpj0_irrig = [] ;
end


%% Loop though GGCM emulators

e2p__child





