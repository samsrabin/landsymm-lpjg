%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate emulator-projected yields as deltas from LPJ-GUESS baseline %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

which_file = 'yield' ;
% which_file = 'gsirrigation' ;

% Development vs. production
save_interp_figs = false ;
save_out_figs = true ;

% Behaviors
excl_lowBL_agmerra = true ;
excl_lowBL_emu = true ;
interp_infs = true ;
remove_outliers = true ;

% Run info
gcm = 'IPSL-CM5A-MR_r1i1p1' ;
ggcm_list = {'LPJ-GUESS', 'LPJmL', 'pDSSAT', 'EPIC-TAMU'} ;
% ggcm_list = {'LPJ-GUESS'} ;
% ggcm_list = {'EPIC-TAMU', 'LPJmL', 'pDSSAT'} ;
rcp = 'rcp45' ;
thisVer = '20200310' ;

baseline_y1 = 2001 ;
baseline_yN = 2010 ;
future_ts = 10 ; % Number of years in future time step
future_yN_lpj = 2100 ;
future_yN_emu = 2099 ;


%% Setup 

cd '/Users/sam/Documents/Dropbox/GGCMI2PLUM_DB/emulation/matlab/emu2plum'

topDir_lpj = sprintf( ...
    '/Users/Shared/GGCMI2PLUM_sh/lpj-guess_runs/GGCMIPLUM_2001-2100_remap6p7_forPotYields_%s', ...
    rcp) ;
topDir_emu_bl = sprintf('/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcropsBaseline_%s', thisVer) ;
topDir_emu_fu = sprintf('/Users/Shared/GGCMI2PLUM_sh/emulation/outputs/outputs_GGCMIcrops_%s', thisVer) ;
outDir = sprintf('/Users/Shared/GGCMI2PLUM_sh/send_to_plum/%s_%s_v%s', gcm, rcp, thisVer) ;
if remove_outliers
    outDir = [outDir '_rmol'] ;
end

if ~contains(topDir_lpj, rcp)
    error('~contains(topDir_lpj, rcp)')
elseif ~exist(topDir_lpj, 'dir')
    error('topDir_lpj does not exist:\n %s', topDir_lpj)
end

getbasename = @(x) regexprep(x,'i?\d\d\d$','') ;
getbasenamei = @(x) regexprep(x,'\d\d\d$','') ;
getN = @(x) x(end-2:end) ;
get_unneeded = @(x) cellfun(@isempty,regexp(x, '.*\d\d+')) ;

outDir_lpj = sprintf('%s/sim_LPJ-GUESS', outDir) ;

outDir_interp_figs = sprintf('%s/interp_figs', outDir) ;
outDir_yield_figs = sprintf('%s/yield_figs', outDir) ;
outDir_irrig_figs = sprintf('%s/irrig_figs', outDir) ;

if ~exist(outDir, 'dir')
    mkdir(outDir) ;
end
if ~exist(outDir_lpj, 'dir')
    mkdir(outDir_lpj) ;
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


%% Loop though GGCM emulators

e2p__child





