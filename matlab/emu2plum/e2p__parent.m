%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate emulator-projected yields %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

whichfile_list = {'yield', 'gsirrigation'} ;

% Development vs. production
figure_visibility = 'off' ; % 'off' or 'on'. Determines whether figures are shown on screen
figure_extension = 'png' ; % fig or png
save_excl_figs = false ;
save_interp_figs = false ;
save_out_figs = false ;
which_out_figs = {'max'} ; % {'max', 'first', 'first0', '4th', '4th0'}
save_txt_files = false ;
load_existing_file = false ;

% Behaviors
excl_lowBL_agmerra = true ;
excl_lowBL_emu = true ;
interp_infs = true ;
when_remove_outliers = 'end' ; % end, before_interp, off
fake1k = true ;
overwrite_existing_txt = true ;
overwrite_existing_figs = true ;

% Run info
gcm_list = {'UKESM1-0-LL'} ;
ggcm_list = {'GEPIC', 'EPIC-TAMU', 'pDSSAT'} ;
ssp_list = {'ssp126', 'ssp585'} ;
thisVer = '20210226' ;
emuVer = 'v2.5' ;
adaptation = 1 ;

future_y1 = 2015 ;
baseline_y1 = 2001 ;
baseline_yN = 2010 ;
future_ts = 10 ; % Number of years in future time step
future_yN_emu = 2084 ;


%% Setup

current_dir = pwd ;
if strcmp(current_dir(1:6), '/Users') || strcmp(current_dir(1:6), '/Volum')
    which_system = 'mymac' ;
    topdir_db = '/Users/sam/Documents/Dropbox/2016_KIT/GGCMI/GGCMI2PLUM_DB' ;
    topDir_agmipout = '/Volumes/Reacher/GGCMI/AgMIP.output' ;
    topDir_emu = sprintf('%s/CMIP_emulated', topDir_agmipout) ;
elseif strcmp(current_dir(1:3), '/pd')
    which_system = 'keal' ;
    topdir_db = '/pd/data/lpj/sam/ggcmi2plum' ;
    topDir_agmipout = sprintf('%s/AgMIP.output', topdir_db) ;
    topDir_emu = '/pd/data/lpj/sam/ggcmi2plum/CMIP_emulated' ;
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

% Check outlier removal setting
if ~any(strcmp({'off', '', 'end', 'before_interp'}, when_remove_outliers))
    error('when_remove_outliers not recognized: %s', when_remove_outliers)
end
remove_outliers = ~strcmp(when_remove_outliers, 'off') ...
    & ~isempty(when_remove_outliers) ;

% Check that directories exist
if ~exist(topDir_emu, 'dir')
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

% Set up N lists
Nlist_lpj = [10 60 200 1000] ;
Nlist_emu = [10 60 200] ;
Nlist_out = [10 60 1000] ;

irrList_in = {'rf', 'ir'} ;
irrList_out = {'', 'i'} ;

% Figure out timesteps
Nyears_ts = future_ts ;
tsN_y1 = floor(future_yN_emu/Nyears_ts)*Nyears_ts ;
ts1_list = future_y1:Nyears_ts:tsN_y1 ;
tsN_list = ts1_list + Nyears_ts - 1 ;
Ntpers = length(ts1_list) ;

% What is the last year of the last time period that should be read of the
% LPJ-GUESS run?
future_yN_lpj = max(tsN_list) ;

if save_out_figs && isempty(which_out_figs)
    warning('save_out_figs is true but which_out_figs is empty. Will not make figures.')
    save_out_figs = false ;
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
% cropIrrNlist_out = {} ;
for c = 1:length(cropIrrList_out)
    thisCropIrr = cropIrrList_out{c} ;
    for n = 1:length(Nlist_out)
        thisN = pad(num2str(Nlist_out(n)), 3, 'left', '0') ;
        thisCropIrrN = [thisCropIrr thisN] ;
%         cropIrrNlist_out{end+1} = thisCropIrrN ; %#ok<SAGROW>
        header_out = [header_out '\t' thisCropIrrN] ; %#ok<AGROW>
        format_out = [format_out '\t%0.3f'] ; %#ok<AGROW>
    end
end
header_out = [header_out '\n'] ;
format_out = [format_out '\n'] ;


%% Loop though GGCM emulators

e2p__child





