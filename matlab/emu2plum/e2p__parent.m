%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate emulator-projected yields %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

whichfile_list = {'yield', 'gsirrigation'} ;

% Development vs. production
figure_visibility = 'on' ; % 'off' or 'on'. Determines whether figures are shown on screen
figure_extension = 'png' ; % fig or png
save_excl_figs = false ;
save_interp_figs = false ;
save_out_figs = false ;
which_out_figs = {'max'} ; % {'max', 'first', 'first0', '4th', '4th0'}
save_txt_files_emu = false ;
save_txt_files_lpjg = false ;
load_existing_file = false ;
overwrite_existing_txt = false ;
overwrite_existing_figs = false ;

% Behavior for combining crops
%  {:,1} = Output crop
%  {:,2} = "Mid" crops (i.e., the crops in the calibration factor set)
% % % combineCrops = { ...
% % %     'CerealsC3', {'CerealsC3s', 'CerealsC3w'} ;
% % %     'Oilcrops', {'OilNfix', 'OilOther'} ;
% % %     'Sugar', {'Sugarbeet', 'Sugarcane'} ;
% % %     } ;
combineCrops = { ...
    'CerealsC3', {'CerealsC3s', 'CerealsC3w'} ;
    } ;

% Other behaviors
excl_lowBL_agmerra = true ;
excl_lowBL_emu = true ;
interp_infs = true ;
when_remove_outliers = 'end' ; % end, before_interp, off
fake1k = true ;
scale_200to1000 = true ;
force_consistent_baseline = true ;
emulated_baseline = true ;

% Run info
gcm_list = {'UKESM1-0-LL'} ;
% ggcm_list = {'GEPIC', 'EPIC-TAMU', 'pDSSAT'} ;
ggcm_list = {'EPIC-TAMU'} ;
ssp_list = {'ssp126'} ;
thisVer = '20210702' ;
emuVer = 'v2.5' ;
adaptation = 1 ;

future_y1 = 2005 ;
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

% Define location of calibration factor files:
cfDir = sprintf('%s/emulation/calibration_factors/calibration_factors_20210526', ...
    topdir_db) ;

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
% getN_char = @(x) regexprep(regexprep(x, 'CerealsC[34]', ''), '^[a-zA-Z_]+', '') ;
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

% Make sure burnIn sources are sorted alphabetically
for c = 1:size(combineCrops,1)
    combineCrops{c,2} = sort(combineCrops{c,2}) ;
end; clear c


%% Set up crop lists

% The crops present in the GGCMI outputs
cropList_in = {'spring_wheat', 'winter_wheat', 'maize', 'soy', 'rice'} ;

% The crops present in the calibration factors
cropList_mid = {'CerealsC3', 'CerealsC4', 'Rice', 'Pulses', 'StarchyRoots', ...
    'OilNfix', 'OilOther', 'Sugarbeet', 'Sugarcane', 'FruitAndVeg'} ;

% The crops present in the outputs of this code
% % % cropList_out = {'CerealsC3', 'CerealsC4', 'Rice', 'Pulses', 'StarchyRoots', ...
% % %     'Oilcrops', 'Sugar', 'FruitAndVeg'} ;
cropList_out = cropList_mid ;

% Sanity checks for combineCrops
if isempty(combineCrops)
    if ~isequal(shiftdim(sort(cropList_mid)), shiftdim(sort(cropList_out)))
        error('If cropList_mid and cropList_out aren''t the same, you must specify a transformation key using combineCrops')
    end
else
    
    % Make sure "source" types in combineCrops are in cropList_mid
    % OR the destination type is
    for c = 1:size(combineCrops, 1)
        thisDest = combineCrops{c,1} ;
        theseSources = combineCrops{c,2} ;
        D = setdiff(theseSources, cropList_mid) ;
        if ~isempty(D) && ~any(strcmp(cropList_mid, thisDest))
            error('Neither %s nor all of its sources were found in cropList_mid')
        end
    end; clear c D thisDest theseSources IA
    
    % Make sure "destination" types in combineCrops are in cropList_out
    combineCrops_dest = combineCrops(:,1) ;
    if length(intersect(cropList_out, combineCrops_dest)) ~= length(combineCrops_dest)
        error('Missing %d expected burnIn destination crops from cropList_out', ...
            length(combineCrops_dest) - length(intersect(cropList_out, combineCrops_dest)))
    end
        
end

Ncrops_in = length(cropList_in) ;
Ncrops_out = length(cropList_out) ;

% Get output file header and other cropList-derived vars
[cropIrrList_mid, ~, ~, ...
    varNames_mid, varNames_mid_allN, varNames_mid_overlapN] = ...
    process_cropLists(cropList_mid, irrList_out, ...
    Nlist_out, Nlist_emu, Nlist_lpj) ;
[cropIrrList_out, header_out, format_out, ...
    varNames_out, varNames_out_allN, varNames_out_overlapN] = ...
    process_cropLists(cropList_out, irrList_out, ...
    Nlist_out, Nlist_emu, Nlist_lpj) ;


%% Loop though GGCM emulators

e2p__child


%% FUNCTIONS

function [cropIrrList_this, header_this, format_this, ...
    varNames_this, varNames_this_allN, varNames_this_overlapN] = ...
    process_cropLists(cropList_this, irrList_this, ...
    Nlist_this, Nlist_emu, Nlist_lpj)

cropIrrList_this = [ ...
    strcat(cropList_this, irrList_this{1}) ...
    strcat(cropList_this, irrList_this{2})] ;
header_this = 'Lon\tLat' ;
format_this = '%0.2f\t%0.2f' ;
varNames_this = {} ;
varNames_this_allN = {} ;
varNames_this_overlapN = {} ;
for c = 1:length(cropIrrList_this)
    thisCropIrr = cropIrrList_this{c} ;
    for n = 1:length(Nlist_this)
        thisN_num = Nlist_this(n) ;
        thisN = pad(num2str(thisN_num), 4, 'left', '0') ;
        thisCropIrrN = [thisCropIrr thisN] ;
        varNames_this{end+1} = thisCropIrrN ; %#ok<AGROW>
        if any(Nlist_emu == thisN_num)
            varNames_this_overlapN{end+1} = thisCropIrrN ; %#ok<AGROW>
        end
        header_this = [header_this '\t' thisCropIrrN] ; %#ok<AGROW>
        format_this = [format_this '\t%0.3f'] ; %#ok<AGROW>
    end
    for n = 1:length(Nlist_lpj)
        thisN_num = Nlist_lpj(n) ;
        thisN = pad(num2str(thisN_num), 4, 'left', '0') ;
        thisCropIrrN = [thisCropIrr thisN] ;
        varNames_this_allN{end+1} = thisCropIrrN ; %#ok<AGROW>
    end
end
header_this = [header_this '\n'] ;
format_this = [format_this '\n'] ;

end

