%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process C-table outputs from LPJ-GUESS for PLUM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/Users/Shared/lpj-guess_git-svn_20190828/ssr/landsymm-dev-forestry/landsymm-for-02_bySource' ;

% What year is the last in the time horizon of the project?
runInfo.dest_yN = 2010 ;

% How long is the project time horizon?
runInfo.dest_duration = 30 ;

% How many years does LU stay as source before transition?
% This should always be the same as dest_duration, to ensure that source
% forestry is cut immediately before transition (as long as you have
% manageforest_dec31 enabled).
runInfo.source_duration = runInfo.dest_duration ;


%% Setup

if ~exist(topDir, 'dir')
    error('topDir does not exist: %s', topDir)
end

% Get list of stands (last 4 characters of output dir name)
get_standname_from_dirname = @(x) x(end-3:end) ;
resultDirs = dir(sprintf('%s/out_*', topDir)) ;
resultDirs = {resultDirs.name} ;
standList = cellfun(@(x) x(end-3:end), resultDirs, ...
    'UniformOutput', false) ;
standList(strcmp(standList, 'PLUM')) = [] ;
Nstands = length(standList) ;

% Set up directory for outputs from this script
outDir = sprintf('%s/out_PLUM', topDir) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end


%% Import

runInfo.target = {} ;
runInfo.standList = standList ;
runInfo.standList_files = {} ;
runInfo.yearList = [] ;

disp('Importing from single-source runs:')
for s = 1:Nstands
    thisStand = standList{s} ;
    fprintf('    %s (%d/%d)...\n', thisStand, s, Nstands)
    
    if strcmp(thisStand, 'ntrl')
        thisDir = sprintf('%s/out_spinup_%s', topDir, thisStand) ;
    else
        thisDir = sprintf('%s/out_%s', topDir, thisStand) ;
    end
    
    % Import necessary files
    if s==1
        [tmp_slowprod, runInfo] = import_file('cmass_wood_harv_toprod', thisDir, runInfo) ;
        gridlist.lonlats = runInfo.target{1} ;
        gridlist.list2map = runInfo.target{2} ;
        Ncells = length(gridlist.list2map) ;
    else
        tmp_slowprod = import_file('cmass_wood_harv_toprod', thisDir, runInfo) ;
    end
    tmp_cflux_noslowh = import_file('cflux_noslowh', thisDir, runInfo) ;
    tmp_cflux_slowh = import_file('cflux_slowh', thisDir, runInfo) ;
    
    % Save to big arrays
    if s==1
        % Set up structures
        c_slowprod_luc = rmfield(tmp_slowprod, 'garr_xdZ') ;
        c_slowprod_rot = rmfield(tmp_slowprod, 'garr_xdZ') ;
        c_flux = rmfield(tmp_cflux_noslowh, 'garr_xdZ') ;
        % Set up arrays: gridcell x destination x source
        c_slowprod_luc.garr_xds = nan([Ncells Nstands Nstands]) ;
        c_slowprod_rot.garr_xds = nan([Ncells Nstands Nstands]) ;
        c_flux.garr_xds = nan([Ncells Nstands Nstands]) ;
    end
    
    % Process products
    c_slowprod_luc.garr_xds(:,:,s) = tmp_slowprod.garr_xdZ(:,:,1) ;
    c_slowprod_rot.garr_xds(:,:,s) = tmp_slowprod.garr_xdZ(:,:,2) ;
    clear tmp_slowprod
    
    % Process fluxes
    tmp_cflux_xd = tmp_cflux_noslowh.garr_xdZ + tmp_cflux_slowh.garr_xdZ ;
    tmp_cflux_xd = tmp_cflux_xd - repmat(tmp_cflux_xd(:,s), [1 Nstands]) ;
    c_flux.garr_xds(:,:,s) = tmp_cflux_xd ;
    clear tmp_cflux_noslowh tmp_cflux_slowh
    
    clear tmp*
    
end
disp('Done.')


%% Save

save_file(sprintf('%s/wood_luc_%d-%d.out', ...
    outDir, runInfo.yearList(1), runInfo.yearList(end)), ...
    c_slowprod_luc, gridlist)
save_file(sprintf('%s/wood_rot_%d-%d.out', ...
    outDir, runInfo.yearList(1), runInfo.yearList(end)), ...
    c_slowprod_rot, gridlist)
save_file(sprintf('%s/cflux_net_%d-%d.out', ...
    outDir, runInfo.yearList(1), runInfo.yearList(end)), ...
    c_flux, gridlist)




%% FUNCTIONS

function save_file(file_out, S, gridlist)

S = reshape_struct(S) ;

out_xv = cat(2, gridlist.lonlats, S.garr_xv) ;

% Get output header and data-row format
header_out = 'Lon\tLat' ;
Nyears = length(S.varNames) ;
for v = 1:Nyears
    header_out = [header_out '\t' S.varNames{v}] ; %#ok<AGROW>
end
header_out = [header_out '\n'] ;
format_out = ['%0.2f\t%0.2f' repmat('\t%0.4f', [1 Nyears]) '\n'] ;

% Open file and save header
fid = fopen(file_out,'w') ;
fprintf(fid, header_out) ;

% Write data
fprintf(fid, ...
    format_out, ...
    out_xv' ...
    ) ;
fclose(fid) ;

end


function S = reshape_struct(S)

garr_xds = S.garr_xds ;
S = rmfield(S, 'garr_xds') ;

Ncells = size(garr_xds, 1) ;
standList_in = shiftdim(S.varNames) ;
Nstands = length(standList_in) ;

tmp = repmat(standList_in, [1 5])' ;
standList_out = tmp(:) ;
standList_out = strcat(standList_out, '_to_') ;
try
    standList_out = strcat(standList_out, repmat(standList_in, [Nstands 1])) ;
catch ME
    keyboard
end


S.varNames = standList_out ;
S.garr_xv = reshape(garr_xds, [Ncells Nstands^2]) ;

end


function garr_xdZ = perform_operation(garr_xdy, thisVar)

if contains(thisVar, 'cflux')
    operation = 'sum' ;
elseif strcmp(thisVar, 'cmass_wood_harv_toprod')
    operation = 'first_thenSumRest' ;
else
    error('thisVar %s not recognized', thisVar)
end

if strcmp(operation, 'sum')
    garr_xdZ = sum(garr_xdy, 3) ;
elseif strcmp(operation, 'first_thenSumRest')
    garr_xdZ = cat(3, garr_xdy(:,:,1), sum(garr_xdy(:,:,2:end),3)) ;
else
    error('Operation %s not recognized', operation)
end

end


function [S_out, runInfo] = import_file(thisVar, thisDir, runInfo)

thisFile = sprintf('%s/%s_sts.out', thisDir, thisVar) ;
S_out = lpjgu_matlab_read2geoArray(thisFile, ...
    'force_mat_save', false, ...
    'force_mat_nosave', true, ...
    'verboseIfNoMat', false, ...
    'xres', 0.5, 'yres', 0.5, ...
    'target', runInfo.target) ;

if isempty(runInfo.target)
    runInfo.target = {S_out.lonlats, S_out.list2map} ;
end
if isempty(runInfo.standList_files)
    runInfo.standList_files = S_out.varNames ;
end

% Perform operation across years of interest
if isempty(runInfo.yearList)
    y1 = runInfo.dest_yN - runInfo.dest_duration ;
    yN = y1 + runInfo.dest_duration - 1 ;
    runInfo.yearList = y1:yN ;
end
[~, IA] = intersect(S_out.yearList, runInfo.yearList) ;
if length(IA) ~= length(runInfo.yearList)
    error('Expected to find %d matching years in S_out.yearList; found %d', ...
        runInfo.dest_duration, length(IA))
end
garr_xdy = S_out.garr_xvy(:,:,IA) ;
S_out = rmfield(S_out, 'garr_xvy') ;
S_out.yearList = runInfo.yearList ;
S_out.garr_xdZ = perform_operation(garr_xdy, thisVar) ;
% Rearrange stand list, if needed
if ~isequal(runInfo.standList, S_out.varNames)
    [~, ~, IB] = intersect(runInfo.standList, S_out.varNames, 'stable') ;
    if length(IB) ~= length(S_out.varNames)
        error('varNames mismatch')
    end
    S_out.garr_xdZ = S_out.garr_xdZ(:,IB,:) ;
    S_out.varNames = S_out.varNames(IB) ;
end




end