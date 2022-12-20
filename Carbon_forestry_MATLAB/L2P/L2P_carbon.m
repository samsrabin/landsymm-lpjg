%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process LPJ-GUESS carbon outputs for PLUM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yearList_proj = 1980:2009 ;

inDir = '/Users/Shared/lpj-guess_git-svn_20190828/ssr/landsymm-dev-forestry/landsymm-for-01/out.test.20210115' ;


%% Setup

% Get output directory
outDir = sprintf('%s/forPLUM', inDir) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

py1 = min(yearList_proj) ;
py2 = max(yearList_proj) ;


%% Import

% C mass
ctotal = lpjgu_matlab_read2geoArray(sprintf('%s/ctotal_sts.out', inDir), ...
    'xres', 0.5, 'yres', 0.5, ...
    'force_mat_save', false, 'force_mat_nosave', true) ;

% Harvest to long-term pools
c2harvslow = lpjgu_matlab_read2geoArray(sprintf('%s/cmass_wood_harv_toprod_sts.out', inDir), ...
    'xres', 0.5, 'yres', 0.5, ...
    'force_mat_save', false, 'force_mat_nosave', true) ;

% Get stand info
stands = rmfield(ctotal, 'garr_xvy') ;

% Check year coverage
if length(intersect(stands.yearList, [py1 py2+1])) ~= 2
    error('Did not find single exact matches for beginning and end of project period')
end

% % Land use, assuming GLOBAL_YEARLY
% file_lu = get_ins_value(sprintf('%s/test.ins', inDir), 'file_lu', true) ;
% lu = readtable(file_lu) ;


%% Process

% Get indices for different LUs
In = get_indices(stands.varNames, 'Natural', [], 1) ;
Ic = get_indices(stands.varNames, 'TeWW', [], 1) ;
If = get_indices(stands.varNames, [], [In;Ic]) ;
Nf = length(find(If)) ;

% Get forest C change
ctotal_chg_f_xv = ...
    ctotal.garr_xvy(:,If,stands.yearList==py2+1) ...
    - ctotal.garr_xvy(:,If,stands.yearList==py1) ;

% Get crop C change
ctotal_chg_c_x = ...
    ctotal.garr_xvy(:,Ic,stands.yearList==py2+1) ...
    - ctotal.garr_xvy(:,Ic,stands.yearList==py1) ;

% Get forest C change MINUS crop C change
ctotal_chg2_f_xv = ...
    ctotal_chg_f_xv ...
    - repmat(ctotal_chg_c_x, [1 Nf 1])


%% Prepare to save

% Get output filename
file_out = sprintf('%s/forPLUM.out', outDir) ;

% Get output header and data-row format
header_out = 'Lon\tLat' ;
for f = 1:Nf
    header_out = [header_out '\t' stands.varNames{If(f)}] ; %#ok<AGROW>
end
header_out = [header_out '\n'] ;
format_out = ['%0.2f\t%0.2f' repmat('\t%0.3f', [1 Nf]) '\n'] ;


%% Save Ctotal

% Open file and save header
fid = fopen(file_out,'w') ;
fprintf(fid, header_out) ;

data_YV = [stands.lonlats ctotal_chg2_f_xv] ;

% Write data
fprintf(fid, ...
    format_out, ...
    data_YV' ...
    ) ;
fclose(fid) ;
disp('Done saving')


%% FUNCTIONS

function I = get_indices(varList, inclList, exclList, varargin)

Nexpected = NaN ;
if ~isempty(varargin)
    Nexpected = varargin{1} ;
    if length(varargin) > 1
        error('get_indices accepts at most 1 option argument (Nexpected)')
    end
end

% In case inclList is empty
inclFound = varList ;
IA = 1:length(inclFound) ;

% Get included variables
if ~isempty(inclList)
    if ischar(inclList) || iscellstr(inclList) %#ok<ISCLSTR>
        [inclFound, IA] = intersect(varList, inclList, 'stable') ;
    elseif all(isint(inclList)) && min(ndims(inclList))==1
        inclFound = varList(inclList) ;
        [~,IA] = intersect(varList, inclFound, 'stable') ;
    else
        error('Format of inclList not recognized')
    end
end

% Remove excluded variables
if ~isempty(exclList)
    if ischar(exclList) || iscellstr(exclList) %#ok<ISCLSTR>
        inclFound = setdiff(inclFound, exclList, 'stable') ;
    elseif all(isint(exclList)) && min(size(exclList))==1
        [~,IA] = setdiff(IA, exclList, 'stable') ;
        inclFound = varList(IA) ;
    else
        error('Format of exclList not recognized')
    end
end

[~, I] = intersect(varList, inclFound, 'stable') ;

if ~isnan(Nexpected) && length(I)~=Nexpected
    error('Found %d, expected %d', length(I), Nexpected)
end

end