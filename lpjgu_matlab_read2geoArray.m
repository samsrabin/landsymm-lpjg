function out_struct = lpjgu_matlab_read2geoArray(in_file,varargin)

% Set up & parse input arguments
is_2elem_cell = @(x) iscell(x) && numel(x)==2 ;
p = inputParser ;
addRequired(p,'in_file',@ischar) ;
addOptional(p,'xres',NaN,@isnumeric) ;
addOptional(p,'yres',NaN,@isnumeric) ;
addOptional(p,'lat_orient','',@isstr) ;
addOptional(p,'lon_orient','',@isstr) ;
addOptional(p,'verbose',false,@islogical) ;
addOptional(p,'verboseIfNoMat',true,@islogical) ;
addOptional(p,'force_mat_save',true,@islogical) ;
addOptional(p,'force_mat_nosave',false,@islogical) ;
addOptional(p,'dataType','double',@isstr) ;
addOptional(p,'in_prec',2,@isint) ;
% addOptional(p,'lonlats_target',[],is_Nx2_array) ;
% addOptional(p,'list2map_target',[]) ;
addOptional(p,'target',{},is_2elem_cell) ;
parse(p,in_file,varargin{:});

pFields = fieldnames(p.Results) ;
Nfields = length(pFields) ;
for f = 1:Nfields
    thisField = pFields{f} ;
    if ~exist(thisField,'var')
        eval(['global ' thisField ';']) ;
        eval([thisField ' = p.Results.' thisField ' ;']) ;
    end
    clear thisField
end ; clear f
clear p pFields

% Check inputs
if ~isempty(lat_orient) && ~(strcmp(lat_orient,'lower') || strcmp(lat_orient,'center')  || strcmp(lat_orient,'upper'))
    error(['If providing lat_orient, it must be either lower, center, or upper. (' lat_orient ')'])
end
if ~isempty(lon_orient) && ~(strcmp(lon_orient,'left') || strcmp(lon_orient,'center')  || strcmp(lon_orient,'right'))
    error(['If providing lon_orient, it must be either left, center, or right. (' lon_orient ')'])
end
if force_mat_nosave && force_mat_save
    warning('force_mat_save and force_mat_nosave can''t both be true. MATfile will not be saved.')
    force_mat_save = false ;
end
if verboseIfNoMat && verbose %#ok<*NODEF>
%     warning('Because verboseIfNoMat==true, setting verbose=false.')
    verbose = false ;
end

% Process target gridlist
lonlats_target = [] ;
list2map_target = [] ;
if ~isempty(target) %#ok<USENS>
    lonlats_target = target{1} ; %#ok<IDISVAR>
    if ~(ismatrix(lonlats_target) && size(lonlats_target,2)==2)
        error('lonlats_target is malformed (must be Nx2 array)')
    end
    list2map_target = target{2} ;
    if ~isvector(list2map_target)
        error('list2map_target is malformed (must be vector)')
    end
end

% Find file
[NAME, EXT] = process_filename(in_file) ;

% Import
if contains(in_file, '.garr.mat') && strcmp('.garr.mat', in_file(end-(length('.garr.mat')-1):end))
    in_matfile_garr = in_file ;
else
    in_matfile_garr = [in_file '.garr.mat'] ;
end
if exist(in_matfile_garr,'file')
    if verbose
        disp([NAME EXT ':'])
        disp('   Loading geoArray MAT-file...')
    end
    try
        load(in_matfile_garr) ; %#ok<LOAD>
    catch ME
        warning('Problem loading MAT file. Will try again (once) in 10 minutes.\nIntercepted error message follows:\n%s\n', ...
            ME.message) ;
        pause(600)
    end
else
    if verboseIfNoMat || verbose
        disp([NAME EXT ':'])
    end
    
    % Read table
    table_in = lpjgu_matlab_readTable(in_file,...
        'verbose',verbose,...
        'verboseIfNoMat',verboseIfNoMat,...
        'dont_save_MAT',true,...
        'do_save_MAT',false) ;
    is_gridlist = size(table_in,2)==2 ;
    if is_gridlist
        error('This file appears to be a gridlist. Will not make a geoArray.')
    end
    
    if verboseIfNoMat || verbose
        disp('    Getting metadata...')
    end
    
    % Get lat/lons and map indices
    lonlats_in = unique(table2array(table_in(:,1:2)), ...
        'rows', 'stable') ;
    out_struct.lonlats = lonlats_in ;
    Ncells = length(lonlats_in) ;
    if ~isempty(lonlats_target)
        if ~isequal(lonlats_target, lonlats_in)
            if ~isequal(sortrows(lonlats_target), sortrows(lonlats_in))
                error('Gridlist mismatch, but just a sort issue: Add code here to sort to match target')
            else
                error('Gridlist mismatch')
            end
        end
        out_struct.list2map = list2map_target ;
    else
        out_struct.list2map = get_indices(lonlats_in,xres,yres,list2map_target,lat_orient,lon_orient,verboseIfNoMat,verbose,in_prec) ;
    end
    
    % Get variable names
    varNames = setdiff(table_in.Properties.VariableNames, ...
        {'Lon','Lat','Year'}, ...
        'stable') ;
    Nvars = length(varNames) ;
    out_struct.varNames = varNames ;
    
    % Get years (if necessary)
    has_years = any(strcmp(table_in.Properties.VariableNames,'Year')) ;
    if has_years
        yearList = unique(table_in.Year) ;
        Nyears = length(yearList) ;
        if length(table_in.Lon) ~= Ncells*Nyears
            error('length(table_in.Lon) ~= Ncells*Nyears')
        end
        out_struct.yearList = yearList ;
    end
    
    % Reshape to array
    if verboseIfNoMat || verbose
        disp('   Reshaping...')
    end
    if ~has_years
        garr_xv = table2array(table_in(:,3:end)) ;
        if ~strcmp(dataType, 'double')
            eval(sprintf('garr_xv = %s(garr_xv) ;', dataType)) ;
        end
        out_struct.garr_xv = garr_xv ;
    else
        garr_yxv = lpjgu_matlab_table2array(table_in(:,4:end), [Nyears Ncells Nvars]) ;
        if ~strcmp(dataType, 'double')
            eval(sprintf('garr_yxv = %s(garr_yxv) ;', dataType)) ;
        end
        out_struct.garr_yxv = garr_yxv ;
    end
    
    % Save to MAT-file
    if ~force_mat_nosave
        lpjgu_matlab_save_to_matfile(out_struct,in_matfile_garr,force_mat_save,verboseIfNoMat,verbose) ;
    end
    if verboseIfNoMat || verbose
        disp('   Done.')
    end
end

end


function [NAME, EXT] = process_filename(in_file)

if strcmp(in_file(end-2:end),'.gz')
    in_file = in_file(1:end-3) ;
elseif strcmp(in_file(end-3:end),'.mat')
    in_file = in_file(1:end-4) ;
end

% If in_file is symlink, replace it with its target
[s,w] = unix(['[[ -L ' in_file ' ]] && echo true']) ;
if s==0 && contains(w,'true') % is symlink
    disp('Symlink; pointing to target instead.')
    [~,w] = unix(['stat -f "%Y" ' in_file]) ;
    in_file = regexprep(w,'[\n\r]+','') ; % Remove extraneous newline
end

% If in_file has wildcard, expand into full filename (fail if not exactly 1
% match)
if contains(in_file,'*')
    filelist = dir(in_file) ;
    if isempty(filelist)
        error('No match found for %s', in_file)
    elseif length(filelist) > 1
        keyboard
        error('More than one match found for %s', in_file)
    else
        tmp = sprintf('%s/%s', filelist(1).folder, filelist(1).name) ;
        fprintf('Resolving\n%s\ninto\n%s\n', in_file, tmp)
        in_file =  tmp ;
    end
end

% Get info
[~,NAME,EXT] = fileparts(in_file) ;

end


function list_to_map = get_indices(lonlats_in,xres,yres,list2map_target,lat_orient,lon_orient,verboseIfNoMat,verbose,in_prec)

% Get table info
in_lons = lonlats_in(:,1) ;
in_lats = lonlats_in(:,2) ;

% Sort out map resolution
[xres,yres] = lpjgu_process_resolution(xres,yres,in_lons,in_lats,verboseIfNoMat,verbose) ;

% Get ready for mapping
[lons_map,lats_map] = lpjgu_set_up_maps(xres,yres,in_lons,in_lats,lat_orient,lon_orient,verboseIfNoMat,verbose) ;

% Get indices for mapping
if isempty(list2map_target)
    list_to_map = lpjgu_get_map_indices(in_lons,in_lats,lons_map,lats_map,verboseIfNoMat,verbose,in_prec) ;
else
    if length(in_lons) ~= length(list2map_target)
        warning('length(in_lons) ~= length(list2map_target)! Ignoring list2map_target.')
        list_to_map = lpjgu_get_map_indices(in_lons,in_lats,lons_map,lats_map,verboseIfNoMat,verbose,in_prec) ;
    else
        list_to_map = list2map_target ;
    end
end

end






