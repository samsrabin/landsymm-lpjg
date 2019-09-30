function out_struct = lpjgu_matlab_readTable_then2map(in_file,varargin)

% Set up & parse input arguments
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
addOptional(p,'list_to_map_in',[]) ;
addOptional(p,'dataType','double',@isstr) ;
addOptional(p,'in_prec',2,@isint) ;
parse(p,in_file,varargin{:});

xres = p.Results.xres ;
yres = p.Results.yres ;
verbose = p.Results.verbose ;
verboseIfNoMat = p.Results.verboseIfNoMat ;
force_mat_save = p.Results.force_mat_save ;
force_mat_nosave = p.Results.force_mat_nosave ;
list_to_map_in = p.Results.list_to_map_in ;
lat_orient = p.Results.lat_orient ;
lon_orient = p.Results.lon_orient ;
dataType = p.Results.dataType ;
in_prec = p.Results.in_prec ;

if ~isempty(lat_orient) && ~(strcmp(lat_orient,'lower') || strcmp(lat_orient,'center')  || strcmp(lat_orient,'upper'))
    error(['If providing lat_orient, it must be either lower, center, or upper. (' lat_orient ')'])
end
if ~isempty(lon_orient) && ~(strcmp(lon_orient,'left') || strcmp(lon_orient,'center')  || strcmp(lon_orient,'right'))
    error(['If providing lon_orient, it must be either left, center, or right. (' lon_orient ')'])
end

% pFields = fieldnames(p.Results) ;
% Nfields = length(pFields) ;
% for f = 1:Nfields
%     thisField = pFields{f} ;
%     if ~exist(thisField,'var')
%         eval(['global ' thisField ';']) ;
%         eval([thisField ' = p.Results.' thisField ' ;']) ;
%     end
%     clear thisField
% end ; clear f
clear p

if force_mat_nosave && force_mat_save
    warning('force_mat_save and force_mat_nosave can''t both be true. MATfile will not be saved.')
    force_mat_save = false ;
end
if verboseIfNoMat && verbose
    %     warning('Because verboseIfNoMat==true, setting verbose=false.')
    verbose = false ;
end

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


% Display info
[~,NAME,EXT] = fileparts(in_file) ;

in_matfile_maps = [in_file '.maps.mat'] ;
if exist(in_matfile_maps,'file')
    if verbose
        disp([NAME EXT ':'])
        disp('   Loading maps MAT-file...')
    end
    try
        load(in_matfile_maps) ;
    catch ME
        warning('Problem loading MAT file. Will try again (once) in 10 minutes.\nIntercepted error message follows:\n%s\n', ...
            ME.message) ;
        pause(600)
    end
    
    % Fix gridlist, if somehow screwed up
    if isfield(out_struct,'list_to_map') && (isfield(out_struct,'maps_YXvy') || isfield(out_struct,'maps_YXv'))
        if isfield(out_struct,'maps_YXvy')
            theseMaps = out_struct.maps_YXvy ;
        else
            theseMaps = out_struct.maps_YXv ;
        end
        gltestmap = mean(mean(theseMaps,4),3) ;
        possiblenewgl = find(~isnan(gltestmap)) ;
        if ~isequal(sort(possiblenewgl),sort(out_struct.list_to_map))
            consistentnantestmap = sum(sum(isnan(theseMaps),4),3) ;
            if ~isequal(unique(consistentnantestmap(:)),[0;size(theseMaps,3)*size(theseMaps,4)])
                warning(sprintf(['out_struct.list_to_map (N=' num2str(length(out_struct.list_to_map)) ') does not match gltestmap (N=' num2str(length(sort(possiblenewgl))) ').'...
                    '\n...but NaNs are not consistent across variables (and/or years, if included), so not changing.']))
            else
                warning(sprintf(['out_struct.list_to_map (N=' num2str(length(out_struct.list_to_map)) ') does not match gltestmap (N=' num2str(length(sort(possiblenewgl))) ').\n'...
                    '...Fixing.']))
                out_struct.list_to_map = possiblenewgl ;
            end
        end
    end
else
    
    if verboseIfNoMat || verbose
        disp([NAME EXT ':'])
    end
    
    % Read table
    in_table = lpjgu_matlab_readTable(in_file,...
        'verbose',verbose,...
        'verboseIfNoMat',verboseIfNoMat,...
        'dont_save_MAT',force_mat_nosave,...
        'do_save_MAT',force_mat_save) ;
    is_gridlist = size(in_table,2)==2 ;
    
    % Make maps
    if verboseIfNoMat || verbose
        disp('   Making maps...')
    end
    [yearList,multi_yrs,varNames,list_to_map,xres,yres,found] = get_indices(in_table,xres,yres,list_to_map_in,lat_orient,lon_orient,verboseIfNoMat,verbose,in_prec) ;
    if is_gridlist
        mask_YX = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,dataType,verboseIfNoMat,verbose) ;
    else
        if multi_yrs
            maps_YXvy = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,dataType,verboseIfNoMat,verbose) ;
        else
            maps_YXv = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,dataType,verboseIfNoMat,verbose) ;
        end
    end
    
    % Make output structure
    out_struct.list_to_map = list_to_map ;
    if is_gridlist
        out_struct.mask_YX = mask_YX ;
    else
        out_struct.varNames = varNames ;
        if multi_yrs
            out_struct.maps_YXvy = maps_YXvy ;
            out_struct.yearList = yearList ;
        else
            out_struct.maps_YXv = maps_YXv ;
        end
    end
    
    % Save to MAT-file
    if ~force_mat_nosave
        lpjgu_matlab_save_to_matfile(out_struct,in_matfile_maps,force_mat_save,verboseIfNoMat,verbose) ;
    end
    if verboseIfNoMat || verbose
        disp('   Done.')
    end
end

end






function [yearList,multi_yrs,varNames,list_to_map,xres,yres,found] = get_indices(in_table,xres,yres,list_to_map_in,lat_orient,lon_orient,verboseIfNoMat,verbose,in_prec)

% Get table info
in_lons = in_table.Lon ;
in_lats = in_table.Lat ;
in_Ncells = length(in_lats) ;

% Extract one-year lats/lons, if necessary
[~,varNames,multi_yrs,yearList] = get_names(in_table) ;
if multi_yrs
    in_years = in_table.Year ;
    in_lons = in_lons(in_years==min(in_years)) ;
    in_lats = in_lats(in_years==min(in_years)) ;
    in_Ncells = length(in_lats) ;
end

% Sort out map resolution
[xres,yres] = lpjgu_process_resolution(xres,yres,in_lons,in_lats,verboseIfNoMat,verbose) ;

% Get ready for mapping
[lons_map,lats_map] = lpjgu_set_up_maps(xres,yres,in_lons,in_lats,lat_orient,lon_orient,verboseIfNoMat,verbose) ;

% Get indices for mapping
if isempty(list_to_map_in)
    [list_to_map,found] = lpjgu_get_map_indices(in_lons,in_lats,lons_map,lats_map,verboseIfNoMat,verbose,in_prec) ;
else
    if length(in_lons) ~= length(list_to_map_in)
        warning('length(in_lons) ~= length(list_to_map_in)! Ignoring list_to_map_in.')
        [list_to_map,found] = lpjgu_get_map_indices(in_lons,in_lats,lons_map,lats_map,verboseIfNoMat,verbose,in_prec) ;
    else
        list_to_map = list_to_map_in ;
        found = true(size(list_to_map)) ;
    end
end

end


function [colNames,varNames,multi_yrs,yearList] = get_names(in_table)

% Get names of variables other than Lat, Lon, and Year
colNames = in_table.Properties.VariableNames ;
varNames = {} ;
multi_yrs = false ;
for c = 1:length(colNames)
    thisName = colNames{c} ;
    if ~(strcmp(thisName,'Lat') || strcmp(thisName,'Lon') || strcmp(thisName,'Year'))
        varNames{end+1} = thisName ;
    elseif strcmp(thisName,'Year')
        if length(unique(in_table.Year))>1
            multi_yrs = true ;
            yearList = unique(in_table.Year) ;
        end
    end
end
if ~multi_yrs
    yearList = [] ;
end

end


function out_maps = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,...
    dataType, verboseIfNoMat, verbose)

Nlon = 360 / xres ;
Nlat = 180 / yres ;

if is_gridlist
    if verboseIfNoMat || verbose
        disp('      Making mask...')
    end
    out_maps = false(Nlat,Nlon) ;
    out_maps(list_to_map) = true ;
else
    if verboseIfNoMat || verbose
        disp('      Making maps...')
    end
    Nvars = length(varNames) ;
    if multi_yrs
        Nyears = length(yearList) ;
        out_maps = nan(Nlat,Nlon,Nvars,Nyears,dataType) ;
    else
        out_maps = nan(Nlat,Nlon,Nvars,dataType) ;
    end
    for v = 1:Nvars
        thisVar = varNames{v} ;
        if verboseIfNoMat || verbose
            disp(['         ' thisVar ' (' num2str(v) ' of ' num2str(Nvars) ')...'])
        end
        if multi_yrs
            for y = 1:Nyears
                thisYear = yearList(y) ;
                tmp = nan(Nlat,Nlon,dataType) ;
                try
                    tmp(list_to_map) = table2array(in_table(in_table.Year==thisYear,thisVar)) ;
                catch ME
                    x=1;
                end
                out_maps(:,:,v,y) = tmp ;
            end
        else
            tmp = nan(Nlat,Nlon,dataType) ;
            tmp(list_to_map) = table2array(in_table(found,thisVar)) ;
            out_maps(:,:,v) = tmp ;
        end
    end
end


end