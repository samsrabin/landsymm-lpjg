function out_struct = lpjgu_matlab_readTable_then2map(in_file,varargin)

% Set up & parse input arguments
p = inputParser ;
addRequired(p,'in_file',@ischar) ;
addOptional(p,'xres',NaN,@isnumeric) ;
addOptional(p,'yres',NaN,@isnumeric) ;
addOptional(p,'verbose',true,@islogical) ;
addOptional(p,'verboseIfNoMat',false,@islogical) ;
addOptional(p,'force_mat_save',false,@islogical) ;
addOptional(p,'force_mat_nosave',false,@islogical) ;
addOptional(p,'list_to_map_in',[]) ;
parse(p,in_file,varargin{:});

xres = p.Results.xres ;
yres = p.Results.yres ;
verbose = p.Results.verbose ;
verboseIfNoMat = p.Results.verboseIfNoMat ;
force_mat_save = p.Results.force_mat_save ;
force_mat_nosave = p.Results.force_mat_nosave ;
list_to_map_in = p.Results.list_to_map_in ;

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
    error('force_mat_save and force_mat_nosave can''t both be true.')
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

% Display info
[~,NAME,EXT] = fileparts(in_file) ;

in_matfile_maps = [in_file '.maps.mat'] ;
if exist(in_matfile_maps,'file')
    if verbose
        disp([NAME EXT ':'])
        disp('   Loading maps MAT-file...')
    end
    load(in_matfile_maps)
    
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
    [yearList,multi_yrs,varNames,list_to_map,xres,yres,found] = get_indices(in_table,xres,yres,list_to_map_in,verboseIfNoMat,verbose) ;
    if is_gridlist
        mask_YX = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,verboseIfNoMat,verbose) ;
    else
        if multi_yrs
            maps_YXvy = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,verboseIfNoMat,verbose) ;
        else
            maps_YXv = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,verboseIfNoMat,verbose) ;
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
        save_to_matfile(out_struct,in_matfile_maps,force_mat_save,verboseIfNoMat,verbose) ;
    end
    if verboseIfNoMat || verbose
        disp('   Done.')
    end
end

end






function [yearList,multi_yrs,varNames,list_to_map,xres,yres,found] = get_indices(in_table,xres,yres,list_to_map_in,verboseIfNoMat,verbose)

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
[xres,yres] = process_resolution(xres,yres,in_lons,in_lats,verboseIfNoMat,verbose) ;

% Get ready for mapping
[lons_map,lats_map] = set_up_maps(xres,yres,in_lons,in_lats) ;

% Get indices for mapping
if isempty(list_to_map_in)
    [list_to_map,found] = get_map_indices(in_lons,in_lats,lons_map,lats_map,verboseIfNoMat,verbose) ;
else
    if length(in_lons) ~= length(list_to_map_in)
        error('length(in_lons) ~= length(list_to_map_in)')
    end
    list_to_map = list_to_map_in ;
    found = true(size(list_to_map)) ;
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


function [xres_out,yres_out] = process_resolution(xres_in,yres_in,...
    in_lons,in_lats,verboseIfNoMat,verbose)

if xres_in>0 && yres_in>0
    xres_out = xres_in ;
    yres_out = yres_in ;
elseif xres_in>0 && ~(yres_in>0)
    % If only xres provided, set yres to same value
    xres_out = xres_in ;
    yres_out = xres_in ;
elseif ~(xres_in>0) && yres_in>0
    % If only yres provided, set xres to same value
    xres_out = yres_in ;
    yres_out = yres_in ;
else
    % Determine X and Y resolution
    if verboseIfNoMat || verbose
        disp('      Determining X and Y resolution...')
    end
%     minX = Inf ;
%     minY = Inf ;
%     Ncells = length(in_lons) ;
%     progress = 0 ;
%     progress_step_pct = 25 ;
%     tic ;
%     for c = 1:Ncells
%         this_lon = in_lons(c) ;
%         this_lat = in_lats(c) ;
%         if ~(yres_in>0)
%             those_lats = in_lats(in_lons==this_lon & in_lats~=this_lat) ;
%             those_lats_diff = abs(those_lats - this_lat) ;
%             minY = min([minY min(those_lats_diff)]) ;
%         end
%         if ~(xres_in>0)
%             those_lons = in_lons(in_lats==this_lat & in_lons~=this_lon) ;
%             those_lons_diff = abs(those_lons - this_lon) ;
%             minX = min([minX min(those_lons_diff)]) ;
%         end
%         
%         % Update progress
%         if rem(c,ceil(Ncells*progress_step_pct/100))==0
%             progress = progress + progress_step_pct ;
%             disp(['         ' num2str(progress) '% complete (' toc_hms(toc) ')'])
%         end
%     end

    if ~(xres_in>0)
        unique_lats = unique(in_lats) ;
        minX = Inf ;
        for L = 1:length(unique_lats)
            this_lat = unique_lats(L) ;
            those_lons = sort(in_lons(in_lats==this_lat)) ;
            those_lons_diff = abs(those_lons(1:end-1) - those_lons(2:end)) ;
            minX = min([minX min(those_lons_diff)]) ;
        end
        xres_out = minX ;
    end
    if ~(yres_in>0)
        unique_lons = unique(in_lons) ;
        minY = Inf ;
        for L = 1:length(unique_lons)
            this_lon = unique_lons(L) ;
            those_lats = sort(in_lats(in_lons==this_lon)) ;
            those_lats_diff = abs(those_lats(1:end-1) - those_lats(2:end)) ;
            minY = min([minY min(those_lats_diff)]) ;
        end
        yres_out = minY ;
    end
    if verboseIfNoMat || verbose
        disp(['      Assuming X res. = ' num2str(minX) ', Y res. = ' num2str(minY)])
    end
    
end



end


function [out_lons_map,out_lats_map] = set_up_maps(xres,yres,in_lons,in_lats)
% Set up maps
lon_min = -180 ;
lat_min = -90 ;
lon_max = 180-xres ;
lat_max = 90-yres ;
if rem(in_lons(1)*100,50)~=0
    if xres==3 && yres==3 && rem(max(in_lons-0.25),1.5)==0
        % Specifically set up to deal with some weird gridlist
        lon_min = -178.25 ;
        lon_max = 178.75 ;
    elseif 100*xres/2 ~= abs(rem(in_lons(1)*100,50))
        error('How do I deal with this orientation (lon.)?')
    else
        lon_min = lon_min + xres/2 ;
        lon_max = lon_max + xres/2 ;
    end
end
if rem(in_lats(1)*100,50)~=0
    if xres==3 && yres==3 && rem(max(in_lats-0.25),1.5)==0
        % Specifically set up to deal with some weird gridlist
        lat_min = -88.25 ;
        lat_max = 88.75 ;
    elseif 100*yres/2 ~= abs(rem(in_lats(1)*100,50))
        error('How do I deal with this orientation (lat.)?')
    else
        lat_min = lat_min + yres/2 ;
        lat_max = lat_max + yres/2 ;
    end
end
lons = lon_min:xres:lon_max ;
lats = lat_min:yres:lat_max ;
out_lons_map = repmat(lons,[length(lats) 1]) ;
out_lats_map = repmat(lats',[1 length(lons)]) ;
end


function [list_to_map,found] = get_map_indices(in_lons,in_lats,in_lonsMap,in_latsMap,verboseIfNoMat,verbose)

if verboseIfNoMat || verbose
    disp('      Getting indices to convert list to map...')
end
Ncells = length(in_lons) ;
list_to_map = nan(Ncells,1) ;
found = true(size(list_to_map)) ;
progress = 0 ;
progress_step_pct = 25 ;
tic ;
for c = 1:Ncells
    thisIndex = find(in_latsMap==in_lats(c) & in_lonsMap==in_lons(c)) ;
    if ~isempty(thisIndex)
        list_to_map(c) = thisIndex ;
    else
        found(c) = false ;
        if length(find(~found))==1
            warning(['Some cells from input are getting ignored '...
                'at output resolution ' num2str(360/size(in_lonsMap,2)) ' lon. x'...
                num2str(180/size(in_lonsMap,1)) ' lat.'])
        end
    end
    
    % Update progress
    if verboseIfNoMat || verbose
        if rem(c,ceil(Ncells*progress_step_pct/100))==0
            progress = progress + progress_step_pct ;
            disp(['         ' num2str(progress) '% complete (' toc_hms(toc) ')'])
        end
    end
end

% There should be at least one cell mapped!
if ~any(~isnan(list_to_map))
    error('Something went badly wrong.')
end

% Deal with skipped cells
list_to_map(~found) = [] ;

% Sanity checks
if any(isnan(list_to_map))
    error('Somehow list_to_map contains NaN.')
end
if length(list_to_map) ~= Ncells
    error('length(list_to_map) ~= Ncells')
end
end


function out_maps = make_maps(xres,yres,multi_yrs,yearList,varNames,in_table,list_to_map,is_gridlist,found,...
                              verboseIfNoMat,verbose)

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
        out_maps = nan(Nlat,Nlon,Nvars,Nyears,'single') ;
    else
        out_maps = nan(Nlat,Nlon,Nvars,'single') ;
    end
    for v = 1:Nvars
        thisVar = varNames{v} ;
        if verboseIfNoMat || verbose
            disp(['         ' thisVar ' (' num2str(v) ' of ' num2str(Nvars) ')...'])
        end
        if multi_yrs
            for y = 1:Nyears
                thisYear = yearList(y) ;
                tmp = nan(Nlat,Nlon,'single') ;
                try
                    tmp(list_to_map) = table2array(in_table(in_table.Year==thisYear,thisVar)) ;
                catch ME
                    x=1;
                end
                out_maps(:,:,v,y) = tmp ;
            end
        else
            tmp = nan(Nlat,Nlon,'single') ;
            tmp(list_to_map) = table2array(in_table(found,thisVar)) ;
            out_maps(:,:,v) = tmp ;
        end
    end
end

end


function save_to_matfile(out_struct,in_matfile_maps,force_mat_save,verboseIfNoMat,verbose)
if ~isstruct(out_struct)
    error('out_struct must be a struct.')
end
% Ask to save to MAT-file
if force_mat_save
    if verboseIfNoMat || verbose
        disp('         Saving MAT-file...')
    end
    save(in_matfile_maps,'out_struct','-v7.3') ;
else
    ok = false ;
    while ~ok
        if exist(in_matfile_maps,'file')
            disp('      Save, overwriting existing MAT-file? Y or [N]. 10 seconds...')
            default_save = false ;
        else
            disp('      Save to MAT-file? [Y] or N. 10 seconds...')
            default_save = true ;
        end
        dbl = getkeywait_ssr(10) ;
        if (dbl==-1 && default_save) || strcmp(char(dbl),'y') || strcmp(char(dbl),'Y')
            ok = true ;
            disp('         Saving MAT-file...')
            save(in_matfile_maps,'out_struct','-v7.3') ;
        elseif (dbl==-1 && ~default_save) || strcmp(char(dbl),'n') || strcmp(char(dbl),'N')
            ok = true ;
        else
            warning(['Input (' char(dbl) ') not recognized.'])
        end
    end ; clear ok
end

end