function [out_XYm, yearList_out] = get_cmip5_temp(...
    in_dir, thisParam, filepattern, ncvar, yearList_thisPeriod)

% Find file
[status, result] = unix(sprintf( ...
    'grep -he "^param \\"%s" %s/main_plum_lpjg.ins', ...
    thisParam, in_dir)) ;
if status~=0
    error('Error in unix() call')
end
tmp = strsplit(result,'"') ;
if length(thisParam)>=4 && strcmp(thisParam(1:4), 'file')
    warning('Reading NetCDF file instead of the actual binary file!')
    thisParam_file = tmp{4} ;
    thisParam_file = strrep(thisParam_file,'/home/fh1-project-lpjgpi/lr8247','/Users/Shared/lpj-guess') ;
    if ~exist(thisParam_file, 'file')
        error('thisParam_file (%s=%s) not found!', thisParam, thisParam_file)
    end
    nc_dir = dir(thisParam_file) ;
    nc_dir = nc_dir.folder ;
elseif length(thisParam)>=4 && strcmp(thisParam(1:4), 'path')
    nc_dir = tmp{4} ;
    nc_dir = removeslashifneeded(strrep(nc_dir,'/home/fh1-project-lpjgpi/lr8247','/Users/Shared/lpj-guess')) ;
    if ~exist(nc_dir, 'dir')
        error('nc_dir (%s) not found!', nc_dir)
    end
else
    error('thisParam (%s) not recognized', thisParam)
end
file_nc = dir(sprintf('%s/%s', nc_dir, filepattern)) ;
if isempty(file_nc)
    error('No files found matching %s/%s', nc_dir, filepattern)
elseif length(file_nc) > 1
    error('Multiple files matching %s/%s', nc_dir, filepattern)
end
file_nc = sprintf('%s/%s', nc_dir, file_nc.name) ;

% What years are in this file?
[yearList_temp_bl, temp_size] = get_netcdf_yearList(file_nc, ncvar) ;
yearKey_temp_bl_hist = repmat(yearList_temp_bl,[12 1]) ;
yearKey_temp_bl_hist = yearKey_temp_bl_hist(:) ;
ok_time = find(yearKey_temp_bl_hist>=min(yearList_thisPeriod) & yearKey_temp_bl_hist<=max(yearList_thisPeriod)) ;
yearList_out = unique(yearKey_temp_bl_hist(ok_time)) ;
if length(ok_time) ~= length(yearList_out)*12
    error('length(ok_time) ~= length(yearList_out)*12')
end

% Read this file (only years we care about right now)
start = [1 1 1] ;
count = inf(1,3) ;
lonDim = find(temp_size==720) ;
if length(lonDim) ~= 1
    error('length(lonDim) ~= 1')
end
latDim = find(temp_size==360) ;
if length(latDim) ~= 1
    error('length(latDim) ~= 1')
end
timeDim = find(temp_size~=720 & temp_size~=360) ;
start(timeDim) = min(ok_time) ;
count(timeDim) = length(ok_time) ;
out_XYm = ncread(file_nc, ncvar, start, count) ;

% There's probably a way to generalize, but:
if lonDim==2 && latDim==1
    out_XYm = permute(out_XYm, [2 1 3]) ;
end


end



function [yearList, temp_size] = get_netcdf_yearList(file_nc, ncvar)


ncid = netcdf.open(file_nc) ;
timeID = netcdf.inqVarID(ncid,'time') ;
time_units = netcdf.getAtt(ncid, timeID, 'units') ;
if strcmp(time_units, 'day as %Y%m%d.%f')
    nc_time = ncread(file_nc, 'time') ;
    nc_years = floor(nc_time*1e-4) ;
    yearList = unique(nc_years) ;
    Nmonths = length(nc_years) ;
    if Nmonths ~= length(yearList)*12
        error('Nmonths ~= length(yearList)*12')
    end
else
    if length(time_units)~=30
        error('length(time_units)~=30')
    elseif ~strcmp(time_units(1:11), 'days since ')
        error('~strcmp(time_units(1:11), ''days since '')')
    elseif ~strcmp(time_units(16:end), '-01-01 00:00:00')
        error('~strcmp(time_units(16:end), ''-01-01 00:00:00'')')
    end
    time_baseYear = str2double(time_units(12:15)) ;
    time_bounds = ncread(file_nc, 'time_bnds') ;
    if ~isequal(time_bounds(:,1), [0;31])
        error('~isequal(time_bounds(:,1), [0;31])')
    end
    Nmonths = size(time_bounds,2) ;
    yearList = time_baseYear:(time_baseYear + Nmonths/12 - 1) ;
end

finfo = ncinfo(file_nc) ;
temp_size = finfo.Variables(1+netcdf.inqVarID(ncid, ncvar)).Size ;

netcdf.close(ncid) ;

end