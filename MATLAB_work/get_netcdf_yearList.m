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
    try
        time_bounds = ncread(file_nc, 'time_bnds') ;
        if ~isequal(time_bounds(:,1), [0;31])
            error('~isequal(time_bounds(:,1), [0;31])')
        end
    catch ME
        if ~strcmp(ME.message,'Could not find variable or group ''time_bnds'' in file.')
            error(ME.message)
        end
    end

    
    Nmonths = size(time_bounds,2) ;
    yearList = time_baseYear:(time_baseYear + Nmonths/12 - 1) ;
end

finfo = ncinfo(file_nc) ;
temp_size = finfo.Variables(1+netcdf.inqVarID(ncid, ncvar)).Size ;

netcdf.close(ncid) ;

end