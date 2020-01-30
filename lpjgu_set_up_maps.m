function [out_lons_map,out_lats_map] = lpjgu_set_up_maps(xres,yres,in_lons,in_lats,lat_orient,lon_orient,verboseIfNoMat,verbose)

% Determine orientation, if needed
if isempty(lat_orient)
    if any(in_lats==-90)
        lat_orient = 'lower' ;
    elseif any(in_lats==90)
        lat_orient = 'upper' ;
    elseif any(in_lats-yres/2==-90 | in_lats+yres/2==90)
        lat_orient = 'center' ;
    else
        if any(rem(in_lats,yres)==0)
            lat_orient = 'lower' ;
        else
            lat_orient = 'center' ;
        end
    end
    if verboseIfNoMat || verbose
        disp(['      Assuming lat_orient = ' lat_orient '.'])
    end
end
if isempty(lon_orient)
    if any(in_lons==-180)
        lon_orient = 'left' ;
    elseif any(in_lons==180)
        lon_orient = 'right' ;
    elseif any(in_lons-xres/2==-180 | in_lons+xres/2==180)
        lon_orient = 'center' ;
    else
        if any(rem(in_lons,xres)==0)
            lon_orient = 'left' ;
        else
            lon_orient = 'center' ;
        end
    end
    if verboseIfNoMat || verbose
        disp(['      Assuming lon_orient = ' lon_orient '.'])
    end
end

% Set up maps
if strcmp(lat_orient,'lower')
    lat_min = -90 ;
    lat_max = 90-yres ;
elseif strcmp(lat_orient,'upper')
    lat_min = -90+yres ;
    lat_max = 90 ;
elseif strcmp(lat_orient,'center')
    lat_min = -90+yres/2 ;
    lat_max = 90-yres/2 ;
end
if strcmp(lon_orient,'left')
    lon_min = -180 ;
    lon_max = 180-xres ;
elseif strcmp(lon_orient,'right')
    lon_min = -180+xres ;
    lon_max = 180 ;
elseif strcmp(lon_orient,'center')
    lon_min = -180+xres/2 ;
    lon_max = 180-xres/2 ;
end
lons = lon_min:xres:lon_max ;
lats = lat_min:yres:lat_max ;
out_lons_map = repmat(lons,[length(lats) 1]) ;
out_lats_map = repmat(lats',[1 length(lons)]) ;
end
