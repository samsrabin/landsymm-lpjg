function [out_array, list2Map] = lpjgu_matlab_maps2xvy(in_maps)

if ndims(in_maps)==3 || ndims(in_maps)==4
    ok_YX = ~isnan(mean(mean(in_maps,4),3)) ;
    list2Map = find(ok_YX) ;
    in_maps_size = size(in_maps) ;
    Nvars = size(in_maps,3) ;
    Nyears = size(in_maps,4) ;
    out_array = reshape(in_maps,[prod(in_maps_size(1:2)) Nvars Nyears]) ;
    out_array = out_array(list2Map,:,:) ;
else
    error('in_maps must have either 3 or 4 dimensions!')
end


end