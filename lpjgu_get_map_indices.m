function [list_to_map,found] = lpjgu_get_map_indices(in_lons,in_lats,in_lonsMap,in_latsMap,verboseIfNoMat,verbose,in_prec)

if verboseIfNoMat || verbose
    disp('      Getting indices to convert list to map...')
end
Ncells = length(in_lons) ;
lons_vec = in_lonsMap(1,:) ;
lats_vec = in_latsMap(:,1) ;

in_lons = round(in_lons,in_prec) ;
in_lats = round(in_lats,in_prec) ;
lons_vec = round(lons_vec,in_prec) ;
lats_vec = round(lats_vec,in_prec) ;

[~,lons_inds] = ismember(in_lons,lons_vec) ;
[~,lats_inds] = ismember(in_lats,lats_vec) ;
found = ~(lons_inds==0 | lats_inds==0) ;
if any(~found)
    warning([num2str(length(find(~found))) ' cells being ignored.'])
    lons_inds(~found) = [] ;
    lats_inds(~found) = [] ;
end
list_to_map = sub2ind(size(in_lonsMap),lats_inds,lons_inds) ;

% Sanity checks
if any(isnan(list_to_map))
    error('Somehow list_to_map contains NaN.')
end
if length(list_to_map) ~= Ncells-length(find(~found))
    error('length(list_to_map) ~= Ncells-length(find(~found))')
end
end
