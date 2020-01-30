function map_YXz = lpjgu_xz_to_YXz(in_matrix, map_size, list2map)

if ~ismatrix(in_matrix) || min(size(in_matrix))==1
    error('in_matrix must be a matrix')
elseif size(in_matrix,1) ~= length(list2map)
    error('in_matrix (size1 %d) must have same first-dim size as list2map (size1 %d)', size(in_matrix,1), length(list2map))
end

Nreps = size(in_matrix,2) ;
list2map_all = lpjgu_get_list2map_all(list2map, map_size, Nreps) ;

map_YXz = nan([map_size Nreps]) ;
map_YXz(list2map_all) = in_matrix(:) ;


end