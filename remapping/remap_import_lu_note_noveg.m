function remap_import_lu_note_noveg(out_lu)
% Note cells with no vegetated land according to LU dataset

[~, IA] = intersect(out_lu.varNames, {'BARREN'}) ;
if length(IA) ~= 1
    error('Expected 1 match of ''BARREN'' in out_lu.varNames; found %d', ...
        length(IA))
end
bad_x = sum(out_lu.garr_xvy(:,IA,1), 2)==1 ;
Nbad = length(find(bad_x)) ;
if Nbad
    fprintf('%d cells not vegetated\n', Nbad)
end

end