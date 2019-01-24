function struct_out = import_croppy(file_in, yearList)

struct_out = lpjgu_matlab_readTable_then2map(file_in) ;
if ~isfield(struct_out,'maps_YXvy')
    [struct_out.maps_YXv, struct_out.varNames] = trim_miscanthus(struct_out.maps_YXv, struct_out.varNames) ;
    tmp_YXv = struct_out.maps_YXv ;
    struct_out = rmfield(struct_out,'maps_YXv') ;
    struct_out.yearList = yearList ;
    struct_out.maps_YXvy = repmat(tmp_YXv,[1 1 1 length(struct_out.yearList)]) ;
    clear tmp
else
    struct_out = trim_years(struct_out, yearList) ;
    [struct_out.maps_YXvy, struct_out.varNames] = trim_miscanthus(struct_out.maps_YXvy, struct_out.varNames) ;
end

end
