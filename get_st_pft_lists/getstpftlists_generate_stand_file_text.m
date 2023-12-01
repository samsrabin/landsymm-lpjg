function getstpftlists_generate_stand_file_text(thisCrop_st, thisCrop_pft, isIrr, fid, n, ...
    Nformat_appfert, thisN)

% Generate stand file text
fprintf(fid, 'st "%s" (\n', thisCrop_st) ;
fprintf(fid, '\tcrop_stand\n\tstinclude 1\n') ;
fprintf(fid, '\tpft "%s"\n', thisCrop_pft) ;
if isIrr
    fprintf(fid, '\thydrology "irrigated_sat"') ;
    fprintf(fid, '\n') ;
end
if n > 0
    fprintf(fid, '\tisforpotyield 1\n') ;
    if ~isempty(Nformat_appfert)
        fprintf(fid, ['\tN_appfert_mt ' Nformat_appfert '\n'], thisN*1e-4) ;
    end
end
fprintf(fid, ')\n\n') ;


end