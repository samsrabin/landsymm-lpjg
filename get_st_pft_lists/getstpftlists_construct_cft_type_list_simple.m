function getstpftlists_construct_cft_type_list_simple( ...
    thisCrop, thisCFT, include_cropphencol, cropList, fid)

% Generate CFT file text
fprintf(fid, 'pft "%s" (\n', thisCrop) ;
fprintf(fid, '\t%s_nlim\n', thisCFT) ;
fprintf(fid, '\tinclude 1\n') ;
if include_cropphencol
    getstpftlists_cropphen_col(thisCrop, cropList, fid)
end
fprintf(fid, ')\n\n') ;


end