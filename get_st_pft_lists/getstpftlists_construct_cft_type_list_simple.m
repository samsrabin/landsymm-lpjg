function getstpftlists_construct_cft_type_list_simple( ...
    thisCrop, thisCFT, include_cropphencol, cropList, fid)

% Generate CFT file text
fprintf(fid, 'pft "%s" (\n', thisCrop) ;
fprintf(fid, '\t%s_nlim\n', thisCFT) ;
fprintf(fid, '\tinclude 1\n') ;
if include_cropphencol
    this_cropphencol = thisCrop ;
    if strcmp(thisCrop, 'ExtraCrop')
        if any(strcmp(cropList, 'CerealsC3s'))
            this_cropphencol = 'CerealsC3s' ;
        elseif any(strcmp(cropList, 'CerealsC3'))
            this_cropphencol = 'CerealsC3' ;
        else
            error('Which crop should I use for cropphen_col for ExtraCrop?')
        end
    end
    fprintf(fid, '\tcropphen_col "%s"\n', this_cropphencol) ;
end
fprintf(fid, ')\n\n') ;


end