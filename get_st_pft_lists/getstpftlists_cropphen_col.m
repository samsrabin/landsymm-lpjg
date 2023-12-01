function getstpftlists_cropphen_col(thisCrop, cropList, fid)

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
