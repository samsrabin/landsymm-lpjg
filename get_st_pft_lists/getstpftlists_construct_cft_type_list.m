function getstpftlists_construct_cft_type_list( ...
    irrList, thisCrop, thisCFT, Nformat_stand, Nformat_appfert, include_cropphencol, ...
    Nlist, cropList, fid)

Nirr = length(irrList) ;
Nn = length(Nlist) ;

for ii = 1:Nirr
    thisIrr = irrList{ii} ;
    isIrr = strcmp(thisIrr, 'i') ;
    if isIrr && strcmp(thisCrop, 'ExtraCrop')
        continue
    end
    thisCropIrr = [thisCrop thisIrr] ;
    for n = 0:Nn
        if n > 0
            if strcmp(thisCrop, 'ExtraCrop')
                break
            end
            thisN = Nlist(n) ;
            thisCropIrrN = sprintf(['%s' Nformat_stand], thisCropIrr, thisN) ;
        else
            thisCropIrrN = thisCropIrr ;
        end

        % Generate CFT file text
        fprintf(fid, 'pft "%s" (\n', thisCropIrrN) ;
        fprintf(fid, '\t%s_nlim\n', thisCFT) ;
        fprintf(fid, '\tinclude 1\n') ;
        if n > 0
            fprintf(fid, ['\tN_appfert ' Nformat_appfert '\n'], thisN*1e-4) ;
            fprintf(fid, '\tisforpotyield 1\n') ;
        end
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
end


end
