function get_st_pft_lists(thisVer, remapVer, include_cropphencol, out_dir)
% Produce LPJ-GUESS ins-files with crop stands and PFTs
%
% ARGUMENTS:
%     thisVer: String given to get_remapv2_keys() in order to retrieve the mapping of crop
%              types between LandSyMM and MIRCA (source of crop fraction data).
%              E.g., 'WithFruitVeg_sepSugar_sepOil'.
%     remapVer: String to be included in output filenames. E.g.,
%                   sprintf('crop_n_stlist.remap%s.ins', remapVer)
%     include_cropphencol: If true, include CFT attribute cropphen_col. Set to true if you
%                          want to use crop calendar/phenology forcing files that don't
%                          have one column for every CFT. (It's rare that you will want to
%                          use such forcing files, so usually leave this as false.)
%     out_dir: Output directory.

%%%%%%%%%%%%%
%%% Setup %%%
%%%%%%%%%%%%%

cropList = get_remapv2_keys(thisVer) ;
cropList{end+1} = 'ExtraCrop' ;
Ncrops = length(cropList) ;

irrList = {'', 'i'} ;
Nirr = length(irrList) ;

Nlist = [10 60 200 1000] ;
Nwidth = ceil(log10(max(Nlist))) + 1 ;
if ~all(isint(Nlist))
    error('Rework function to allow non-integer N levels')
end
Nformat_stand = ['%0' num2str(Nwidth) 'd'] ;
Nformat_appfert = ['%0.' num2str(max(-log10(1e-4*Nlist))) 'f'] ;
Nn = length(Nlist) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save stand type list %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_file = sprintf('%s/crop_n_stlist.remap%s.ins', out_dir, remapVer) ;
fid = fopen(out_file, 'w') ;

out_file_simple = sprintf('%s/crop_n_stlist.simplePFT.remap%s.ins', out_dir, remapVer) ;
fids = fopen(out_file_simple, 'w') ;

out_file_simple_noForPotYield = sprintf('%s/crop_n_stlist.simplePFT.noForPotYield.remap%s.ins', out_dir, remapVer) ;
fids_noForPotYield = fopen(out_file_simple_noForPotYield, 'w') ;

for c = 1:Ncrops
    thisCrop = cropList{c} ;
    for ii = 1:Nirr
        thisIrr = irrList{ii} ;
        isIrr = strcmp(thisIrr, 'i') ;
        if isIrr && strcmp(thisCrop, 'ExtraCrop')
            continue
        end
        thisCropIrr = [thisCrop thisIrr] ;

        % Generate stand file text WITHOUT isforpotyield crops
        fprintf(fids_noForPotYield, 'st "%s" (\n', thisCropIrr) ;
        fprintf(fids_noForPotYield, '\tcrop_stand\n\tstinclude 1\n') ;
        fprintf(fids_noForPotYield, '\tpft "%s"\n', thisCrop) ;
        if isIrr
            fprintf(fids_noForPotYield, '\thydrology "irrigated_sat"') ;
            fprintf(fids_noForPotYield, '\n') ;
        end
        fprintf(fids_noForPotYield, ')\n\n') ;

        for n = 0:Nn
            if n > 0
                if strcmp(thisCrop, 'ExtraCrop')
                    break
                end
                thisN = Nlist(n) ;
                thisCropIrrN = sprintf(['%s' Nformat_stand], thisCropIrr, thisN) ;
            else
                thisCropIrrN = thisCropIrr ;
                thisN = [] ;
            end

            % Complex
            getstpftlists_generate_stand_file_text(thisCropIrrN, thisCropIrrN, isIrr, fid, n, ...
                [], [])

            % Simple
            getstpftlists_generate_stand_file_text(thisCropIrrN, thisCrop, isIrr, fids, n, ...
                Nformat_appfert, thisN)
        end
    end
end

fclose(fid) ;
fclose(fids) ;
fclose(fids_noForPotYield) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Construct CFT type list %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_file = sprintf('%s/crop_n_pftlist.remap%s.ins', out_dir, remapVer) ;
fid = fopen(out_file, 'w') ;
out_file_simple = sprintf('%s/crop_n_pftlist.simplePFT.remap%s.ins', out_dir, remapVer) ;
fids = fopen(out_file_simple, 'w') ;

for c = 1:Ncrops
    thisCrop = cropList{c} ;
    thisCFT = getstpftlists_get_cft_from_crop(thisCrop) ;
    getstpftlists_construct_cft_type_list(irrList, thisCrop, thisCFT, Nformat_stand, ...
        Nformat_appfert, include_cropphencol, Nlist, cropList, fid) ;
    getstpftlists_construct_cft_type_list_simple(thisCrop, thisCFT, include_cropphencol, cropList, fids)
end

fclose(fid) ;



end