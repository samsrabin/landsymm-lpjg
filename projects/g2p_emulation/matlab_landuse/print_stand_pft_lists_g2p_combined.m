%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make stand and pft list for ins-files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version for crop mappings
% thisVer = 'WithFruitVegSugar_b' ; remapVer = '7b' ;
thisVer = 'WithFruitVeg_sepSugar_sepOil' ; remapVer = '10_g2p' ;
% thisVer = 'WithFruitVeg_sepSugar_sepOil_sepC3' ; remapVer = '11_g2p' ;
% thisVer = 'ggcmi5' ; remapVer = '12_g2p' ;
% thisVer = 'ggcmi5_preBNF' ; remapVer = '13_g2p' ;

include_cropphencol = true ;


%% Set up

addpath(genpath(landsymm_lpjg_path()))

% out_dir = sprintf('/Users/Shared/G2P/inputs/LU/remaps_v%s/',remapVer) ;
out_dir = '/Users/sam/Downloads/create-ins/ins-all' ;
if ~exist(out_dir,'dir')
    mkdir(out_dir) ;
end

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
Nlevels_token = regexprep(num2str(Nlist), '\s*', '-') ;
Nn = length(Nlist) ;


%% Save stand type list

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
            generate_stand_file_text(thisCropIrrN, thisCropIrrN, isIrr, fid, n, ...
                [], [])

            % Simple
            generate_stand_file_text(thisCropIrrN, thisCrop, isIrr, fids, n, ...
                Nformat_appfert, thisN)
        end
    end
end

fclose(fid) ;
fclose(fids) ;
fclose(fids_noForPotYield) ;


%% Construct CFT type list

out_file_simple = sprintf('%s/crop_n_pftlist.remap%s.ins', out_dir, remapVer) ;
fid = fopen(out_file_simple, 'w') ;
out_file_simple = sprintf('%s/crop_n_pftlist.simplePFT.remap%s.ins', out_dir, remapVer) ;
fids = fopen(out_file_simple, 'w') ;

for c = 1:Ncrops
    thisCrop = cropList{c} ;
    thisCFT = get_cft_from_crop(thisCrop) ;
    construct_cft_type_list(irrList, thisCrop, thisCFT, Nformat_stand, ...
        Nformat_appfert, include_cropphencol, Nlist, cropList, fid) ;
    construct_cft_type_list_simple(thisCrop, thisCFT, include_cropphencol, cropList, fids)
end

fclose(fid) ;
 

%% FUNCTIONS

function generate_stand_file_text(thisCrop_st, thisCrop_pft, isIrr, fid, n, ...
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

function construct_cft_type_list( ...
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


function construct_cft_type_list_simple( ...
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


function cft = get_cft_from_crop(crop)

switch crop
    case {'CerealsC3', 'CerealsC3w'}
        cft = 'TeWW' ;
    case {'CerealsC3s', 'StarchyRoots', 'FruitAndVeg', 'Sugar', ...
            'Sugarbeet', 'OilOther', 'ExtraCrop'}
        cft = 'TeSW' ;
    case {'CerealsC4', 'Sugarcane'}
        cft = 'TeCo' ;
    case {'Rice'}
        cft = 'TrRi' ;
    case {'Oilcrops', 'Pulses', 'OilNfix'}
        cft = 'TeSo' ;
    otherwise
        error('Crop %s not recognized in get_cft_from_crop()', crop)
end


end