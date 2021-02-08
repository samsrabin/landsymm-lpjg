%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make stand and pft list for ins-files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version for crop mappings
% thisVer = 'WithFruitVegSugar_b' ; remapVer = '7b' ;
% thisVer = 'WithFruitVeg_sepSugar_sepOil' ; remapVer = '10_g2p' ;
thisVer = 'WithFruitVeg_sepSugar_sepOil_sepC3' ; remapVer = '11_g2p' ;


%% Set up

cd '/Users/sam/Documents/git_repos/g2p_emulation/matlab_landuse'
addpath(genpath(pwd))

out_dir = sprintf('/Volumes/Reacher/G2P/inputs/LU/remaps_v%s/',remapVer) ;
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

for c = 1:Ncrops
    thisCrop = cropList{c} ;
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
            
            % Generate stand file text
            fprintf(fid, 'st "%s" (\n', thisCropIrrN) ;
            fprintf(fid, '\tcrop_stand\n\tstinclude 1\n') ;
%             fprintf(fid, '\t!rotation 1\n') ;
            fprintf(fid, '\tpft "%s"\n', thisCropIrrN) ;
%             fprintf(fid, ['\tN_appfert ' Nformat_appfert '\n'], thisN*1e-4) ;
            if isIrr
                fprintf(fid, '\thydrology "irrigated_sat"') ;
                fprintf(fid, '\n') ;
            end
            if n > 0
                fprintf(fid, '\tisforpotyield 1\n') ;
            end
            fprintf(fid, ')\n\n') ;
        end
    end
end

fclose(fid) ;


%% Construct CFT type list

out_file = sprintf('%s/crop_n_pftlist.remap%s.ins', out_dir, remapVer) ;
fid = fopen(out_file, 'w') ;

for c = 1:Ncrops
    thisCrop = cropList{c} ;
    thisCFT = get_cft_from_crop(thisCrop) ;
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
            fprintf(fid, ')\n\n') ;
        end
    end
end

fclose(fid) ;
 

%% FUNCTIONS

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