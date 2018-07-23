function total_fao_Ccy = ...
            FAO_to_Ccy2_totalOnly(fao,in2out_key,listCountries_map_present,...
            listCrops_fa2i,listCrops_fa2o,listYears_fao,...
            verbose,ignoreInfYield,ignoreNoData,faoCommBalElement)
            

% Extract to separate tables
total_fao = fao(exact_string_in_cellarray(fao.ElementName,faoCommBalElement,true,true),:) ;

% Create country-crop-year arrays
Ncountries = length(listCountries_map_present) ;
Nyears_fao = length(listYears_fao) ;
Ncrops_in = length(listCrops_fa2i) ;
Ncrops_out = length(listCrops_fa2o) ;
total_fao_Ccy = nan(Ncountries,Ncrops_out,Nyears_fao) ;
found_in_fao = false(size(listCountries_map_present)) ;
% found_crops = false(
for c = 1:Ncrops_in
    thisCrop_FAO = listCrops_fa2i{c} ;
    %%%
    disp(['thisCrop_FAO = ' thisCrop_FAO])
    if verbose ; disp(thisCrop_FAO) ; end
    
    % Is this crop type in the in2out_key?
    is_in_key = false ;
    for i = 1:Ncrops_out
        tmp = find(strcmp(in2out_key{i},thisCrop_FAO),1) ;
        if ~isempty(tmp)
            is_in_key = true ;
            break   % Loop thru LPJG crop types
            % i gets preserved
        end
    end         % Loop thru LPJG crop types
    clear tmp
    
    % If this crop type not found in in2out_key, skip it
    if ~is_in_key
        disp([thisCrop_FAO ' not found in in2out_key. Skipping.'])
        continue
    end
    
    % Otherwise, get FAO data for this crop.
    total_thisCrop = total_fao(strcmp(total_fao.ItemName,thisCrop_FAO),:) ;
    
    % Sanity check
    if isempty(total_thisCrop)
        if ignoreNoData
            warning(['isempty(total_thisCrop): ' thisCrop_FAO '. Skipping.'])
            continue
        else
            error(['isempty(total_thisCrop): ' thisCrop_FAO])
        end
    end
    
    % Put into _Ccy table
    for C = 1:Ncountries
        thisCountry = listCountries_map_present{C} ;
        
        % Total production
        thisCountry_total_indices = find(strcmp(total_thisCrop.AreaName,thisCountry)) ;
        if any(thisCountry_total_indices)
             % Get this country's production of this crop
            tmp_thisCountry = total_thisCrop(thisCountry_total_indices,:) ;
            % Get this country's existing production in this category
            tmp = squeeze(total_fao_Ccy(C,i,ismember(listYears_fao,tmp_thisCountry.Year))) ;
            tmp(isnan(tmp)) = 0 ;   % If this country-PLUMcrop-year combo has no data yet, convert from NaN to zero for addition
            % Add this crop to existing production
            tmp = tmp + tmp_thisCountry.Value ;
            total_fao_Ccy(C,i,ismember(listYears_fao,tmp_thisCountry.Year)) = tmp ;
            clear tmp*
            found_in_fao(C) = true ;
        end
                        
        clear thisCountry*
    end ; clear C   % Country loop
    clear i *thisCrop*
end ; clear c   % Crop loop

for C = 1:Ncountries
    if ~found_in_fao(C)
        thisCountry = listCountries_map_present{C} ;
        warning([thisCountry ' not found in FAO data.'])
    end
end ; clear C

for c = 1:Ncrops_out
    if ~any(~isnan(total_fao_Ccy(:,c,:)))
        error(['total_fao_Ccy(:,' num2str(c) ',:) has no non-NaN values.'])
    end
end ; clear c