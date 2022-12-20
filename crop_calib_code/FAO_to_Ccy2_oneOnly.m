function out_array = ...
            FAO_to_Ccy2_oneOnly(fao, in2out_key, ...
            listCrops_fa2i, listCrops_fa2o, listYears_fao, ...
            area_or_prod, ignoreNoData, faoCommBalElement, ...
            varargin)

need_countries = false ;
if ~isempty(varargin)
    if length(varargin)~=1
        error('0 or 1 optional argument required: listCountries_map_present')
    end
    listCountries_map_present = varargin{1} ;
    need_countries = true ;
end

% Change listCrops_fa2i and in2out_key if necessary
if strcmp(area_or_prod, 'area')
    thisElem_subtable = fao(exact_string_in_cellarray(fao.ElementName,'Area harvested',true,true),:) ;
    thisElem_subtable_uniqueitems = unique(thisElem_subtable.ItemName) ;
    if length(find(contains(listCrops_fa2i,'Rice, paddy')))==1 ...
            && ~any(contains(thisElem_subtable_uniqueitems,'Rice, paddy')) ...
            && length(find(contains(thisElem_subtable_uniqueitems,'Rice paddy')))==1
        listCrops_fa2i = strrep(listCrops_fa2i,'Rice, paddy','Rice paddy') ;
    end
    for c = 1:length(in2out_key)
        thisKey = in2out_key{c} ;
        if length(find(contains(thisKey,'Rice, paddy')))==1 ...
                && ~any(contains(thisElem_subtable_uniqueitems,'Rice, paddy')) ...
                && length(find(contains(thisElem_subtable_uniqueitems,'Rice paddy')))==1
            thisKey = strrep(thisKey,'Rice, paddy','Rice paddy') ;
            in2out_key{c} = thisKey ;
        end
    end
elseif strcmp(area_or_prod, 'prod')
    thisElem_subtable = fao(exact_string_in_cellarray(fao.ElementName,faoCommBalElement,true,true),:) ;
else
    error('area_or_prod not recognized: %s', area_or_prod)
end

% Create country-crop-year arrays
Nyears_fao = length(listYears_fao) ;
Ncrops_in = length(listCrops_fa2i) ;
Ncrops_out = length(listCrops_fa2o) ;
if need_countries
    Ncountries = length(listCountries_map_present) ;
    found_in_fao = false(size(listCountries_map_present)) ;
    out_Ccy = nan(Ncountries,Ncrops_out,Nyears_fao) ;
else
    out_cy = nan(Ncrops_out,Nyears_fao) ;
end
for c = 1:Ncrops_in
    thisCrop_FAO = listCrops_fa2i{c} ;
    
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
    fprintf('%s: %s <- %s\n', area_or_prod, listCrops_fa2o{i}, thisCrop_FAO)
    
    % Otherwise, get FAO data for this crop.
    thisCrop_subtable = thisElem_subtable(strcmp(thisElem_subtable.ItemName,thisCrop_FAO),:) ;
    
    % Sanity check
    if isempty(thisCrop_subtable)
        if ignoreNoData
            warning(['isempty(thisCrop_subtable): ' thisCrop_FAO '. Skipping.'])
            continue
        else
            error(['isempty(thisCrop_subtable): ' thisCrop_FAO])
        end
    end
    
    % Put into _(C)cy table
    if need_countries
        [out_Ccy, found_in_fao] = fill_Ccy( ...
            out_Ccy, found_in_fao, ...
            listCountries_map_present, thisCrop_subtable, i, ...
            listYears_fao) ;
    else
%         keyboard
        for y = 1:Nyears_fao
            thisYear = listYears_fao(y) ;
            out_cy(i,y) = max(0,out_cy(i,y)) + sum(thisCrop_subtable.Value(thisCrop_subtable.Year==thisYear)) ;
        end
    end
    
    clear i *thisCrop*
end ; clear c   % Crop loop

% Sanity check(s) and output assignment
if need_countries
    for C = 1:Ncountries
        if ~found_in_fao(C)
            thisCountry = listCountries_map_present{C} ;
            warning([thisCountry ' not found in FAO data.'])
        end
    end ; clear C
    for c = 1:Ncrops_out
        if ~any(~isnan(out_Ccy(:,c,:)))
            error(['out_Ccy(:,' num2str(c) ',:) has no non-NaN values.'])
        end
    end ; clear c
    out_array = out_Ccy ;
else
    for c = 1:Ncrops_out
        if ~any(~isnan(out_cy(c,:)))
            error(['out_cy(' num2str(c) ',:) has no non-NaN values.'])
        end
    end ; clear c
    out_array = out_cy ;
end

end


function [out_Ccy, found_in_fao] = fill_Ccy( ...
    out_Ccy, found_in_fao, ...
    listCountries_map_present, thisCrop_subtable, i, ...
    listYears_fao)

Ncountries = length(listCountries_map_present) ;
for C = 1:Ncountries
    thisCountry = listCountries_map_present{C} ;
    
    % Total {area, production}
    thisCountry_metric_indices = find(strcmp(thisCrop_subtable.AreaName,thisCountry)) ;
    if any(thisCountry_metric_indices)
        % Get this country's {area, production} of this crop
        tmp_thisCountry = thisCrop_subtable(thisCountry_metric_indices,:) ;
        % Get this country's existing {area, production} in this category
        tmp = squeeze(out_Ccy(C,i,ismember(listYears_fao,tmp_thisCountry.Year))) ;
        tmp(isnan(tmp)) = 0 ;   % If this country-PLUMcrop-year combo has no data yet, convert from NaN to zero for addition
        % Add this crop to existing {area, production}
        tmp = tmp + tmp_thisCountry.Value ;
        out_Ccy(C,i,ismember(listYears_fao,tmp_thisCountry.Year)) = tmp ;
        clear tmp*
        found_in_fao(C) = true ;
    end
    
    clear thisCountry*
end

end




