function [croparea_fao_Ccy,total_fao_Ccy,yield2_fao_Ccy] = ...
            FAO_to_Ccy2(fao,in2out_key,listCountries_map_present,...
            listCrops_fa2i,listCrops_fa2o,listYears_fao,...
            verbose,ignoreInfYield,ignoreNoData)
            

% Extract to separate tables
croparea_fao = fao(exact_string_in_cellarray(fao.ElementName,'Area harvested',true,true),:) ;
total_fao = fao(exact_string_in_cellarray(fao.ElementName,'Production',true,true),:) ;

% Create country-crop-year arrays
Ncountries = length(listCountries_map_present) ;
Nyears_fao = length(listYears_fao) ;
Ncrops_in = length(listCrops_fa2i) ;
Ncrops_out = length(listCrops_fa2o) ;
croparea_fao_Ccy = nan(Ncountries,Ncrops_out,Nyears_fao) ;
total_fao_Ccy = nan(Ncountries,Ncrops_out,Nyears_fao) ;
found_in_fao_Cc = false(length(listCountries_map_present), Ncrops_in) ;
% found_crops = false(
for c = 1:Ncrops_in
    thisCrop_FAO = listCrops_fa2i{c} ;
    %%%
    if verbose ; disp(['thisCrop_FAO = ' thisCrop_FAO]) ; end
    
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
    croparea_thisCrop = croparea_fao(strcmp(croparea_fao.ItemName,thisCrop_FAO),:) ;
    total_thisCrop = total_fao(strcmp(total_fao.ItemName,thisCrop_FAO),:) ;

    % Sanity checks
    if isempty(croparea_thisCrop)
        if ignoreNoData
            warning(['isempty(croparea_thisCrop): ' thisCrop_FAO '. Skipping.'])
            continue
        else
            error(['isempty(croparea_thisCrop): ' thisCrop_FAO])
        end
    end
    if isempty(total_thisCrop)
        if ignoreNoData
            warning(['isempty(total_thisCrop): ' thisCrop_FAO '. Skipping.'])
            continue
        else
            error(['isempty(total_thisCrop): ' thisCrop_FAO])
        end
    end

    % Only use country-years for which we have both area and production.
    tmpArea = get_table_for_intersect_check(croparea_thisCrop) ;
    tmpProd = get_table_for_intersect_check(total_thisCrop) ;
    [~, Iarea, Iprod] = intersect(tmpArea, tmpProd, 'rows') ;
    if isempty(Iarea)
        if ignoreNoData
            warning(['No country-years have both area and production records for ' thisCrop_FAO '. Skipping.'])
            continue
        else
            error(['No country-years have both area and production records for ' thisCrop_FAO])
        end
    end
    if length(Iarea) ~= size(croparea_thisCrop, 1)
%         warning('%d country-years for %s dropped from FAOSTAT area because they''re not in FAOSTAT production', ...
%             size(croparea_thisCrop, 1) - length(Iarea), thisCrop_FAO)
        croparea_thisCrop = croparea_thisCrop(Iarea,:) ;
    end
    if length(Iprod) ~= size(total_thisCrop, 1)
%         warning('%d country-years for %s dropped from FAOSTAT production because they''re not in FAOSTAT area', ...
%             size(total_thisCrop, 1) - length(Iprod), thisCrop_FAO)
        total_thisCrop = total_thisCrop(Iprod,:) ;
    end

    % FAO data shouldn't have production where there is no area
    is_inf_yield = find(croparea_thisCrop.Value == 0 & total_thisCrop.Value>0) ;
    if ~isempty(is_inf_yield)
        warning('FAOSTAT %s has %d country-years with positive production but zero area. Dropping.', ...
            thisCrop_FAO, length(is_inf_yield))
        croparea_thisCrop(is_inf_yield,:) = [] ;
        total_thisCrop(is_inf_yield,:) = [] ;
    end
    
    % Put into _Ccy table
    [croparea_fao_Ccy, total_fao_Ccy, found_in_fao_Cc] = fill_Ccy( ...
        Ncountries, listCountries_map_present, i, ...
        croparea_fao_Ccy, total_fao_Ccy, croparea_thisCrop, total_thisCrop, ...
        found_in_fao_Cc, listYears_fao, c) ;

    if ~isempty(find(total_fao_Ccy>0 & croparea_fao_Ccy==0, 1))
        error('Some country-year(s) in FAOSTAT %s has positive production but zero area!', ...
            thisCrop_FAO)
    end
    
    clear i *thisCrop*
end ; clear c   % Crop loop

for C = 1:Ncountries
    thisCountry = listCountries_map_present{C} ;
    if ~any(found_in_fao_Cc(C,:))
        if any(strcmp(unique(fao.AreaName),thisCountry))
            warning('%s not found in FAO data for any included crops (but is in FAO data for something else).', thisCountry )
        else
            warning('%s not found in FAO data.', thisCountry )
        end
%     elseif ~all(found_in_fao_Cc(C,:))
%         warning('%s not found in FAO data for %d/%d crops.', ...
%             thisCountry, Ncrops_in - length(find(found_in_fao_Cc(C,:))), Ncrops_in) ;
    end
end ; clear C

for c = 1:Ncrops_out
    if ~any(~isnan(total_fao_Ccy(:,c,:)))
        warning(['total_fao_Ccy(:,' num2str(c) ',:) has no non-NaN values.'])
    elseif ~any(~isnan(croparea_fao_Ccy(:,c,:)))
        warning(['croparea_fao_Ccy(:,' num2str(c) ',:) has no non-NaN values.'])
    end
end ; clear c

% Reconcile FAO NaNs
any_fao_nan = (isnan(total_fao_Ccy) | isnan(croparea_fao_Ccy)) ;
total_fao_Ccy(any_fao_nan) = NaN ;
croparea_fao_Ccy(any_fao_nan) = NaN ;

% There shouldn't be any FAO production if no FAO data
if any(total_fao_Ccy>0 & croparea_fao_Ccy==0)
    error('At least one cell has production without crop area.')
end

% Calculate yield (tDM/ha)
yield2_fao_Ccy = total_fao_Ccy ./ croparea_fao_Ccy ;

% Sanity checks
if any(isinf(yield2_fao_Ccy(:)))
    if ~ignoreInfYield
        error('At least one member of yield2_fao_Ccy is Inf!')
    else
        warning('At least one member of yield2_fao_Ccy is Inf! Max production = %0.1g. IGNORING', ...
            max(total_fao_Ccy(isinf(yield2_fao_Ccy))))
    end
end


end



function [croparea_fao_Ccy, total_fao_Ccy, found_in_fao_Cc] = fill_Ccy( ...
    Ncountries, listCountries_map_present, i, ...
    croparea_fao_Ccy, total_fao_Ccy, croparea_thisCrop, total_thisCrop, ...
    found_in_fao_Cc, listYears_fao, cropIndex)

for C = 1:Ncountries
    thisCountry = listCountries_map_present{C} ;
    
    % Crop area
    thisCountry_croparea_indices = find(strcmp(croparea_thisCrop.AreaName,thisCountry)) ;
    if any(thisCountry_croparea_indices)
        % Get this country's area of this crop
        tmp_thisCountry = croparea_thisCrop(thisCountry_croparea_indices,:) ;
        % Get this country's existing area in this category
        tmp = squeeze(croparea_fao_Ccy(C,i,ismember(listYears_fao,tmp_thisCountry.Year))) ;
        tmp(isnan(tmp)) = 0 ;   % If this country-PLUMcrop-year combo has no data yet, convert from NaN to zero for addition
        % Add this crop to existing area
        tmp = tmp + tmp_thisCountry.Value ;
        croparea_fao_Ccy(C,i,ismember(listYears_fao,tmp_thisCountry.Year)) = tmp ;
        clear tmp*
        found_in_fao_Cc(C,cropIndex) = true ;
    end
    
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
        found_in_fao_Cc(C,cropIndex) = true ;
    end
    
    clear thisCountry*
end ; clear C   % Country loop


end


function T = get_table_for_intersect_check(T)

colNames = {'AreaName', 'ItemName', 'Year'} ;
[C, ~, IB] = intersect(colNames, T.Properties.VariableNames, 'stable') ;
if ~isequal(C, colNames)
    error('AreaName, ItemName, and Year not all found in FAO table')
end
T = T(:,IB) ;

end















