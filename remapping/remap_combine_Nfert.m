function in_x = remap_combine_Nfert(thisCrop_in, croparea_in, nfert_dir, gridlist)
% Calculate area-weighted average Nfert (mass/area) for a given crop (thisCrop_in).

% This is where we define what crops from the N fertilization data will be combined. Note
% that there must be one crop name specified from the frac/area dataset for each crop to
% be used from the N fertilization dataset.
if strcmp(thisCrop_in, 'combined_sugars')
    croparea_names = {'Sugarbeet', 'Sugarcane'} ;
    Nfert_names    = {'sugarbeet', 'sugarcane'} ;
else
    error('Not recognized: %s', thisCrop_in)
end

Ncrops_to_combine = length(croparea_names) ;
if length(Nfert_names) ~= Ncrops_to_combine
    error('Different number of members in Nfert_names and croparea_names')
end

Ncells = size(croparea_in.garr_xv, 1) ;
area_xv = nan(Ncells, Ncrops_to_combine) ;
for c = 1:Ncrops_to_combine
    thisCrop_area = croparea_names{c} ;
    is_this_crop = contains(croparea_in.varNames, thisCrop_area) ;
    if sum(is_this_crop) ~= 2
        error('Expected 2 croparea_in.varNames containing %s; found %d', ...
            thisCrop_area, sum(is_this_crop))
    end
    area_xv(:,c) = sum(croparea_in.garr_xv(:,is_this_crop), 2) ;
end

% Get total area, ensuring no NaNs
area_total_x = sum(area_xv, 2) ;
if any(isnan(area_total_x))
    error('Unexpected NaN in crop area')
end

% Avoid dividing by 0 in next step; will result in weights = 0 there
area_total_x(area_total_x==0) = 1 ;

% Calculate weights
weights_xv = area_xv ./ repmat(area_total_x, [1 Ncrops_to_combine]) ;

% Calculate area-weighted average Nfert
in_x = zeros(Ncells, 1) ;
for c = 1:Ncrops_to_combine
    thisCrop_Nfert = Nfert_names{c} ;
    file_N = fullfile(nfert_dir, ...
        sprintf('agmip_%s_apprate_fill_NPK_0.5.nc4', ...
        thisCrop_Nfert)) ;
    map_YX = flipud(transpose(ncread(file_N, 'Napprate'))) ;
    in_x = in_x + weights_xv(:,c).*map_YX(gridlist.list2map) ;
end

end