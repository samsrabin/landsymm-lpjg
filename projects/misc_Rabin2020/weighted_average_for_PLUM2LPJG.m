function array_out = weighted_average_for_PLUM2LPJG(...
    thisMap, colNames, array_in_Xv, ...
    colNames_cropfracs, cropfracs_in_Xv)

Ncrops_thisMap = length(thisMap) ;

[~,IA_this] = intersect(colNames,strcat(thisMap,'i')) ;
if length(IA_this) ~= Ncrops_thisMap
    error('length(IA_this) ~= Ncrops_thisMap')
end

[~,IA_cropfracs] = intersect(colNames_cropfracs,strcat(thisMap,'i')) ;
if length(IA_cropfracs) ~= Ncrops_thisMap
    error('length(IA_cropfracs) ~= Ncrops_thisMap')
end

% Get weights
cropfrac_total = sum(cropfracs_in_Xv(:,IA_cropfracs),2) ;
cropfrac_total_Xv = repmat(cropfrac_total,[1 Ncrops_thisMap]) ;
weights = cropfracs_in_Xv(:,IA_cropfracs) ./ cropfrac_total_Xv ;
weights(cropfrac_total_Xv==0) = 0 ;
if any(isnan(weights))
    error('any(isnan(weights))')
elseif any(isinf(weights))
    error('any(isinf(weights))')
end

% Get weighted average
theseCols = array_in_Xv(:,IA_this) ;
array_out = sum(theseCols .* weights,2) ;
array_out(cropfrac_total==0) = 0 ;
if any(isnan(array_out))
    error('any(isnan(array_out))')
elseif any(isinf(array_out))
    error('any(isinf(array_out))')
end

end