function [landcover, cropland] = set_cropfracs_zero(landcover, cropland, cropfracs)


cropfracs_in_sum = sum(table2array(cropfracs(:,3:end)),2) ;
if any(cropfracs_in_sum==0 & cropland>0)
    arebad = find(cropfracs_in_sum==0 & cropland>0) ;
    nbad = length(arebad) ;
    max_cropFrac_inBad = max(cropland(arebad)) ;
    warning(['cropfracs sums to zero in ' num2str(nbad) ' crop-containing row(s)! Setting to zero. (Max crop frac = ' num2str(max_cropFrac_inBad) ')'])
    if max_cropFrac_inBad>1e-4
        error('max_cropFrac_inBad is too large for me to be comfortable with this!')
    end
    landcover.CROPLAND(arebad) = 0 ;
    cropland = landcover.CROPLAND ;
    clear arebad nbad max_cropFrac_inBad
end


end