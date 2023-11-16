function S = PLUMharm_harmonize_vegdbare(S, notBare, bareFrac_y0_YX, landArea_YX)

bareFrac_y1_YX = sum(S.maps_YXv(:,:,~notBare), 3) ;
if ~isequaln(bareFrac_y0_YX, bareFrac_y1_YX)
    diff_YX = (bareFrac_y1_YX - bareFrac_y0_YX) .* landArea_YX*1e-6 ;
    total_bare_y0 = nansum(nansum(bareFrac_y0_YX .* landArea_YX*1e-6)) ;
    netDiff = nansum(nansum(diff_YX)) ;
    grossDiff = nansum(nansum(abs(diff_YX))) ;
    netDiff_pct = 100 * netDiff / total_bare_y0 ;
    grossDiff_pct = 100 * grossDiff / total_bare_y0 ;
    disp('    Harmonizing vegd+bare fractions')
    fprintf('        Sum bare  diff  = %g km2 (%0.2f%% of y0)\n', netDiff, netDiff_pct)
    fprintf('        Sum bare |diff| = %g km2 (%0.2f%% of y0)\n', grossDiff, grossDiff_pct)
    vegdFrac_y1_YX = sum(S.maps_YXv(:,:,notBare),3) ;
    vegdFrac_y1_YXrep = repmat(vegdFrac_y1_YX, [1 1 sum(notBare)]) ;
    vegdFrac_y1_YXv = S.maps_YXv(:,:,notBare) ./ vegdFrac_y1_YXrep ;
    vegdFrac_y1_YXv(vegdFrac_y1_YXrep==0) = 0 ;
    S.maps_YXv(:,:,notBare) = vegdFrac_y1_YXv .* (1 - bareFrac_y0_YX) ;
    S.maps_YXv(:,:,~notBare) = bareFrac_y0_YX ;
    tol = 1e-12 ;
    maxdiff = nanmax(nanmax(abs(sum(S.maps_YXv, 3) - 1))) ;
    if maxdiff > tol
        error('Land use fractions don''t sum to 1 within tolerance %g; max abs. difference %g', ...
            tol, maxdiff)
    end
end

end