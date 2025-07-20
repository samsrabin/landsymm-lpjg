function S = PLUMharm_harmonize_vegdbare(S, notBare, bareFrac_y0_YX, landArea_YX, ...
    allow_unveg, debug)

if debug
    disp('=== PLUMharm_harmonize_vegdbare ===')
end

bareFrac_y1_YX = sum(S.maps_YXv(:,:,~notBare), 3) ;
if ~isequaln(bareFrac_y0_YX, bareFrac_y1_YX)
    % Print info
    if debug
        debug_S_area(S, landArea_YX, 'before');
        debug_print(bareFrac_y1_YX, bareFrac_y0_YX, landArea_YX);
    end
    
    vegdFrac_y1_YX = sum(S.maps_YXv(:,:,notBare),3) ;
    
    if ~allow_unveg && any(any(vegdFrac_y1_YX > 0 & landArea_YX == 0))
        error('PLUM has vegetated fraction where baseline LU has either no land or all unvegetated land')
    end
    
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
    % Print info
    if debug
        debug_S_area(S, landArea_YX, 'after');
    end
end

if debug
    disp('===================================')
end


end

function debug_print(bareFrac_y1_YX,bareFrac_y0_YX,landArea_YX)

diff_YX = (bareFrac_y1_YX - bareFrac_y0_YX) .* landArea_YX*1e-6 ;
total_bare_y0 = nansum(nansum(bareFrac_y0_YX .* landArea_YX*1e-6)) ;
netDiff = nansum(nansum(diff_YX)) ;
grossDiff = nansum(nansum(abs(diff_YX))) ;
netDiff_pct = 100 * netDiff / total_bare_y0 ;
grossDiff_pct = 100 * grossDiff / total_bare_y0 ;
disp('    Harmonizing vegd+bare fractions')
fprintf('        Sum bare  diff  = %g km2 (%0.2f%% of y0)\n', netDiff, netDiff_pct)
fprintf('        Sum bare |diff| = %g km2 (%0.2f%% of y0)\n', grossDiff, grossDiff_pct)

end

function debug_S_area(S, landArea_YX, txt)
fprintf('S land area %s: %0.4e m2\n', txt, nansum(nansum(sum(S.maps_YXv, 3).*landArea_YX)))
end