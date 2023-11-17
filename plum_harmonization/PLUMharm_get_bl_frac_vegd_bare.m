function [base_vegd_YX, base_bare_YX, base_vegdFrac_YX, base_bareFrac_YX] = ...
    PLUMharm_get_bl_frac_vegd_bare(base, landArea_YX, notBare, LUnames, doHarm)

% Get baseline fraction that is vegetated, barren
disp('    Get baseline fraction that is vegetated, barren')
if doHarm
    base_vegd_YX = sum(base.maps_YXv(:,:,notBare),3) ;
    base_bare_YX = base.maps_YXv(:,:,strcmp(LUnames,'BARREN')) ;
else
    base_vegd_YX = sum(base.maps_YXvy(:,:,notBare,1),3) ;
    base_bare_YX = base.maps_YXvy(:,:,strcmp(LUnames,'BARREN'),1) ;
end
base_vegdFrac_YX = base_vegd_YX ./ landArea_YX ;
base_bareFrac_YX = base_bare_YX ./ landArea_YX ;

end
