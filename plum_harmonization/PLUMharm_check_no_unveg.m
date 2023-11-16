function PLUMharm_check_no_unveg(base, landArea_YX, which_baseline)

if isfield(base, 'maps_YXv')
    bad_YX = sum(base.maps_YXv(:,:,~strcmp(base.varNames,'BARREN')),3)==0 ;
    
else
    bad_YX1y = sum(base.maps_YXvy(:,:,~strcmp(base.varNames,'BARREN'),:),3)==0 ;
    bad_YX = any(bad_YX1y, 4) ;
end
bad_YX = bad_YX & landArea_YX>0 ;

Nbad = length(find(bad_YX)) ;
if Nbad
    keyboard
    error('%s baseline has %d non-vegetated gridcells where landArea>0!', ...
        which_baseline, Nbad)
end

end