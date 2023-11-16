function max_diff = remap_get_xvy_sumLU_max_diff(garr_xvy, target)

max_diff = max(max(abs(target - sum(garr_xvy, 2)))) ;

end