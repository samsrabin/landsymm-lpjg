function PLUMharm_check_preserved_global_deltas( ...
    out_y0_a_YXv, out_y1_a_YXv, ...
    out_y0_b_YXv, out_y1_b_YXv, ...
    conserv_tol_pct, varNames, check_name, areaType)
% Possibly redundant with PLUMharm_checkCons* functions, but this provides an alternative
% method for checking conservation of deltas.

delta_2deg = sum(sum(sum(out_y1_a_YXv - out_y0_a_YXv))) ;
delta = sum(sum(sum(out_y1_b_YXv - out_y0_b_YXv))) ;
delta_diff = delta - delta_2deg ;
delta_diff_pct = 100 * delta_diff / abs(delta_2deg) ;
if abs(delta_diff_pct) > conserv_tol_pct
    warning('%s changes ∆ %s area: Diff %g (%0.2f%%, tolerance %0.2f%%)\n', ...
        check_name, areaType, delta_diff, delta_diff_pct, conserv_tol_pct)
end
for v = 1:size(out_y1_b_YXv, 3)
    delta_2deg = sum(sum(out_y1_a_YXv(:,:,v) - out_y0_a_YXv(:,:,v))) ;
    delta = sum(sum(out_y1_b_YXv(:,:,v) - out_y0_b_YXv(:,:,v))) ;
    delta_diff = delta - delta_2deg ;
    delta_diff_pct = 100 * delta_diff / abs(delta_2deg) ;
    if abs(delta_diff_pct) > conserv_tol_pct
        warning('%s changes ∆ %s area: Diff %g (%0.2f%%, tolerance %0.2f%%)\n', ...
            check_name, varNames{v}, delta_diff, delta_diff_pct, conserv_tol_pct)
    end
end

end