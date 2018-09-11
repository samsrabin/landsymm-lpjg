function bad = PLUMharm_check_conservation(...
    out_y0_mgmt_YXv, out_y0_area_YXv, ...
    mid_y1_mgmt_YXv, mid_y1_area_YXv, ...
    in_y0_mgmt_YXv, in_y0_area_YXv, ...
    in_y1_mgmt_YXv, in_y1_area_YXv, ...
    total_unmet_mgmt_YXv, LPJGcrops, conserv_tol_pct, ...
    check_name)

bad = 0 ;

Ncrops_lpjg = length(LPJGcrops) ;

out_y0_mgmt_YXvTot_YXv = out_y0_mgmt_YXv .* out_y0_area_YXv ;
mid_y1_mgmtTot_YXv = mid_y1_mgmt_YXv .* mid_y1_area_YXv ;
in_y0_mgmt_YXvTot_YXv = in_y0_mgmt_YXv .* in_y0_area_YXv ;
in_y1_mgmt_YXvTot_YXv = in_y1_mgmt_YXv .* in_y1_area_YXv ;
mgmt_d_YXv = in_y1_mgmt_YXvTot_YXv - in_y0_mgmt_YXvTot_YXv ;
for i = 1:Ncrops_lpjg
    mgmt_d2_YX = mid_y1_mgmtTot_YXv(:,:,i) - out_y0_mgmt_YXvTot_YXv(:,:,i) + total_unmet_mgmt_YXv(:,:,i) ;

    mgmt_d_YX = mgmt_d_YXv(:,:,i) ;
    mgmt_d_glob_1 = sum(mgmt_d_YX(:)) ;
    mgmt_d_glob_2 = sum(mgmt_d2_YX(:)) ;
    if isnan(mgmt_d_glob_1)
        error(['NaN produced in check ' check_name ', mgmt_d_glob_1'])
    elseif isnan(mgmt_d_glob_2)
        error(['NaN produced in check ' check_name ', mgmt_d_glob_2'])
    end
    if abs((mgmt_d_glob_2-mgmt_d_glob_1)/mgmt_d_glob_1*100) > conserv_tol_pct
        if sum(sum(mid_y1_mgmtTot_YXv(:,:,i)))==0
            bad = -1 ;
            warning('Global %s mgmt changes are not conserved to within %0.2f percent because they hit zero (check %s)\n', ...
                LPJGcrops{i}, conserv_tol_pct, check_name)
        else
            bad = 1 ;
            warning('Global %s mgmt changes are not conserved to within %0.2f percent! (check %s; err = %0.2f)\n', ...
                LPJGcrops{i}, conserv_tol_pct, check_name, (mgmt_d_glob_2-mgmt_d_glob_1)/mgmt_d_glob_1*100)
        end
    end
end


end