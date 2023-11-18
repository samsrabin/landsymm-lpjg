function bad = PLUMharm_checkCons_mgmt(...
    out_y0_mgmt_YXv, out_y0_area_YXv, ...
    out_y1_mgmt_YXv, out_y1_area_YXv, ...
    in_y0_mgmt_YXv, in_y0_area_YXv, ...
    in_y1_mgmt_YXv, in_y1_area_YXv, ...
    total_unmet_mgmt_YXv, LPJGcrops, conserv_tol_pct, notEnough, ...
    check_name, do_warn)

in_y0_mgmtTot_YXv = in_y0_mgmt_YXv .* in_y0_area_YXv ;
in_y1_mgmtTot_YXv = in_y1_mgmt_YXv .* in_y1_area_YXv ;
out_y0_mgmtTot_YXv = out_y0_mgmt_YXv .* out_y0_area_YXv ;
out_y1_mgmtTot_YXv = out_y1_mgmt_YXv .* out_y1_area_YXv ;

mgmt_d_YXv = in_y1_mgmtTot_YXv - in_y0_mgmtTot_YXv ;

bad = 0 ;
for i = 1:length(LPJGcrops)

    mgmt_d2_YX = out_y1_mgmtTot_YXv(:,:,i) - out_y0_mgmtTot_YXv(:,:,i) + total_unmet_mgmt_YXv(:,:,i) ;

    mgmt_d_YX = mgmt_d_YXv(:,:,i) ;
    mgmt_d_glob_1 = sum(mgmt_d_YX(:)) ;
    mgmt_d_glob_2 = sum(mgmt_d2_YX(:)) ;
    
    if isnan(mgmt_d_glob_1)
        error(['NaN produced in check ' check_name ', mgmt_d_glob_1'])
    elseif isnan(mgmt_d_glob_2)
        error(['NaN produced in check ' check_name ', mgmt_d_glob_2'])
    end
    pct_error = (mgmt_d_glob_2-mgmt_d_glob_1)/mgmt_d_glob_1*100 ;
    if abs(pct_error) > conserv_tol_pct
        if notEnough(i)
            bad = -1 ;
            if do_warn && sum(sum(out_y1_mgmtTot_YXv(:,:,i)))==0
                warning('Global %s mgmt changes are not conserved to within %0.2f percent because they hit zero (check %s)', ...
                    LPJGcrops{i}, conserv_tol_pct, check_name)
            elseif do_warn && notEnough(i)
                warning('Global %s mgmt changes are not conserved to within %0.2f percent because not enough available headroom (check %s)', ...
                    LPJGcrops{i}, conserv_tol_pct, check_name)
            end
        else
            bad = 1 ;
            error('Global %s mgmt changes are not conserved to within %0.2f%%! (check %s; err = %0.2f%%)', ...
                LPJGcrops{i}, conserv_tol_pct, check_name, pct_error)
        end
    end
end


end