function bad = PLUMharm_checkCons_area(...
    out_y0_area_YXv, out_y1_area_YXv, ...
    in_y0_area_YXv, in_y1_area_YXv, ...
    total_unmet_area_YXv, LUnames, ...
    conserv_tol_pct, conserv_tol_area, ...
    check_name)

area_d_YXv = in_y1_area_YXv - in_y0_area_YXv ;

bad = 0 ;
for i = 1:length(LUnames)

    area_d2_YX = (out_y1_area_YXv(:,:,i) + total_unmet_area_YXv(:,:,i)) - out_y0_area_YXv(:,:,i)  ;

    area_d_YX = area_d_YXv(:,:,i) ;
    area_d_glob_1 = sum(area_d_YX(:)) ;
    area_d_glob_2 = sum(area_d2_YX(:)) ;
    
    if isnan(area_d_glob_1)
        error(['NaN produced in check ' check_name ', area_d_glob_1'])
    elseif isnan(area_d_glob_2)
        error(['NaN produced in check ' check_name ', area_d_glob_2'])
    end
    
    err = area_d_glob_2-area_d_glob_1 ;
    err_pct = (err)/area_d_glob_1 * 100 ;
    if area_d_glob_1 == 0
        % Could use actual conserv_tol_area (m2), but let's be strict here.
        conserv_tol_area = 1 ;
        is_bad = abs(err) > conserv_tol_area ;
    else
        is_bad = abs(err_pct) > conserv_tol_pct ;
    end
    if is_bad && abs(err_pct) > conserv_tol_area
        if sum(sum(out_y1_area_YXv(:,:,i)))==0
            bad = -1 ;
            warning('Global %s area changes are not conserved to within %g percent because they hit zero (check %s)\n', ...
                LUnames{i}, conserv_tol_pct, check_name)
        else
            bad = 1 ;
            if area_d_glob_1 == 0
                warning([...
                    'Global %s area changes are not conserved to within %g m2!\n'...
                    '(check %s; err = %0.2f)\n'], ...
                    LUnames{i}, conserv_tol_area, check_name, err)

            else
                warning([...
                    'Global %s area changes are not conserved to within %0.2f percent!\n'...
                    '(check %s; err = %0.2f)\n'], ...
                    LUnames{i}, conserv_tol_pct, check_name, err_pct)
            end
        end
        x=1;
    end
end


end