function e2p_save_outlier_info(outlier_info, outDir, which_file, out_cols)

if size(outlier_info.thresh_vt,2) == 1
    error('This assumes >1 time period')
end

T = make_out_table(outlier_info.thresh_vt, shiftdim(outlier_info.varNames), out_cols) ;
writetable(T, sprintf('%s/outliers_%s_thresh.csv', outDir, which_file))

T = make_out_table(outlier_info.count_vt, shiftdim(outlier_info.varNames), out_cols) ;
writetable(T, sprintf('%s/outliers_%s_count.csv', outDir, which_file))

T = make_out_table(outlier_info.max_vt, shiftdim(outlier_info.varNames), out_cols) ;
writetable(T, sprintf('%s/outliers_%s_max.csv', outDir, which_file))


end


function T = make_out_table(array_vt, varNames, out_cols)

T1 = table(varNames) ;
T2 = array2table(array_vt) ;
T = [T1 T2] ;
T.Properties.VariableNames = ["Crop" ; out_cols] ;

end