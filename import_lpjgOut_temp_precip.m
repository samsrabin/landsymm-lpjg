% Historical
disp('Importing historical...')
file_in_hist = sprintf('%s/seasonality.out', baselineDir) ;
tmp = lpjgu_matlab_read2geoArray(file_in_hist, 'verboseIfNoMat', false) ;
tempH_yx = transpose(squeeze(tmp.garr_xvy(:,strcmp(tmp.varNames,'temp_mean'),:))) ;
precH_yx = transpose(squeeze(tmp.garr_xvy(:,strcmp(tmp.varNames,'prec'),:))) ;
yearListH_tmp = tmp.yearList ;
list2mapH_tmp = tmp.list2map ;
clear tmp file_in_hist

% Future
for r = 1:Nruns
    file_in = sprintf('%s/seasonality.out', runDirs{r}) ;
    fprintf('Importing future run %d of %d...\n', r, Nruns) ;
    tmp = lpjgu_matlab_read2geoArray(file_in, 'verboseIfNoMat', false) ;
    if r == 1
        yearListF = tmp.yearList ;
        NyearsF = length(yearListF) ;
        lonlats = tmp.lonlats ;
        list2mapF_tmp = tmp.list2map ;
        if ~isequal(list2mapH_tmp, list2mapF_tmp)
            error('isequal(list2mapH, list2mapF)')
        end
        Ncells_tmp = length(list2mapF_tmp) ;
        tempF_yxr = nan(NyearsF, Ncells_tmp, Nruns) ;
        precF_yxr = nan(NyearsF, Ncells_tmp, Nruns) ;
    end
    tempF_yxr(:,:,r) = transpose(squeeze(tmp.garr_xvy(:,strcmp(tmp.varNames,'temp_mean'),:))) ;
    precF_yxr(:,:,r) = transpose(squeeze(tmp.garr_xvy(:,strcmp(tmp.varNames,'prec'),:))) ;
    clear tmp file_in
end
clear Ncells_tmp

% Convert deg C to K (as expected by make_big_bar_graph.m
tempH_yx = tempH_yx + 273.15 ;
tempF_yxr = tempF_yxr + 273.15 ;

disp('Finishing...')
land_area_x = land_area_YX(list2mapH_tmp) ;
land_area_weights_YX = land_area_YX ./ nansum(nansum(land_area_YX)) ;
land_area_weights_x = land_area_weights_YX(list2mapH_tmp) ;

ts_temp_bl = sum(tempH_yx .* repmat(permute(land_area_weights_x, [2 1]), [Nyears_bl 1]), 2) ;
ts_temp_yr = squeeze(sum(tempF_yxr .* repmat(permute(land_area_weights_x, [2 1]), [Nyears_fu 1 Nruns]), 2)) ;
ts_prec_bl = sum(precH_yx .* repmat(permute(land_area_weights_x, [2 1]), [Nyears_bl 1]), 2) ;
ts_prec_yr = squeeze(sum(precF_yxr .* repmat(permute(land_area_weights_x, [2 1]), [Nyears_fu 1 Nruns]), 2)) ;

% Convert mm to m3
ts_prec_bl = ts_prec_bl*1e-3 * sum(land_area_x) ;
ts_prec_yr = ts_prec_yr*1e-3 * sum(land_area_x) ;

clear tempH_yxr tempF_yxr precH_yx precF_yxr list2map*_tmp yearList*_tmp

disp('Done importing.')