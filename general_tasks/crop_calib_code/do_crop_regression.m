function [calib_factors_u,calib_factors_w] = ...
            do_crop_regression(yield_in_obs_Cyc,... %yield_fa2_4cal_Cyc
                               yield_in_sim_Cyc,... %yield_lpj_4cal_Cyc
                               ignore_obs_Cc,... % ignore_fa2_Cc
                               ignore_sim_Cc,... % ignore_lpj_Cc
                               weights_regr_Cyc,... %weights_fa2_4cal_Cyc
                               weights_plot_Cyc,... %weights4pts_Cyc
                               listCrops_plot,... % listCrops_fa2o
                               listCrops_data,... % listCrops_4cal
                               plot2data_key,... % FA2_to_PLUM_key
                               scatter_style,...
                               varargin)

% Set up input arguments
p = inputParser ;
addRequired(p,'yield_in_obs_Cyc',@(x) isnumeric(x) & (ismatrix(x)|ndims(x)==3)) ;
addRequired(p,'yield_in_sim_Cyc',@(x) isnumeric(x) & (ismatrix(x)|ndims(x)==3)) ;
addRequired(p,'ignore_obs_Cc',@(x) islogical(x) & ismatrix(x)) ;
addRequired(p,'ignore_sim_Cc',@(x) islogical(x) & ismatrix(x)) ;
addRequired(p,'weights_regr_Cyc',@(x) (isnumeric(x) & ndims(x)==3) | isempty(x)) ;
addRequired(p,'weights_plot_Cyc',@(x) (isnumeric(x) & ndims(x)==3) | isempty(x)) ;
addRequired(p,'listCrops_plot',@iscellstr) ;
addRequired(p,'listCrops_data',@iscellstr) ;
addRequired(p,'plot2data_key',@iscell) ;
addRequired(p,'scatter_style',@ischar) ;
addParameter(p,'plot_title_prefix','',@ischar) ;
addParameter(p,'max_prctile',100,@isnumeric) ;
addParameter(p,'avg_over_yrs',false,@islogical) ;
addParameter(p,'reg_line_width',3,@isnumeric) ;
addParameter(p,'fig_font_size',14,@isnumeric) ;
addParameter(p,'miscanthus_file','',@ischar) ;
addParameter(p,'miscanthus_x_os',[],@isnumeric) ;
addParameter(p,'slope_symbol','',@ischar) ;
addParameter(p,'marker_size',10,@isscalar) ;
addParameter(p,'separate_figs',false,@islogical) ;
addParameter(p,'regression_type','',@ischar) ;
addParameter(p,'dir_outFigs','',@ischar) ;
addParameter(p,'listCountries_map_present',{},@iscellstr) ;
addParameter(p,'yearList',[],@isnumeric) ;
isOK_outlierThresh = @(x) numel(x)==1 & x>0 ;
addParameter(p,'outlier_thresh',Inf,isOK_outlierThresh) ;
parse(p,yield_in_obs_Cyc,yield_in_sim_Cyc,...
        ignore_obs_Cc, ignore_sim_Cc, ...
        weights_regr_Cyc,weights_plot_Cyc, ...
        listCrops_plot,listCrops_data,...
        plot2data_key,scatter_style,varargin{:});
pr = p.Results ;
do_wtd_reg = ~isempty(pr.weights_regr_Cyc) ;
do_wtd_plot = strcmp('scatter_style','size_weighted') ;

% Set up for output saving
dir_outData = '' ;
if ~isempty(pr.dir_outFigs)
    if any(ignore_obs_Cc(:)) || any(ignore_sim_Cc(:))
        error('How do you want to save tmpO and tmpS where you''re ignoring some from the get-go?')
    end

    dir_outData = fullfile(pr.dir_outFigs, 'scatter_data_CSVs') ;
    if ~exist(dir_outData, 'dir')
        mkdir(dir_outData)
    end

    if isempty(pr.listCountries_map_present)
        error('If specifying dir_outData, you must also specify listCountries_map_present')
    elseif length(pr.listCountries_map_present) ~= size(yield_in_obs_Cyc, 1)
        error('Size mismatch between listCountries_map_present and yield_in_obs_Cyc')
    end
    countries_Cy = repmat(pr.listCountries_map_present(:), [1 size(yield_in_obs_Cyc, 2)]) ;

    if isempty(pr.yearList)
        error('If specifying dir_outData, you must also specify yearList')
    elseif length(pr.yearList) ~= size(yield_in_obs_Cyc, 2)
        error('Size mismatch between yearList and yield_in_obs_Cyc')
    end
    years_Cy = repmat(shiftdim(pr.yearList)', [length(pr.listCountries_map_present) 1]) ;
end

% Handle Miscanthus inputs
do_miscanthus = ~isempty(pr.miscanthus_file) || ~isempty(pr.miscanthus_x_os) ;
if do_miscanthus && (do_wtd_reg || do_wtd_plot)
    warning('No way programmed to do weighting for Miscanthus data/plot; will do unweighted only.')
end
if ~isempty(pr.miscanthus_file) && ~isempty(pr.miscanthus_x_os)
    error('Specify only one of miscanthus_file and miscanthus_x_os')
end
if ~isempty(pr.miscanthus_x_os) && size(pr.miscanthus_x_os, 2) ~= 2
    error('Expected 2 columns in miscanthus_x_os; got %d', size(pr.miscanthus_x_os, 2))
end

Ncrops_plot = length(listCrops_plot) ;
if do_miscanthus
    Ncrops_plot = Ncrops_plot + 1 ;
end
Ncrops_data = length(listCrops_data) ;

if Ncrops_plot==5 || Ncrops_plot==6
    nsubp_x = 3 ;
    nsubp_y = 2 ;
elseif Ncrops_plot==7 || Ncrops_plot==8
    nsubp_x = 4 ;
    nsubp_y = 2 ;
elseif Ncrops_plot==9
    nsubp_x = 3 ;
    nsubp_y = 3 ;
elseif Ncrops_plot>=10 && Ncrops_plot<=12
    nsubp_x = 4 ;
    nsubp_y = 3 ;
elseif Ncrops_plot==2
    nsubp_x = 2 ;
    nsubp_y = 1 ;
elseif Ncrops_plot==1
    nsubp_x = 1 ;
    nsubp_y = 1 ;
elseif Ncrops_plot==3 || Ncrops_plot==4
    nsubp_x = 2 ;
    nsubp_y = 2 ;
else
    error(['Ncrops_plot = ' num2str(Ncrops_plot) '. Set up subplot layout.'])
end

if strcmp(pr.regression_type, 'slope-only')
    calib_factors_u = nan(Ncrops_plot,1) ;
    calib_factors_w = nan(Ncrops_plot,1) ;
elseif strcmp(pr.regression_type, 'full linear but 0->0')
    % Columns: Y-intercept, slope
    calib_factors_u = nan(Ncrops_plot,2) ;
    calib_factors_w = nan(Ncrops_plot,2) ;
else
    error('regression_type ''%s'' not recognized', pr.regression_type)
end


if ~pr.separate_figs
    figure('Color', 'w', 'Position', get(0,'ScreenSize')) ;
end

for c_plot = 1:Ncrops_plot
    if pr.separate_figs
        figure('Color','w') ;
    end
    is_miscanthus = do_miscanthus && c_plot==Ncrops_plot ;
    if is_miscanthus
        if ~isempty(pr.miscanthus_file)
            miscanthus_data = load(pr.miscanthus_file) ;
            tmpO = miscanthus_data.yields_obs ;
            tmpS = miscanthus_data.yields_mod ;
        elseif ~isempty(pr.miscanthus_x_os)
            tmpO = pr.miscanthus_x_os(:,1) ;
            tmpS = pr.miscanthus_x_os(:,2) ;
        else
            error('Not using miscanthus_file OR miscanthus_x_os??')
        end
        thisCrop_plot = 'Miscanthus' ;
        disp('thisCrop = Miscanthus')
    else
        thisCrop_plot = listCrops_plot{c_plot} ;
        if isempty(plot2data_key)
            c_data = c_plot ;
        else
            thisCrop_found = false ;
            for c_data = 1:Ncrops_data
                if any(find(strcmp(plot2data_key{c_data},thisCrop_plot)))
                    thisCrop_found = true ;
                    break
                end
            end
            if ~thisCrop_found
                error([thisCrop_plot ' not found in plot2data_key!'])
            end
        end
        thisCrop_data = listCrops_data{c_data} ;
        disp(['thisCrop_plot = ' thisCrop_plot  ' (index ' num2str(c_plot) ')'])
        disp(['thisCrop_data = ' thisCrop_data ' (index ' num2str(c_data) ')'])

% % %         tmpO = yield_in_obs_Cyc(:,:,c_plot) ;
% % %         tmpS = yield_in_sim_Cyc(:,:,c_data) ;
% % %         includeThese = true(size(ignore_obs_Cc,1),1) ;
        ok_obs = ~ignore_obs_Cc(:,c_plot) ;
        if ~any(ok_obs)
            warning('%s: None included (obs)!', thisCrop_plot)
        end
        ok_sim = ~ignore_sim_Cc(:,c_data) ;
        if ~any(ok_sim)
            warning('%s: None included (sim)!', thisCrop_plot)
        end
        includeThese = ok_obs & ok_sim ;
        if ~any(includeThese)
            warning('%s: None included (combination)!', thisCrop_plot)
        end
        disp(['Including ' num2str(length(find(includeThese))) '/' num2str(length(includeThese)) ' countries.'])
        tmpO = yield_in_obs_Cyc(includeThese,:,c_plot) ;
        if ~isequal(size(yield_in_obs_Cyc),size(yield_in_sim_Cyc))
            yield_in_sim_Cyc = permute(yield_in_sim_Cyc,[2 1 3]) ;
            if ~isequal(size(yield_in_obs_Cyc),size(yield_in_sim_Cyc))
                error('yield_in_obs_Cyc and yield_in_sim_Cyc are different sizes, even after attempted fix!')
            end
        end
        tmpS = yield_in_sim_Cyc(includeThese,:,c_data) ;
    end

    if do_wtd_reg && ~is_miscanthus
%             tmpW = weights_regr_Cyc(:,:,c_plot) ;
%             tmpWp = weights_plot_Cyc(:,:,c_plot) ;
        tmpW = weights_regr_Cyc(includeThese,:,c_plot) ;
        tmpWp = weights_plot_Cyc(includeThese,:,c_plot) ;
    else
        tmpW = [] ;
        tmpWp = [] ;
    end

    if pr.avg_over_yrs && ~is_miscanthus
        tmpO = nanmean(tmpO,2) ;
        tmpS = nanmean(tmpS,2) ;
        if do_wtd_reg
            tmpW = nanmean(tmpW,2) ;
            tmpWp = nanmean(tmpWp,2) ;
        end
    end

    % Trim where NaN
    cond1 = isnan(tmpO) ;
    cond2 = isnan(tmpS) ;
    if all(cond2)
        warning(['No data found for ' thisCrop_plot '! Skipping.'])
        continue
    end
    bad = cond1 | cond2 ;
    tmpO = tmpO(~bad) ;
    tmpS = tmpS(~bad) ;

    % Save data to CSV
    if ~isempty(dir_outData)
        country = countries_Cy(~bad) ;
        year = years_Cy(~bad) ;
        obs = tmpO ;
        sim = tmpS ;
        T = table(country, year, obs, sim) ;
        filename = fullfile(dir_outData, [thisCrop_plot '.csv']) ;
        writetable(T, filename)
    end

    % Trim extremes
    if do_wtd_reg && ~is_miscanthus
        tmpW = tmpW(~bad) ;
        isoutlier_tmpO = isoutlier(tmpO(:), 'quartiles',...
            'ThresholdFactor',pr.outlier_thresh) ;
        isoutlier_tmpO = reshape(isoutlier_tmpO, size(tmpO)) ;
        bad = isnan(tmpW) | tmpO>prctile(tmpO(:),pr.max_prctile) | isoutlier_tmpO ;
    else
%             bad = isnan(tmpO) | isnan(tmpS) | tmpO>prctile(tmpO(:),pr.max_prctile) ;
        cond3 = tmpO>prctile(tmpO(:),pr.max_prctile) ;
        outliersO = get_outliers(tmpO, false(size(tmpO)), false(size(tmpO)), cond3, pr) ;
        outliersS = get_outliers(tmpS, false(size(tmpO)), false(size(tmpO)), cond3, pr) ;
        bad = cond3 | outliersO.cond4 | outliersS.cond4 ;
    end
    if any(isinf(tmpO))
        error('any(isinf(tmpO))')
    end
    
    
    % Troubleshooting
    if ~do_wtd_reg || is_miscanthus
        disp('EXCLUDING:')
        if ~is_miscanthus
            disp(['   NaN in obs and sim: ' num2str(length(find(cond1 & cond2)))])
            disp(['   NaN in obs not sim: ' num2str(length(find(cond1 & ~cond2)))])
            disp(['   NaN in sim not obs: ' num2str(length(find(~cond1 & cond2)))])
        end
        disp(['   Low percentile (obs). :     ' num2str(length(find(cond3)))])
        if ~isinf(pr.outlier_thresh)
            fprintf('   Outliers (ThresholdFactor: %.9g):\n', pr.outlier_thresh)
            fprintf('      Obs: %d\n', length(find(outliersO.cond4)))
            if any(outliersO.cond4(:))
                print_outlier_info(outliersO)
            end
            fprintf('      Sim: %d\n', length(find(outliersS.cond4)))
            if any(outliersS.cond4(:))
                print_outlier_info(outliersS)
            end
        end
    end
    
    tmpO = tmpO(~bad) ;
    tmpS = tmpS(~bad) ;
    disp(['N = ' num2str(length(tmpS))])
    if do_wtd_reg && ~is_miscanthus
        tmpW = tmpW(~bad) ;
        tmpWp = tmpWp(~bad) ;
    end
    
%     TF = isoutlier(tmpO(:),'quartiles','ThresholdFactor',4) ;
%     Noutliers = length(find(TF)) ;
%     fprintf('# outliers: %d\n', Noutliers) ;
%     if Noutliers>0
%         ajbrejrbe =1 ;
%     end

    
    % Get calibration factors and lines
    regLine_u_x = [min(tmpO) max(tmpO)] ;
    if strcmp(pr.regression_type, 'slope-only')
        [calib_factors_u, calib_factors_w, ...
            regLine_u_y, regLine_w_y, legend_arg2] = doReg_slopeOnly( ...
            calib_factors_u, calib_factors_w, c_plot, tmpS, tmpO, tmpW, ...
            pr, do_wtd_reg && ~is_miscanthus, regLine_u_x) ;
    elseif strcmp(pr.regression_type, 'full linear but 0->0')
        [calib_factors_u, calib_factors_w, ...
            regLine_u_y, regLine_w_y, legend_arg2] = doReg_linBut00( ...
            calib_factors_u, calib_factors_w, c_plot, tmpS, tmpO, tmpW, ...
            do_wtd_reg && ~is_miscanthus, regLine_u_x) ;
    else
        error('regression_type ''%s'' not recognized', pr.regression_type)
    end
     
    % Make plot
    if ~pr.separate_figs
        subplot(nsubp_y,nsubp_x,c_plot)
    end
    if strcmp(scatter_style,'size_uniform')
        plot(tmpO,tmpS,'.b','MarkerSize',pr.marker_size)
    elseif strcmp(scatter_style,'size_weighted')
        thisColor = 'b' ;
        hscatter = scatter(tmpO,tmpS,1000*tmpWp, 'filled', ...
            'MarkerFaceColor', thisColor, ...
            'MarkerEdgeColor', thisColor) ;
        alpha(hscatter, 0.1)
    else
        error(['scatter_style "' scatter_style '" not recognized.'])
    end
    % Make identical X and Y limits
    xlims = get(gca,'XLim') ;
    ylims = get(gca,'YLim') ;
    newlims = [min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)]) ] ;
    set(gca,'XLim',newlims,'YLim',newlims) ;
    % Plot 1:1 line
    hold on
    plot([min(newlims) max(newlims)],[min(newlims) max(newlims)],'--k')
    % Plot regression lines
    regLine_u = plot(regLine_u_x,regLine_u_y,'-m','LineWidth',pr.reg_line_width) ;
    if do_wtd_reg && ~is_miscanthus
        regLine_w = plot(regLine_u_x,regLine_w_y,'-c','LineWidth',pr.reg_line_width) ;
    end
    hold off
    % Labels
    ht = title([pr.plot_title_prefix thisCrop_plot]) ;
    if strcmpi(thisCrop_plot,'Starchy Roots')
        xlabel('Observed avg. yield (t/ha)') % http://www.fao.org/economic/ess/ess-standards/commodity/comm-chapters/en/
    else
        xlabel('Observed avg. yield (tDM/ha)')
    end
    ylabel('Simulated avg. yield (tDM/ha)')
    set(gca,'FontSize',pr.fig_font_size)
    
    % Legend
    if do_wtd_reg && ~is_miscanthus
        hl = legend([regLine_u,regLine_w], legend_arg2, 'Location','north') ;
    else
        hl = legend(regLine_u, legend_arg2, 'Location','north') ;
    end
    
    if strcmp(pr.regression_type, 'full linear but 0->0') && calib_factors_u(c_plot, 2)<0
        warning('Negative slope!')
        hl.TextColor = 'r' ;
        ht.Color = 'r' ;
    end
    disp(' ')
            
end


end


function [calib_factors_u, calib_factors_w, ...
    regLine_u_y, regLine_w_y, legend_arg2] = doReg_slopeOnly( ...
    calib_factors_u, calib_factors_w, c_plot, tmpS, tmpO, tmpW, ...
    pr, do_wtd_reg, regLine_u_x)

% Get calibration factor(s)
calib_factors_u(c_plot) = lscov(tmpS(:),tmpO(:)) ;
disp(['Calibration factor, unweighted = ' num2str(calib_factors_u(c_plot))])
if do_wtd_reg
    calib_factors_w(c_plot) = lscov(tmpS,tmpO,tmpW) ;
    disp(['Calibration factor, weighted   = ' num2str(calib_factors_w(c_plot))])
end

if calib_factors_u(c_plot) == 0
    error('Calibration factor (unweighted) is 0!')
elseif calib_factors_w(c_plot) == 0
    error('Calibration factor (weighted) is 0!')
end

% Get line(s) and legend text
regLine_u_y = regLine_u_x / calib_factors_u(c_plot) ;

if do_wtd_reg
    regLine_w_y = regLine_u_x / calib_factors_w(c_plot) ;
    legend_arg2 = ...
        {[pr.slope_symbol '_u = ' num2str(round(calib_factors_u(c_plot),3))],...
        [pr.slope_symbol '_w = ' num2str(round(calib_factors_w(c_plot),3))]} ;
else
    regLine_w_y = [] ;
    legend_arg2 = [pr.slope_symbol ' = ' num2str(round(calib_factors_u(c_plot),3))] ;
end


end


function [calib_factors_u, calib_factors_w, ...
    regLine_u_y, regLine_w_y, legend_arg2] = doReg_linBut00( ...
    calib_factors_u, calib_factors_w, c_plot, tmpS, tmpO, tmpW, ...
    do_wtd_reg, regLine_u_x)

% Get calibration factor(s)
y = tmpO(tmpS>0) ;
if ~isempty(tmpW)
    w = tmpW(tmpS>0) ;
end
x = tmpS(tmpS>0) ;
X = [ones(length(x),1) x];
result = X \ y ;
bu = result(1) ;
mu = result(2) ;
calib_factors_u(c_plot, 1) = bu ; % Y-intercept
calib_factors_u(c_plot, 2) = mu ; % Slope
legend_arg2 = sprintf('y = %0.3g + %0.3gx', bu, mu) ;
disp(legend_arg2) ;
if do_wtd_reg
    result = ([w w].*X) \ (w .* y) ;
    bw = result(1) ;
    mw = result(2) ;
    calib_factors_w(c_plot, 1) = bw ; % Y-intercept
    calib_factors_w(c_plot, 2) = mw ; % Slope
    legend_arg2 = {legend_arg2} ;
    legend_arg2{2} = sprintf('Wtd: y = %0.3g + %0.3gx', bu, mu) ;
    disp(legend_arg2{2}) ;
end

% Get line(s)
regLine_u_y = (regLine_u_x - bu) / mu ;
regLine_w_y = [] ;
if do_wtd_reg
    regLine_w_y = (regLine_w_x - bw) / mw ;
end


end



function S = get_outliers(tmpX, cond1, cond2, cond3, pr)

S.tmp = tmpX(:) ;
S.tmp(cond1 | cond2 | cond3) = NaN ;
[S.TF, S.LO, S.UP] = isoutlier(S.tmp, ...
    'quartiles', ...
    'ThresholdFactor', pr.outlier_thresh) ;
S.cond4 = reshape(S.TF, size(tmpX)) ;


end


function print_outlier_info(S)

if any(S.tmp<S.LO)
    Nlo = length(find(S.tmp<S.LO)) ;
    loMin = min(S.tmp(S.tmp<S.LO)) ;
    if Nlo > 1
        loMax = max(S.tmp(S.tmp<S.LO)) ;
        fprintf('        Below %0.2f: %d (%0.2f to %0.2f)\n', S.LO, Nlo, loMin, loMax) ;
    else
        fprintf('        Below %0.2f: %d (%0.2f)\n', S.LO, Nlo, loMin) ;
    end
end
if any(S.tmp>S.UP)
    Nup = length(find(S.tmp>S.UP)) ;
    upMin = min(S.tmp(S.tmp>S.UP)) ;
    if Nup > 1
        upMax = max(S.tmp(S.tmp>S.UP)) ;
        fprintf('        Above %0.2f: %d (%0.2f to %0.2f)\n', S.UP, Nup, upMin, upMax) ;
    else
        fprintf('        Above %0.2f: %d (%0.2f)\n', S.UP, Nup, upMin) ;
    end
end

end

