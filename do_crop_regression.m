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
addRequired(p,'yield_in_obs_Cyc',@(x) isnumeric(x) & ndims(x)==3) ;
addRequired(p,'yield_in_sim_Cyc',@(x) isnumeric(x) & ndims(x)==3) ;
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
addParameter(p,'slope_symbol','',@ischar) ;
addParameter(p,'marker_size',10,@isscalar) ;
addParameter(p,'separate_figs',false,@islogical) ;
parse(p,yield_in_obs_Cyc,yield_in_sim_Cyc,...
        ignore_obs_Cc, ignore_sim_Cc, ...
        weights_regr_Cyc,weights_plot_Cyc, ...
        listCrops_plot,listCrops_data,...
        plot2data_key,scatter_style,varargin{:});
pr = p.Results ;
do_wtd_reg = ~isempty(pr.weights_regr_Cyc) ;
do_wtd_plot = strcmp('scatter_style','size_weighted') ;
do_miscanthus = ~isempty(pr.miscanthus_file) ;
if do_miscanthus && (do_wtd_reg || do_wtd_plot)
    error('No way programmed to do weighting for Miscanthus data/plot!')
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
elseif Ncrops_plot==10
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

calib_factors_u = nan(Ncrops_plot,1) ;
calib_factors_w = nan(Ncrops_plot,1) ;

if ~pr.separate_figs
    figure('Color','w') ;
end

for c_plot = 1:Ncrops_plot
    if pr.separate_figs
        figure('Color','w') ;
    end
    if do_miscanthus && c_plot==Ncrops_plot
        miscanthus_data = load(pr.miscanthus_file) ;
        tmpO = miscanthus_data.yields_obs ;
        tmpS = miscanthus_data.yields_mod ;
        thisCrop_plot = 'Miscanthus' ;
    else
        thisCrop_plot = listCrops_plot{c_plot} ;
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
        thisCrop_data = listCrops_data{c_data} ;
        disp(['thisCrop_plot = ' thisCrop_plot  ' (index ' num2str(c_plot) ')'])
        disp(['thisCrop_data = ' thisCrop_data ' (index ' num2str(c_data) ')'])

% % %         tmpO = yield_in_obs_Cyc(:,:,c_plot) ;
% % %         tmpS = yield_in_sim_Cyc(:,:,c_data) ;
% % %         includeThese = true(size(ignore_obs_Cc,1),1) ;
        includeThese = ~ignore_obs_Cc(:,c_plot) & ~ignore_sim_Cc(:,c_data) ;
        disp(['Including ' num2str(length(find(includeThese))) '/' num2str(length(includeThese)) ' countries.'])
        tmpO = yield_in_obs_Cyc(includeThese,:,c_plot) ;
        if ~isequal(size(yield_in_obs_Cyc),size(yield_in_sim_Cyc))
            yield_in_sim_Cyc = permute(yield_in_sim_Cyc,[2 1 3]) ;
            if ~isequal(size(yield_in_obs_Cyc),size(yield_in_sim_Cyc))
                error('yield_in_obs_Cyc and yield_in_sim_Cyc are different sizes, even after attempted fix!')
            end
        end
        tmpS = yield_in_sim_Cyc(includeThese,:,c_data) ;
        if do_wtd_reg
%             tmpW = weights_regr_Cyc(:,:,c_plot) ;
%             tmpWp = weights_plot_Cyc(:,:,c_plot) ;
            tmpW = weights_regr_Cyc(includeThese,:,c_plot) ;
            tmpWp = weights_plot_Cyc(includeThese,:,c_plot) ;
        else
            tmpW = [] ;
            tmpWp = [] ;
        end

        if pr.avg_over_yrs
            tmpO = nanmean(tmpO,2) ;
            tmpS = nanmean(tmpS,2) ;
            if do_wtd_reg
                tmpW = nanmean(tmpW,2) ;
                tmpWp = nanmean(tmpWp,2) ;
            end
        end

        if do_wtd_reg
            bad = isnan(tmpO) | isnan(tmpS) | isnan(tmpW) | tmpO>prctile(tmpO(:),pr.max_prctile) ;
        else
%             bad = isnan(tmpO) | isnan(tmpS) | tmpO>prctile(tmpO(:),pr.max_prctile) ;
            cond1 = isnan(tmpO) ;
            cond2 = isnan(tmpS) ;
            cond3 = tmpO>prctile(tmpO(:),pr.max_prctile) ;
            bad = cond1 | cond2 | cond3 ;
        end
        if ~any(find(~bad))
%             error('No cells found!')
            if all(cond2)
                warning(['No data found for ' thisCrop_plot '! Skipping.'])
                continue
            end
        end
        if any(isinf(tmpO))
            error('any(isinf(tmpO))')
        end
        
        
        % Troubleshooting
        disp('EXCLUDING:')
        disp(['   NaN in obs and sim: ' num2str(length(find(cond1 & cond2)))])
        disp(['   NaN in obs not sim: ' num2str(length(find(cond1 & ~cond2)))])
        disp(['   NaN in sim not obs: ' num2str(length(find(~cond1 & cond2)))])
        disp(['   Low percentile:     ' num2str(length(find(cond3)))])
        
        tmpO = tmpO(~bad) ;
        tmpS = tmpS(~bad) ;
        disp(['N = ' num2str(length(tmpS))])
        if do_wtd_reg
            tmpW = tmpW(~bad) ;
            tmpWp = tmpWp(~bad) ;
        end
    end
        
    calib_factors_u(c_plot) = lscov(tmpS(:),tmpO(:)) ;
    disp(['Calibration factor, unweighted = ' num2str(calib_factors_u(c_plot))])
    if do_wtd_reg
        calib_factors_w(c_plot) = lscov(tmpS,tmpO,tmpW) ;
        disp(['Calibration factor, weighted   = ' num2str(calib_factors_w(c_plot))])
    end
    
    disp(' ')
    
    if ~pr.separate_figs
        subplot(nsubp_y,nsubp_x,c_plot)
    end
    if strcmp(scatter_style,'size_uniform')
        plot(tmpO,tmpS,'.b','MarkerSize',pr.marker_size)
    elseif strcmp(scatter_style,'size_weighted')
        scatter(tmpO,tmpS,1000*tmpWp,'.b')
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
    hold off
    % Plot regression lines
    regLine_u_x = [min(tmpO) max(tmpO)] ;
    regLine_u_y = regLine_u_x / calib_factors_u(c_plot) ;
    
    hold on
    regLine_u = plot(regLine_u_x,regLine_u_y,'-m','LineWidth',pr.reg_line_width) ;
    if do_wtd_reg
        regLine_w_y = regLine_u_x / calib_factors_w(c_plot) ;
        regLine_w = plot(regLine_u_x,regLine_w_y,'-c','LineWidth',pr.reg_line_width) ;
    end
    hold off
    % Labels
    title([pr.plot_title_prefix thisCrop_plot]) ;
    xlabel('Observed avg. yield (tDM/ha)')
    ylabel('Simulated avg. yield (tDM/ha)')
    set(gca,'FontSize',pr.fig_font_size)
    % Legend
    if do_wtd_reg
        legend([regLine_u,regLine_w],...
            {[pr.slope_symbol '_u = ' num2str(round(calib_factors_u(c_plot),3))],...
            [pr.slope_symbol '_w = ' num2str(round(calib_factors_w(c_plot),3))]},...
            'Location','north')
    else
        legend(regLine_u,...
            [pr.slope_symbol ' = ' num2str(round(calib_factors_u(c_plot),3))],...
            'Location','north')
    end
end



end