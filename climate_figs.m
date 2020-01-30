%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make climate figures for a set of runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_in_hist = '/Users/Shared/PLUM/trunk_runs/LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851/seasonality.out.gz' ;

top_dir = '/Users/Shared/PLUM/trunk_runs' ;
subfiles = { ...
    'LPJGPLUM_2011-2100_harm3_SSP1_RCP45/output-2019-02-27-103914/seasonality.out.gz', ...
    'LPJGPLUM_2011-2100_harm3_SSP4_RCP60/output-2019-02-27-093259/seasonality.out.gz', ...
    'LPJGPLUM_2011-2100_harm3_SSP5_RCP85/output-2019-02-27-104120/seasonality.out.gz'} ;
plot_labels = strcat('RCP', { ...
    '4.5', '6.0', '8.5'}) ;


%% Import

% Historical
disp('Importing historical...')
tmp = lpjgu_matlab_read2geoArray(file_in_hist, 'verboseIfNoMat', false) ;
tempH_yx = tmp.garr_yxv(:,:,strcmp(tmp.varNames,'temp_mean')) ;
precH_yx = tmp.garr_yxv(:,:,strcmp(tmp.varNames,'prec')) ;
yearListH = tmp.yearList ;
list2mapH = tmp.list2map ;
clear tmp

% Future
Nruns = length(subfiles) ;
for r = 1:Nruns
    file_in = sprintf('%s/%s', top_dir, subfiles{r}) ;
    fprintf('Importing future run %d of %d...\n', r, Nruns) ;
    tmp = lpjgu_matlab_read2geoArray(file_in, 'verboseIfNoMat', false) ;
    if r == 1
        yearListF = tmp.yearList ;
        NyearsF = length(yearListF) ;
        lonlats = tmp.lonlats ;
        list2mapF = tmp.list2map ;
        if ~isequal(list2mapH, list2mapF)
            error('isequal(list2mapH, list2mapF)')
        end
        Ncells = length(list2mapF) ;
        tempF_yxr = nan(NyearsF, Ncells, Nruns) ;
        precF_yxr = nan(NyearsF, Ncells, Nruns) ;
    end
    tempF_yxr(:,:,r) = tmp.garr_yxv(:,:,strcmp(tmp.varNames,'temp_mean')) ;
    precF_yxr(:,:,r) = tmp.garr_yxv(:,:,strcmp(tmp.varNames,'prec')) ;
    clear tmp
end
disp('Done importing.')

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear tmp land_frac_YXqd land_area_YXqd
land_area_x = land_area_YX(list2mapH) ;
land_area_weights_x = land_area_x / sum(land_area_x) ;

%% Plot temperature change map

plotYearsH = 2001:2010 ;
plotYearsF = 2091:2100 ;

o.spacing = [0.05, 0.03] ;   % v, h
o.ny = 3 ;
o.nx = 1 ;
o.thisPos = [1    34   720   771] ;
o.thisVar = 'temp.' ;
o.fontSize = 14 ;
o.diff_as_pct = false ;
o.units = 'degrees C' ;
o.noTitle = false ;

figure('Position', o.thisPos, 'Color', 'w') ;
make_map(tempH_yx, tempF_yxr, ...
    yearListH, yearListF, plotYearsH, plotYearsF, ...
    plot_labels, list2mapH, o, land_area_weights_x)


%% Plot precip change map

plotYearsH = 1971:2000 ;
plotYearsF = 2071:2100 ;

o.spacing = [0.05, 0.03] ;   % v, h
o.ny = 3 ;
o.nx = 1 ;
o.thisPos = [1    34   720   771] ;
o.thisVar = 'precip.' ;
o.fontSize = 14 ;
o.diff_as_pct = true ;
o.units = 'mm yr^{-1}' ;
o.noTitle = false ;

figure('Position', o.thisPos, 'Color', 'w') ;
make_map(precH_yx, precF_yxr, ...
    yearListH, yearListF, plotYearsH, plotYearsF, ...
    plot_labels, list2mapH, o, land_area_weights_x)


%% Plot both together

o.spacing = [0.05, 0.03] ;   % v, h
o.ny = 3 ;
o.nx = 2 ;
o.thisPos = [1          33        1243         772] ;
o.fontSize = 14 ;
o.noTitle = true ;

figure('Position', o.thisPos, 'Color', 'w') ;
hs = {} ;

o.thisVar = 'temp.' ;
o.diff_as_pct = false ;
o.units = 'degrees C' ;
plotYearsH_temp = 2001:2010 ;
plotYearsF_temp = 2091:2100 ;
hs = {} ;
hs = make_map(tempH_yx, tempF_yxr, ...
    yearListH, yearListF, plotYearsH_temp, plotYearsF_temp, ...
    plot_labels, list2mapH, [1 3 5], hs, o, land_area_weights_x) ;

o.thisVar = 'precip.' ;
o.diff_as_pct = true ;
o.units = 'mm yr^{-1}' ;
plotYearsH_prec = 1971:2000 ;
plotYearsF_prec = 2071:2100 ;
hs = make_map(precH_yx, precF_yxr, ...
    yearListH, yearListF, plotYearsH_prec, plotYearsF_prec, ...
    plot_labels, list2mapH, [2 4 6], hs, o, land_area_weights_x) ;

label_mag = 1.25 ;
if o.noTitle
    % Set up invisible axes for plot labels
    a = axes;
    a.Position = [0 0 1 0.95] ;
    axis off
    % Add row labels (runs)
    for r = 1:Nruns
        plot_ind = (r-1)*o.nx + 1 ;
        t1 = text( ...
            hs{plot_ind}.Position(1), ...
            hs{plot_ind}.Position(2)+hs{plot_ind}.Position(4)*0.65, ...
            plot_labels{r}) ;
        t1.HorizontalAlignment = 'center' ;
        t1.VerticalAlignment = 'bottom' ;
        t1.Rotation = 90 ;
        t1.FontSize = o.fontSize*label_mag ;
        t1.FontWeight = 'bold' ;
    end
    % Add column labels (meteorological variable)
    var_labels = strcat('\Delta', { ...
        sprintf(' Temperature (%d-%d to %d-%d)', ...
            min(plotYearsH_temp), max(plotYearsH_temp), ...
            min(plotYearsF_temp), max(plotYearsF_temp)), ...
        sprintf(' Precipitation (%d-%d to %d-%d)', ...
            min(plotYearsH_prec), max(plotYearsH_prec), ...
            min(plotYearsF_prec), max(plotYearsF_prec))}) ;
    for v = 1:o.nx
        t1 = text( ...
            hs{v}.Position(1)+hs{v}.Position(3)*0.5, ...
            1, ...
            var_labels{v} ) ;
        t1.HorizontalAlignment = 'center' ;
        t1.VerticalAlignment = 'bottom' ;
        t1.FontSize = o.fontSize*label_mag ;
        t1.FontWeight = 'bold' ;
    end
end
stop
export_fig('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/writing/LaTeX_v6.0.1/figures_supp_results/climate_maps.png', '-r300') ;


%% FUNCTIONS

function hs = make_map(dataH_yx, dataF_yxr, ...
    yearListH, yearListF, plotYearsH, plotYearsF, ...
    plot_labels, list2map, plot_inds, hs, o, land_area_weights_x)

is_yearsH = yearListH>=plotYearsH(1) & yearListH<=plotYearsH(end) ;
is_yearsF = yearListF>=plotYearsF(1) & yearListF<=plotYearsF(end) ;

Nruns = size(dataF_yxr, 3) ;
hist_xr = repmat(shiftdim(mean(dataH_yx(is_yearsH,:), 1)), [1 Nruns]) ;
diff_xr = squeeze(mean(dataF_yxr(is_yearsF,:,:),1)) - hist_xr ;
% mean_diff_r = sum(diff_xr .* repmat(land_area_weights_x, [1 Nruns]), 1)
if o.diff_as_pct
    diff_xr = 100 * diff_xr./hist_xr ;
end

if strcmp(o.thisVar, 'temp.')
    if min(min(min(diff_xr)))<0
        error('Some temperature change <0!')
    else
        new_clim = minmax_ssr(diff_xr) ;
    end
elseif strcmp(o.thisVar, 'precip.')
    if o.diff_as_pct
        new_clim = 100*[-1 1] ;
    else
        new_clim = max(max(max(abs(diff_xr))))*[-1 1] ;
    end
else
    error('thisVar (%s) not recognized', o.thisVar) ;
end

for r = 1:Nruns
    this_plot_ind = plot_inds(r) ;
    hold on
    hsp = subplot_tight(o.ny, o.nx, this_plot_ind, o.spacing) ;
    hs{this_plot_ind} = hsp ;
    hold off
    map_YX = nan(360,720) ;
    map_YX(list2map) = diff_xr(:,r) ;
    map_YX = map_YX(60:end,:) ;
    pcolor(hsp, map_YX) ; shading flat; axis tight off
%     pbaspect([2 1 1])
    axis equal
    caxis(new_clim) ;
    if strcmp(o.thisVar, 'temp.')
        colormap(gca,brewermap(64,'YlOrRd'))
    elseif strcmp(o.thisVar, 'precip.')
        colormap(gca,brewermap(64,'BrBG_ssr'))
    else
        error('thisVar (%s) not recognized', o.thisVar) ;
    end
    hcb = colorbar(hsp) ;

    if ~o.noTitle
        title(sprintf('%s mean %s change: %d-%d to %d-%d', ...
            plot_labels{r}, o.thisVar, min(plotYearsH), max(plotYearsH), min(plotYearsF), max(plotYearsF))) ;
    end
    
    set(hsp, 'FontSize', o.fontSize) ;
    hcb.TickLabels(hcb.Ticks>=0) = strcat('+', hcb.TickLabels(hcb.Ticks>=0)) ;
    if o.diff_as_pct
        hcb.TickLabels = strcat(hcb.TickLabels, '%') ;
        if strcmp(o.thisVar, 'precip.')
            hcb.TickLabels{end} = ['>=' hcb.TickLabels{end}] ;
        end
    else
        ylabel(hcb, o.units) ;
    end
    
end

end









