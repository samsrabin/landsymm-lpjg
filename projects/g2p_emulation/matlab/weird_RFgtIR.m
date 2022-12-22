
ggcm = 'EPIC-TAMU' ;
% thisN = 10 ; crop = 'soy' ; crop_short = 'soy' ;
thisN = 60 ; crop = 'spring_wheat' ; crop_short = 'swh' ;

adapt = 1 ;


%% Setup

% Timesteps
future_y1 = 2005 ;
baseline_y1 = 2001 ;
baseline_yN = 2010 ;
future_ts = 10 ; % Number of years in future time step
future_yN_emu = 2084 ;
Nyears_ts = future_ts ;
tsN_y1 = floor(future_yN_emu/Nyears_ts)*Nyears_ts ;
ts1_list = future_y1:Nyears_ts:tsN_y1 ;
tsN_list = ts1_list + Nyears_ts - 1 ;
Ntpers = length(ts1_list) ;

topDir = '/Volumes/Reacher/GGCMI/AgMIP.output/CMIP_emulated/yields/CMIP6' ;
outDir = sprintf('%s/troubleshooting/irr_reduces_yield/%s_%s_A%d_N%d', ...
    topDir, ggcm, crop_short, adapt, thisN) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

yearList_fu = 2011:2084 ;
ts1_list(1) = 2011 ;

% Add space for two baselines
ts1_list = [1980 1980 ts1_list] ;
tsN_list = [2010 2010 tsN_list] ;
Ntpers = length(ts1_list) ;


%% Import

irrList = {'rf', 'ir'} ;
Wlist = {'0', 'inf'} ;

maps_YXit = nan(360, 720, 2, Ntpers) ;
filenames_it = cell(2, Ntpers) ;
for ii = 1:2
    thisIrr = irrList{ii} ;
    thisW = Wlist{ii} ;
    varName = sprintf('yield_%s_%s', thisIrr, crop_short) ;
    fprintf('Importing %s...\n', thisIrr)
    
    % Baseline simulated
    t = 1 ;
    filepattern = ['/Volumes/Reacher/GGCMI/AgMIP.output/%s/phase2/%s/A%d/yield/' ...
        '%s_agmerra_fullharm_yield_%s_global_annual_1980_2010_C360_T0_W%s_N%d_A%d.nc4' ] ;
    if strcmp(ggcm, 'EPIC-TAMU') && (strcmp(crop, 'soy') || strcmp(crop, 'spring_wheat'))
        adapt_ph2 = 0 ;
    else
        error('Baseline simulation file not specified for %s %s', ...
            ggcm, crop)
    end
    filename = sprintf(filepattern, ...
        ggcm, crop, adapt_ph2, lower(ggcm), crop_short, thisW, thisN, adapt_ph2) ;
    filenames_it{ii,t} = filename ;
    varName_ph2 = sprintf('yield_%s', crop_short) ;
    if exist(filename, 'file')
        tmp_YXy = flip(permute(ncread(filename, varName_ph2), [2 1 3]), 1) ;
        maps_YXit(:,:,ii,t) = mean(tmp_YXy(:,:,1:end-1), 3) ;
    else
        warning('%s not found', filename)
    end
    
    % Baseline emulated
    t = 2 ;
    filename = sprintf(['%s/A%d_N%d/ssp126/%s/' ...
        'cmip6_ssp126_UKESM1-0-LL_r1i1p1_%s_%s_A%d_N%d_' ...
        'emulated_yield_baseline_1980_2010_average_v2.5.nc4'], ...
        topDir, adapt, thisN, ggcm, crop, ggcm, adapt, thisN) ;
    filenames_it{ii,t} = filename ;
    maps_YXit(:,:,ii,t) = flip(permute(ncread(filename, varName), [2 1 3]), 1) ;
    
    % Future
    filename = sprintf(['%s/A1_N%d/ssp126/%s/' ...
        'cmip6_ssp126_UKESM1-0-LL_r1i1p1_%s_%s_A%d_N%d_' ...
        'emulated_yield_movingwindow_2011_2084_v2.5.nc4'], ...
        topDir, thisN, ggcm, crop, ggcm, adapt, thisN) ;
    tmp_YXy = flip(permute(ncread(filename, varName), [2 1 3]), 1) ;
    for t = 3:Ntpers
        filenames_it{ii,3} = filename ;
        yearList_thisTper = ts1_list(t):tsN_list(t) ;
        [~,I] = intersect(yearList_thisTper, yearList_fu) ;
        maps_YXit(:,:,ii,t) = mean(tmp_YXy(:,:,I), 3) ;
    end
end
disp('Done importing.')


%% Plot

%%%%%%%%%%%%%%%%%%%%%%%%
thisPos = figurePos ;
fontSize = 14 ;
yRange = 65:360 ;
cbarLoc = 'SouthOutside' ;
cmap_negpos = 'PiYG' ;
spacing = [0.085 0.05] ; % v h
%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:Ntpers
    
    figure('Color', 'w', 'Position', thisPos)
    
    map1_YX = maps_YXit(:,:,1,t) ;
    map2_YX = maps_YXit(:,:,2,t) ;
    
    nanmask_YX = map1_YX==0 & map2_YX==0 ;
    map1_YX(nanmask_YX) = NaN ;
    map2_YX(nanmask_YX) = NaN ;
    
    skip1 = all(all(isnan(map1_YX))) ;
    skip2 = all(all(isnan(map2_YX))) ;
    
    new_caxis = [0 max(max(max(cat(3, map1_YX, map2_YX))))] ;
    
    subplot_tight(2, 2, 1, spacing)
    if skip1
        textstr = sprintf('%s\nnot found', ...
            get_basename(filenames_int{ii, n, t})) ;
        textstr = strrep(textstr, '_', '\_') ;
        text(0.5, 0.5, textstr, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle')
        axis off
    else
        ha1 = pcolor(map1_YX(yRange,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colormap(gca, 'jet'); hcb = colorbar('Location', cbarLoc) ;
        ylabel(hcb, 'tons/ha')
        set(gca, 'FontSize', fontSize)
    end
    title(irrList{1}, ...
        'FontSize', fontSize)
    
    subplot_tight(2, 2, 2, spacing)
    if skip2
        textstr = sprintf('%s\nnot found', ...
            get_basename(filenames_int{ii, n, t})) ;
        textstr = strrep(textstr, '_', '\_') ;
        text(0.5, 0.5, textstr, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle')
        axis off
    else
        ha2 = pcolor(map2_YX(yRange,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colormap(gca, 'jet'); hcb = colorbar('Location', cbarLoc) ;
        ylabel(hcb, 'tons/ha')
        set(gca, 'FontSize', fontSize)
        
        mapD_YX = map2_YX - map1_YX ;
        mapD_max = max(max(abs(mapD_YX))) ;
        if any(any(mapD_YX < 0))
            cmap = cmap_negpos ;
            new_caxis = [-1 1] * mapD_max ;
        else
            cmap = 'Greens' ;
            new_caxis = [0 1] * mapD_max ;
        end
        new_cmap = brewermap(256, cmap) ;
    end
    title(irrList{2}, ...
        'FontSize', fontSize)
    
    subplot_tight(2, 1, 2, spacing)
    if skip1 || skip2
        axis off
    else
        ha3 = pcolor(mapD_YX(yRange,:)); shading flat; axis equal tight off
        caxis(new_caxis)
        colormap(gca, new_cmap); hcb = colorbar('Location', cbarLoc) ;
        ylabel(hcb, 'tons/ha')
        title(strrep(sprintf('%s - %s', ...
            irrList{2}, irrList{1}), ...
            '_', '\_'))
        set(gca, 'FontSize', fontSize)
    end
    
    y1 = ts1_list(t) ;
    yN = tsN_list(t) ;
    if t == 1
        hsgt = sgtitle(strrep(sprintf('%s %s N%d A%d (sim): %d-%d', ...
            ggcm, crop, thisN, adapt_ph2, y1, yN-1), ...
            '_', '\_')) ;
        filename = sprintf('%s/%s_%s_A%d_rf-ir_%d-%d_sim.png', ...
            outDir, ggcm, crop_short, adapt, y1, yN-1) ;
    else
        hsgt = sgtitle(strrep(sprintf('%s %s N%d A%d (emu): %d-%d', ...
            ggcm, crop, thisN, adapt, y1, yN), ...
            '_', '\_')) ;
        filename = sprintf('%s/%s_%s_A%d_rf-ir_%d-%d_emu.png', ...
            outDir, ggcm, crop_short, adapt, y1, yN) ;
    end
    
    set(hsgt, 'FontSize', fontSize+4, 'FontWeight', 'bold')
    
    export_fig(filename, '-r200') ;
    close
    
end


%% FUNCTIONS

function basename = get_basename(filename)

fileparts = strsplit(filename, '/') ;
basename = fileparts{end} ;

end





