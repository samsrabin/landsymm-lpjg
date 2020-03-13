function e2p_save_out_figs(data_fu_lpj, data_fu_emu, data_fu_out, ...
    ggcm, getbasename, getbasenamei, getN, outDir_figs, which_file)

this_colormap = 'parula' ;
% this_colormap = 'jet' ;

yield_fu_lpj_max_xv = max(data_fu_lpj.garr_xvt,[],3) ;
yield_fu_emu_max_xv = max(data_fu_emu.garr_xvt,[],3) ;
yield_fu_out_max_xv = max(data_fu_out.garr_xvt,[],3) ;

% Apply calibration factors used in ES paper
cropList_lpj = unique(getbasename(data_fu_lpj.varNames)) ;
Ncrops = length(cropList_lpj) ;
if strcmp(which_file, 'yield')
    cf_lpj = nan(Ncrops,1) ;
    cf_lpj(strcmp(cropList_lpj, 'CerealsC3')) = 1.056 ;
    cf_lpj(strcmp(cropList_lpj, 'CerealsC4')) = 0.738 ;
    cf_lpj(strcmp(cropList_lpj, 'Rice')) = 1.052 ;
    cf_lpj(strcmp(cropList_lpj, 'Oilcrops')) = 0.687 ;
    cf_lpj(strcmp(cropList_lpj, 'Pulses')) = 0.865 ;
    cf_lpj(strcmp(cropList_lpj, 'StarchyRoots')) = 5.443 ;
    for c = 1:Ncrops
        thisCrop = cropList_lpj{c} ;
%         fprintf('%s lpj: %0.3f\n', thisCrop, cf_lpj(c)) ;
        
        % Apply to LPJ-GUESS sim
        isThisCrop = contains(data_fu_lpj.varNames, thisCrop) ;
        if length(find(isThisCrop)) ~= Ncrops
            error('length(find(isThisCrop)) ~= Ncrops')
        end
        yield_fu_lpj_max_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield_fu_lpj_max_xv(:,isThisCrop) ;
        
        % Apply to final outputs
        isThisCrop = contains(data_fu_out.varNames, thisCrop) ;
        if length(find(isThisCrop)) ~= Ncrops
            error('length(find(isThisCrop)) ~= Ncrops')
        end
        yield_fu_out_max_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield_fu_out_max_xv(:,isThisCrop) ;
    end
end

if strcmp(which_file, 'yield')
    tmp_which_file = 'Yield' ;
    units = 'tons/ha' ;
    % Convert kg/m2 to tons/ha
    yield_fu_lpj_max_xv = 10 * yield_fu_lpj_max_xv ;
    yield_fu_emu_max_xv = 10 * yield_fu_emu_max_xv ;
    yield_fu_out_max_xv = 10 * yield_fu_out_max_xv ;
elseif strcmp(which_file, 'gsirrigation')
    tmp_which_file = 'Irrig' ;
    units = 'mm' ;
else
    error('which_file (%s) not recognized', which_file) ;
end

title_left = 'LPJ-GUESS sim' ;
title_right = sprintf('LPJ-GUESS sim bl + %s emu %ss', ggcm, '\Delta') ;

for v = 1:length(data_fu_out.varNames)
    
    % What crop?
    thisVar_out = data_fu_out.varNames{v} ;
    thisCrop_out = getbasename(thisVar_out) ;
    thisCropi_out = getbasenamei(thisVar_out) ;
    switch thisCrop_out
        case 'CerealsC3'
            thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'max_wheat') ;
        case 'CerealsC4'
            thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'maize') ;
        case 'StarchyRoots'
            thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'spring_wheat') ;
        case 'Rice'
            switch ggcm
                case 'LPJ-GUESS'
                    thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'spring_wheat') ;
                otherwise
                    thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'rice') ;
            end
        case {'Oilcrops', 'Pulses'}
            switch ggcm
                case 'LPJ-GUESS'
                    thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'spring_wheat') ;
                otherwise
                    thisCropi_emu = strrep(thisCropi_out, thisCrop_out, 'soy') ;
            end
        otherwise
            error('thisCrop_out not recognized')
    end
    title_center = sprintf('Raw %s emu (%s)', ggcm, strrep(thisCropi_emu, '_', '\_')) ;
    
    % Get maps
    yield_fu_lpj_max_YX = lpjgu_vector2map( ...
        yield_fu_lpj_max_xv(:,strcmp(data_fu_lpj.varNames, thisVar_out)), ...
        [360 720], data_fu_lpj.list2map) ;
    yield_fu_out_max_YX = lpjgu_vector2map( ...
        yield_fu_out_max_xv(:,strcmp(data_fu_out.varNames, thisVar_out)), ...
        [360 720], data_fu_out.list2map) ;
    yield_fu_out_max_YX( ...
        isnan(yield_fu_out_max_YX) ...
        & ~isnan(yield_fu_lpj_max_YX)) ...
        = 0 ;
    thisVar_emu = sprintf('%s%s', thisCropi_emu, getN(thisVar_out)) ;
    if any(strcmp(data_fu_emu.varNames, thisVar_emu))
        yield_fu_emu_max_YX = lpjgu_vector2map( ...
            yield_fu_emu_max_xv(:,strcmp(data_fu_emu.varNames, thisVar_emu)), ...
            [360 720], data_fu_emu.list2map) ;
    else
        warning('No matching variable found in emu: %s', thisVar_emu)
        yield_fu_emu_max_YX = [] ;
    end
    
    % Set up figure
    spacing = [0.025 0.025] ; % v h
    yrange = 65:360 ;
    fontSize = 14 ;
    thisPos = figurePos ;
%     thisPos = [1         112        1440         693] ;
    figure('Color','w','Position',thisPos) ;
    
    % Maps with caxis "fit"
    subplot_tight(2,3,1,spacing) ;
    pcolor(yield_fu_lpj_max_YX(yrange,:)); shading flat; axis equal tight off
    hcb = colorbar('Location','SouthOutside') ;
    colormap(gca, this_colormap) ;
    xlabel(hcb, units)
    title(title_left)
    set(gca,'FontSize',fontSize)
    
    if ~isempty(yield_fu_emu_max_YX)
        subplot_tight(2,3,2,spacing) ;
        pcolor(yield_fu_emu_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        title(title_center)
        set(gca,'FontSize',fontSize)
    end
    
    subplot_tight(2,3,3,spacing) ;
    pcolor(yield_fu_out_max_YX(yrange,:)); shading flat; axis equal tight off
    hcb = colorbar('Location','SouthOutside') ;
    colormap(gca, this_colormap) ;
    xlabel(hcb, units)
    title(title_right)
    set(gca,'FontSize',fontSize)
    
    % Get limited caxis
    max_raw = max(yield_fu_lpj_max_YX(:)) ;
    if ~isempty(yield_fu_emu_max_YX)
        max_raw = max( ...
            max(yield_fu_emu_max_YX(:)), ...
            max_raw ...
            ) ;
    end
    max_out = max(yield_fu_out_max_YX(:)) ;
    new_caxis = [0 min(max_raw*1.25, max_out)] ;
%     new_caxis = [0 max_raw] ;

    % Maps with caxis limited
    subplot_tight(2,3,4,spacing) ;
    pcolor(yield_fu_lpj_max_YX(yrange,:)); shading flat; axis equal tight off
    hcb = colorbar('Location','SouthOutside') ;
    colormap(gca, this_colormap) ;
    xlabel(hcb, units)
    caxis(new_caxis)
    title(title_left)
    set(gca,'FontSize',fontSize)
    
    if ~isempty(yield_fu_emu_max_YX)
        subplot_tight(2,3,5,spacing) ;
        pcolor(yield_fu_emu_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        caxis(new_caxis)
        title(title_center)
        set(gca,'FontSize',fontSize)
    end
    
    subplot_tight(2,3,6,spacing) ;
    pcolor(yield_fu_out_max_YX(yrange,:)); shading flat; axis equal tight off
    hcb = colorbar('Location','SouthOutside') ;
    colormap(gca, this_colormap) ;
    xlabel(hcb, units)
    caxis(new_caxis)
    title(title_right)
    set(gca,'FontSize',fontSize)
    
    % Add overall title
    hold on; ax2 = axes(); hold off
    ax2.Position = [0 0 1 1] ;
    ht = text(ax2, 0.5, 0.97, ...
        sprintf('Max yield over future: %s', thisVar_out), ...
        'FontSize', fontSize*1.5, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center') ;
    ht.Units = 'normalized' ;
    ht.Position(1) = 0.5 ;
    
    % Add info about top vs. bottom rows
    hold on
    plot([0.05 0.95], 0.5*[1 1], '-k')
    ax2.Visible = 'off' ;
    text(ax2, 0.4, 0.53, ...
        '\uparrow Color axis: Fit', ...
        'FontSize', fontSize) ;
    text(ax2, 0.4, 0.47, ...
        '\downarrow Color axis: Limited', ...
        'FontSize', fontSize) ;
    hold off
    
    % Save
    filename = sprintf('%s/max%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisCropi_out) ;
    if ~exist(outDir_figs, 'dir')
        mkdir(outDir_figs)
    end
    export_fig(filename, '-r150')
    close
    
end
    
    
end