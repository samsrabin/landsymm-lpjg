function e2p_save_out_figs(data_fu_lpj, ...
    data_fu_emu, data_fu_out, ...
    ggcm, getN, outDir_figs, ...
    which_file, cropList_lpj_asEmu, figure_visibility, ...
    figure_extension, which_out_figs, overwrite_existing_figs, renderer)

if ~exist(outDir_figs, 'dir')
    mkdir(outDir_figs)
end

do_max = any(strcmp(which_out_figs, 'max')) ;
do_first = strcmp(which_out_figs, 'first') ;
do_4th = strcmp(which_out_figs, '4th') ;

this_colormap = 'parula' ;
% this_colormap = 'jet' ;

% Convert kg/m2 to tons/ha
if strcmp(which_file, 'yield')
    data_fu_lpj.garr_xvt = 10 * data_fu_lpj.garr_xvt ;
    data_fu_emu.garr_xvt = 10 * data_fu_emu.garr_xvt ;
    data_fu_out.garr_xvt = 10 * data_fu_out.garr_xvt ;
end

% Get arrays for "max" figures
yield_fu_lpj_max_xv = max(data_fu_lpj.garr_xvt,[],3) ;
yield_fu_emu_max_xv = max(data_fu_emu.garr_xvt,[],3) ;
yield_fu_out_max_xv = max(data_fu_out.garr_xvt,[],3) ;

% Get arrays for "first" figures
yield1_lpj_xv = data_fu_lpj.garr_xvt(:,:,1) ;
yield1_out_xv = data_fu_out.garr_xvt(:,:,1) ;

% Get arrays for "4th" figures
yield4_lpj_xv = data_fu_lpj.garr_xvt(:,:,4) ;
yield4_out_xv = data_fu_out.garr_xvt(:,:,4) ;

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
        yield_fu_lpj_max_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield_fu_lpj_max_xv(:,isThisCrop) ;
        yield1_lpj_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield1_lpj_xv(:,isThisCrop) ;
        yield4_lpj_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield4_lpj_xv(:,isThisCrop) ;
        
        % Apply to final outputs
        isThisCrop = contains(data_fu_out.varNames, thisCrop) ;
        yield_fu_out_max_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield_fu_out_max_xv(:,isThisCrop) ;
        yield1_out_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield1_out_xv(:,isThisCrop) ;
        yield4_out_xv(:,isThisCrop) = ...
            cf_lpj(c) * yield4_out_xv(:,isThisCrop) ;
    end
end

if strcmp(which_file, 'yield')
    tmp_which_file = 'Yield' ;
    units = 'tons/ha' ;
elseif strcmp(which_file, 'gsirrigation')
    tmp_which_file = 'Irrig' ;
    units = 'mm' ;
else
    error('which_file (%s) not recognized', which_file) ;
end

% Set up titles for "max over future" figure
title_left_max = 'LPJ-GUESS sim' ;
title_right_max_orig = sprintf('Phase 2 %s sim bl + emu %ss', ...
    ggcm, '\Delta') ;

% Set up titles for "Nth timestep" figures
title_left_Nth = 'LPJ-GUESS sim' ;

% Set up figure properties
spacing_max = [0.025 0.025] ; % v h
yrange = 65:360 ;
fontSize = 14 ;
thisPos = figurePos ;

for v = 1:length(data_fu_out.varNames)
    
    if isfield(data_fu_out, 'actually_emu_char') ...
    && ~strcmp('sim', data_fu_out.actually_emu_char{v})
        title_right_max = strrep(title_right_max_orig, ...
            'sim bl', [data_fu_out.actually_emu_char{v} ' bl']) ;
    else
        title_right_max = title_right_max_orig ;
    end
    
    % What crop?
    thisVar_out = data_fu_out.varNames{v} ;
    thisCrop_out = getbasename(thisVar_out) ;
    thisCropi_out = getbasenamei(thisVar_out) ;
    thisCropi_emu = strrep(thisCropi_out, thisCrop_out, ...
        cropList_lpj_asEmu{strcmp(cropList_lpj, thisCrop_out)}) ;
    thisVar_emu = sprintf('%s%s', thisCropi_emu, getN(thisVar_out)) ;
    
    % Skip if looking at irrigation of a rainfed crop
    if strcmp(which_file, 'gsirrigation') && strcmp(thisCrop_out, thisCropi_out)
        continue
    end
    
    % Get filenames
    filename_max = sprintf('%s/max%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisVar_out) ;
    if strcmp(figure_extension, 'fig')
        filename_max = strrep(filename_max, '.png', '.fig') ;
    end
    filename_first = sprintf('%s/first%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisVar_out) ;
    if strcmp(figure_extension, 'fig')
        filename_first = strrep(filename_first, '.png', '.fig') ;
    end
    filename_4th = sprintf('%s/fourth%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisVar_out) ;
    if strcmp(figure_extension, 'fig')
        filename_4th = strrep(filename_4th, '.png', '.fig') ;
    end
    
    % Skip if figure files already exist
    if ((do_max && exist(filename_max, 'file')) || ~do_max) ...
    && ((do_first && exist(filename_first, 'file')) || ~do_first) ...
    && ((do_4th && exist(filename_4th, 'file')) || ~do_4th) ...
    && ~overwrite_existing_figs
        fprintf('        %s %s figure skipped (exists)\n', thisVar_out, which_file)
        continue
    end
    
    fprintf('        %s %s...\n', thisVar_out, which_file)
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Maximum yield over future %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_max && (overwrite_existing_figs || ~exist(filename_max, 'file'))
    
        title_center_max = sprintf('Raw %s emu (%s)', ggcm, strrep(thisCropi_emu, '_', '\_')) ;

        % Get maps
        yield_fu_lpj_max_YX = lpjgu_vector2map( ...
            yield_fu_lpj_max_xv(:,strcmp(data_fu_lpj.varNames, thisVar_out)), ...
            [360 720], data_fu_lpj.list2map) ;
        yield_fu_out_max_YX = lpjgu_vector2map( ...
            yield_fu_out_max_xv(:,strcmp(data_fu_out.varNames, thisVar_out)), ...
            [360 720], data_fu_out.list2map) ;
        if any(strcmp(data_fu_emu.varNames, thisVar_emu))
            yield_fu_emu_max_YX = lpjgu_vector2map( ...
                yield_fu_emu_max_xv(:,strcmp(data_fu_emu.varNames, thisVar_emu)), ...
                [360 720], data_fu_emu.list2map) ;
        else
            warning('No matching variable found in emu: %s', thisVar_emu)
            yield_fu_emu_max_YX = [] ;
        end

        % Set up figure
        figure('Color','w','Position',thisPos, 'Visible', figure_visibility) ;

        % Maps with caxis "fit"
        subplot_tight(2,3,1,spacing_max) ;
        pcolor(yield_fu_lpj_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        title(title_left_max)
        set(gca,'FontSize',fontSize)

        if ~isempty(yield_fu_emu_max_YX)
            subplot_tight(2,3,2,spacing_max) ;
            pcolor(yield_fu_emu_max_YX(yrange,:)); shading flat; axis equal tight off
            hcb = colorbar('Location','SouthOutside') ;
            colormap(gca, this_colormap) ;
            xlabel(hcb, units)
            title(title_center_max)
            set(gca,'FontSize',fontSize)
        end

        subplot_tight(2,3,3,spacing_max) ;
        pcolor(yield_fu_out_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        title(title_right_max)
        set(gca,'FontSize',fontSize)

        % Get limited caxis
        max_raw = max(yield_fu_lpj_max_YX(:)) ;
    %     if ~isempty(yield_fu_emu_max_YX)
    %         max_raw = max( ...
    %             max(yield_fu_emu_max_YX(:)), ...
    %             max_raw ...
    %             ) ;
    %     end
        max_out = max(yield_fu_out_max_YX(:)) ;
        if max_out >= max_raw
            new_caxis = [0 min(max_raw*1.25, max_out)] ;
        else
            new_caxis = [0 max_raw] ;
        end

        % Maps with caxis limited
        subplot_tight(2,3,4,spacing_max) ;
        pcolor(yield_fu_lpj_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        caxis(new_caxis)
        title(title_left_max)
        set(gca,'FontSize',fontSize)

        if ~isempty(yield_fu_emu_max_YX)
            subplot_tight(2,3,5,spacing_max) ;
            pcolor(yield_fu_emu_max_YX(yrange,:)); shading flat; axis equal tight off
            hcb = colorbar('Location','SouthOutside') ;
            colormap(gca, this_colormap) ;
            xlabel(hcb, units)
            caxis(new_caxis)
            title(title_center_max)
            set(gca,'FontSize',fontSize)
        end

        subplot_tight(2,3,6,spacing_max) ;
        pcolor(yield_fu_out_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        caxis(new_caxis)
        title(title_right_max)
        set(gca,'FontSize',fontSize)

        % Add overall title
        hold on; ax2 = axes(); hold off
        ax2.Position = [0 0 1 1] ;
        ht = text(ax2, 0.5, 0.97, ...
            sprintf('Max %s over future: %s', lower(tmp_which_file), thisVar_out), ...
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

        %    fprintf('%s... ', toc_hms(toc))

        % Save
        tic
        if strcmp(figure_extension, 'png')
            export_fig(filename_max, '-r100', renderer)
        elseif strcmp(figure_extension, 'fig')
            savefig(filename_max)
        else
            error('figure_extension %s not recognized', figure_extension)
        end
        close
        %    fprintf('%s.\n ', toc_hms(toc))
    
    end


    %%%%%%%%%%%%%%%%%%%%%%
    %%% First timestep %%%
    %%%%%%%%%%%%%%%%%%%%%%
    if do_first && (overwrite_existing_figs || ~exist(filename_first, 'file'))
        
        Ny = 1 ;
        thisPos = [1 325 1440 480] ;
        spacing_first = [0.025 0.025] ; % v h

        % Get maps
        yield1_lpj_YX = lpjgu_vector2map( ...
            yield1_lpj_xv(:,strcmp(data_fu_lpj.varNames, thisVar_out)), ...
            [360 720], data_fu_lpj.list2map) ;
        yield1_out_YX = lpjgu_vector2map( ...
            yield1_out_xv(:,strcmp(data_fu_out.varNames, thisVar_out)), ...
            [360 720], data_fu_out.list2map) ;

        % Set up figure
        figure('Color', 'w', ...
            'Position', thisPos, ...
            'Visible', figure_visibility) ;
        new_caxis = [0 max(max(max(yield1_lpj_YX)), max(max(yield1_out_YX)))] ;
        title_right_first = sprintf('Emulated %s (%s) after processing', ...
            ggcm, strrep(thisCropi_emu, '_', '\_')) ;
        
        % LPJ-GUESS simulation
        subplot_tight(Ny, 2, 1, spacing_first) ;
        pcolor(yield1_lpj_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        caxis(new_caxis) ;
        xlabel(hcb, units)
        title(title_left_Nth)
        set(gca, 'FontSize', fontSize)
        
        % Emulator
        subplot_tight(Ny, 2, 2, spacing_first) ;
        pcolor(yield1_out_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        caxis(new_caxis) ;
        xlabel(hcb, units)
        title(title_right_first)
        set(gca, 'FontSize', fontSize)

        % Add overall title
        hold on; ax2 = axes(); hold off
        ax2.Position = [0 0 1 1] ;
        ht = text(ax2, 0.5, 0.97, ...
            sprintf('First-timestep %s: %s', lower(tmp_which_file), thisVar_out), ...
            'FontSize', fontSize*1.5, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center') ;
        ht.Units = 'normalized' ;
        ht.Position(1) = 0.5 ;
        ax2.Visible = 'off' ;
        
        % Save
        tic
        if strcmp(figure_extension, 'png')
            export_fig(filename_first, '-r100', renderer)
        elseif strcmp(figure_extension, 'fig')
            savefig(filename_first)
        else
            error('figure_extension %s not recognized', figure_extension)
        end
        close
        %    fprintf('%s.\n ', toc_hms(toc))
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Fourth timestep %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    if do_4th && (overwrite_existing_figs || ~exist(filename_4th, 'file'))
        
        Ny = 1 ;
        thisPos = [1 325 1440 480] ;
        spacing_first = [0.025 0.025] ; % v h

        % Get maps
        yield4_lpj_YX = lpjgu_vector2map( ...
            yield4_lpj_xv(:,strcmp(data_fu_lpj.varNames, thisVar_out)), ...
            [360 720], data_fu_lpj.list2map) ;
        yield4_out_YX = lpjgu_vector2map( ...
            yield4_out_xv(:,strcmp(data_fu_out.varNames, thisVar_out)), ...
            [360 720], data_fu_out.list2map) ;

        % Set up figure
        figure('Color', 'w', ...
            'Position', thisPos, ...
            'Visible', figure_visibility) ;
        new_caxis = [0 max(max(max(yield4_lpj_YX)), max(max(yield4_out_YX)))] ;
        title_right_first = sprintf('Emulated %s (%s) after processing', ...
            ggcm, strrep(thisCropi_emu, '_', '\_')) ;
        
        % LPJ-GUESS simulation
        subplot_tight(Ny, 2, 1, spacing_first) ;
        pcolor(yield4_lpj_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        caxis(new_caxis) ;
        xlabel(hcb, units)
        title(title_left_Nth)
        set(gca, 'FontSize', fontSize)
        
        % Emulator
        subplot_tight(Ny, 2, 2, spacing_first) ;
        pcolor(yield4_out_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        caxis(new_caxis) ;
        xlabel(hcb, units)
        title(title_right_first)
        set(gca, 'FontSize', fontSize)

        % Add overall title
        hold on; ax2 = axes(); hold off
        ax2.Position = [0 0 1 1] ;
        ht = text(ax2, 0.5, 0.97, ...
            sprintf('Fourth-timestep %s: %s', lower(tmp_which_file), thisVar_out), ...
            'FontSize', fontSize*1.5, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center') ;
        ht.Units = 'normalized' ;
        ht.Position(1) = 0.5 ;
        ax2.Visible = 'off' ;
        
        % Save
        tic
        if strcmp(figure_extension, 'png')
            export_fig(filename_4th, '-r100', renderer)
        elseif strcmp(figure_extension, 'fig')
            savefig(filename_4th)
        else
            error('figure_extension %s not recognized', figure_extension)
        end
        close
        %    fprintf('%s.\n ', toc_hms(toc))
    
    end


end
    
    
end
