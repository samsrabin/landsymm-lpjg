function e2p_save_out_figs(data_fu_lpj, data_fu_lpj0, ...
    data_fu_emu, data_fu_out, ...
    ggcm, getN, outDir_figs, ...
    which_file, cropList_lpj_asEmu, figure_visibility, ...
    figure_extension, which_out_figs, overwrite_existing_figs, renderer)

if ~exist(outDir_figs, 'dir')
    mkdir(outDir_figs)
end

do_max = any(strcmp(which_out_figs, 'max')) ;
do_first = any(strcmp(which_out_figs, 'first') | strcmp(which_out_figs, 'first0')) ;
do_first0 = any(strcmp(which_out_figs, 'first0')) ;
do_4th = any(strcmp(which_out_figs, '4th') | strcmp(which_out_figs, '4th0')) ;
do_4th0 = any(strcmp(which_out_figs, '4th0')) ;
do_Nth0 = do_first0 | do_4th0 ;

if do_Nth0
    if isempty(data_fu_lpj0)
        error('do_Nth0 but isempty(data_fu_lpj0)')
    end
    Nlist_str = unique(cellfun(getN, data_fu_lpj.varNames, 'UniformOutput', false)) ;
    Nlist0_str = unique(cellfun(getN, data_fu_lpj0.varNames, 'UniformOutput', false)) ;
    Nlist_num = str2double(Nlist_str) ;
    Nlist0_num = str2double(Nlist0_str) ;
end

this_colormap = 'parula' ;
% this_colormap = 'jet' ;

% Convert kg/m2 to tons/ha
if strcmp(which_file, 'yield')
    data_fu_lpj.garr_xvt = 10 * data_fu_lpj.garr_xvt ;
    if ~isempty(data_fu_lpj0)
        data_fu_lpj0.garr_xvt = 10 * data_fu_lpj0.garr_xvt ;
    end
    data_fu_emu.garr_xvt = 10 * data_fu_emu.garr_xvt ;
    data_fu_out.garr_xvt = 10 * data_fu_out.garr_xvt ;
end

% Get arrays for "max" figures
yield_fu_lpj_max_xv = max(data_fu_lpj.garr_xvt,[],3) ;
if do_Nth0
    yield_fu_lpj0_max_xv = max(data_fu_lpj0.garr_xvt,[],3) ;
end
yield_fu_emu_max_xv = max(data_fu_emu.garr_xvt,[],3) ;
yield_fu_out_max_xv = max(data_fu_out.garr_xvt,[],3) ;

% Get arrays for "first" figures
yield1_lpj_xv = data_fu_lpj.garr_xvt(:,:,1) ;
if do_first0
    yield1_lpj0_xv = data_fu_lpj0.garr_xvt(:,:,1) ;
end
yield1_out_xv = data_fu_out.garr_xvt(:,:,1) ;

% Get arrays for "4th" figures
yield4_lpj_xv = data_fu_lpj.garr_xvt(:,:,4) ;
if do_4th0
    yield4_lpj0_xv = data_fu_lpj0.garr_xvt(:,:,4) ;
end
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
        
        % Apply to ORIGINAL LPJ-GUESS sim
        if do_first0 || do_4th0
            isThisCrop = contains(data_fu_lpj0.varNames, thisCrop) ;
            if length(find(isThisCrop)) ~= Ncrops
                error('length(find(isThisCrop)) ~= Ncrops')
            end
            yield_fu_lpj0_max_xv(:,isThisCrop) = ...
                cf_lpj(c) * yield_fu_lpj0_max_xv(:,isThisCrop) ;
            if do_first0
                yield1_lpj0_xv(:,isThisCrop) = ...
                    cf_lpj(c) * yield1_lpj0_xv(:,isThisCrop) ;
            end
            if do_4th0
                yield4_lpj0_xv(:,isThisCrop) = ...
                    cf_lpj(c) * yield4_lpj0_xv(:,isThisCrop) ;
            end
        end
        
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
title_right_max = sprintf('LPJ-GUESS sim bl + %s emu %ss', ggcm, '\Delta') ;

% Set up titles for "Nth timestep" figures
title_left_Nth = 'LPJ-GUESS sim' ;
title_left_Nth0 = 'LPJ-GUESS sim (original N%d)' ;

% Set up figure properties
spacing_max = [0.025 0.025] ; % v h
yrange = 65:360 ;
fontSize = 14 ;
thisPos = figurePos ;

for v = 1:length(data_fu_out.varNames)
    
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
    if do_max && ~exist(filename_max, 'file')
    
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
    if do_first && ~exist(filename_first, 'file')
        
        if do_first0
            Ny = 2 ;
            thisPos = figurePos ;
            spacing_first = [0.1 0.025] ; % v h
        else
            Ny = 1 ;
            thisPos = [1 325 1440 480] ;
            spacing_first = [0.025 0.025] ; % v h
        end

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
        
        % ORIGINAL LPJ-GUESS simulation
        if do_first0
            N_out_num = str2double(getN(thisVar_out)) ;
            N_match = any(Nlist0_num==N_out_num) ;
            if N_match || N_out_num==min(Nlist_num)
                if N_out_num==min(Nlist_num)
                    N_orig_num = min(Nlist0_num) ;
                else
                    N_orig_num = N_out_num ;
                end
                N_orig_ind = find(Nlist0_num==N_orig_num) ;
                N_orig_str = Nlist0_str{N_orig_ind} ; %#ok<FNDSB>
                thisVar_orig = [thisCropi_out N_orig_str] ;
                yield1_lpj0_YX = lpjgu_vector2map( ...
                    yield1_lpj0_xv(:,strcmp(data_fu_lpj0.varNames, thisVar_orig)), ...
                    [360 720], data_fu_lpj0.list2map) ;
                new_caxis = [0 max(new_caxis(end), max(max(yield1_lpj0_YX)))] ;
                    
                subplot_tight(Ny, 2, 3, spacing_first) ;
                pcolor(yield1_lpj0_YX(yrange,:)); shading flat; axis equal tight off
                hcb = colorbar('Location','SouthOutside') ;
                colormap(gca, this_colormap) ;
                caxis(new_caxis) ;
                xlabel(hcb, units)
                title(sprintf(title_left_Nth0, N_orig_num))
                set(gca, 'FontSize', fontSize)
            else
                subplot_tight(Ny, 2, 3, spacing_first) ;
                pcolor(NaN*yield1_lpj_YX(yrange,:)); shading flat; axis equal tight off
                hcb = colorbar('Location','SouthOutside') ;
                colormap(gca, this_colormap) ;
                caxis(new_caxis) ;
                xlabel(hcb, units)
                title('No original LPJ-GUESS sim with similar N')
                set(gca, 'FontSize', fontSize)
            end
        end

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
    if do_4th && ~exist(filename_4th, 'file')
        
        if do_4th0
            Ny = 2 ;
            thisPos = figurePos ;
            spacing_first = [0.1 0.025] ; % v h
        else
            Ny = 1 ;
            thisPos = [1 325 1440 480] ;
            spacing_first = [0.025 0.025] ; % v h
        end

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
        
        % ORIGINAL LPJ-GUESS simulation
        if do_4th0
            N_out_num = str2double(getN(thisVar_out)) ;
            N_match = any(Nlist0_num==N_out_num) ;
            if N_match || N_out_num==min(Nlist_num)
                if N_out_num==min(Nlist_num)
                    N_orig_num = min(Nlist0_num) ;
                else
                    N_orig_num = N_out_num ;
                end
                N_orig_ind = find(Nlist0_num==N_orig_num) ;
                N_orig_str = Nlist0_str{N_orig_ind} ; %#ok<FNDSB>
                thisVar_orig = [thisCropi_out N_orig_str] ;
                yield4_lpj0_YX = lpjgu_vector2map( ...
                    yield4_lpj0_xv(:,strcmp(data_fu_lpj0.varNames, thisVar_orig)), ...
                    [360 720], data_fu_lpj0.list2map) ;
                new_caxis = [0 max(new_caxis(end), max(max(yield4_lpj0_YX)))] ;
                    
                subplot_tight(Ny, 2, 3, spacing_first) ;
                pcolor(yield4_lpj0_YX(yrange,:)); shading flat; axis equal tight off
                hcb = colorbar('Location','SouthOutside') ;
                colormap(gca, this_colormap) ;
                caxis(new_caxis) ;
                xlabel(hcb, units)
                title(sprintf(title_left_Nth0, N_orig_num))
                set(gca, 'FontSize', fontSize)
            else
                subplot_tight(Ny, 2, 3, spacing_first) ;
                pcolor(NaN*yield4_lpj_YX(yrange,:)); shading flat; axis equal tight off
                hcb = colorbar('Location','SouthOutside') ;
                colormap(gca, this_colormap) ;
                caxis(new_caxis) ;
                xlabel(hcb, units)
                title('No original LPJ-GUESS sim with similar N')
                set(gca, 'FontSize', fontSize)
            end
        end

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
