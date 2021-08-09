function e2p_save_out_figs(data_fu_lpj, ...
    data_fu_emuA, data_fu_emuB, ...
    ggcm, outDir_figs, ...
    which_file, figure_visibility, ...
    figure_extension, which_out_figs, overwrite_existing_figs, renderer, ...
    cf_lpj, cf_emu)

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
    data_fu_emuA.garr_xvt = 10 * data_fu_emuA.garr_xvt ;
    data_fu_emuB.garr_xvt = 10 * data_fu_emuB.garr_xvt ;
end

% Get max_wheat(i) for emuA, if needed
cropList_emuA = unique(getbasename(data_fu_emuA.varNames)) ;
combineCrops = {'max_wheat', {'spring_wheat', 'winter_wheat'}} ;
if ~any(strcmp(cropList_emuA, 'max_wheat'))
    data_fu_emuA = e2p_split_combine_burn(data_fu_emuA, combineCrops, ...
        cropList_emuA, [cropList_emuA {'max_wheat'}], []) ;
    cropList_emuA = unique(getbasename(data_fu_emuA.varNames)) ;
end

% Get arrays for "max" figures
yield_fu_lpj_max_xv = max(data_fu_lpj.garr_xvt,[],3) ;
yield_fu_emuA_max_xv = max(data_fu_emuA.garr_xvt,[],3) ;
yield_fu_emuB_max_xv = max(data_fu_emuB.garr_xvt,[],3) ;

% Get arrays for "first" figures
yield1_lpj_xv = data_fu_lpj.garr_xvt(:,:,1) ;
yield1_emuB_xv = data_fu_emuB.garr_xvt(:,:,1) ;

% Get arrays for "4th" figures
yield4_lpj_xv = data_fu_lpj.garr_xvt(:,:,4) ;
yield4_emuB_xv = data_fu_emuB.garr_xvt(:,:,4) ;

% Apply calibration factors
cropList_lpj = unique(getbasename(data_fu_lpj.varNames)) ;
Ncrops = length(cropList_lpj) ;
if strcmp(which_file, 'yield')
%     cf_lpj = get_CFs_Rabin2020(cropList_lpj) ;
%     cf_lpj = get_CFs_e2p(cropList_lpj) ;
    if isempty(cf_lpj) || isempty(cf_emu)
        return
    end
    for c = 1:Ncrops
        thisCrop = cropList_lpj{c} ;
        thisCF_lpj = cf_lpj(c) ;
        if isnan(thisCF_lpj)
            error('LPJ-GUESS calibration factor not defined for %s', ...
                thisCrop)
        end
%         fprintf('%s lpj: %0.3f\n', thisCrop, thisCF_lpj) ;
        thisCF_emu = cf_emu(c) ;
        if isnan(thisCF_emu)
            error('%s calibration factor not defined for %s', ...
                ggcm, thisCrop)
        end
%         fprintf('%s emu: %0.3f\n', thisCrop, thisCF_emu) ;
        
        % Apply to LPJ-GUESS sim
        isThisCrop = contains(data_fu_lpj.varNames, thisCrop) ;
        yield_fu_lpj_max_xv(:,isThisCrop) = ...
            thisCF_lpj * yield_fu_lpj_max_xv(:,isThisCrop) ;
        yield1_lpj_xv(:,isThisCrop) = ...
            thisCF_lpj * yield1_lpj_xv(:,isThisCrop) ;
        yield4_lpj_xv(:,isThisCrop) = ...
            thisCF_lpj * yield4_lpj_xv(:,isThisCrop) ;
        
        % Apply to final outputs
        isThisCrop = contains(data_fu_emuB.varNames, thisCrop) ;
        yield_fu_emuB_max_xv(:,isThisCrop) = ...
            thisCF_emu * yield_fu_emuB_max_xv(:,isThisCrop) ;
        yield1_emuB_xv(:,isThisCrop) = ...
            thisCF_emu * yield1_emuB_xv(:,isThisCrop) ;
        yield4_emuB_xv(:,isThisCrop) = ...
            thisCF_emu * yield4_emuB_xv(:,isThisCrop) ;
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

for v = 1:length(data_fu_emuB.varNames)
    
    if isfield(data_fu_emuB, 'actually_emuBL_char') ...
    && ~strcmp('sim', data_fu_emuB.actually_emuBL_char{v})
        title_right_max = strrep(title_right_max_orig, ...
            'sim bl', [data_fu_emuB.actually_emuBL_char{v} ' bl']) ;
    else
        title_right_max = title_right_max_orig ;
    end
    
    % What crop?
    thisVar_emuB = data_fu_emuB.varNames{v} ;
    thisCrop_emuB = getbasename(thisVar_emuB) ;
    thisCropi_emuB = getbasenamei(thisVar_emuB) ;
    if any(strcmp(cropList_emuA, thisCrop_emuB))
        thisCrop_emuA = thisCrop_emuB ;
    else
        thisCrop_emuA = e2p_translate_crops_2emu({thisCrop_emuB}, ...
            cropList_emuA, '', false) ;
        thisCrop_emuA = thisCrop_emuA{1} ;
    end
    thisCropi_emuA = strrep(thisCropi_emuB, thisCrop_emuB, thisCrop_emuA) ;
    thisVar_emuA = sprintf('%s%s', thisCropi_emuA, getN_char(thisVar_emuB)) ;
    
    % Skip if looking at irrigation of a rainfed crop
    if strcmp(which_file, 'gsirrigation') && strcmp(thisCrop_emuB, thisCropi_emuB)
        continue
    end
    
    % Get filenames
    filename_max = sprintf('%s/max%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisVar_emuB) ;
    if strcmp(figure_extension, 'fig')
        filename_max = strrep(filename_max, '.png', '.fig') ;
    end
    filename_first = sprintf('%s/first%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisVar_emuB) ;
    if strcmp(figure_extension, 'fig')
        filename_first = strrep(filename_first, '.png', '.fig') ;
    end
    filename_4th = sprintf('%s/fourth%ss_%s_%s.png', outDir_figs, tmp_which_file, ggcm, thisVar_emuB) ;
    if strcmp(figure_extension, 'fig')
        filename_4th = strrep(filename_4th, '.png', '.fig') ;
    end
    
    % Skip if figure files already exist
    if ((do_max && exist(filename_max, 'file')) || ~do_max) ...
    && ((do_first && exist(filename_first, 'file')) || ~do_first) ...
    && ((do_4th && exist(filename_4th, 'file')) || ~do_4th) ...
    && ~overwrite_existing_figs
        fprintf('        %s %s figure skipped (exists)\n', thisVar_emuB, which_file)
        continue
    end
    
    fprintf('        %s %s...\n', thisVar_emuB, which_file)
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Maximum yield over future %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_max && (overwrite_existing_figs || ~exist(filename_max, 'file'))
    
        title_center_max = sprintf('Raw %s emu (%s)', ggcm, strrep(thisCropi_emuA, '_', '\_')) ;
        
        % Get maps
        yield_fu_lpj_max_YX = lpjgu_vector2map( ...
            yield_fu_lpj_max_xv(:,strcmp(data_fu_lpj.varNames, thisVar_emuB)), ...
            [360 720], data_fu_lpj.list2map) ;
        yield_fu_emuB_max_YX = lpjgu_vector2map( ...
            yield_fu_emuB_max_xv(:,strcmp(data_fu_emuB.varNames, thisVar_emuB)), ...
            [360 720], data_fu_emuB.list2map) ;
        if any(strcmp(data_fu_emuA.varNames, thisVar_emuA))
            yield_fu_emuA_max_YX = lpjgu_vector2map( ...
                yield_fu_emuA_max_xv(:,strcmp(data_fu_emuA.varNames, thisVar_emuA)), ...
                [360 720], data_fu_emuA.list2map) ;
        else
            warning('No matching variable found in emu: %s', thisVar_emuA)
            yield_fu_emuA_max_YX = [] ;
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

        if ~isempty(yield_fu_emuA_max_YX)
            subplot_tight(2,3,2,spacing_max) ;
            pcolor(yield_fu_emuA_max_YX(yrange,:)); shading flat; axis equal tight off
            hcb = colorbar('Location','SouthOutside') ;
            colormap(gca, this_colormap) ;
            xlabel(hcb, units)
            title(title_center_max)
            set(gca,'FontSize',fontSize)
        end

        subplot_tight(2,3,3,spacing_max) ;
        pcolor(yield_fu_emuB_max_YX(yrange,:)); shading flat; axis equal tight off
        hcb = colorbar('Location','SouthOutside') ;
        colormap(gca, this_colormap) ;
        xlabel(hcb, units)
        title(title_right_max)
        set(gca,'FontSize',fontSize)

        % Get limited caxis
        max_raw = max(yield_fu_lpj_max_YX(:)) ;
    %     if ~isempty(yield_fu_emuA_max_YX)
    %         max_raw = max( ...
    %             max(yield_fu_emuA_max_YX(:)), ...
    %             max_raw ...
    %             ) ;
    %     end
        max_emuB = max(yield_fu_emuB_max_YX(:)) ;
        if isnan(max_raw) || isnan(max_emuB)
            error('yield_fu_[lpj and/or emuB]_max_YX all NaN!')
        elseif max_emuB >= max_raw
            new_caxis = [0 min(max_raw*1.25, max_emuB)] ;
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

        if ~isempty(yield_fu_emuA_max_YX)
            subplot_tight(2,3,5,spacing_max) ;
            pcolor(yield_fu_emuA_max_YX(yrange,:)); shading flat; axis equal tight off
            hcb = colorbar('Location','SouthOutside') ;
            colormap(gca, this_colormap) ;
            xlabel(hcb, units)
            caxis(new_caxis)
            title(title_center_max)
            set(gca,'FontSize',fontSize)
        end

        subplot_tight(2,3,6,spacing_max) ;
        pcolor(yield_fu_emuB_max_YX(yrange,:)); shading flat; axis equal tight off
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
            sprintf('Max %s over future: %s', lower(tmp_which_file), thisVar_emuB), ...
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
            yield1_lpj_xv(:,strcmp(data_fu_lpj.varNames, thisVar_emuB)), ...
            [360 720], data_fu_lpj.list2map) ;
        yield1_emuB_YX = lpjgu_vector2map( ...
            yield1_emuB_xv(:,strcmp(data_fu_emuB.varNames, thisVar_emuB)), ...
            [360 720], data_fu_emuB.list2map) ;

        % Set up figure
        figure('Color', 'w', ...
            'Position', thisPos, ...
            'Visible', figure_visibility) ;
        new_caxis = [0 max(max(max(yield1_lpj_YX)), max(max(yield1_emuB_YX)))] ;
        title_right_first = sprintf('Emulated %s (%s) after processing', ...
            ggcm, strrep(thisCropi_emuA, '_', '\_')) ;
        
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
        pcolor(yield1_emuB_YX(yrange,:)); shading flat; axis equal tight off
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
            sprintf('First-timestep %s: %s', lower(tmp_which_file), thisVar_emuB), ...
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
            yield4_lpj_xv(:,strcmp(data_fu_lpj.varNames, thisVar_emuB)), ...
            [360 720], data_fu_lpj.list2map) ;
        yield4_emuB_YX = lpjgu_vector2map( ...
            yield4_emuB_xv(:,strcmp(data_fu_emuB.varNames, thisVar_emuB)), ...
            [360 720], data_fu_emuB.list2map) ;

        % Set up figure
        figure('Color', 'w', ...
            'Position', thisPos, ...
            'Visible', figure_visibility) ;
        new_caxis = [0 max(max(max(yield4_lpj_YX)), max(max(yield4_emuB_YX)))] ;
        title_right_first = sprintf('Emulated %s (%s) after processing', ...
            ggcm, strrep(thisCropi_emuA, '_', '\_')) ;
        
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
        pcolor(yield4_emuB_YX(yrange,:)); shading flat; axis equal tight off
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
            sprintf('Fourth-timestep %s: %s', lower(tmp_which_file), thisVar_emuB), ...
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


function cf_lpj = get_CFs_Rabin2020(cropList_lpj)

cf_lpj = nan(length(cropList_lpj),1) ;
cf_lpj(strcmp(cropList_lpj, 'CerealsC3')) = 1.204 ;
cf_lpj(strcmp(cropList_lpj, 'CerealsC4')) = 0.738 ;
cf_lpj(strcmp(cropList_lpj, 'Rice')) = 1.052 ;
cf_lpj(strcmp(cropList_lpj, 'Oilcrops')) = 0.687 ;
cf_lpj(strcmp(cropList_lpj, 'Pulses')) = 0.865 ;
cf_lpj(strcmp(cropList_lpj, 'StarchyRoots')) = 5.443 ;

end


function cf_lpj = get_CFs_e2p(cropList_lpj)
% /Volumes/Reacher/G2P/outputs_LPJG/remap12_2016/calibration_Ks3_Oilcrops_ownDates/output-2021-05-03-230307

cf_lpj = nan(length(cropList_lpj),1) ;
cf_lpj(strcmp(cropList_lpj, 'CerealsC3')) = 1.056 ;
cf_lpj(strcmp(cropList_lpj, 'CerealsC4')) = 0.626 ;
cf_lpj(strcmp(cropList_lpj, 'Rice')) = 1.576 ;
cf_lpj(strcmp(cropList_lpj, 'OilNfix')) = 0.906 ;
cf_lpj(strcmp(cropList_lpj, 'OilOther')) = 0.686 ;
cf_lpj(strcmp(cropList_lpj, 'Pulses')) = 0.693 ;
cf_lpj(strcmp(cropList_lpj, 'StarchyRoots')) = 7.537 ;
cf_lpj(strcmp(cropList_lpj, 'Sugarbeet')) = 16.082 ;
cf_lpj(strcmp(cropList_lpj, 'Sugarcane')) = 11.489 ;
cf_lpj(strcmp(cropList_lpj, 'FruitAndVeg')) = 10.433 ;

end

