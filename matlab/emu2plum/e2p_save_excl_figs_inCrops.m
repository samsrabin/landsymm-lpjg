function e2p_save_excl_figs_inCrops( ...
    ggcm, which_file, gridlist, excl_vecs, cropList_emu_basei, ...
    figure_visibility, figure_extension, ...
    outDir_excl_figs_inCrops, overwrite_existing_figs, renderer)

% Unpack
missing_emu_xc = excl_vecs{1} ;
if strcmp(which_file, 'yield')
    missing_agmerra_xc = excl_vecs{2} ;
    isexcl_lowBL_emu_xc = excl_vecs{3} ;
    isexcl_lowBL_agmerra_xc = excl_vecs{4} ;
elseif strcmp(which_file,'gsirrigation')
    missing_yield_xc = excl_vecs{2} ;
    isexcl_lowBLyield_xc = excl_vecs{3} ;
    isexcl_all0_emu_c = excl_vecs{4} ;
else
    error('which_file (%s) not recognized', which_file)
end

% Set up for figures
if ~exist(outDir_excl_figs_inCrops, 'dir')
    mkdir(outDir_excl_figs_inCrops) ;
end
spacing = [0.025 0.025] ; % v h
yrange = 65:360 ;
fontSize = 14 ;
if strcmp(which_file, 'yield')
    blue = [61 45 165]/255 ;
    green = [38 190 183]/255 ;
    yellow = [249 248 59]/255 ;
    gray = 0.5*[1 1 1] ;
    this_colormap = [yellow ; blue ; green ; gray];
    x_bgy = 0 ;
    y_bgy_top = 0.3 ;
    down_bgy = 0.055 ;
    y_b = y_bgy_top ;
    y_g = y_b - down_bgy ;
    y_y = y_g - down_bgy ;
    y_gray = y_y - down_bgy ;
elseif strcmp(which_file, 'gsirrigation')
    this_colormap = 'parula' ;
    x_bgy = 0 ;
    y_bgy_top = 0.3 ;
    down_bgy = 0.055 ;
    y_g = y_bgy_top ;
    y_b = y_g - down_bgy ;
    y_y = y_b - down_bgy ;
else
    error('which_file (%s) not recognized', which_file)
end
fontSize_bgy = fontSize - 2 ;
fontSize_big = fontSize + 4 ;

% Set up for filename
if strcmp(which_file,'yield')
    thisFile = 'yield' ;
elseif strcmp(which_file,'gsirrigation')
    thisFile = 'irrig' ;
else
    error('which_file (%s) not recognized', which_file)
end

for c = 1:length(cropList_emu_basei)
    thisCrop = cropList_emu_basei{c} ;
    thisCrop_title = strrep(thisCrop, '_', '\_') ;
    
    % If irrigation, skip rainfed crops
    if strcmp(which_file,'gsirrigation') && ~strcmp(thisCrop(end), 'i')
        continue
    end
    
    filename = sprintf('%s/exclusions_%s_%s_%s.png', ...
        outDir_excl_figs_inCrops, ggcm, thisCrop, thisFile) ;
    if strcmp(figure_extension, 'fig')
        filename = strrep(filename, '.png', '.fig') ;
    elseif ~strcmp(figure_extension, 'png')
        error('figure_extension %s not recognized', figure_extension)
    end
    
    if exist(filename, 'file') && ~overwrite_existing_figs
        fprintf('        %s %s exclusions figure skipped (exists)\n', thisCrop, which_file)
        continue
    end
        
    if strcmp(which_file,'yield')
        
        figure('Color', 'w', 'Position', figurePos, ...
            'Visible', figure_visibility) ;
        legend_text = { ...
            'Blue: Excluded in a previous step', ...
            'Green: Excluded in a previous step AND this step', ...
            'Yellow: Excluded in this step', ...
            'Gray: Still included after this step'} ;
        
        % Where is emulator yield missing?
        h1 = subplot_tight(2,2,1,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        missing_emu_YX = lpjgu_vector2map( ...
            missing_emu_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isyellow_YX = gridlist.mask_YX & missing_emu_YX==1 ;
        map_YX(gridlist.mask_YX & missing_emu_YX==0) = 3 ;
        map_YX(isyellow_YX) = 0 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 3])
        ok_sofar_YX = gridlist.mask_YX & missing_emu_YX==0 ;
        bad_sofar_YX = missing_emu_YX==1 ;
        title('1: Missing yield: Emulator')
        set(gca, 'FontSize', fontSize)
        nyellow = length(find(isyellow_YX)) ;
        nblue = length(find(map_YX==1)) ;
        ngreen = length(find(map_YX==2)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_y = sprintf('Yellow: %d', nyellow) ;
        text_g = sprintf('Green: %d', ngreen) ;
        text_gray = sprintf('Gray: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_gray, text_gray, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Where there is emulator yield, where is AgMERRA yield missing?
        h2 = subplot_tight(2,2,2,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        map_YX(ok_sofar_YX) = 3 ;
        map_YX(bad_sofar_YX) = 1 ;
        missing_agmerra_YX = lpjgu_vector2map( ...
            missing_agmerra_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isyellow_YX = ok_sofar_YX & missing_agmerra_YX==1 ;
        map_YX(isyellow_YX) = 0 ;
        map_YX(bad_sofar_YX & missing_agmerra_YX==1) = 2 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 3])
        ok_sofar_YX = ok_sofar_YX & missing_agmerra_YX==0 ;
        bad_sofar_YX = bad_sofar_YX | missing_agmerra_YX==1 ;
        title('2: Missing yield: Phase 2 baseline sim')
        set(gca, 'FontSize', fontSize)
        nyellow = length(find(isyellow_YX)) ;
        nblue = length(find(map_YX==1)) ;
        ngreen = length(find(map_YX==2)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_y = sprintf('Yellow: %d', nyellow) ;
        text_g = sprintf('Green: %d', ngreen) ;
        text_gray = sprintf('Gray: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_gray, text_gray, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Where there is emulator and AgMERRA yield, where is there
        % exclusion based on too-low yield in emulator?
        h3 = subplot_tight(2,2,3,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        map_YX(ok_sofar_YX) = 3 ;
        map_YX(bad_sofar_YX) = 1 ;
        excl_emu_YX = lpjgu_vector2map( ...
            isexcl_lowBL_emu_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isyellow_YX = ok_sofar_YX & excl_emu_YX==1 ;
        map_YX(isyellow_YX) = 0 ;
        map_YX(bad_sofar_YX & excl_emu_YX==1) = 2 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 3])
        ok_sofar_YX = ok_sofar_YX & excl_emu_YX==0 ;
        bad_sofar_YX = bad_sofar_YX | excl_emu_YX==1 ;
        title('3: Yield too low: Emulator')
        set(gca, 'FontSize', fontSize)
        nyellow = length(find(isyellow_YX)) ;
        nblue = length(find(map_YX==1)) ;
        ngreen = length(find(map_YX==2)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_y = sprintf('Yellow: %d', nyellow) ;
        text_g = sprintf('Green: %d', ngreen) ;
        text_gray = sprintf('Gray: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_gray, text_gray, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Where there is emulator and AgMERRA yield, and emulator yield is
        % not too low, where is there exclusion based on too-low yield in
        % AgMERRA?
        h4 = subplot_tight(2,2,4,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        map_YX(ok_sofar_YX) = 3 ;
        map_YX(bad_sofar_YX) = 1 ;
        excl_agmerra_YX = lpjgu_vector2map( ...
            isexcl_lowBL_agmerra_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isyellow_YX = ok_sofar_YX & excl_agmerra_YX==1 ;
        map_YX(isyellow_YX) = 0 ;
        map_YX(bad_sofar_YX & excl_agmerra_YX==1) = 2 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 3])
        ok_sofar_YX = ok_sofar_YX & excl_agmerra_YX==0 ;
        bad_sofar_YX = bad_sofar_YX | excl_agmerra_YX==1 ; %#ok<NASGU>
        title('4: Yield too low: Phase 2 baseline sim')
        set(gca, 'FontSize', fontSize)
        nyellow = length(find(isyellow_YX)) ;
        nblue = length(find(map_YX==1)) ;
        ngreen = length(find(map_YX==2)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_y = sprintf('Yellow: %d', nyellow) ;
        text_g = sprintf('Green: %d', ngreen) ;
        text_gray = sprintf('Gray: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_gray, text_gray, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Add overall title
        hst = sgtitle( ...
            sprintf('Exclusions: %s %s %s', ggcm, thisCrop_title, thisFile), ...
            'FontWeight', 'bold') ;
        hst.FontSize = fontSize_big ;
        
        % Shift subplot positions
        h1.Position(2) = h1.Position(2) - 0.115 ;
        h2.Position(2) = h2.Position(2) - 0.115 ;
        h3.Position(2) = h3.Position(2) - 0.05 ;
        h4.Position(2) = h4.Position(2) - 0.05 ;
        
        % Add legend
        ia = axes('Position', [0 0 1 1], 'Visible', 'off') ;
        hpos = 0.13 ;
        topy = 0.93 ;
        down = 0.02 ;
        text(ia, hpos, topy, legend_text{1}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
        text(ia, hpos, topy-down, legend_text{2}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
        text(ia, hpos, topy-2*down, legend_text{3}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
        text(ia, hpos, topy-3*down, legend_text{4}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
        
    elseif strcmp(which_file,'gsirrigation')
        
        figure('Color', 'w', 'Position', [470    42   971   481], ...
            'Visible', figure_visibility) ;
        
        do_excl_check = c <= length(isexcl_all0_emu_c) ;
        if do_excl_check && isexcl_all0_emu_c(c)
            error('Add some kind of text in this figure saying there is no irrigation')
        else
        
            legend_text = { ...
                'Blue: Excluded because missing irrigation emulator', ...
                'Green: Excluded during yield exclusions', ...
                'Yellow: Included in yield and irrigation'} ;

            % Where there was yield, where is irrigation missing?
            missing_yield_YX = lpjgu_vector2map( ...
                missing_yield_xc(:,c) | isexcl_lowBLyield_xc(:,c), size(gridlist.mask_YX), ...
                gridlist.list_to_map) ;
            ok_sofar_YX = gridlist.mask_YX & missing_yield_YX==0 ;
            bad_sofar_YX = gridlist.mask_YX & missing_yield_YX==1 ;
            ngreen = length(find(bad_sofar_YX)) ;
            map_YX = nan(size(gridlist.mask_YX)) ;
            map_YX(ok_sofar_YX) = 2 ;
            map_YX(bad_sofar_YX) = 1 ;
            missing_emu_YX = lpjgu_vector2map( ...
                missing_emu_xc(:,c), size(gridlist.mask_YX), ...
                gridlist.list_to_map) ;
            isblue_YX = ok_sofar_YX & missing_emu_YX==1 ;
            map_YX(isblue_YX) = 0 ;
            pcolor(map_YX(yrange,:))
            shading flat; axis equal tight off
            colormap(gca, this_colormap)
            caxis([0 2])
            ok_sofar_YX = ok_sofar_YX & missing_emu_YX==0 ;
            nblue = length(find(isblue_YX)) ;
            text_b = sprintf('Blue: %d', nblue) ;
            text_g = sprintf('Green: %d', ngreen) ;
            text_y = sprintf('Yellow: %d', length(find(ok_sofar_YX))) ;
            text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
            text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
            text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)

            % Add overall title
            hst = sgtitle( ...
                sprintf('Exclusions: %s %s %s', ggcm, thisCrop_title, thisFile), ...
                'FontWeight', 'bold') ;
            hst.FontSize = fontSize_big ;

            % Add legend
            ia = axes('Position', [0 0 1 1], 'Visible', 'off') ;
            hpos = 0.13 ;
            topy = 0.91 ;
            down = 0.03 ;
            text(ia, hpos, topy, legend_text{1}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
            text(ia, hpos, topy-down, legend_text{2}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
            text(ia, hpos, topy-2*down, legend_text{3}, 'HorizontalAlignment', 'left', 'FontSize', fontSize)
        end
                
    else
        error('which_file (%s) not recognized', which_file)
    end
    
    if strcmp(figure_extension, 'png')
        export_fig(filename, '-r100', renderer)
    elseif strcmp(figure_extension, 'fig')
        savefig(filename)
    else
        error('figure_extension %s not recognized', figure_extension)
    end
    close
    
end


end