function e2p_save_excl_figs_outCrops( ...
    ggcm, which_file, gridlist, excl_vecs, ...
    cropList_lpj, cropList_lpj_basei, ...
    cropList_emu, cropList_emu_basei, ...
    figure_visibility, figure_extension, ...
    outDir_excl_figs_outCrops)

% Unpack
missing_emu_xc = excl_vecs{1} ;
if strcmp(which_file, 'yield')
    missing_agmerra_xc = excl_vecs{2} ;
    isexcl_lowBL_emu_xc = excl_vecs{3} ;
    isexcl_lowBL_agmerra_xc = excl_vecs{4} ;
elseif strcmp(which_file,'gsirrigation')
    missing_yield_xc = excl_vecs{2} ;
    isexcl_lowBLyield_xc = excl_vecs{3} ;
else
    error('which_file (%s) not recognized', which_file)
end

% Add max wheats (rainfed)
wheat_inds = [find(strcmp(cropList_emu_basei, 'winter_wheat')) ...
    find(strcmp(cropList_emu_basei, 'spring_wheat'))] ;
if length(wheat_inds) ~= 2
    error('Error finding wheat_inds: %d found', length(wheat_inds))
end
missing_emu_xc(:,end+1) = ...
    missing_emu_xc(:,wheat_inds(1)) ...
    & missing_emu_xc(:,wheat_inds(2)) ;
if strcmp(which_file, 'yield')
    missing_agmerra_xc(:,end+1) = ...
        missing_agmerra_xc(:,wheat_inds(1)) ...
        & missing_agmerra_xc(:,wheat_inds(2)) ;
    isexcl_lowBL_emu_xc(:,end+1) = ...
        isexcl_lowBL_emu_xc(:,wheat_inds(1)) ...
        & isexcl_lowBL_emu_xc(:,wheat_inds(2)) ;
    isexcl_lowBL_agmerra_xc(:,end+1) = ...
        isexcl_lowBL_agmerra_xc(:,wheat_inds(1)) ...
        & isexcl_lowBL_agmerra_xc(:,wheat_inds(2)) ;
elseif strcmp(which_file,'gsirrigation')
    missing_yield_xc(:,end+1) = ...
        missing_yield_xc(:,wheat_inds(1)) ...
        & missing_yield_xc(:,wheat_inds(2)) ;
    isexcl_lowBLyield_xc(:,end+1) = ...
        isexcl_lowBLyield_xc(:,wheat_inds(1)) ...
        & isexcl_lowBLyield_xc(:,wheat_inds(2)) ;
else
    error('which_file (%s) not recognized', which_file)
end
cropList_emu_basei{end+1} = 'max_wheat' ;

% Add max wheats (irrigated)
wheat_inds = [find(strcmp(cropList_emu_basei, 'winter_wheati')) ...
    find(strcmp(cropList_emu_basei, 'spring_wheati'))] ;
if length(wheat_inds) ~= 2
    error('Error finding wheat_inds: %d found', wheat_inds)
end
missing_emu_xc(:,end+1) = ...
    missing_emu_xc(:,wheat_inds(1)) ...
    & missing_emu_xc(:,wheat_inds(2)) ;
if strcmp(which_file, 'yield')
    missing_agmerra_xc(:,end+1) = ...
        missing_agmerra_xc(:,wheat_inds(1)) ...
        & missing_agmerra_xc(:,wheat_inds(2)) ;
    isexcl_lowBL_emu_xc(:,end+1) = ...
        isexcl_lowBL_emu_xc(:,wheat_inds(1)) ...
        & isexcl_lowBL_emu_xc(:,wheat_inds(2)) ;
    isexcl_lowBL_agmerra_xc(:,end+1) = ...
        isexcl_lowBL_agmerra_xc(:,wheat_inds(1)) ...
        & isexcl_lowBL_agmerra_xc(:,wheat_inds(2)) ;
elseif strcmp(which_file,'gsirrigation')
    missing_yield_xc(:,end+1) = ...
        missing_yield_xc(:,wheat_inds(1)) ...
        & missing_yield_xc(:,wheat_inds(2)) ;
    isexcl_lowBLyield_xc(:,end+1) = ...
        isexcl_lowBLyield_xc(:,wheat_inds(1)) ...
        & isexcl_lowBLyield_xc(:,wheat_inds(2)) ;
else
    error('which_file (%s) not recognized', which_file)
end
cropList_emu_basei{end+1} = 'max_wheati' ;
cropList_emu{end+1} = 'max_wheat' ;

% Set up for figures
if ~exist(outDir_excl_figs_outCrops, 'dir')
    mkdir(outDir_excl_figs_outCrops) ;
end
spacing = [0.025 0.025] ; % v h
yrange = 65:360 ;
fontSize = 14 ;
this_colormap = 'parula' ;
x_bgy = 0 ;
y_bgy_top = 0.3 ;
down_bgy = 0.055 ;
y_g = y_bgy_top ;
y_b = y_g - down_bgy ;
y_y = y_b - down_bgy ;
fontSize_bgy = fontSize - 2 ;
fontSize_big = fontSize + 4 ;

% Translate crop names
verbose = false ;
cropList_lpj_asEmu = e2p_translate_crops( ...
    cropList_lpj, cropList_emu, verbose) ;

% Set up for filename
if strcmp(which_file,'yield')
    thisFile = 'yield' ;
elseif strcmp(which_file,'gsirrigation')
    thisFile = 'irrig' ;
else
    error('which_file (%s) not recognized', which_file)
end

for ci_lpj = 1:length(cropList_lpj_basei)
    thisCropi = cropList_lpj_basei{ci_lpj} ;
    thisCrop = thisCropi ;
    isIrr = strcmp(thisCropi(end), 'i') ;
    if isIrr
        thisCrop = thisCrop(1:end-1) ;
    end
    thisCropi_title = strrep(thisCropi, '_', '\_') ;
    
    % If irrigation, skip rainfed crops
    if strcmp(which_file,'gsirrigation') && ~isIrr
        continue
    end
    
    % Get index of equivalent GGCMI crop
    c_lpj = find(strcmp(cropList_lpj,thisCrop)) ;
    if length(c_lpj) ~= 1
        error('Error finding c_lpj: %d found', length(c_lpj))
    end
    thisCrop_emu = cropList_lpj_asEmu{c_lpj} ;
    thisCropi_emu = thisCrop_emu ;
    if isIrr
        thisCropi_emu = [thisCrop_emu 'i'] ;
    end
    c = find(strcmp(cropList_emu_basei, thisCropi_emu)) ;
    if length(c) ~= 1
        error('Error finding index of thisCropi_emu %s', thisCropi_emu)
    end
%     keyboard
    
    filename = sprintf('%s/exclusions_%s_%s_%s.png', ...
        outDir_excl_figs_outCrops, ggcm, thisCropi, thisFile) ;
    if strcmp(figure_extension, 'fig')
        filename = strrep(filename, '.png', '.fig') ;
    elseif ~strcmp(figure_extension, 'png')
        error('figure_extension %s not recognized', figure_extension)
    end
%     if exist(filename, 'file')
%         continue
%     end
    
    ngreen = 0 ;
    
    if strcmp(which_file,'yield')
        
        figure('Color', 'w', 'Position', figurePos, ...
            'Visible', figure_visibility) ;
        legend_text = { ...
            'Blue: Excluded in this step', ...
            'Green: Excluded in a previous step', ...
            'Yellow: Still included after this step'} ;
        
        % Where is emulator yield missing?
        h1 = subplot_tight(2,2,1,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        missing_emu_YX = lpjgu_vector2map( ...
            missing_emu_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isblue_YX = gridlist.mask_YX & missing_emu_YX==1 ;
        map_YX(gridlist.mask_YX & missing_emu_YX==0) = 2 ;
        map_YX(isblue_YX) = 0 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 2])
        ok_sofar_YX = gridlist.mask_YX & missing_emu_YX==0 ;
        bad_sofar_YX = missing_emu_YX==1 ;
        title('1: Missing yield: Emulator')
        set(gca, 'FontSize', fontSize)
        nblue = length(find(isblue_YX)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_g = sprintf('Green: %d', ngreen) ;
        ngreen = ngreen + nblue ;
        text_y = sprintf('Yellow: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Where there is emulator yield, where is AgMERRA yield missing?
        h2 = subplot_tight(2,2,2,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        map_YX(ok_sofar_YX) = 2 ;
        map_YX(bad_sofar_YX) = 1 ;
        missing_agmerra_YX = lpjgu_vector2map( ...
            missing_agmerra_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isblue_YX = ok_sofar_YX & missing_agmerra_YX==1 ;
        map_YX(isblue_YX) = 0 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 2])
        ok_sofar_YX = ok_sofar_YX & missing_agmerra_YX==0 ;
        bad_sofar_YX = bad_sofar_YX | missing_agmerra_YX==1 ;
        title('2: Missing yield: Phase 2 baseline sim')
        set(gca, 'FontSize', fontSize)
        nblue = length(find(isblue_YX)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_g = sprintf('Green: %d', ngreen) ;
        ngreen = ngreen + nblue ;
        text_y = sprintf('Yellow: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Where there is emulator and AgMERRA yield, where is there
        % exclusion based on too-low yield in emulator?
        h3 = subplot_tight(2,2,3,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        map_YX(ok_sofar_YX) = 2 ;
        map_YX(bad_sofar_YX) = 1 ;
        excl_emu_YX = lpjgu_vector2map( ...
            isexcl_lowBL_emu_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isblue_YX = ok_sofar_YX & excl_emu_YX==1 ;
        map_YX(isblue_YX) = 0 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 2])
        ok_sofar_YX = ok_sofar_YX & excl_emu_YX==0 ;
        bad_sofar_YX = bad_sofar_YX | excl_emu_YX==1 ;
        title('3: Yield too low: Emulator')
        set(gca, 'FontSize', fontSize)
        nblue = length(find(isblue_YX)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_g = sprintf('Green: %d', ngreen) ;
        ngreen = ngreen + nblue ;
        text_y = sprintf('Yellow: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Where there is emulator and AgMERRA yield, and emulator yield is
        % not too low, where is there exclusion based on too-low yield in
        % AgMERRA?
        h4 = subplot_tight(2,2,4,spacing) ;
        map_YX = nan(size(gridlist.mask_YX)) ;
        map_YX(ok_sofar_YX) = 2 ;
        map_YX(bad_sofar_YX) = 1 ;
        excl_agmerra_YX = lpjgu_vector2map( ...
            isexcl_lowBL_agmerra_xc(:,c), size(gridlist.mask_YX), ...
            gridlist.list_to_map) ;
        isblue_YX = ok_sofar_YX & excl_agmerra_YX==1 ;
        map_YX(isblue_YX) = 0 ;
        pcolor(map_YX(yrange,:))
        shading flat; axis equal tight off
        colormap(gca, this_colormap)
        caxis([0 2])
        ok_sofar_YX = ok_sofar_YX & excl_agmerra_YX==0 ;
        bad_sofar_YX = bad_sofar_YX | excl_agmerra_YX==1 ; %#ok<NASGU>
        title('4: Yield too low: Phase 2 baseline sim')
        set(gca, 'FontSize', fontSize)
        nblue = length(find(isblue_YX)) ;
        text_b = sprintf('Blue: %d', nblue) ;
        text_g = sprintf('Green: %d', ngreen) ;
        ngreen = ngreen + nblue ; %#ok<NASGU>
        text_y = sprintf('Yellow: %d', length(find(ok_sofar_YX))) ;
        text(x_bgy, y_b, text_b, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_g, text_g, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        text(x_bgy, y_y, text_y, 'Units', 'normalized', 'FontSize', fontSize_bgy)
        
        % Add overall title
        hst = sgtitle( ...
            sprintf('Exclusions: %s %s %s', ggcm, thisCropi_title, thisFile), ...
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
        
    elseif strcmp(which_file,'gsirrigation')
        
        figure('Color', 'w', 'Position', [470    42   971   481], ...
            'Visible', figure_visibility) ;
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
            sprintf('Exclusions: %s %s %s', ggcm, thisCropi_title, thisFile), ...
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
                
    else
        error('which_file (%s) not recognized', which_file)
    end
    
    if strcmp(figure_extension, 'png')
        export_fig(filename, '-r100')
    elseif strcmp(figure_extension, 'fig')
        savefig(filename)
    else
        error('figure_extension %s not recognized', figure_extension)
    end
    close
    
end


end