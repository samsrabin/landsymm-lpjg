function e2p_save_combineCrops_figs(combineCrops_in, gridlist, ...
    thisModel, figure_visibility, figure_extension, ...
    outDir_combineCrops_figs, overwrite_existing_figs, renderer)

Ncombines = length(combineCrops_in) ;

if ~exist(outDir_combineCrops_figs, 'dir')
    mkdir(outDir_combineCrops_figs)
end

for c = 1:Ncombines
    
    % Get info
    destCrop = combineCrops_in(c).destCrop ;
    sourceCrops = combineCrops_in(c).sourceCrops_cf ;
    Nsources = length(sourceCrops) ;
    whichmax_xvt = combineCrops_in(c).whichmax_xvt ;
    maxzero_xvt = combineCrops_in(c).maxzero_xvt ;
    varNames_sourceA = combineCrops_in(c).varNames_sourceA ;
    Ncells = size(whichmax_xvt, 1) ;
    Nvars = length(varNames_sourceA) ;
    Ntpers = size(whichmax_xvt, 3) ;
    
    whichmax_xvt(maxzero_xvt) = NaN ;
    
    choseThis_xvts = nan([size(whichmax_xvt) Nsources], 'single') ;
    for s = 1:length(sourceCrops)
        choseThis_xvt = single(whichmax_xvt == s) ;
        choseThis_xvt(maxzero_xvt) = NaN ;
        choseThis_xvts(:,:,:,s) = choseThis_xvt ;
    end
    choseThis_xsvt = permute(choseThis_xvts, [1 4 2 3]) ;
    clear choseThis_xvts
        
    % Overall map
    %    Would be simpler to just do mean(mean(choseThis_xsvt, 4), 3),
    %    but that's not robust to NaNs. Weighting would need to happen
    %    since we're averaging across multiple dimensions.
    choseThis_xsZ = reshape(choseThis_xsvt, [Ncells Nsources Nvars*Ntpers]) ;
    choseThis_xs = nanmean(choseThis_xsZ, 3) ;
    filename = sprintf('%s/combinedCrops_allTx_allYr_%s_%s.%s', ...
        outDir_combineCrops_figs, thisModel, destCrop, figure_extension) ;
    if overwrite_existing_figs || ~exist(filename, 'file')
        make_maps(choseThis_xs, sourceCrops, destCrop, ...
            gridlist, 'treatments and time periods', thisModel, ...
            figure_visibility, figure_extension, filename, ...
            renderer)
    end
    
    % Map for each N level
    
    % Map for each irrigation level
    
    % Map for each variable
    
end

end


function make_maps(choseThis_xs, sourceCrops, destCrop, ...
    gridlist, acrossAllWhat, thisModel, figure_visibility, ...
    figure_extension, filename, renderer)

% Figure options
spacing = [0.025 0.025] ;
fontSize = 14 ;

if ~ismatrix(choseThis_xs) || isvector(choseThis_xs)
    error('make_maps() requires 2-d choseThis_xs')
end
map_size = size(gridlist.mask_YX) ;
list2map = gridlist.list_to_map ;

Nsources = length(sourceCrops) ;
if Nsources == 2
    ny = 1 ;
    nx = 2 ;
    thisPos = [1 357 1440 448] ;
else
    error('Set up figure arrangement for Nsources==%d', ...
        Nsources)
end

figure('Color', 'w', 'Position', thisPos, ...
       'Visible', figure_visibility) ;

for s = 1:Nsources
    subplot_tight(ny, nx, s, spacing) ;
    
    map_YX = lpjgu_vector2map(choseThis_xs(:,s), ...
        map_size, list2map) ;
    pcolor(map_YX); shading flat; axis equal tight off
    
%     title(sprintf('Fraction %s was chosen for %s across all %s', ...
%         sourceCrops{s}, destCrop, acrossAllWhat))
    title(strrep(sourceCrops{s}, '_', '\_'))
    caxis([0 1]); colorbar
    set(gca, 'FontSize', fontSize)
end
sgtitle(sprintf('%s: Across all %s, fraction of %s from...', ...
    thisModel, acrossAllWhat, strrep(destCrop, '_', '\_')), ...
    'FontSize', fontSize + 4, ...
    'FontWeight', 'bold') ;

if strcmp(figure_extension, 'png')
    export_fig(filename, '-r100', renderer)
elseif strcmp(figure_extension, 'fig')
    savefig(filename)
else
    error('figure_extension %s not recognized', figure_extension)
end
close


end