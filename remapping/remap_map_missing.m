function isbad_YX = remap_map_missing(garr, gridlist, thisTitle, fig_dir)

isbad_x = isnan(garr) ;
while length(find(size(isbad_x)>1)) > 1
    isbad_x = squeeze(any(isbad_x, 2)) ;
end

isbad_YX = lpjgu_vector2map(isbad_x, size(gridlist.mask_YX), gridlist.list2map) ;

fig_outfile = fullfile(fig_dir, ...
    sprintf('cells_missing_from_%s.png', thisTitle)) ;
Nmissing = length(find(isbad_x)) ;
thisTitle = sprintf('%d cells missing from %s', ...
    Nmissing, thisTitle) ;
remap_shademap(isbad_YX, thisTitle, fig_outfile) ;

end
