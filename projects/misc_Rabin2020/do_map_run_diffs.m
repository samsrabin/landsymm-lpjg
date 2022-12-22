function do_map_run_diffs(do_save, maps_d1, maps_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim)

do_pct = false ;
map_run_diffs(maps_d1, maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;
if do_save
    filename = [filename_base '_diff'] ;
    if do_pct ; filename = [filename 'Pct'] ; end
    export_fig([filename '_2010s-2090s.png'],['-r' num2str(pngres)])
    close
end

do_pct = true ;
map_run_diffs(maps_d1, maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, pct_clim) ;
if do_save
    filename = [filename_base '_diff'] ;
    if do_pct ; filename = [filename 'Pct'] ; end
    export_fig([filename '_2010s-2090s.png'],['-r' num2str(pngres)])
    close
end


end