function do_map_run_diffs_fromEndHist(do_save, maps_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim)

do_pct = false ;
map_run_diffs_fromEndHist(maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim) ;
if do_save
    filename = [filename_base '_diff'] ;
    if do_pct 
        filename = [filename 'Pct'] ;
    elseif ~isempty(prctile_clim)
        if isint(prctile_clim)
            filename = sprintf('%s_limPrctile%d', filename, prctile_clim) ;
        else
            filename = sprintf('%s_limPrctile%0.1f', filename, prctile_clim) ;
        end
    end
    filename = [filename '_2000s-2090s.png'] ;
    export_fig(filename,['-r' num2str(pngres)])
    close
end

do_pct = true ;
map_run_diffs_fromEndHist(maps_d9, title_text, sumvars, ...
    do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim) ;
if do_save
    filename = [filename_base '_diff'] ;
    if do_pct 
        filename = [filename 'Pct'] ;
    elseif ~isempty(prctile_clim)
        if isint(prctile_clim)
            filename = sprintf('%s_limPrctile%df', filename, prctile_clim) ;
        else
            filename = sprintf('%s_limPrctile%0.1f', filename, prctile_clim) ;
        end
    end
    filename = [filename '_2000s-2090s.png'] ;
    export_fig(filename,['-r' num2str(pngres)])
    close
end


end