function do_map_run_diffs_fromEndHist(do_save, data_d9, sumvars, title_text, filename_base, ...
    equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, pngres, ...
    thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
    conv_fact_map, units_map, conv_fact_total, units_total, ...
    pct_clim, prctile_clim, R, gtif_missing, varargin)

already_maps = true ;
if ~isempty(varargin)
    if length(varargin)==2
        map_size = varargin{1} ;
        list2map = varargin{2} ;
        already_maps = false ;
    else
        error('do_map_run_diffs_fromEndHist() accepts 0 or 2 optional arguments: map_size and list2map')
    end
end

do_pct = false ;
if already_maps
    maps_YXr = map_run_diffs_fromEndHist(data_d9, title_text, sumvars, ...
        do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
        thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
        conv_fact_map, units_map, conv_fact_total, units_total, ...
        pct_clim, prctile_clim) ;
else
    maps_YXr = map_run_diffs_fromEndHist(data_d9, title_text, sumvars, ...
        do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
        thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
        conv_fact_map, units_map, conv_fact_total, units_total, ...
        pct_clim, prctile_clim, map_size, list2map) ;
end
if do_save
    actually_save(filename_base, do_pct, prctile_clim, pngres, ...
        maps_YXr, runList, R, gtif_missing)
end


do_pct = true ;
if already_maps
    maps_YXr = map_run_diffs_fromEndHist(data_d9, title_text, sumvars, ...
        do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
        thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
        conv_fact_map, units_map, conv_fact_total, units_total, ...
        pct_clim, prctile_clim) ;
else
    maps_YXr = map_run_diffs_fromEndHist(data_d9, title_text, sumvars, ...
        do_pct, equalize_cbars, fontSize, spacing, textX, textY_1, textY_2, ...
        thisPos, nx, ny, colorBarLoc, runList, do_caps, this_land_area_map, ...
        conv_fact_map, units_map, conv_fact_total, units_total, ...
        pct_clim, prctile_clim, map_size, list2map) ;
end
if do_save
    actually_save(filename_base, do_pct, prctile_clim, pngres, ...
        maps_YXr, runList, R, gtif_missing)
end


end


function actually_save(filename_base, do_pct, prctile_clim, pngres, ...
    maps_YXr, runList, R, gtif_missing)

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
save_geotiffs(maps_YXr, filename, runList, R, gtif_missing)
close

end
