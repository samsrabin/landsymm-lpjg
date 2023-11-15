function remap_interp_soil(files_soil, gridlist, out_dir, out_dir_figs, ...
    inpaint_method, remapVer, outWidth, delimiter, overwrite, fancy)
% Interpolate and save soil

Ncells = length(gridlist.list2map) ;

for f = 1:length(files_soil)
    file_soil = files_soil{f} ;
    [~, soil_filename, soil_ext] = fileparts(file_soil) ;
    disp([soil_filename soil_ext ':'])
 
    soil = lpjgu_matlab_read2geoArray(file_soil, ...
        'verboseIfNoMat', false, 'force_mat_save', false, 'force_mat_nosave', true, ...
        'target', gridlist, 'fill_missing_cells_with_nan', true) ;

    remap_map_missing(soil.garr_xv, gridlist, ['soil_' soil_filename], out_dir_figs) ;
    out_soil.varNames = soil.varNames ;
    out_soil.lonlats = gridlist.lonlats ;
    out_soil.list2map = gridlist.list2map ;
    Nvars = length(soil.varNames) ;
    out_soil.garr_xv = nan(Ncells, Nvars) ;
    for v = 1:Nvars
        fprintf('    Interpolating soil %s...\n', soil.varNames{v})
        tmp0_YX = lpjgu_xz_to_YXz(soil.garr_xv(:,v), size(gridlist.mask_YX), soil.list2map) ;
        tmp1_YX = inpaint_nans(tmp0_YX, inpaint_method) ;
        out_soil.garr_xv(:,v) = tmp1_YX(gridlist.list2map) ;
    end
    
    if any(any(isnan(out_soil.garr_xv)))
        error('NaN in out_soil.garr_xv')
    elseif any(any(out_soil.garr_xv < 0))
        error('Negative in out_soil.garr_xv')
    end
    
    out_soil_header_cell = [ ...
        {'lon', 'lat'}, ...
        out_soil.varNames] ;
    
    out_file_soil = fullfile(out_dir, sprintf('%s.remapv%s%s', ...
        soil_filename, remapVer, soil_ext)) ;
    
    disp('    Saving...')
    lpjgu_matlab_saveTable(out_soil_header_cell, out_soil, out_file_soil,...
        'outPrec', 3, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', 20, ...
        'verbose', false, ...
        'gzip', false) ;
end

disp('Done with soil.')

end