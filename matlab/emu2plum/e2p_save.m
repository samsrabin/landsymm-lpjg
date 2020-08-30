function e2p_save(outDir, y1, yN, out_header_cell, lonlats, garr_xv, which_file, ...
    interp_infs, remove_outliers, overwrite)

outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
fancy = false ;

outDir_thisT = sprintf('%s/%d-%d', outDir, y1, yN) ;
if ~exist(outDir_thisT, 'dir')
    s = unix(sprintf('mkdir -p %s', outDir_thisT)) ;
    if s~=0
        error('Error in mkdir -p')
    end
end

% Get output file name
outfile = sprintf('%s/%s', outDir_thisT, which_file) ;
if interp_infs
    outfile = [outfile '_intpinfs'] ;
end
if remove_outliers
    outfile = [outfile '_rmol'] ;
end
outfile = [outfile '.out'] ;

% Save file
out_array = cat(2, lonlats, garr_xv) ;
lpjgu_matlab_saveTable(out_header_cell, out_array, outfile,...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy, ...
    'progress_step_pct', 20, ...
    'verbose', false) ;
clear out_array

% Compress file
out_gz = [outfile '.gz'] ;
if exist(out_gz, 'file')
    s = unix(sprintf('rm %s', out_gz)) ;
    if s~=0
        error('Error using rm')
    end
end
s = unix(sprintf('gzip %s', outfile)) ;
if s~=0
    error('Error using gzip')
end


end
