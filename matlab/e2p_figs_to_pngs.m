%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert MATLAB fig files to PNG %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topdir = '/Users/Shared/GGCMI2PLUM_sh/CMIP_emulated_work' ;
var_list = {'yield', 'irrig'} ;

dir_list = dir(sprintf('%s/*ssp*', topdir)) ;
for d = 1:length(dir_list)
    thisDir_top = sprintf('%s/%s', topdir, dir_list(d).name) ;
    fprintf('%s (%d of %d)...\n', dir_list(d).name, d, length(dir_list)) ;
    for v = 1:length(var_list)
        thisVar = var_list{v} ;
        fprintf('    %s (%d of %d)...\n', thisVar, v, length(var_list)) ;
        thisDir = sprintf('%s/%s_figs', thisDir_top, thisVar) ;
        file_list = dir(thisDir) ;
        for f = 1:length(file_list)
            fprintf('        %s (%d of %d)...\n', file_list(f).name, f, length(file_list)) ;
            thisFile_in = sprintf('%s/%s_figs', thisDir, file_list(f).name) ;
            thisFile_out = strrep(thisFile_in, '.fig', '.png') ;
            openfig(thisFile_in, 'visible') ;
            export_fig(thisFile_out, '-r100')
            close
        end
    end
end
disp('Done')

%%
file_list = dir(sprintf('%s/**/*.fig', topdir)) ;
for f = 63:length(file_list)
    fprintf('        %s (%d of %d)...\n', file_list(f).name, f, length(file_list)) ;
    thisFile_in = sprintf('%s/%s', file_list(f).folder, file_list(f).name) ;
    thisFile_out = strrep(thisFile_in, '.fig', '.png') ;
    openfig(thisFile_in, 'visible') ;
    export_fig(thisFile_out, '-r100')
    close
end

disp('Done')