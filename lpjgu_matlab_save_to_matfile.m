function lpjgu_matlab_save_to_matfile(out_struct,in_matfile,force_mat_save,verboseIfNoMat,verbose)

if ~isstruct(out_struct)
    error('out_struct must be a struct.')
end
% Ask to save to MAT-file
if force_mat_save
    if verboseIfNoMat || verbose
        disp('         Saving MAT-file...')
    end
    save(in_matfile,'out_struct','-v7.3') ;
else
    ok = false ;
    while ~ok
        if exist(in_matfile,'file')
            disp('      Save, overwriting existing MAT-file? Y or [N]. 10 seconds...')
            default_save = false ;
        else
            disp('      Save to MAT-file? [Y] or N. 10 seconds...')
            default_save = true ;
        end
        dbl = getkeywait_ssr(10) ;
        if (dbl==-1 && default_save) || strcmp(char(dbl),'y') || strcmp(char(dbl),'Y')
            ok = true ;
            disp('         Saving MAT-file...')
            save(in_matfile,'out_struct','-v7.3') ;
        elseif (dbl==-1 && ~default_save) || strcmp(char(dbl),'n') || strcmp(char(dbl),'N')
            ok = true ;
        else
            warning(['Input (' char(dbl) ') not recognized.'])
        end
    end ; clear ok
end


end