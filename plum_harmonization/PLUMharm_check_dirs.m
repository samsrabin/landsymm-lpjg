function dirList = PLUMharm_check_dirs(dirList, rw)

Ndirs = length(dirList) ;

for d = 1:Ndirs
    thisDir = dirList{d} ;

    if ~exist(thisDir, 'dir')
        if ~contains(rw, 'w')
            error('%s not found', ...
                thisDir)
        end
        mkdir(thisDir)
    end

    [~, fa] = fileattrib(thisDir) ;

    % Get absolute path
    thisDir = fa.Name ;

    % Check
    if contains(rw, 'r') && ~fa.UserRead
        error('%s is not readable!', thisDir)
    end
    if contains(rw, 'w') && ~fa.UserWrite
        error('%s is not writeable!', thisDir)
    end
    if ~fa.directory
        error('%s is not a directory!', thisDir)
    end

    % Finish up
    thisDir = addslashifneeded(thisDir) ;
    dirList{d} = thisDir ;
end

end