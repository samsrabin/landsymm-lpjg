function dirList = PLUMharm_check_dirs(dirList, rw)

Ndirs = length(dirList) ;

for d = 1:Ndirs
    thisDir = dirList{d} ;
    thisDir = removeslashifneeded(thisDir) ;

    if ~exist(thisDir, 'dir')
        if ~contains(rw, 'w')
            error('%s not found', ...
                thisDir)
        end
        try
            mkdir(thisDir)
        catch ME
            error('Trying to make %s:\n%s', ...
                thisDir, ME.message)
        end
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
    dirList{d} = thisDir ;
end

end
