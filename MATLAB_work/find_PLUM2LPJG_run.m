function dir_out = find_PLUM2LPJG_run(dir_in)

if exist(dir_in,'dir')
    dir_out = addslashifneeded(dir_in) ;
else
    topDir = addslashifneeded('/Users/Shared/PLUM/trunk_runs') ;
    dir_out = addslashifneeded([topDir dir_in]) ;
    if ~exist(dir_out,'dir')
        dir_out = '' ;
        externalDirs = dir('/Volumes') ;
        externalDirs = {externalDirs.name} ;
        externalDirs(contains(externalDirs,{'.','..','Macintosh HD','MobileBackups'})) = [] ;
        for d = 1:length(externalDirs)
            thisDir = addslashifneeded(['/Volumes/' externalDirs{d}]) ;
            topDir = addslashifneeded([thisDir 'PLUM/trunk_runs/' dir_in]) ;
            if exist(topDir,'dir')
                dir_out = topDir ;
                break
            end
        end
        
        if isempty(dir_out)
            error(['Input directory not found! (' dir_in ')'])
        end
    end
end

% If needed, get latest output directory
outdirs = dir([dir_out 'output*']) ;
if ~isempty(outdirs)
    dir_out = addslashifneeded(sprintf('%s/%s', outdirs(end).folder, outdirs(end).name)) ;
end


end
