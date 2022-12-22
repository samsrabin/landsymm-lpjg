function dir_out = find_PLUM2LPJG_inputs(dir_in,varargin)

error_if_not_found = true ;
if ~isempty(varargin)
    error_if_not_found = varargin{1} ;
end

topDir = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;
dir_out = addslashifneeded([topDir dir_in]) ;
if ~exist(dir_out,'dir')
% %     topDir = addslashifneeded('/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG') ;
%     topDir = addslashifneeded('/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG') ;
%     dir_out = addslashifneeded([topDir dir_in]) ;
%     if ~exist(dir_out,'dir')
%         if error_if_not_found
%             error(['Input directory not found! (' dir_in ')'])
%         else
%             dir_out = [] ;
%         end
%     end
    dir_out = '' ;
    externalDirs = dir('/Volumes') ;
    externalDirs = {externalDirs.name} ;
    externalDirs(contains(externalDirs,{'.','..','Macintosh HD','MobileBackups'})) = [] ;
    for d = 1:length(externalDirs)
        thisDir = addslashifneeded(['/Volumes/' externalDirs{d}]) ;
        topDir = addslashifneeded([thisDir 'PLUM/PLUM_outputs_for_LPJG/' dir_in]) ;
        if exist(topDir,'dir')
            dir_out = topDir ;
            break
        end
    end
    
    if isempty(dir_out)
        testDir = sprintf('/Volumes/Reacher/LandSyMM/inputs/LU/%s/', dir_in) ;
        if exist(testDir, 'dir')
            dir_out = testDir ;
        end
    end

    if isempty(dir_out)
        error(['Input directory not found! (' dir_in ')'])
    end

end


end