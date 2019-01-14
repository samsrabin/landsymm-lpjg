function file_out = find_PLUM2LPJG_input_file(file_in,varargin)

error_if_not_found = true ;
if ~isempty(varargin)
    error_if_not_found = varargin{1} ;
end

topDir = addslashifneeded('/Users/Shared/PLUM/PLUM_outputs_for_LPJG') ;
file_out = [topDir file_in] ;
% if ~exist(file_out,'file') && ~exist([file_out '.gz'],'file')
%     topDir = addslashifneeded('/Volumes/WDMPP_Storage/Shared/PLUM/PLUM_outputs_for_LPJG') ;
%     file_out = [topDir file_in] ;
%     if ~exist(file_out,'file') && ~exist([file_out '.gz'],'file')
%         if error_if_not_found
%             error(['Input file not found! (' file_in ')'])
%         else
%             file_out = [] ;
%         end
%     end
% end
if ~exist(file_out,'file') && ~exist([file_out '.gz'],'file')
    file_out = '' ;
    externalDirs = dir('/Volumes') ;
    externalDirs = {externalDirs.name} ;
    externalDirs(contains(externalDirs,{'.','..','Macintosh HD','MobileBackups'})) = [] ;
    for d = 1:length(externalDirs)
        topDir = addslashifneeded(['/Volumes/' externalDirs{d} '/PLUM/PLUM_outputs_for_LPJG/']) ;
        file_try = [topDir file_in] ;
        if exist(file_try,'file') || exist([file_try '.gz'],'file')
            file_out = file_try ;
            break
        end
    end
    
    if isempty(file_out)
        error(['Input file not found! (' file_in ')'])
    end
end


end