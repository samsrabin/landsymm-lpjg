function [thisDir, yN] = e2p_get_thisDir(topDir, y1, end_yN)

% Get this directory
thisDir_name = sprintf('%s/*%d-*', topDir, y1) ;
thisDir = dir(thisDir_name) ;
% disp(thisDir)
if length(thisDir) ~= 1
    error('Error finding directory: %d found (%s)', length(thisDir), thisDir_name)
end
thisDir_years = strsplit(thisDir.name,'-') ;

% Check that you're not outside the bounds of the specified period
yN = str2double(thisDir_years{2}) ;
if yN > end_yN
    error('yN (%d) > end_yN (%d)', yN, end_yN) ;
end
% fprintf('%d-%d\n', y1, yN)

end