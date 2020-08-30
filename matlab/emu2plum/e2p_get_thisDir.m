function [thisDir, yN] = e2p_get_thisDir(topDir, y1, end_yN, which_file)

% Get this directory
thisDir_name = sprintf('%s/*%d-*', topDir, y1) ;
thisDir = dir(thisDir_name) ;

% If multiple found, eliminate any that don't have a which_file* file
if length(thisDir) > 1
	test = false(length(thisDir),1) ;
	for d = 1:length(thisDir)
		thisFile = sprintf('%s/%s/%s*',thisDir(d).folder,thisDir(d).name, which_file) ;
		if ~isempty(dir(thisFile))
			test(d) = true ;
		end
	end
	thisDir(~test) = [] ;
end

% If you still have multiple, throw an error
if length(thisDir) ~= 1
    error('Error finding directory: %d found (%s)', length(thisDir), thisDir_name)
end

% Check that you're not outside the bounds of the specified period
thisDir_years = strsplit(thisDir.name,'-') ;
yN = str2double(thisDir_years{2}) ;
if yN > end_yN
    error('yN (%d) > end_yN (%d)', yN, end_yN) ;
end

end
