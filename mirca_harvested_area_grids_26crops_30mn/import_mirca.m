function out_YX = import_mirca(filename)

out_YX = dlmread(filename,' ',6,0) ;
out_YX = flipud(out_YX) ;
out_YX(out_YX==-9999) = NaN ;
out_YX = out_YX(:,1:720) ;   % Extra column because files have spaces at the end of lines

% Sanity checks
if ~isempty(find(out_YX<0,1))
    warning('Map has negative value(s).')
end

end