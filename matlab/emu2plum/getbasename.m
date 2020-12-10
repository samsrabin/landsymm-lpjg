function out = getbasename(in)

mid = in ;

% Make sure actual crop names have no numbers in them
crop_names_w_nums_in = {'CerealsC3', 'CerealsC4'} ;
crop_names_w_nums_mid = {'CerealsCthree', 'CerealsCfour'} ;
if length(crop_names_w_nums_in) ~= length(crop_names_w_nums_mid)
    error('length(crop_names_w_nums_in) ~= length(crop_names_w_nums_mid)')
end
for c = 1:length(crop_names_w_nums_in)
    mid = strrep(mid, crop_names_w_nums_in{c}, crop_names_w_nums_mid{c}) ;
end

% Remove (i)NNN...
out = regexprep(mid, 'i?\d+$', '') ;

% Convert back to actual crop names
for c = 1:length(crop_names_w_nums_in)
    out = strrep(out, crop_names_w_nums_mid{c}, crop_names_w_nums_in{c}) ;
end

end