function out = getN_num(x) 

mid = getN_char(x) ;
out = str2double(mid) ;
if any(isnan(out))
    I = find(isnan(out), 1) ;
    error('Error translating %s to double', mid(I))
end

end