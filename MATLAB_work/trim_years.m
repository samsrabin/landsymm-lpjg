function struct_in = trim_years(struct_in, yearList)

Nyears = length(yearList) ;
[~,IA] = intersect(struct_in.yearList,yearList) ;
if length(IA) ~= Nyears
    error('length(IA) ~= Nyears')
end
struct_in.maps_YXvy = struct_in.maps_YXvy(:,:,:,IA) ;
struct_in.yearList = yearList ;

end