function save_last_decade(struct_in, struct_name, filename)

is_matfile = @(x) ischar(x) & length(x)>=5 & strcmp(x(end-3:end),'.mat') ;

p = inputParser ;
addRequired(p,'struct_in',@isstruct) ;
addRequired(p,'struct_name',@ischar) ;
addRequired(p,'filename',is_matfile) ;
parse(p, struct_in, struct_name, filename);

struct_in.maps_YXvy = struct_in.maps_YXvy(:,:,:,end-9:end) ;
struct_in.yearList = struct_in.yearList(end-9:end) ;
if ~exist(filename,'file')
    save(filename,struct_name) ;
else
    save(filename,struct_name,'-append') ;
end

end