function title_out = get_commods_title(name_in)

title_out = name_in ;
title_out(1) = upper(title_out(1)) ;
if strcmp(title_out,'StarchyRoots')
    title_out = 'Starchy roots' ;
elseif strcmp(title_out,'Crops')
    title_out = 'Total crops' ;
elseif strcmp(title_out,'Livestock')
    title_out = 'Total livestock' ;
end