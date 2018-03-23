function [strip_fao_nans, fix_cotedivoire, combine_sudans,...
    combine_subChinas, combine_subChinas_map, combine_serbmont] ...
    = get_FAOread_options(calib_ver)

latest_calib_ver = 17 ;

% Strip NaNs from FAO data?
if calib_ver>=5 && calib_ver<=latest_calib_ver
    strip_fao_nans = true ;
elseif calib_ver<=4
    strip_fao_nans = false ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Fix Cote d'Ivoire?
if calib_ver>=5 && calib_ver<=latest_calib_ver
    fix_cotedivoire = true ;
elseif calib_ver<=4
    fix_cotedivoire = false ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Combine Sudan & South Sudan? (ONLY IF YEARN<2011)
if calib_ver>=5 && calib_ver<=latest_calib_ver
    combine_sudans = true ;
elseif calib_ver<=4
    combine_sudans = false ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Combine sub-Chinas in FAO data?
if calib_ver>=5 && calib_ver<=latest_calib_ver
    combine_subChinas = true ;
elseif calib_ver<=4
    combine_subChinas = false ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Combine sub-Chinas in map?
if calib_ver<=8
    combine_subChinas_map = false ;
elseif calib_ver<=latest_calib_ver
    combine_subChinas_map = true ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Combine Serbia & Montenegro?
if calib_ver>=5 && calib_ver<=latest_calib_ver
    combine_serbmont = true ;
elseif calib_ver<=4
    combine_serbmont = false ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end


end