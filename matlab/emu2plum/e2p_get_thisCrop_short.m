function thisCrop_short = e2p_get_thisCrop_short(thisCrop)

switch thisCrop
    case 'spring_wheat'
        thisCrop_short = 'swh' ;
    case 'winter_wheat'
        thisCrop_short = 'wwh' ;
    case 'maize'
        thisCrop_short = 'mai' ;
    case 'soy'
        thisCrop_short = 'soy' ;
    case 'rice'
        thisCrop_short = 'ric' ;
    otherwise
        error('thisCrop (%s) not recognized', thisCrop)
end

end