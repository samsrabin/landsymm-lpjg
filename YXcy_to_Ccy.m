function [croparea_out_Ccy,total_out_Ccy,yield_out_Ccy] = ...
    YXcy_to_Ccy(total_in_YXcy,croparea_in_YXcy,...
                countries_YX,countries_key,listCountries_map_present,...
                verbose)

% Process empties
do_croparea = true ;
if isempty(croparea_in_YXcy)
    do_croparea = false ;
    croparea_out_Ccy = [] ;
    yield_out_Ccy = [] ;
end

% Set up empty _Ccy arrays
Ncountries = length(listCountries_map_present) ;
Ncrops = size(total_in_YXcy,3) ;
Nyears = size(total_in_YXcy,4) ;
total_out_Ccy = nan(Ncountries,Ncrops,Nyears) ;
if do_croparea
    croparea_out_Ccy = nan(Ncountries,Ncrops,Nyears) ;
end

% Fill in _Ccy arrays
for c = 1:Ncountries
    thisCountry_name = listCountries_map_present{c} ;
    thisCountry = countries_key.numCode(strcmp(countries_key.Country,thisCountry_name)) ;
    inThisCountry = countries_YX==thisCountry ;
    if any(inThisCountry(:))
        if verbose ; disp([thisCountry_name '...']) ; end
        N_inThisCountry = length(find(inThisCountry)) ;
        
        % Total production
        tmp_x = total_in_YXcy(repmat(inThisCountry,[1 1 Ncrops Nyears])) ;
        tmp_xvy = reshape(tmp_x,[N_inThisCountry Ncrops Nyears]) ;
        total_out_Ccy(c,:,:) = sum(tmp_xvy,1) ;
        clear tmp*
        
        % Crop area
        if do_croparea
            tmp_x = croparea_in_YXcy(repmat(inThisCountry,[1 1 Ncrops Nyears])) ;
            tmp_xvy = reshape(tmp_x,[N_inThisCountry Ncrops Nyears]) ;
            croparea_out_Ccy(c,:,:) = sum(tmp_xvy,1) ;
            clear tmp*
        end
    elseif verbose
        warning(['No cells found in ' thisCountry_name '.'])
    end
end ; clear c
if verbose ; disp('Done.') ; end

% Calculate average per-country yield (tDM/ha)
if do_croparea
    yield_out_Ccy = total_out_Ccy ./ croparea_out_Ccy ;
end



end