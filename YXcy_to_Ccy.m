function [croparea_out_Ccy,cropareaR_out_Ccy,R_out_Ccy,total_out_Ccy,yield_out_Ccy] = ...
    YXcy_to_Ccy(total_in_YXcy,croparea_in_YXcy,...
                cropareaR_in_YXcy, ...
                countries_YX,countries_key,listCountries_map_present,...
                verbose, varargin)
            
ctry_excluded_area_thresh = Inf ;
if ~isempty(varargin)
    ctry_excluded_area_thresh = varargin{1} ;
    varargin(1) = [] ;
    if ~isempty(varargin)
        error('Too many optional arguments included')
    end
end

% Process empties
do_croparea = true ;
if isempty(croparea_in_YXcy)
    do_croparea = false ;
    croparea_out_Ccy = [] ;
    yield_out_Ccy = [] ;
    cropareaR_out_Ccy = [] ;
    R_out_Ccy = [] ;
end
removed_area_dueto_NaNsim = ~isempty(cropareaR_in_YXcy) ;
if ~removed_area_dueto_NaNsim
    cropareaR_out_Ccy = [] ;
    R_out_Ccy = [] ;
elseif ~do_croparea
    error('removed_area_dueto_NaNsim requires croparea_in_YXcy')
end

% Set up empty _Ccy arrays
Ncountries = length(listCountries_map_present) ;
Ncrops = size(total_in_YXcy,3) ;
Nyears = size(total_in_YXcy,4) ;
total_out_Ccy = nan(Ncountries,Ncrops,Nyears) ;
if do_croparea
    croparea_out_Ccy = nan(Ncountries,Ncrops,Nyears) ;
    if removed_area_dueto_NaNsim
        cropareaR_out_Ccy = nan(Ncountries,Ncrops,Nyears) ;
    end
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
            if removed_area_dueto_NaNsim
                tmpR_x = cropareaR_in_YXcy(repmat(inThisCountry,[1 1 Ncrops Nyears])) ;
                tmpR_xvy = reshape(tmpR_x,[N_inThisCountry Ncrops Nyears]) ;
                cropareaR_out_Ccy(c,:,:) = sum(tmpR_xvy,1) ;
            end
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
    if removed_area_dueto_NaNsim
        cropareaR_frac_Ccy = cropareaR_out_Ccy ...
            ./ (cropareaR_out_Ccy + croparea_out_Ccy) ;
        R_out_Ccy = cropareaR_frac_Ccy > ctry_excluded_area_thresh ;
        yield_out_Ccy(R_out_Ccy) = NaN ;
        total_out_Ccy(R_out_Ccy) = NaN ;
    end
end



end