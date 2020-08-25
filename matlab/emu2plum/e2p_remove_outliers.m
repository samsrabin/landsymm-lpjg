function [Sdata, outlier_info] = e2p_remove_outliers(Sdata, which_file)

smad_mult = get_smad_mult(which_file) ;

if isfield(Sdata, 'garr_xvt')
    data_xvt = Sdata.garr_xvt ;
else
    data_xvt = Sdata.garr_xv ;
end

info_size = size(data_xvt) ;
info_size = info_size(2:end) ;

Nvars = info_size(1) ;
if length(info_size)==2
    Ntpers = info_size(2) ;
elseif length(info_size)==1
    Ntpers = 1 ;
else
    error('???')
end

outlier_info.thresh_vt = nan(info_size) ;
outlier_info.count_vt = nan(info_size) ;
outlier_info.max_vt = nan(info_size) ;
outlier_info.varNames = Sdata.varNames ;

any_found = false ;
for v = 1:Nvars
    for t = 1:Ntpers
        
        % Calculate outlier threshold: > median + smad_mult*scaled median absolute
        % difference, of cells with any yield. Round up to next 10^-6.
        y = data_xvt(:,v,t) ;
        ygt0 = y(y>0) ;
        scaled_MAD = (-1/(sqrt(2)*erfcinv(3/2)))*median(abs(ygt0-median(ygt0))) ;
        upper_outl_thresh = median(ygt0) + smad_mult*scaled_MAD ;
        upper_outl_thresh = ceil(upper_outl_thresh*1e6) * 1e-6 ;
        
        % Find outliers
        is_upper_outlier = y > upper_outl_thresh ;
        
        % Save info
        outlier_info.thresh_vt(v,t) = upper_outl_thresh ;
        outlier_info.count_vt(v,t) = length(find(is_upper_outlier)) ;
        if any(is_upper_outlier)
            outlier_info.max_vt(v,t) = max(y(is_upper_outlier)) ;
            if outlier_info.count_vt(v,t)==0
                error('???')
            end
        end
        
        % Set outliers to upper_outl_thresh
        if any(is_upper_outlier)
            any_found = true ;
            y(is_upper_outlier) = upper_outl_thresh ;
            data_xvt(:,v,t) = y ;
        end
        
    end
end

% Round info "max" array to nearest 10^-6
outlier_info.max_vt = round(outlier_info.max_vt, 6) ;

if any_found
    if isfield(Sdata, 'garr_xvt')
        Sdata.garr_xvt = data_xvt ;
    else
        Sdata.garr_xv = data_xvt ;
    end
end



end


function smad_mult = get_smad_mult(which_file)

if strcmp(which_file,'yield')
    smad_mult = 30 ;
elseif strcmp(which_file,'gsirrigation')
    warning('Might have to change gsirrigation smad_mult once properly pre-thresholding')
    smad_mult = 3 ;
else
    error('which_file (%s) not recognized', which_file)
end

end