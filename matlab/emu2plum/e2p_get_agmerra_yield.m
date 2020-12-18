function yield_agmerraBL_xv = e2p_get_agmerra_yield(...
    varNames_emu, topDir_phase2, ggcm, list2map, getN, adaptation)

warning('on','all')

% Get info
Nvars_emu = length(varNames_emu) ;

% Set up output array
yield_agmerraBL_xv = nan(length(list2map),Nvars_emu) ;

% Get AgMERRA-forced yield
for v = 1:Nvars_emu
    
    % Get info about this crop
    thisCropIrr = getbasenamei(varNames_emu{v}) ;
    thisN = str2double(getN(varNames_emu{v})) ;
    isirrig = strcmp(thisCropIrr(end),'i') ;
    thisCrop = thisCropIrr ;
    if isirrig
        thisCrop = thisCrop(1:end-1) ;
    end
    
    % Get filename
    switch thisCrop
        case 'maize'; thisCrop_short = 'mai' ;
        case 'soy'; thisCrop_short = 'soy' ;
        case 'rice'; thisCrop_short = 'ric' ;
        case 'spring_wheat'; thisCrop_short = 'swh' ;
        case 'winter_wheat'; thisCrop_short = 'wwh' ;
        otherwise; error('thisCrop (%s) not recognized', thisCrop)
    end
    thisFile = sprintf('%s/%s/A%d/yield/%s_agmerra_fullharm_yield_%s_global_annual_1980_2010_C360_T0_W0_N%d_A%d.nc4', ...
        topDir_phase2, thisCrop, adaptation, lower(ggcm), thisCrop_short, thisN, ...
        adaptation) ;
    if isirrig
        thisFile = strrep(thisFile, 'W0', 'Winf') ;
    end
    nofilefound = false ;
    if ~exist(thisFile, 'file')
        thisAN = sprintf('A%d_N%d', adaptation, thisN) ;
        msg = sprintf('Phase 2 AgMERRA %s %s not found, nor were any substitutes', ...
            thisCropIrr, thisAN) ;
        msg = sprintf('%s (A%d_N%d, A%d_NNA, or A%d_NNA)', ...
            msg, ~adaptation, thisN, adaptation, ~adaptation) ;
        % Try A[other]
        tmp_thisAdapt = ~adaptation ;
        tmp_thisN = thisN ;
        thisFile = sprintf('%s/%s/A%d/yield/%s_agmerra_fullharm_yield_%s_global_annual_1980_2010_C360_T0_W0_N%d_A%d.nc4', ...
            topDir_phase2, thisCrop, tmp_thisAdapt, lower(ggcm), thisCrop_short, thisN, ...
            tmp_thisAdapt) ;
        if ~exist(thisFile, 'file')
            % Try NNA
            tmp_thisAdapt = adaptation ;
            tmp_thisN = NaN ;
            thisFile = sprintf('%s/%s/A%d/yield/%s_agmerra_fullharm_yield_%s_global_annual_1980_2010_C360_T0_W0_NNA_A%d.nc4', ...
                topDir_phase2, thisCrop, tmp_thisAdapt, lower(ggcm), thisCrop_short, ...
                tmp_thisAdapt) ;
            if ~exist(thisFile, 'file')
                % Try A[other] NNA
                tmp_thisAdapt = ~adaptation ;
                thisFile = sprintf('%s/%s/A%d/yield/%s_agmerra_fullharm_yield_%s_global_annual_1980_2010_C360_T0_W0_NNA_A%d.nc4', ...
                    topDir_phase2, thisCrop, tmp_thisAdapt, lower(ggcm), thisCrop_short, ...
                    tmp_thisAdapt) ;
                if ~exist(thisFile, 'file')
                    nofilefound = true ;
                    warning('%s\n  Will not check for missing or low baseline AgMERRA yield!', ...
                        msg) ;
                end
            end
        end
        
        % Warn if using a substitute file
        if ~nofilefound
            if isnan(tmp_thisN)
                replacementAN = sprintf('A%d_NNA', tmp_thisAdapt) ;
            else
                replacementAN = sprintf('A%d_N%d', tmp_thisAdapt, thisN) ;
            end
            warning('Phase 2 AgMERRA %s %s not found; using %s instead', ...
                thisCropIrr, thisAN, replacementAN) ;
        end
    end
    
    % Import (exclude last timestep to avoid incomplete final seasons)
    if nofilefound
        yield_agmerraBL_xv(:,v) = Inf ;
    else
        tmp_XYt = ncread(thisFile,['yield_' thisCrop_short]) ;
        tmp_YX = flipud(transpose(nanmean(tmp_XYt(:,:,1:end-1), 3))) ;
        tmp = tmp_YX(list2map) ;
        if ~any(~isnan(tmp))
            warning('%s is all NaN', thisCrop)
        end
        yield_agmerraBL_xv(:,v) = tmp ;
    end
    
end

end
