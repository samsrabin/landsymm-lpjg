%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Doing LPJ-GUESS crop calibration for PLUM runs %%%
%%%                Realized yields                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

script_setup_cropCalibration


%% Read LPJ-GUESS data and finish countries processing

if calib_ver<=4
    script_import_lpj_yields
elseif calib_ver>=5 && calib_ver<=21
    script_import_lpj_yields_noCCy
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

script_adjust_countries


% %% DIVIDE YIELDS BY 2
% 
% yield_lpj.maps_YXvy = yield_lpj.maps_YXvy / 2 ;
% yield_lpj_comb.maps_YXvy = yield_lpj_comb.maps_YXvy / 2 ;
% total_lpj_YXcy_comb = total_lpj_YXcy_comb / 2 ;


%% Read FAO data

% UNITS FOR FAO INPUTS
%    Area harvested: ha
%    Production:     metric tons
%    Yield:          Hg/ha

faoCommBalElement = 'Production' ;
% faoCommBalElement = 'Feed' ;
% faoCommBalElement = 'Food' ;
% faoCommBalElement = 'Waste' ;
% faoCommBalElement = 'Processing' ;
% faoCommBalElement = 'Other uses' ;

if ~strcmp(faoCommBalElement,'Production')
    error(['faoCommBalElement is set to something other than Production. '...
           'If this is intentional, continue manually.'])
end

[total_fa2_Ccy, croparea_fa2_Ccy, yield_fa2_Ccy, ...
    fao_itemNames, ...
    listCrops_fa2o, Ncrops_fa2o, ...
    listCountries_map_present_all, ...
    is_tropical, is_xtratrop, ...
    calib_ver_used, twofiles] ...
    = get_fao_data(year1,yearN,calib_ver,...
    Ncountries, listCountries_map_present, countries_YX, countries_key, ...
    faoCommBalElement) ;

verbose = false ;
getPi = @(x) find(strcmp(listCrops_lpj_comb,x)) ;


%% (OPTIONAL) Save FAO data for 1961-2010

% tmp_fao_yearList = 1961:2010 ;
% [tmp_total_fa2_Ccy, tmp_croparea_fa2_Ccy, tmp_yield_fa2_Ccy, ...
%     ~, listCrops_fa2o, ~, ...
%     listCountries_map_present_all, ~, ~, ~, ~, ...
%     yieldWasInf_fa2_Ccy] ...
%     = get_fao_data(tmp_fao_yearList(1),tmp_fao_yearList(end),calib_ver,...
%     Ncountries, listCountries_map_present, countries_YX, countries_key, ...
%     faoCommBalElement) ;
% %
% tmp_outFile = ['/Users/Shared/PLUM/crop_calib_data/fao/FAOdata_' num2str(tmp_fao_yearList(1)) '-' num2str(tmp_fao_yearList(end)) '_calibVer' num2str(calib_ver) '_' faoCommBalElement '.mat'] ;
% if ~exist(tmp_outFile,'file')
%     save(tmp_outFile,'tmp_total_fa2_Ccy','tmp_croparea_fa2_Ccy',...
%         'tmp_yield_fa2_Ccy','listCrops_fa2o','tmp_fao_yearList',...
%         'listCountries_map_present_all','yieldWasInf_fa2_Ccy') ;
% else
%     error('tmp_outFile already exists!')
% end
% 
% clear tmp*
% disp('Done saving.')


%% Make Ccy arrays from LPJ data

if ~twofiles
    listCountries_map_present_all = listCountries_map_present ;
end

disp('Making LPJ _Ccy datasets... ')
verbose = false ;
total_lpj_YXcy_comb2 = total_lpj_YXcy_comb(:,:,:,yield_lpj_comb.yearList>=year1 & yield_lpj_comb.yearList<=yearN) ;
croparea_lpj_YXcy_comb2 = croparea_lpj_YXcy_comb(:,:,:,yield_lpj_comb.yearList>=year1 & yield_lpj_comb.yearList<=yearN) ;
[croparea_lpj_Ccy,total_lpj_Ccy,yield_lpj_Ccy] = ...
    YXcy_to_Ccy(total_lpj_YXcy_comb2,croparea_lpj_YXcy_comb2,...
    countries_YX,countries_key,listCountries_map_present_all,...
    verbose) ;
disp('Done.')



%% Get calibration factors via slope-only regression; also plot Miscanthus calibration

% Options %%%%%%%%%

% Threshold multiple of IQR above/below median to consider something an
% outlier. To find no outliers, set to Inf. (MATLAB default is 1.5.)
outlier_thresh = Inf ;

%%% Base regression weights on FA2 total harvest or crop area? (Empty for no
%%% weighting)
% regWeight_basedOn = 'total' ;
% regWeight_basedOn = 'croparea' ;
regWeight_basedOn = '' ;

% Weight point size based on regression weights?
scatter_style = 'size_uniform' ;
% scatter_style = 'size_weighted' ;

% Restrict where observed values <= Xth percentile
max_prctile = 100 ;
% max_prctile = 90 ;

% Average over years?
avg_over_yrs = false ;
% avg_over_yrs = true ;

% Regression line width
reg_line_width = 3 ;

% Font size in figure
fig_font_size = 16 ;

% Miscanthus file
if calib_ver==11 || calib_ver==17 || is_ggcmi || calib_ver==21
    miscanthus_file = '' ;
elseif calib_ver<=16 || (calib_ver>=18 && calib_ver<=20)
    warning('Using horrible Miscanthus kludge from miscanthus_calibration_kludge.m!')
    miscanthus_file = 'Miscanthus_yields_for_plot.mat' ;
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized in "Miscanthus file"'])
end

% Symbol for slope in legend
% slope_symbol = '\beta' ;
slope_symbol = 'CF' ;

%%%%%%%%%%%%%%%%%%%

if size(yield_lpj_Ccy,3) ~= size(yield_fa2_Ccy,3)
    % Restrict years
    yI_lpj = yield_lpj_comb.yearList>=year1 & yield_lpj_comb.yearList<=yearN ;
    yI_fao = listYears_fao>=year1 & listYears_fao<=yearN ;
    % Extract data for these years
    if size(yield_lpj_Ccy,3) == length(yI_lpj)
        yield_lpj_4cal_Ccy = yield_lpj_Ccy(:,:,yI_lpj) ;
        croparea_lpj_4cal_Ccy = croparea_lpj_Ccy(:,:,yI_lpj) ;
        total_lpj_4cal_Ccy = total_lpj_Ccy(:,:,yI_lpj) ;
    else
        yield_lpj_4cal_Ccy = yield_lpj_Ccy ;
        croparea_lpj_4cal_Ccy = croparea_lpj_Ccy ;
        total_lpj_4cal_Ccy = total_lpj_Ccy ;
    end
    if size(yield_fa2_Ccy,3) == length(yI_fao)
        yield_fa2_4cal_Ccy = yield_fa2_Ccy(:,:,yI_fao) ;
        croparea_fa2_4cal_Ccy = croparea_fa2_Ccy(:,:,yI_fao) ;
        total_fa2_4cal_Ccy = total_fa2_Ccy(:,:,yI_fao) ;
    else
        yield_fa2_4cal_Ccy = yield_fa2_Ccy ;
        croparea_fa2_4cal_Ccy = croparea_fa2_Ccy ;
        total_fa2_4cal_Ccy = total_fa2_Ccy ;
    end
else
    % Restrict years
    yI_lpj = ones(size(croparea_lpj_Ccy,3)) ;
    yI_fao = ones(size(croparea_lpj_Ccy,3)) ;
    % Extract data for these years
    yield_lpj_4cal_Ccy = yield_lpj_Ccy ;
    yield_fa2_4cal_Ccy = yield_fa2_Ccy ;
    croparea_lpj_4cal_Ccy = croparea_lpj_Ccy ;
    croparea_fa2_4cal_Ccy = croparea_fa2_Ccy ;
    total_fa2_4cal_Ccy = total_fa2_Ccy ;
    total_lpj_4cal_Ccy = total_lpj_Ccy ;
end

% Convert to form for regression
if calib_ver==11 || calib_ver==21
    listCrops_4cal = {'TrSo','GlyM','FaBe'} ;
    Ncrops_4cal = length(listCrops_4cal) ;
    yield_lpj_4cal_tmp.TrSo_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TrSo'),:)) ;
    yield_lpj_4cal_tmp.GlyM_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('GlyM'),:)) ;
    yield_lpj_4cal_tmp.FaBe_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('FaBe'),:)) ;
    yield_lpj_4cal_Cyc = nan([size(yield_lpj_4cal_tmp.TrSo_Cy) ...
                              length(listCrops_4cal)]) ;
    for c = 1:length(listCrops_4cal)
        thisCrop = listCrops_4cal{c} ;
        eval(['yield_lpj_4cal_Cyc(:,:,c) = yield_lpj_4cal_tmp.' thisCrop '_Cy ;']) ;
    end
    clear yield_lpj_4cal_tmp
elseif (calib_ver>=12 && calib_ver<=16) || (calib_ver>=18 && calib_ver<=20)
    listCrops_4cal = listCrops_lpj_comb ;
    Ncrops_4cal = length(listCrops_4cal) ;
    for c = 1:Ncrops_4cal
        thisCrop = listCrops_4cal{c} ;
        eval(['yield_lpj_4cal_tmp.' thisCrop '_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi(thisCrop),:)) ;']) ;
    end
%     yield_lpj_4cal_Cyc = nan([size(yield_lpj_4cal_tmp.Rice_Cy) length(listCrops_4cal)]) ;
    tmp = size(total_lpj_Ccy) ;
    yield_lpj_4cal_Cyc = nan(tmp([1 3 2])) ;
    for c = 1:length(listCrops_4cal)
        thisCrop = listCrops_4cal{c} ;
        eval(['yield_lpj_4cal_Cyc(:,:,c) = yield_lpj_4cal_tmp.' thisCrop '_Cy ;']) ;
    end
    clear yield_lpj_4cal_tmp
elseif calib_ver<=10
    listCrops_4cal = {'TeWWorSW','TeSW','TeCo','TrRi'} ;
    Ncrops_4cal = length(listCrops_4cal) ;
    yield_lpj_4cal_tmp.TeWW_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TeWW'),:)) ;
    yield_lpj_4cal_tmp.TeSW_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TeSW'),:)) ;
    yield_lpj_4cal_tmp.TeCo_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TeCo'),:)) ;
    yield_lpj_4cal_tmp.TrRi_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TrRi'),:)) ;
    yield_lpj_4cal_tmp.TeWWorSW_Cy = max(cat(3,yield_lpj_4cal_tmp.TeWW_Cy,yield_lpj_4cal_tmp.TeSW_Cy),[],3) ;
    yield_lpj_4cal_Cyc = nan([size(yield_lpj_4cal_tmp.TeWW_Cy) length(listCrops_4cal)]) ;
    for c = 1:length(listCrops_4cal)
        thisCrop = listCrops_4cal{c} ;
        eval(['yield_lpj_4cal_Cyc(:,:,c) = yield_lpj_4cal_tmp.' thisCrop '_Cy ;']) ;
    end
    clear yield_lpj_4cal_tmp
elseif calib_ver==17
    listCrops_4cal = {'TeWW','TeSW','TeCo','TrRi'} ;
    Ncrops_4cal = length(listCrops_4cal) ;
    yield_lpj_4cal_tmp.TeWW_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TeWW'),:)) ;
    yield_lpj_4cal_tmp.TeSW_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TeSW'),:)) ;
    yield_lpj_4cal_tmp.TeCo_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TeCo'),:)) ;
    yield_lpj_4cal_tmp.TrRi_Cy = squeeze(yield_lpj_4cal_Ccy(:,getPi('TrRi'),:)) ;
    yield_lpj_4cal_Cyc = nan([size(yield_lpj_4cal_tmp.TeWW_Cy) length(listCrops_4cal)]) ;
    for c = 1:length(listCrops_4cal)
        thisCrop = listCrops_4cal{c} ;
        eval(['yield_lpj_4cal_Cyc(:,:,c) = yield_lpj_4cal_tmp.' thisCrop '_Cy ;']) ;
    end
    clear yield_lpj_4cal_tmp
else
    error(['calib_ver ' num2str(calib_ver) ' not recognized in "Convert to form for regression."']) ;
end

% Set up mapping of calibration factors
getPi2 = @(x)find(strcmp(listCrops_4cal,x)) ;
if calib_ver==1 || (calib_ver>=3 && calib_ver<=10)
    FA2_to_PLUM_key{getPi2('TeWWorSW')}  = {'Wheat','Oilcrops'} ;
    FA2_to_PLUM_key{getPi2('TeSW')}  = {'Pulses','Starchy roots'} ;
    FA2_to_PLUM_key{getPi2('TeCo')}  = {'Maize'} ;
    FA2_to_PLUM_key{getPi2('TrRi')}  = {'Rice'} ;
elseif calib_ver==2
    FA2_to_PLUM_key{getPi2('TeWWorSW')}  = {'Wheat','Barley','Oats','Oilcrops'} ;
    FA2_to_PLUM_key{getPi2('TeSW')}  = {'Pulses','Starchy roots'} ;
    FA2_to_PLUM_key{getPi2('TeCo')}  = {'Maize','Millet','Sorghum'} ;
    FA2_to_PLUM_key{getPi2('TrRi')}  = {'Rice'} ;
elseif calib_ver==11 || calib_ver==21
    FA2_to_PLUM_key{getPi2('TrSo')}  = {'Sorghum'} ;
    FA2_to_PLUM_key{getPi2('GlyM')}  = {'Soybean'} ;
    FA2_to_PLUM_key{getPi2('FaBe')}  = {'Faba bean'} ;
elseif calib_ver==12 || calib_ver==15 || calib_ver==16 || calib_ver==18
%     FA2_to_PLUM_key{getPi2('CerealsC3')}  = {'Wheat'} ;
%     FA2_to_PLUM_key{getPi2('CerealsC4')}  = {'Maize'} ;
%     FA2_to_PLUM_key{getPi2('Rice')}  = {'Rice'} ;
%     FA2_to_PLUM_key{getPi2('Oilcrops')}  = {'Oilcrops'} ;
%     FA2_to_PLUM_key{getPi2('StarchyRoots')}  = {'Starchy roots'} ;
%     FA2_to_PLUM_key{getPi2('Pulses')}  = {'Pulses'} ;
    % "If" tests added for compatibility with GGCMI runs
    if any(strcmp(listCrops_4cal,'CerealsC3'))
        FA2_to_PLUM_key{getPi2('CerealsC3')}  = {'Wheat'} ;
    end
    if any(strcmp(listCrops_4cal,'CerealsC4'))
        FA2_to_PLUM_key{getPi2('CerealsC4')}  = {'Maize'} ;
    end
    if any(strcmp(listCrops_4cal,'Rice'))
        FA2_to_PLUM_key{getPi2('Rice')}  = {'Rice'} ;
    end
    if any(strcmp(listCrops_4cal,'Oilcrops'))
        FA2_to_PLUM_key{getPi2('Oilcrops')}  = {'Oilcrops'} ;
    end
    if any(strcmp(listCrops_4cal,'StarchyRoots'))
        FA2_to_PLUM_key{getPi2('StarchyRoots')}  = {'Starchy roots'} ;
    end
    if any(strcmp(listCrops_4cal,'Pulses'))
        FA2_to_PLUM_key{getPi2('Pulses')}  = {'Pulses'} ;
    end
elseif calib_ver==13
    FA2_to_PLUM_key{getPi2('Wheat')}     = {'Wheat'} ;
    FA2_to_PLUM_key{getPi2('Maize')}     = {'Maize'} ;
    FA2_to_PLUM_key{getPi2('Rice')}      = {'Rice'} ;
    FA2_to_PLUM_key{getPi2('Soybeans')}  = {'Soybeans'} ;
    FA2_to_PLUM_key{getPi2('Sorghum')}   = {'Sorghum'} ;
    FA2_to_PLUM_key{getPi2('Pulses')}    = {'Pulses'} ;
elseif calib_ver==14
    FA2_to_PLUM_key{getPi2('Wheat')}     = {'Wheat'} ;
    FA2_to_PLUM_key{getPi2('Maize')}     = {'Maize'} ;
    FA2_to_PLUM_key{getPi2('Rice')}      = {'Rice'} ;
    FA2_to_PLUM_key{getPi2('Sorghum')}   = {'Sorghum'} ;
elseif calib_ver==17
    FA2_to_PLUM_key{getPi2('TeWW')}  = {'CerealsC3'} ;
    FA2_to_PLUM_key{getPi2('TeSW')}  = {'OtherC3'} ;
    FA2_to_PLUM_key{getPi2('TeCo')}  = {'CerealsC4'} ;
    FA2_to_PLUM_key{getPi2('TrRi')}  = {'Rice'} ;
elseif calib_ver==19 || calib_ver==20
    FA2_to_PLUM_key{getPi2('CerealsC3')}     = {'Wheat'} ;
    FA2_to_PLUM_key{getPi2('CerealsC4')}     = {'Maize'} ;
    FA2_to_PLUM_key{getPi2('Rice')}          = {'Rice'} ;
    FA2_to_PLUM_key{getPi2('Oilcrops')}      = {'Oilcrops'} ;
    FA2_to_PLUM_key{getPi2('StarchyRoots')}  = {'Starchy roots'} ;
    FA2_to_PLUM_key{getPi2('Pulses')}        = {'Pulses'} ;
    FA2_to_PLUM_key{getPi2('Sugar')}         = {'Sugar'} ;
    if calib_ver==19
        FA2_to_PLUM_key{getPi2('DateCitGrape')}    = {'DateCitGrape'} ;
    elseif calib_ver==20
        FA2_to_PLUM_key{getPi2('FruitAndVeg')}    = {'FruitAndVeg'} ;
    end
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end

% Remove countries where LPJG and/or FAO have 0 area
if calib_ver==11 || calib_ver==17 || calib_ver==21
    ignore_lpj_Cc = countries2ignore(croparea_lpj_4cal_Ccy) ;
    ignore_fa2_Cc = countries2ignore(croparea_fa2_4cal_Ccy) ;
elseif (calib_ver>=1 && calib_ver<=16) || (calib_ver>=18 && calib_ver<=20)
    ignore_lpj_Cc = false(size(croparea_lpj_4cal_Ccy,1),size(croparea_lpj_4cal_Ccy,2)) ;
    ignore_fa2_Cc = false(size(croparea_fa2_4cal_Ccy,1),size(croparea_fa2_4cal_Ccy,2)) ;
else
    error(['calib_ver (' num2str(calib_ver) ') not recognized! In "Remove countries where LPJG and/or FAO have 0 area"'])
end
if any(sum(ignore_lpj_Cc,1)==0)
    warning('No included countries (LPJ): %s',strjoin(listCrops_lpj_comb(sum(ignore_lpj_Cc,1)==0), ', '))
end
if any(sum(ignore_fa2_Cc,1)==0)
    warning('No included countries (FAO): %s',strjoin(listCrops_fa2o(sum(ignore_fa2_Cc,1)==0), ', '))
end

% Get ancillary arrays
croparea_fa2_4cal_Cy = squeeze(nansum(croparea_fa2_4cal_Ccy,2)) ;
total_fa2_4cal_Cy = squeeze(nansum(total_fa2_4cal_Ccy,2)) ;
croparea_lpj_4cal_Cyc = permute(croparea_lpj_4cal_Ccy,[1 3 2]) ;
croparea_fa2_4cal_Cyc = permute(croparea_fa2_4cal_Ccy,[1 3 2]) ;
total_fa2_4cal_Cyc = permute(total_fa2_4cal_Ccy,[1 3 2]) ;
total_lpj_4cal_Cyc = permute(total_lpj_4cal_Ccy,[1 3 2]) ;
yield_fa2_4cal_Cyc = permute(yield_fa2_4cal_Ccy,[1 3 2]) ;
% yield_lpj_4cal_Cyc = permute(yield_lpj_4cal_Ccy,[1 3 2]) ;

% Get weights for regression
if isempty(regWeight_basedOn)
    weights_fa2_4cal_Cyc = [] ;
else
    if calib_ver==11 || calib_ver==21
        error('Add code to do weights when ignoring countries.')
    elseif ~(calib_ver>=1 || calib_ver<=20)
        error(['calib_ver (' num2str(calib_ver) ') not recognized! In "Get weights for regression"'])
    end
    weights_fa2_4cal_Cyc = nan(size(croparea_fa2_4cal_Ccy,1),size(croparea_fa2_4cal_Ccy,3),size(croparea_fa2_4cal_Ccy,2)) ;
    for c = 1:Ncrops_fa2o
        if strcmp(regWeight_basedOn,'total')
            total_fa2_4cal_CyTHISCROP = total_fa2_4cal_Cyc(:,:,c) ;
            if any(total_fa2_4cal_CyTHISCROP>0 & total_fa2_4cal_Cy==0)
                error(['Somehow crop ' num2str(c) ' has area but total_fa2_4cal_Cy = 0.'])
            end
            total_fa2_4cal_CyTHISCROP(total_fa2_4cal_Cy==0) = NaN ;
            tmp = total_fa2_4cal_CyTHISCROP ./ total_fa2_4cal_Cy ;
            % Sanity check
            if length(find(~isnan(total_fa2_4cal_CyTHISCROP))) ~= length(find(~isnan(tmp)))
                error(['total_fa2_4cal_Cyc(:,:,' num2str(c) ') has ' num2str(length(find(~isnan(total_fa2_4cal_Cyc(:,:,c))))) ' non-NaNs but tmp has ' num2str(length(find(~isnan(tmp))))]) ;
            end
        elseif strcmp(regWeight_basedOn,'croparea')
            croparea_fa2_4cal_CyTHISCROP = croparea_fa2_4cal_Cyc(:,:,c) ;
            if any(croparea_fa2_4cal_CyTHISCROP>0 & croparea_fa2_4cal_Cy==0)
                error(['Somehow crop ' num2str(c) ' has area but croparea_fa2_4cal_Cy = 0.'])
            end
            croparea_fa2_4cal_CyTHISCROP(croparea_fa2_4cal_Cy==0) = NaN ;
            tmp = croparea_fa2_4cal_CyTHISCROP ./ croparea_fa2_4cal_Cy ;
            % Sanity check
            if length(find(~isnan(croparea_fa2_4cal_CyTHISCROP))) ~= length(find(~isnan(tmp)))
                error(['croparea_fa2_4cal_Cyc(:,:,' num2str(c) ') has ' num2str(length(find(~isnan(croparea_fa2_4cal_Cyc(:,:,c))))) ' non-NaNs but tmp has ' num2str(length(find(~isnan(tmp))))]) ;
            end
        else
            error(['regWeight_basedOn not recognized: ' regWeight_basedOn])
        end
        % Save to output array
        weights_fa2_4cal_Cyc(:,:,c) = tmp ;
        clear tmp *THISCROP
    end ; clear c
end

% Get weights for plotting point size
if strcmp(scatter_style,'size_weighted')
    for c = 1:Ncrops_fa2o
        tmp = croparea_fa2_4cal_Cyc(:,:,c) ./ max(max(croparea_fa2_4cal_Cyc(:,:,c))) ;
        tmp = tmp / max(max(tmp)) ;
        if max(max(tmp)) ~= 1
            error(['max(max(tmp)) = ' num2str(max(max(tmp)))])
        end
        weights4pts_Cyc(:,:,c) = tmp ;
        clear tmp
        %     weights4pts_Cyc = weights_fa2_4cal_Cyc(:,:,c)./max(max(weights_fa2_4cal_Cyc(:,:,c))) ;
    end ; clear c
    if any(weights4pts_Cyc<0)
        error('Some value of weights4pts_Cyc is < 0.')
    end
else
%     weights4pts_Cyc = nan(size(croparea_fa2_4cal_Cyc)) ;
    weights4pts_Cyc = [] ;
end

% Rearrange if needed
if Ncrops_4cal==1 && ismatrix(yield_fa2_4cal_Cyc) && ~isvector(yield_fa2_4cal_Cyc)
    yield_fa2_4cal_Cyc = reshape(yield_fa2_4cal_Cyc,[size(yield_fa2_4cal_Cyc,1) 1 size(yield_fa2_4cal_Cyc,2)]) ;
    yield_lpj_4cal_Cyc = reshape(yield_lpj_4cal_Cyc,[size(yield_lpj_4cal_Cyc,1) 1 size(yield_lpj_4cal_Cyc,2)]) ;
    weights_fa2_4cal_Cyc = reshape(weights_fa2_4cal_Cyc,[size(weights_fa2_4cal_Cyc,1) 1 size(weights_fa2_4cal_Cyc,2)]) ;
    weights4pts_Cyc = reshape(weights4pts_Cyc,[size(weights4pts_Cyc,1) 1 size(weights4pts_Cyc,2)]) ;
end

% Do regression: All points

% Stuff for compatibility with GGCMI runs
if ~exist('yield_fa2_4cal_Cyc_ORIG','var')
    yield_fa2_4cal_Cyc_ORIG = yield_fa2_4cal_Cyc ;
else
    yield_fa2_4cal_Cyc = yield_fa2_4cal_Cyc_ORIG ;
end
if ~exist('ignore_fa2_Cc_ORIG','var')
    ignore_fa2_Cc_ORIG = ignore_fa2_Cc ;
else
    ignore_fa2_Cc = ignore_fa2_Cc_ORIG ;
end
if ~exist('weights_fa2_4cal_Cyc_ORIG','var')
    weights_fa2_4cal_Cyc_ORIG = weights_fa2_4cal_Cyc ;
else
    weights_fa2_4cal_Cyc = weights_fa2_4cal_Cyc_ORIG ;
end
if ~exist('weights4pts_Cyc_ORIG','var')
    weights4pts_Cyc_ORIG = weights4pts_Cyc ;
else
    weights4pts_Cyc = weights4pts_Cyc_ORIG ;
end
if ~exist('listCrops_fa2o_ORIG','var')
    listCrops_fa2o_ORIG = listCrops_fa2o ;
else
    listCrops_fa2o = listCrops_fa2o_ORIG ;
end
if ~any(strcmp(listCrops_4cal,'CerealsC3'))
    yield_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Wheat')) = [] ;
    ignore_fa2_Cc(:,strcmp(listCrops_fa2o,'Wheat')) = [] ;
    if ~isempty(weights_fa2_4cal_Cyc)
        weights_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Wheat')) = [] ;
    end
    if ~isempty(weights4pts_Cyc) && any(~isnan(weights4pts_Cyc(:)))
        weights4pts_Cyc(:,:,strcmp(listCrops_fa2o,'Wheat')) = [] ;
    end
    listCrops_fa2o(strcmp(listCrops_fa2o,'Wheat')) = [] ;
end
if ~any(strcmp(listCrops_4cal,'CerealsC4'))
    yield_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Maize')) = [] ;
    ignore_fa2_Cc(:,strcmp(listCrops_fa2o,'Maize')) = [] ;
    if ~isempty(weights_fa2_4cal_Cyc)
        weights_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Maize')) = [] ;
    end
    if ~isempty(weights4pts_Cyc) && any(~isnan(weights4pts_Cyc(:)))
        weights4pts_Cyc(:,:,strcmp(listCrops_fa2o,'Maize')) = [] ;
    end
    listCrops_fa2o(strcmp(listCrops_fa2o,'Maize')) = [] ;
end
if ~any(strcmp(listCrops_4cal,'Rice'))
    yield_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Rice')) = [] ;
    ignore_fa2_Cc(:,strcmp(listCrops_fa2o,'Rice')) = [] ;
    if ~isempty(weights_fa2_4cal_Cyc)
        weights_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Rice')) = [] ;
    end
    if ~isempty(weights4pts_Cyc) && any(~isnan(weights4pts_Cyc(:)))
        weights4pts_Cyc(:,:,strcmp(listCrops_fa2o,'Rice')) = [] ;
    end
    listCrops_fa2o(strcmp(listCrops_fa2o,'Rice')) = [] ;
end
if ~any(strcmp(listCrops_4cal,'Oilcrops'))
    yield_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Oilcrops')) = [] ;
    ignore_fa2_Cc(:,strcmp(listCrops_fa2o,'Oilcrops')) = [] ;
    if ~isempty(weights_fa2_4cal_Cyc)
        weights_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Oilcrops')) = [] ;
    end
    if ~isempty(weights4pts_Cyc) && any(~isnan(weights4pts_Cyc(:)))
        weights4pts_Cyc(:,:,strcmp(listCrops_fa2o,'Oilcrops')) = [] ;
    end
    listCrops_fa2o(strcmp(listCrops_fa2o,'Oilcrops')) = [] ;
end
if ~any(strcmp(listCrops_4cal,'StarchyRoots'))
    yield_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Starchy roots')) = [] ;
    ignore_fa2_Cc(:,strcmp(listCrops_fa2o,'Starchy roots')) = [] ;
    if ~isempty(weights_fa2_4cal_Cyc)
        weights_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Starchy roots')) = [] ;
    end
    if ~isempty(weights4pts_Cyc) && any(~isnan(weights4pts_Cyc(:)))
        weights4pts_Cyc(:,:,strcmp(listCrops_fa2o,'Starchy roots')) = [] ;
    end
    listCrops_fa2o(strcmp(listCrops_fa2o,'Starchy roots')) = [] ;
end
if ~any(strcmp(listCrops_4cal,'Pulses'))
    yield_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Pulses')) = [] ;
    ignore_fa2_Cc(:,strcmp(listCrops_fa2o,'Pulses')) = [] ;
    if ~isempty(weights_fa2_4cal_Cyc)
        weights_fa2_4cal_Cyc(:,:,strcmp(listCrops_fa2o,'Pulses')) = [] ;
    end
    if ~isempty(weights4pts_Cyc) && any(~isnan(weights4pts_Cyc(:)))
        weights4pts_Cyc(:,:,strcmp(listCrops_fa2o,'Pulses')) = [] ;
    end
    listCrops_fa2o(strcmp(listCrops_fa2o,'Pulses')) = [] ;
end


disp('ALL POINTS')
[calib_factors_u,calib_factors_w] = ...
            do_crop_regression(yield_fa2_4cal_Cyc,yield_lpj_4cal_Cyc,...
                               ignore_fa2_Cc,ignore_lpj_Cc,...
                               weights_fa2_4cal_Cyc,weights4pts_Cyc,...
                               listCrops_fa2o,listCrops_4cal,...
                               FA2_to_PLUM_key,scatter_style,...
                               'max_prctile',max_prctile,...
                               'reg_line_width',reg_line_width,...
                               'fig_font_size',fig_font_size,...
                               'miscanthus_file',miscanthus_file,...
                               'slope_symbol',slope_symbol,...
                               'marker_size',25,...
                               'separate_figs',false, ...
                               'outlier_thresh', outlier_thresh) ;

