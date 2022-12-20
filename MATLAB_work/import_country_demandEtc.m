%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Import per-country demand/production/imports (Mt to kg, kcal) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Importing per-country demand/production/imports...')

warning('off','MATLAB:table:ModifiedAndSavedVarnames')
combine_cereals_dom = false ;
for r = 1:Nruns
    thisDir = runDirs_plum{r} ;
    thisTable_dem = readtable(sprintf('%s/countryDemand.txt', thisDir)) ;
    [~,~,sortCols] = intersect({'Country','Commodity','Year'},thisTable_dem.Properties.VariableNames, 'stable') ;
    thisTable_dem = sortrows(thisTable_dem,sortCols) ; % Needed to produce _ymur dimensioned array
    thisfile = sprintf('%s/domestic.txt', thisDir) ;
    try
        thisTable_dom = readtable(thisfile) ;
    catch ME
        if strcmp(ME.identifier, 'MATLAB:readtable:OpenFailed') ...
        && exist([thisfile '.gz'], 'file')
            gunzip([thisfile '.gz']) ;
            thisTable_dom = readtable(thisfile) ;
            gzip(thisfile)
        else
            rethrow(ME)
        end
    end
    
    
    
    
    [~,~,sortCols] = intersect({'Country','Crop','Year'},thisTable_dom.Properties.VariableNames, 'stable') ;
    thisTable_dom = sortrows(thisTable_dom,sortCols) ; % Needed to produce _ymur dimensioned array
    thisTable_dom(contains(thisTable_dom.Crop,{'energycrops','setaside','pasture'}),:) = [] ;
    thisTable_dom.Net_imports = cellfun(@str2num,thisTable_dom.Net_imports) ;
    if r==1
        commods_dem = commods(~strcmp(commods,'crops') & ~strcmp(commods,'livestock')) ;
        Ncommods_dem = length(commods_dem) ;
        if ~isequal(commods_dem,unique(thisTable_dem.Commodity))
            error('Mismatch in commodity list between countryDemand.txt and countryDemand.txt') ;
        elseif ~isequal(yearList_PLUMout,unique(thisTable_dem.Year))
            error('Mismatch in year list between countryDemand.txt and countryDemand.txt') ;
        end
        commods_dom = unique(thisTable_dom.Crop) ;
        Ncommods_dom = length(commods_dom) ;
        if ~isequal(commods_dem,commods_dom)
            commods_dom_tmp = commods_dom ;
            commods_dom_tmp(contains(commods_dom,{'rice','wheat'})) = [] ;
            commods_dom_tmp(strcmp(commods_dom,'maize')) = {'cereals'} ;
            if ~isequal(commods_dem, commods_dom_tmp)
                error('Mismatch in commodity list between domestic.txt and countryDemand.txt; can''t be fixed by combining cereals (at least with current method)')
            end
            clear commods_dom_tmp
            warning('Mismatch in commodity list between domestic.txt and countryDemand.txt. Fixing by combining cereals.') ;
            combine_cereals_dom = true ;
        elseif ~isequal(yearList_PLUMout,unique(thisTable_dom.Year))
            error('Mismatch in year list between domestic.txt and countryDemand.txt') ;
        end
        countryList_dem = unique(thisTable_dem.Country) ;
        Ncountries_dem = length(countryList_dem) ;
        ts_countryDemand_ymur = nan(Nyears_PLUMout, Ncommods_dem, Ncountries_dem, Nruns) ;
        ts_countryProd_ymur = nan(Nyears_PLUMout, Ncommods_dom, Ncountries_dem, Nruns) ;
        ts_countryImps_ymur = nan(Nyears_PLUMout, Ncommods_dom, Ncountries_dem, Nruns) ;
        ts_countryFrum_ymur = nan(Nyears_PLUMout, Ncommods_dom, Ncountries_dem, Nruns) ;
        ts_countryFmon_ymur = nan(Nyears_PLUMout, Ncommods_dom, Ncountries_dem, Nruns) ;
    end
    ts_countryDemand_ymur(:,:,:,r) = reshape(thisTable_dem.Demand,[Nyears_PLUMout Ncommods_dem Ncountries_dem]) ;
    ts_countryProd_ymur(:,:,:,r) = reshape(thisTable_dom.Production,[Nyears_PLUMout Ncommods_dom Ncountries_dem]) ;
    ts_countryImps_ymur(:,:,:,r) = reshape(thisTable_dom.Net_imports,[Nyears_PLUMout Ncommods_dom Ncountries_dem]) ;
    ts_countryFrum_ymur(:,:,:,r) = reshape(thisTable_dom.Rum_feed_amount,[Nyears_PLUMout Ncommods_dom Ncountries_dem]) ;
    ts_countryFmon_ymur(:,:,:,r) = reshape(thisTable_dom.Mon_feed_amount,[Nyears_PLUMout Ncommods_dom Ncountries_dem]) ;
    clear thisDir thisTable* sortCols
end
warning('on','MATLAB:table:ModifiedAndSavedVarnames')

% Get seed & waste rate (Swrt)
ts_countrySwrt_1m = nan(1,Ncommods_dom) ;
for m = 1:Ncommods_dom
    thisCommod = commods_dom{m} ;
    switch thisCommod
        case 'wheat'
            ts_countrySwrt_1m(m) = 9.5;
        case 'maize'
            ts_countrySwrt_1m(m) = 5.1;
        case 'rice'
            ts_countrySwrt_1m(m) = 8.3;
        case 'oilcrops'
            ts_countrySwrt_1m(m) = 4.4;
        case 'pulses'
            ts_countrySwrt_1m(m) = 10.8;
        case 'starchyRoots'
            ts_countrySwrt_1m(m) = 14.3;
        case 'energycrops'
            ts_countrySwrt_1m(m) = 5;
        case 'monogastrics'
            ts_countrySwrt_1m(m) = 3.1;
        case 'ruminants'
            ts_countrySwrt_1m(m) = 2.2;
        otherwise
            error('thisCommod (%s) not recognized for seed/waste rate', thisCommod)
    end
end
ts_countrySwrt_1m = ts_countrySwrt_1m * 0.01 ;
if any(isnan(ts_countrySwrt_1m))
    error('Some commodity is missing seed/waste rate')
end
ts_countrySwrt_ymur = repmat(ts_countrySwrt_1m,[Nyears_PLUMout 1 Ncountries_dem, Nruns]) ;
if combine_cereals_dom
    prod_maizeRiceWheat_ym1r = sum(ts_countryProd_ymur(:,contains(commods_dom,{'maize','rice','wheat'}),:,:),3) ;
    cerealWeights_maizeRiceWheat_ym1r = prod_maizeRiceWheat_ym1r ...
        ./ repmat(sum(prod_maizeRiceWheat_ym1r,2), [1 3 1 1]) ;
    ts_countrySwrt_ymur(:,strcmp(commods_dom,'maize'),:,:) = ...
        sum(ts_countrySwrt_ymur(:,contains(commods_dom,{'maize','rice','wheat'}),:,:) ...
        .* repmat(cerealWeights_maizeRiceWheat_ym1r, [1 1 Ncountries_dem 1]), 2) ;
    ts_countrySwrt_ymur(:,contains(commods_dom,{'rice','wheat'}),:,:) = [] ;
    clear prod_maizeRiceWheat_ym1r cerealWeights_maizeRiceWheat_ym1r
end
[~, i_crop] = setdiff(commods_dem, commods_livestock) ;
prod_ym1r = sum(ts_countryProd_ymur(:,i_crop,:,:),3) ;
weights_ym1r = prod_ym1r ./ repmat(sum(prod_ym1r,2),[1 length(i_crop) 1 1]) ;
ts_countrySwrt_ymur(:,end+1,:,:) = sum(ts_countrySwrt_ymur(:,i_crop,:,:) .* weights_ym1r, 2) ;
[~, i_livestock] = intersect(commods_dem, commods_livestock) ;
prod_ym1r = sum(ts_countryProd_ymur(:,i_livestock,:,:),3) ;
weights_ym1r = prod_ym1r ./ repmat(sum(prod_ym1r,2),[1 length(i_livestock) 1 1]) ;
ts_countrySwrt_ymur(:,end+1,:,:) = sum(ts_countrySwrt_ymur(:,i_livestock,:,:) .* weights_ym1r, 2) ;
clear prod_ym1r weights_ym1r

if combine_cereals_dom
    ts_countryProd_ymur(:,strcmp(commods_dom,'maize'),:,:) = ...
        sum(ts_countryProd_ymur(:,contains(commods_dom,{'maize','rice','wheat'}),:,:),2) ;
    ts_countryProd_ymur(:,contains(commods_dom,{'rice','wheat'}),:,:) = [] ;
    ts_countryImps_ymur(:,strcmp(commods_dom,'maize'),:,:) = ...
        sum(ts_countryImps_ymur(:,contains(commods_dom,{'maize','rice','wheat'}),:,:),2) ;
    ts_countryImps_ymur(:,contains(commods_dom,{'rice','wheat'}),:,:) = [] ;
    ts_countryFrum_ymur(:,strcmp(commods_dom,'maize'),:,:) = ...
        sum(ts_countryFrum_ymur(:,contains(commods_dom,{'maize','rice','wheat'}),:,:),2) ;
    ts_countryFrum_ymur(:,contains(commods_dom,{'rice','wheat'}),:,:) = [] ;
    ts_countryFmon_ymur(:,strcmp(commods_dom,'maize'),:,:) = ...
        sum(ts_countryFmon_ymur(:,contains(commods_dom,{'maize','rice','wheat'}),:,:),2) ;
    ts_countryFmon_ymur(:,contains(commods_dom,{'rice','wheat'}),:,:) = [] ;
end
clear *commods_dom

% Add total crop and livestock demand/production/imports
commods_livestock = {'ruminants','monogastrics'} ;
[~, i_crop] = setdiff(commods_dem, commods_livestock) ;
ts_countryDemand_ymur(:,end+1,:,:) = sum(ts_countryDemand_ymur(:,i_crop,:,:),2) ;
ts_countryProd_ymur(:,end+1,:,:) = sum(ts_countryProd_ymur(:,i_crop,:,:),2) ;
ts_countryImps_ymur(:,end+1,:,:) = sum(ts_countryImps_ymur(:,i_crop,:,:),2) ;
ts_countryFrum_ymur(:,end+1,:,:) = sum(ts_countryFrum_ymur(:,i_crop,:,:),2) ;
ts_countryFmon_ymur(:,end+1,:,:) = sum(ts_countryFmon_ymur(:,i_crop,:,:),2) ;
[~, i_livestock] = intersect(commods_dem, commods_livestock) ;
ts_countryDemand_ymur(:,end+1,:,:) = sum(ts_countryDemand_ymur(:,i_livestock,:,:),2) ;
ts_countryProd_ymur(:,end+1,:,:) = sum(ts_countryProd_ymur(:,i_livestock,:,:),2) ;
ts_countryImps_ymur(:,end+1,:,:) = sum(ts_countryImps_ymur(:,i_livestock,:,:),2) ;
ts_countryFrum_ymur(:,end+1,:,:) = sum(ts_countryFrum_ymur(:,i_livestock,:,:),2) ;
ts_countryFmon_ymur(:,end+1,:,:) = sum(ts_countryFmon_ymur(:,i_livestock,:,:),2) ;

% Convert Mt to kg
ts_countryDemand_ymur = ts_countryDemand_ymur * 1e6*1e3 ;
ts_countryProd_ymur = ts_countryProd_ymur * 1e6*1e3 ;
ts_countryImps_ymur = ts_countryImps_ymur * 1e6*1e3 ;
ts_countryFrum_ymur = ts_countryFrum_ymur * 1e6*1e3 ;
ts_countryFmon_ymur = ts_countryFmon_ymur * 1e6*1e3 ;

% Get domestic production net of seed/waste losses (Prodnet) and fraction
% of true demand satisfied by domestic production (Self)
ts_countryDemandWithFeed_ymur = ts_countryDemand_ymur + ts_countryFrum_ymur + ts_countryFmon_ymur ;
ts_countryProdnet_ymur = ts_countryProd_ymur .* (1 - ts_countrySwrt_ymur) ;
ts_countrySelf_ymur = ts_countryProdnet_ymur ./ ts_countryDemandWithFeed_ymur ;
ts_countrySelf_ymur(ts_countryDemandWithFeed_ymur==0) = NaN ;
ts_countrySelf_ymur(ts_countrySelf_ymur>1) = 1 ;
% Sanity checks
if min(min(min(min(ts_countrySelf_ymur)))) < 0
    error('Minimum value of ts_countrySelf_ymur < 0 (%0.1f)', min(min(min(min(ts_countrySelf_ymur)))))
elseif max(max(max(max(ts_countrySelf_ymur)))) > 1
    error('Maximum value of ts_countrySelf_ymur > 1 (%0.1f)', max(max(max(max(ts_countrySelf_ymur)))))
end

% Get calories (crops only)
ts_countryDemand_kcal_ymur = nan(Nyears_PLUMout, Ncommods, Ncountries_dem, Nruns) ;
% ts_countryProd_kcal_ymur = nan(Nyears_PLUMout, Ncommods, Ncountries_dem, Nruns) ;
% ts_countryImps_kcal_ymur = nan(Nyears_PLUMout, Ncommods, Ncountries_dem, Nruns) ;
for ii = 1:length(i_crop)
    thisCrop = commods{i_crop(ii)} ;
    kcal_per_g = get_kcalDensity(thisCrop) ;
    kcal_per_kg = 1e3 * kcal_per_g ;
    ts_countryDemand_kcal_ymur(:,ii,:,:) = kcal_per_kg*ts_countryDemand_ymur(:,ii,:,:) ;
%     ts_countryProd_kcal_ymur(:,ii,:,:) = kcal_per_kg*ts_countryProd_ymur(:,ii,:,:) ;
%     ts_countryImps_kcal_ymur(:,ii,:,:) = kcal_per_kg*ts_countryImps_ymur(:,ii,:,:) ;
end
ts_countryDemand_kcal_ymur(:,strcmp(commods,'crops'),:,:) = nansum(ts_countryDemand_kcal_ymur,2) ;
% ts_countryProd_kcal_ymur(:,strcmp(commods,'crops'),:,:) = nansum(ts_countryProd_kcal_ymur,2) ;
% ts_countryImps_kcal_ymur(:,strcmp(commods,'crops'),:,:) = nansum(ts_countryImps_kcal_ymur,2) ;

clear *commods_dem
