function [total_fa2_Ccy, croparea_fa2_Ccy, yield_fa2_Ccy, ...
    fao_itemNames, ...
    listCrops_fa2o, Ncrops_fa2o, ...
    listCountries_map_present_all, ...
    is_tropical, is_xtratrop, ...
    calib_ver_used, twofiles, ...
    yieldWasInf_fa2_Ccy] ...
    = get_fao_data(year1,yearN,calib_ver,...
    need_countries, ...
    varargin)

% UNITS
%    Area harvested: ha
%    Production:     metric tons
%    Yield:          Hg/ha

if ~isempty(varargin)
    if length(varargin) ~= 6
        error('0 or 6 optional arguments required: Ncountries, listCountries_map_present, countries_YX, countries_key, faoCommBalElement, is_ggcmi')
    end
    Ncountries = varargin{1} ;
    listCountries_map_present = varargin{2} ;
    countries_YX = varargin{3} ;
    countries_key = varargin{4} ;
    faoCommBalElement = varargin{5} ;
    is_ggcmi = varargin{6} ;
end

yearList = year1:yearN ;

if calib_ver<9
    error('Currently get_fao_data() is only tested for calib_ver>=9.')
end

if need_countries
    [fao, fao1, fao2, twofiles, ...
        listCountries_map_present_all, is_tropical, is_xtratrop] = import_FAO_data(...
        calib_ver, year1, yearN, ...
        need_countries, listCountries_map_present, countries_YX, countries_key) ;
else
    [fao, fao1, fao2, twofiles, ...
        listCountries_map_present_all, is_tropical, is_xtratrop] = import_FAO_data(...
        calib_ver, year1, yearN, ...
        need_countries) ;
end

% Get crop info
[listCrops_fa2o, ...
    listCrops_fa2i, listCrops_fa2i1, listCrops_fa2i2, ...
    FAO_to_FAO_key, FAO_to_FAO_key1, FAO_to_FAO_key2] = ...
    get_FAO_info(calib_ver, twofiles) ;
Ncrops_fa2o = length(listCrops_fa2o) ;

% Get item names
if ~twofiles
    warning('Add code to get itemNames!')
    fao_itemNames = [] ;
else
    fao_itemNames.a = unique(fao1.ItemName(strcmp(fao1.ElementName,'Area harvested'))) ;
    fao_itemNames.p = unique(fao2.ItemName(strcmp(fao2.ElementName,'Production'))) ;
end


% Get _Ccy arrays
ignoreInfYield = true ;
ignoreNoData = true ;
% if ~exist('yield_fa2_Ccy','var') || calib_ver~=calib_ver_used
    verbose = false ;
    if ~twofiles
        disp('Getting FA2 _Ccy arrays...')
        [croparea_fa2_Ccy, total_fa2_Ccy,yield_fa2_Ccy] = ...
            FAO_to_Ccy2(fao,FAO_to_FAO_key,listCountries_map_present,...
            listCrops_fa2i,listCrops_fa2o,yearList,verbose,...
            ignoreInfYield,ignoreNoData) ;
    else
        disp('Getting FA2 _Ccy array: Area...')
        croparea_fa2_Ccy = ...
            FAO_to_Ccy2_areaOnly(fao1,FAO_to_FAO_key1,listCountries_map_present_all,...
            listCrops_fa2i1,listCrops_fa2o,yearList,verbose,...
            ignoreNoData) ;
        disp('Done.')
        disp(' ')
        disp('Getting FA2 _Ccy array: Production...')
        total_fa2_Ccy = ...
            FAO_to_Ccy2_totalOnly(fao2,FAO_to_FAO_key2,listCountries_map_present_all,...
            listCrops_fa2i2,listCrops_fa2o,yearList,verbose,...
            ignoreInfYield,ignoreNoData,faoCommBalElement) ;
        % Reconcile FAO NaNs
        either_fa2_nan = (isnan(total_fa2_Ccy) | isnan(croparea_fa2_Ccy)) ;
        total_fa2_Ccy(either_fa2_nan) = NaN ;
        croparea_fa2_Ccy(either_fa2_nan) = NaN ;
        % There shouldn't be any FAO production if no FAO data
        if any(total_fa2_Ccy>0 & croparea_fa2_Ccy==0)
            error('At least one cell has production without crop area.')
        end
        % Calculate yield (tDM/ha)
        yield_fa2_Ccy = total_fa2_Ccy ./ croparea_fa2_Ccy ;
        yieldWasInf_fa2_Ccy = isinf(yield_fa2_Ccy) ;
        % Sanity check
        if any(yieldWasInf_fa2_Ccy(:))
            if ~ignoreInfYield
                error('At least one member of yield_fa2_Ccy is Inf!')
            else
                warning('At least one member of yield_fa2_Ccy is Inf! IGNORING')
                badCountries = listCountries_map_present_all(squeeze(sum(sum(yieldWasInf_fa2_Ccy,2),3))>0) ;
                badCrops = listCrops_fa2o(squeeze(sum(sum(yieldWasInf_fa2_Ccy,1),3))>0) ;
                badYears = yearList(squeeze(sum(sum(yieldWasInf_fa2_Ccy,2),1))>0) ;
                if length(badCountries)==1 && length(badCrops)==1
                    warning(['This is ' badCrops{1} ' in ' badCountries{1} ' for ' num2str(length(badYears)) ' years.'])
                end
                total_fa2_Ccy(yieldWasInf_fa2_Ccy) = NaN ;
                croparea_fa2_Ccy(yieldWasInf_fa2_Ccy) = NaN ;
                yield_fa2_Ccy(yieldWasInf_fa2_Ccy) = NaN ;
            end
        end
    end
    disp('Done.')
    calib_ver_used = calib_ver ;
% end
% 

if is_ggcmi
    total_fa2_Ccy = nanmean(total_fa2_Ccy,3) ;
    croparea_fa2_Ccy = nanmean(croparea_fa2_Ccy,3) ;
    yield_fa2_Ccy = nanmean(yield_fa2_Ccy,3) ;
end

end




% Process FAO data into country-crop-year arrays
%%%verbose = false ;
% % % % % Define crop map for FAO-->PLUM (based on MIRCA mapping but with names
% % % % % changed to match those used by FAO)
%%%getPi = @(x) find(strcmp(listCrops_lpj_comb,x)) ;
% % % % FAO_to_PLUM_key{getPi('TeWW')}  = {'Wheat','Barley','Rye'} ;
% % % % FAO_to_PLUM_key{getPi('TeSW')}  = {'Pulses,Total','Potatoes','Sugar beet','Cassava','Sunflower seed','Soybeans','Groundnuts, with shell','Rapeseed'} ;
% % % % % FAO_to_PLUM_key{getPi('TeCo')}  = {'Maize','Millet','Sorghum'} ;
% % % % FAO_to_PLUM_key{getPi('TeCo')}  = {'Maize'} ;
% % % % % FAO_to_PLUM_key{getPi('TrRi')}  = {'Rice, paddy'} ;
% % % % FAO_to_PLUM_key{getPi('TrRi')}  = {'Rice, paddy'} ;
% % % % % [croparea_fao_Ccy,total_fao_Ccy,yield_fao_Ccy] = ...
% % % % %             FAO_to_Ccy(fao,FAO_to_PLUM_key,listCountries_map_present,...
% % % % %             listCrops_fao,Ncrops_lpj_comb,listYears_fao,verbose) ;
