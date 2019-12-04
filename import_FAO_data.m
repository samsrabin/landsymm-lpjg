function [fao, fao1, fao2, twofiles, ...
    listCountries_map_present_all, ...
    is_tropical, is_xtratrop] = import_FAO_data(...
    calib_ver, year1, yearN, ...
    countries_YX, countries_key, listCountries_map_present, ...
    strip_fao_nans, fix_cotedivoire, combine_sudans, ...
    combine_subChinas, combine_serbmont)

check_country_names = true ;
Ncountries = length(listCountries_map_present) ;

% Import
if calib_ver <= 4 || calib_ver==6 || calib_ver==7 || calib_ver==11
    fao_filename = 'FAOSTAT_20160928.csv' ;
    fao_filename_trimmed = [fao_filename '.trimmed.txt'] ;
elseif calib_ver==5
    fao_filename_trimmed = 'FAOSTAT_20170410_CommodityBalances_Crops_Production_trimmed_rearr.csv' ;
elseif calib_ver==8
    fao_filename_trimmed = 'FAOSTAT_20170412_Production_Crops_E_All_Data_neg999s_rearr.csv' ;
elseif calib_ver==9 || calib_ver==10 || (calib_ver>=12 && calib_ver<=16) || (calib_ver>=18 && calib_ver<=20)
    fao_filename_trimmed = 'CommodityBalances_Crops_E_All_Data_Norm.csv' ;
elseif calib_ver==17
    fao_filename_trimmed = 'Production_Crops_E_All_Data_(Normalized).csv' ;
elseif calib_ver==21
    fao_filename_trimmed = 'Production_Crops_E_All_Data_(Normalized).csv' ;
else
    error(['calib_ver not recognized: ' num2str(calib_ver)])
end
if ~exist(fao_filename_trimmed,'file')
    twofiles = false ;
    disp('Reading FAO data from CSV...')
    fao = readtable('FAOSTAT_20160928.csv') ;
    
    % Get rid of unnecessary columns
    fao = fao(:,{'AreaName','ElementName','ItemName','Year','Value'}) ;
    
    % Trim metadata row
    fao(cellfun(@isempty,fao.Year),:) = [] ;
    
    % Trim any rows with no data
    fao(cellfun(@isempty,fao.Value),:) = [] ;
    
    % Force numbers to numeric
    fao_orig = fao ;
    disp('Forcing numbers to numeric...')
    tmp = cellfun(@str2num,fao.Year) ;
    fao.Year = tmp ;
    clear tmp
    tmp = cellfun(@str2num,fao.Value) ;
    fao.Value = tmp ;
    clear tmp
    
    % Save
    disp('Saving to TXT file...')
    writetable(fao,fao_filename_trimmed) ;
else
    if calib_ver<=4 || calib_ver==6 || calib_ver==7 || calib_ver==8 || calib_ver==11 || calib_ver==17
        disp('Reading FAO data from TXT file...')
        fao = readtable(fao_filename_trimmed,'Format','%u%s%u%s%u%s%u%u%s%f%s') ;
        twofiles = false ;
    elseif calib_ver==5
        disp('Reading FAO data from TXT file 1...')
        fao1 = readtable('FAOSTAT_20160928.csv.trimmed.txt') ;
        disp('Reading FAO data from TXT file 2...')
        fao2 = readtable('FAOSTAT_20170410_CommodityBalances_Crops_Production_trimmed_rearr.csv') ;
        twofiles = true ;
    elseif calib_ver==9 || calib_ver==10 || (calib_ver>=12 && calib_ver<=16) || (calib_ver>=18 && calib_ver<=20)
        disp('Reading FAO data from TXT file 1...')
        fao1 = readtable('/Users/Shared/PLUM/crop_calib_data/fao/FAOStat-Dec2015_fromPeterAl/Production_Crops_E_All_Data_Norm.csv') ;
        extraneous_columns = {'CountryCode', 'ItemCode', 'ElementCode', ...
            'Unit', 'Flag', 'YearCode'} ;
        fao1(:, contains(fao1.Properties.VariableNames, extraneous_columns)) = [] ;
        fao1 = fao1(:,[1 3 2 4 5]) ;
        fao1.Properties.VariableNames = {'AreaName','ElementName','ItemName','Year','Value'} ;
        disp('Reading FAO data from TXT file 2...')
        fao2 = readtable('/Users/Shared/PLUM/crop_calib_data/fao/FAOStat-Dec2015_fromPeterAl/CommodityBalances_Crops_E_All_Data_Norm.csv') ;
        fao2(:, contains(fao2.Properties.VariableNames, extraneous_columns)) = [] ;
        fao2 = fao2(:,[1 3 2 4 5]) ;
        fao2.Properties.VariableNames = {'AreaName','ElementName','ItemName','Year','Value'} ;
        twofiles = true ;
    elseif calib_ver==21
        disp('Reading FAO data from TXT file...')
        fao = readtable('/Users/Shared/PLUM/crop_calib_data/fao/FAOStat-Dec2015_fromPeterAl/Production_Crops_E_All_Data_Norm.csv') ;
        extraneous_columns = {'CountryCode', 'ItemCode', 'ElementCode', ...
            'Unit', 'Flag', 'YearCode'} ;
        fao(:, contains(fao.Properties.VariableNames, extraneous_columns)) = [] ;
        fao = fao(:,[1 3 2 4 5]) ;
        fao.Properties.VariableNames = {'AreaName','ElementName','ItemName','Year','Value'} ;
        twofiles = false ;
    else
        error(['calib_ver not recognized: ' num2str(calib_ver)])
    end
    
    if calib_ver==17                                             % Block added AB 2018-06-20
        fao(:,contains(fao.Properties.VariableNames,{'Code','Unit','Flag'})) = [] ;
        fao = fao(:,[1 3 2 4 5]) ;
        fao.Properties.VariableNames = {'AreaName','ElementName','ItemName','Year','Value'} ;
        
        % Trim any rows with no data
        fao(cellfun(@isempty,fao.Value),:) = [] ;
        
        % Force numbers to numeric
        fao_orig = fao ;
        disp('Forcing numbers to numeric (year)...')
        tmp = cellfun(@str2num,fao.Year) ;
        fao.Year = tmp ;
        clear tmp
        disp('Forcing numbers to numeric (value)...')
        tmp = cellfun(@str2num,fao.Value) ;    %    'UniformOutput', false added
        fao.Value = tmp ;
        clear tmp
        disp('Done forcing numbers to numeric.')
    end   
    
    % Restrict data to years of interest
    if ~twofiles
        fao = fao(fao.Year>=year1 & fao.Year<=yearN,:) ;
    else
        fao1 = fao1(fao1.Year>=year1 & fao1.Year<=yearN,:) ;
        fao2 = fao2(fao2.Year>=year1 & fao2.Year<=yearN,:) ;
    end
    
    % Remove NaNs
    if strip_fao_nans
        if ~twofiles
            fao(isnan(fao.Value),:) = [] ;
        else
            fao1(isnan(fao1.Value),:) = [] ;
            fao2(isnan(fao2.Value),:) = [] ;
        end
    end
    
    % Remove regions
    regions = {'Africa','Americas','Asia','Australia & New Zealand','Caribbean','Central America','Central Asia','Dominica','Eastern Africa','Eastern Asia','Eastern Europe','Europe','European Union','Land Locked Developing Countries','Least Developed Countries','Low Income Food Deficit Countries','Melanesia','Micronesia','Middle Africa','Net Food Importing Developing Countries','Northern Africa','Northern America','Northern Europe','Oceania','Polynesia','Small Island Developing States','South America','South-Eastern Asia','Southern Africa','Southern Asia','Southern Europe','Western Africa','Western Asia','Western Europe','World'} ;
    for rr = 1:length(regions)
        thisReg = regions{rr} ;
        if ~twofiles
            fao(strcmp(fao.AreaName,thisReg),:) = [] ;
        else
            fao1(strcmp(fao1.AreaName,thisReg),:) = [] ;
            fao2(strcmp(fao2.AreaName,thisReg),:) = [] ;
        end
    end
    
    % Correct country names
    if fix_cotedivoire
        if ~twofiles
            fao.AreaName(strcmp(fao.AreaName,'Côte d''Ivoire')) = {'Cote d''Ivoire'} ;
            fao.AreaName(strcmp(fao.AreaName,'C�te d''Ivoire')) = {'Cote d''Ivoire'} ;
        else
            fao1.AreaName(strcmp(fao1.AreaName,'Côte d''Ivoire')) = {'Cote d''Ivoire'} ;
            fao2.AreaName(strcmp(fao2.AreaName,'Côte d''Ivoire')) = {'Cote d''Ivoire'} ;
        end
    end
    if combine_sudans
        if ~twofiles
            fao.AreaName(strcmp(fao.AreaName,'Sudan (former)')) = {'Sudan'} ;
        else
            fao1.AreaName(strcmp(fao1.AreaName,'Sudan (former)')) = {'Sudan'} ;
            fao2.AreaName(strcmp(fao2.AreaName,'Sudan (former)')) = {'Sudan'} ;
        end
    end
    if ~twofiles
        fao.AreaName(strcmp(fao.AreaName,'Czechia')) = {'Czech Republic'} ;
    else
        fao2.AreaName(strcmp(fao2.AreaName,'Czechia')) = {'Czech Republic'} ;
    end
    
    if combine_subChinas
        if ~twofiles
%             fao = fao(~cellstrfind(fao.AreaName,'China,',true,true),:) ;
            fao = fao(~contains(fao.AreaName,'China,'),:) ;
        else
            if calib_ver==5
                % For PRODUCTION data, get rid of sub-Chinas
                fao1 = fao1(~cellstrfind(fao1.AreaName,'China,',true,true),:) ;
                % For COMMODITY BALANCES data, compile Chinese parts into China
                tmp = fao2(cellstrfind(fao2.AreaName,'China'),:) ;
                china = grpstats(tmp,{'ElementName','ItemName','Year'},'sum','DataVars','Value') ;
                china.Properties.RowNames = {} ;
                china.GroupCount = [] ;
                china = [table(repmat('China',[size(china,1) 1])) china] ;
                china.Properties.VariableNames = fao2.Properties.VariableNames ;
                china.AreaName = cellstr(china.AreaName) ;
                fao2 = fao2(~cellstrfind(fao2.AreaName,'China',true,true),:) ;
                fao2 = vertcat(fao2,china) ;
                clear tmp china
            elseif calib_ver>=9
                % Get rid of sub-Chinas
                fao1 = fao1(~cellstrfind(fao1.AreaName,'China ',true,true),:) ;
                fao2 = fao2(~cellstrfind(fao2.AreaName,'China ',true,true),:) ;
            else
                error(['calib_ver not recognized: ' num2str(calib_ver)])
            end
        end
    end
    
    % For COMMODITY BALANCE data, combine Serbia and Montenegro,
    % because PRODUCTION data don't have them as separate.
    if twofiles
        serbOrMont = strcmp(fao2.AreaName,'Serbia') | ...
            strcmp(fao2.AreaName,'Montenegro') ;
        if combine_serbmont && any(serbOrMont)
            tmp = fao2(serbOrMont,:) ;
            serbmont = grpstats(tmp,{'ElementName','ItemName','Year'},'sum','DataVars','Value') ;
            serbmont.Properties.RowNames = {} ;
            serbmont.GroupCount = [] ;
            serbmont = [table(repmat('Serbia and Montenegro',[size(serbmont,1) 1])) serbmont] ;
            serbmont.Properties.VariableNames = fao2.Properties.VariableNames ;
            serbmont.AreaName = cellstr(serbmont.AreaName) ;
            fao2 = fao2(~serbOrMont,:) ;
            fao2 = vertcat(fao2,serbmont) ;
            clear tmp serbmont
        end
    end
end
disp('Done.')


% What countries, crops, and years are in the dataset?
if ~twofiles
    % listCrops_fao = unique(fao.ItemName) ;
    listCountries_fao = unique(fao.AreaName) ;
    listYears_fao = unique(fao.Year) ;
    % Ncrops_fao = length(listCrops_fao) ;
    Ncountries_fao = length(listCountries_fao) ;
    %%%    Nyears_fao = length(listYears_fao) ;
    if ~isequal(listYears_fao',min(listYears_fao):max(listYears_fao))
        error('Missing year(s)???')
    end
else
    listCountries_fao1 = unique(fao1.AreaName) ;
    listYears_fao1 = unique(fao1.Year) ;
    Ncountries_fao1 = length(listCountries_fao1) ;
    %%%    Nyears_fao1 = length(listYears_fao1) ;
    if ~isequal(listYears_fao1',min(listYears_fao1):max(listYears_fao1))
        error('Missing year(s) from 1???')
    end
    listCountries_fao2 = unique(fao2.AreaName) ;
    listYears_fao2 = unique(fao2.Year) ;
    Ncountries_fao2 = length(listCountries_fao2) ;
    %%%    Nyears_fao2 = length(listYears_fao2) ;
    if ~isequal(listYears_fao2',min(listYears_fao2):max(listYears_fao2))
        error('Missing year(s) from 2???')
    end
end

% Check for country differences
if check_country_names
    if ~twofiles
        listCountries_map_present_all = listCountries_map_present ;
        disp('IN MAP BUT NOT FAO DATA:')
        for c = 1:Ncountries
            thisCountry = listCountries_map_present{c} ;
            if isempty(exact_string_in_cellarray(listCountries_fao,thisCountry,false))
                disp(['   ' thisCountry])
                listCountries_map_present_all(strcmp(listCountries_map_present_all,thisCountry)) = [] ;
            end
        end
        disp('IN FAO DATA BUT NOT MAP:')
        for c = 1:Ncountries_fao
            thisCountry = listCountries_fao{c} ;
            if isempty(exact_string_in_cellarray(listCountries_map_present,thisCountry,false))
                disp(['   ' thisCountry])
            end
        end
    else
        disp('IN MAP BUT NOT FAO DATA1:')
        for c = 1:Ncountries
            thisCountry = listCountries_map_present{c} ;
            %if isempty(exact_string_in_cellarray(listCountries_fao1,thisCountry,false))
            if ~any(strcmp(listCountries_fao1,thisCountry))
                disp(['   ' thisCountry])
            end
        end
        disp('IN FAO DATA1 BUT NOT MAP:')
        for c = 1:Ncountries_fao1
            thisCountry = listCountries_fao1{c} ;
            %if isempty(exact_string_in_cellarray(listCountries_map_present,thisCountry,false))
            if ~any(strcmp(listCountries_map_present,thisCountry))
                disp(['   ' thisCountry])
            end
        end
        disp('IN MAP BUT NOT FAO DATA2:')
        for c = 1:Ncountries
            thisCountry = listCountries_map_present{c} ;
            %if isempty(exact_string_in_cellarray(listCountries_fao2,thisCountry,false))
            if ~any(strcmp(listCountries_fao2,thisCountry))
                disp(['   ' thisCountry])
            end
        end
        disp('IN FAO DATA2 BUT NOT MAP:')
        for c = 1:Ncountries_fao2
            thisCountry = listCountries_fao2{c} ;
            %if isempty(exact_string_in_cellarray(listCountries_map_present,thisCountry,false))
            if ~any(strcmp(listCountries_map_present,thisCountry))
                disp(['   ' thisCountry])
            end
        end
        disp('IN FAO DATA1 BUT NOT FAO DATA2:')
        for c = 1:Ncountries_fao1
            thisCountry = listCountries_fao1{c} ;
            %if isempty(exact_string_in_cellarray(listCountries_fao2,thisCountry,false))
            if ~any(strcmp(listCountries_fao2,thisCountry))
                disp(['   ' thisCountry])
            end
        end
        disp('IN FAO DATA2 BUT NOT FAO DATA1:')
        for c = 1:Ncountries_fao2
            thisCountry = listCountries_fao2{c} ;
            %if isempty(exact_string_in_cellarray(listCountries_fao1,thisCountry,false))
            if ~any(strcmp(listCountries_fao1,thisCountry))
                disp(['   ' thisCountry])
            end
        end
        disp('IN MAP BUT NOT BOTH FAO DATA1 AND DATA2:')
        listCountries_map_present_all = listCountries_map_present ;
        for c = 1:Ncountries
            thisCountry = listCountries_map_present{c} ;
            if ~any(strcmp(listCountries_fao1,thisCountry)) || ~any(strcmp(listCountries_fao2,thisCountry))
                disp(['   ' thisCountry])
                listCountries_map_present_all(strcmp(listCountries_map_present_all,thisCountry)) = [] ;
            end
        end
%         disp('IN FAO DATA1 :')
%         disp(listCountries_fao1)
%         disp('IN FAO DATA2 :')
%         disp(listCountries_fao2)
        
        
    end
    [is_tropical,is_xtratrop] = classify_tropical_countries(...
        listCountries_map_present_all,countries_YX,countries_key) ;
end

if twofiles
    fao = [] ;
else
    fao1 = [] ;
    fao2 = [] ;
end


end