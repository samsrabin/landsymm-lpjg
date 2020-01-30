%% Import

T = readtable('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/livestock_byregion_FAOSTAT_data_12-13-2019.csv') ;

% Remove unnecessary columns
T(:,contains(T.Properties.VariableNames,{'Code', 'Flag', 'Domain'})) = [] ;

% Replace number strings with numbers
T.Year = str2double(T.Year) ;
T.Value = str2double(T.Value) ;

% head(T)

list_countries = unique(T.Country) ;
list_elements = unique(T.Element) ;
list_items = unique(T.Item) ;
list_years = unique(T.Year) ;

% Convert to PLUM feed-equivalent values (Peter Al. email 2019-12-16)
T_feedeq = T(contains(T.Item, {'Poultry Meat', 'Pigmeat', 'Mutton & Goat Meat', 'Bovine Meat'}), :) ;
T_feedeq.Unit = strrep(T_feedeq.Unit, 'tonnes', 'tonnes PLUM feed equiv.') ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Poultry Meat'))       =  3.3 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Poultry Meat')) ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Pigmeat'))            =  6.4 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Pigmeat')) ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Mutton & Goat Meat')) = 15.0 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Mutton & Goat Meat')) ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Bovine Meat'))        = 25.0 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Bovine Meat')) ;


%% Make delta time series by continent

% thisT = T ;
thisT = T_feedeq ;

% thisElement = 'Domestic supply quantity' ;
thisElement = 'Food supply quantity (tonnes)' ;
itemList_tmp = {'M+R', 'Monogastric meat'} ;
line_styles = {'-', '--', ':'} ;
thisPos = [229   147   911   530] ;

base_year = 1990 ;
yearList_tmp = list_years(list_years>=base_year) ;
Nyears_tmp = length(yearList_tmp) ;
T_tmp0 = thisT(thisT.Year>=base_year ...
    & strcmp(thisT.Element, thisElement), :) ;

% continent_list = {'Africa', 'Asia', 'Northern America', 'Lat. America + Carib.', ...
%     'Europe', 'Oceania'} ;
continent_list = {'Africa', 'Asia + Oceania', 'N. Am. + Europe', ...
    'Lat. Am. + Carib.'} ;
Ncontinents = length(continent_list) ;

hf1 = figure('Color', 'w', 'Position', thisPos) ;
hf2 = figure('Color', 'w', 'Position', thisPos) ;

handles_hf1 = [] ;
handles_hf2 = [] ;
ts_ycm = nan(length(yearList_tmp), Ncontinents, length(itemList_tmp)) ;
for ii = 1:length(itemList_tmp)
    thisItem = itemList_tmp{ii} ;
    if strcmp(thisItem, 'Monogastric meat')
        isThisItem = contains(T_tmp0.Item, {'Pigmeat', 'Poultry Meat'}) ;
    elseif strcmp(thisItem, 'Ruminant meat')
        isThisItem = contains(T_tmp0.Item, {'Bovine Meat', 'Mutton & Goat Meat'}) ;
    elseif strcmp(thisItem, 'M+R')
        isThisItem = contains(T_tmp0.Item, {'Pigmeat', 'Poultry Meat', 'Bovine Meat', 'Mutton & Goat Meat'}) ;
    else
        isThisItem = strcmp(T_tmp0.Item, thisItem) ;
    end
    T_tmp = T_tmp0(isThisItem,:) ;
    ts_yc = nan(length(yearList_tmp), Ncontinents) ;
    for y = 1:Nyears_tmp
        thisYear = yearList_tmp(y) ;
        T_tmp_thisYear = T_tmp(T_tmp.Year==thisYear,:) ;
        for c = 1:Ncontinents
            thisCont = continent_list{c} ;
            if strcmp(thisCont, 'Lat. Am. + Carib.')
                isThisCont = contains(T_tmp_thisYear.Country, {'Caribbean', 'Central America', 'South America'}) ;
            elseif strcmp(thisCont, 'Asia + Oceania')
                isThisCont = contains(T_tmp_thisYear.Country, {'Asia', 'Oceania'}) ;
            elseif strcmp(thisCont, 'N. Am. + Europe')
                isThisCont = contains(T_tmp_thisYear.Country, {'Northern America', 'Europe'}) ;
            else
                isThisCont = strcmp(T_tmp_thisYear.Country, thisCont) ;
            end
            ts_yc(y,c) = sum(T_tmp_thisYear.Value(isThisCont)) ;
        end
    end
    
    ts_ycm(:,:,ii) = ts_yc ;
    
    for c = 1:Ncontinents
        thisCont = continent_list{c} ;
        fprintf('%d %s: %0.1e\n', yearList_tmp(end), thisCont, ts_yc(end,c)) ;
    end
    
    set(0,'CurrentFigure',hf1) ;
    h0 = plot(yearList_tmp, ts_yc, line_styles{ii}, 'LineWidth', 3) ;
    handles_hf1 = [handles_hf1 ; h0] ;
    set(gca, 'ColorOrderIndex', 1) ;
    if ii==1
        hl = legend(handles_hf1, continent_list, 'Location', 'Northwest') ;
        hl.AutoUpdate = 'off' ;
        title(sprintf('FAOSTAT: %s', thisElement))
        ylabel(sprintf('Demand (%s)', T_tmp_thisYear.Unit{1}))
        set(gca, 'FontSize', 14)
    end
    if ii==1
        hold on
    end
    
    tsR_yc = ts_yc - repmat(ts_yc(1,:), [Nyears_tmp 1]) ;
    for c = 1:Ncontinents
        thisCont = continent_list{c} ;
        fprintf('%d %s: %0.1e\n', yearList_tmp(end), thisCont, tsR_yc(end,c)) ;
    end
    
    set(0,'CurrentFigure',hf2) ;
    if ii==1
        plot(yearList_tmp, 0*yearList_tmp, '--k')
        set(gca, 'XLim', minmax_ssr(yearList_tmp))
        hold on
    end
    h0 = plot(yearList_tmp, tsR_yc, line_styles{ii}, 'LineWidth', 3) ;
    handles_hf2 = [handles_hf2 ; h0] ;
    set(gca, 'ColorOrderIndex', 1) ;
    if ii==1
        hl = legend(handles_hf2, continent_list, 'Location', 'Northwest') ;
        hl.AutoUpdate = 'off' ;
        title(sprintf('FAOSTAT: %s', thisElement))
        ylabel(sprintf('Change since %d (%s)', base_year, T_tmp_thisYear.Unit{1}))
        set(gca, 'FontSize', 14)
    end
end
plot(yearList_tmp, 0*yearList_tmp, '--k')
set(0,'CurrentFigure',hf1) ;
text(0.25, 0.95, ...
    {'Solid: Monogastrics (poultry and pigs) + ruminants (beef, mutton, and goat)',...
    'Dashed: Monogastrics only'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 12) ;
set(0,'CurrentFigure',hf2) ;
text(0.25, 0.95, ...
    {'Solid: Monogastrics (poultry and pigs) + ruminants (beef, mutton, and goat)',...
    'Dashed: Monogastrics only'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 12) ;

%{'M+R', 'Monogastric meat'}
figure('Color', 'w', 'Position', thisPos) ;
set(gca, 'ColorOrderIndex', 1) ;
ydata = ts_ycm(:,:,strcmp(itemList_tmp,'Monogastric meat')) ...
    ./ ts_ycm(:,:,strcmp(itemList_tmp,'M+R')) ;
plot(yearList_tmp, ydata, 'LineWidth', 3)
title('Monogastric fraction')
set(gca, 'FontSize', 14)
legend(continent_list, 'Location', 'Best')

% clear *tmp
