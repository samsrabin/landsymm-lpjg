include_milk = false ;
include_eggs = false ;


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
included_items = {'Poultry Meat', 'Pigmeat', 'Mutton & Goat Meat', 'Bovine Meat'} ;
if include_milk
    included_items = [included_items {'Milk - Excluding Butter'}] ;
end
if include_eggs
    included_items = [included_items {'Eggs'}] ;
end
T_feedeq = T(contains(T.Item, included_items), :) ;
T_feedeq.Unit = strrep(T_feedeq.Unit, 'tonnes', 'tonnes PLUM feed equiv.') ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Poultry Meat')) ...
    =  3.3 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Poultry Meat')) ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Pigmeat')) ...
    =  6.4 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Pigmeat')) ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Mutton & Goat Meat')) ...
    = 15.0 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Mutton & Goat Meat')) ;
T_feedeq.Value(strcmp(T_feedeq.Item, 'Bovine Meat')) ...
    = 25.0 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Bovine Meat')) ;
if include_milk
    T_feedeq.Value(strcmp(T_feedeq.Item, 'Milk - Excluding Butter')) ...
        = 0.7 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Milk - Excluding Butter')) ;
end
if include_eggs
    T_feedeq.Value(strcmp(T_feedeq.Item, 'Eggs')) ...
        = 2.3 * T_feedeq.Value(strcmp(T_feedeq.Item, 'Eggs')) ;
end


%% Make delta time series by continent

% thisT = T ;
thisT = T_feedeq ;

% thisElement = 'Domestic supply quantity' ;
thisElement = 'Food supply quantity (tonnes)' ;
itemList_tmp = {'M+R', 'Just monogastrics'} ;
line_styles = {'-', '--', ':'} ;
thisPos = [229   147   911   530] ;
fontSize = 18 ;

base_year = 1990 ;
yearList_tmp = list_years(list_years>=base_year) ;
Nyears_tmp = length(yearList_tmp) ;
T_tmp0 = thisT(thisT.Year>=base_year ...
    & strcmp(thisT.Element, thisElement),   :) ;

% % continent_list = {'Africa', 'Asia', 'Northern America', 'Lat. America + Carib.', ...
% %     'Europe', 'Oceania'} ;
% continent_list = {'Africa', 'Asia + Oceania', 'N. Am. + Europe', ...
%     'Lat. Am. + Carib.'} ;
continent_list = {'Africa', 'Asia + Oceania', 'N. Am. + Europe', ...
    'Lat. Am. + Carib.', 'World'} ;
Ncontinents = length(continent_list) ;

hf1 = figure('Color', 'w', 'Position', thisPos) ;
hf2 = figure('Color', 'w', 'Position', thisPos) ;

handles_hf1 = [] ;
handles_hf2 = [] ;
ts_ycm = nan(length(yearList_tmp), Ncontinents, length(itemList_tmp)) ;
for ii = 1:length(itemList_tmp)
    thisItem = itemList_tmp{ii} ;
    if strcmp(thisItem, 'Just monogastrics')
        if include_eggs
            isThisItem = contains(T_tmp0.Item, ...
                {'Pigmeat', 'Poultry Meat', 'Eggs'}) ;
        else
            isThisItem = contains(T_tmp0.Item, ...
                {'Pigmeat', 'Poultry Meat'}) ;
        end
    elseif strcmp(thisItem, 'Just ruminants')
        if include_milk
            isThisItem = contains(T_tmp0.Item, ...
                {'Bovine Meat', 'Mutton & Goat Meat'}) ;
        else
            isThisItem = contains(T_tmp0.Item, ...
                {'Bovine Meat', 'Mutton & Goat Meat', 'Milk - Excluding Butter'}) ;
        end
    elseif strcmp(thisItem, 'M+R')
        isThisItem = contains(T_tmp0.Item, included_items) ;
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
        hl = legend(handles_hf1, continent_list, 'Location', 'West') ;
        hl.AutoUpdate = 'off' ;
        title(sprintf('FAOSTAT: %s', thisElement))
        ylabel(sprintf('Demand (%s)', T_tmp_thisYear.Unit{1}))
        set(gca, 'FontSize', fontSize)
        hold on
    end
    
    % If plotting World, color it black
    if any(contains(continent_list, 'World'))
        kids = get(gca, 'Children') ;
        if ii==1
            w1 = length(kids) - find(strcmp(continent_list, 'World')) + 1 ;
        end
        kids(w1).Color = [0 0 0] ;
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
    set(gca, 'ColorOrderIndex', 1) ;
    h0 = plot(yearList_tmp, tsR_yc, line_styles{ii}, 'LineWidth', 3) ;
    handles_hf2 = [handles_hf2 ; h0] ;
    set(gca, 'ColorOrderIndex', 1) ;
    if ii==1
        hl = legend(handles_hf2, continent_list, 'Location', 'West') ;
        hl.AutoUpdate = 'off' ;
        title(sprintf('FAOSTAT: %s', thisElement))
        ylabel(sprintf('Change since %d (%s)', base_year, T_tmp_thisYear.Unit{1}))
        set(gca, 'FontSize', fontSize)
    end
    % If plotting World, color it black
    if any(contains(continent_list, 'World'))
        kids = get(gca, 'Children') ;
        if ii==1
            w2 = length(kids) - find(strcmp(continent_list, 'World')) ;
        end
        kids(w2).Color = [0 0 0] ;
    end
end
plot(yearList_tmp, 0*yearList_tmp, '--k')
tmp_text = {'Solid: Monogastrics (poultry and pigs) + ruminants (beef, mutton, and goat)',...
    'Dashed: Monogastrics only'} ;
if include_milk
    tmp_text = strrep(tmp_text, 'beef', 'milk not butter, beef') ;
end
if include_eggs
    tmp_text = strrep(tmp_text, 'poultry', 'eggs, poultry,') ;
end
set(0,'CurrentFigure',hf1) ;
text(0.22, 0.95, tmp_text, ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 12) ;
set(0,'CurrentFigure',hf2) ;
text(0.22, 0.95, tmp_text, ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 12) ;

figure('Color', 'w', 'Position', thisPos) ;
set(gca, 'ColorOrderIndex', 1) ;
ydata = ts_ycm(:,:,strcmp(itemList_tmp,'Just monogastrics')) ...
    ./ ts_ycm(:,:,strcmp(itemList_tmp,'M+R')) ;
plot(yearList_tmp, ydata, 'LineWidth', 3)

% If plotting World, color it black
if any(contains(continent_list, 'World'))
    kids = get(gca, 'Children') ;
    if ii==1
        w2 = length(kids) - find(strcmp(continent_list, 'World')) + 1 ;
    end
    kids(w2).Color = [0 0 0] ;
end
    
title('Monogastric fraction')
set(gca, 'FontSize', fontSize, 'XLim', minmax_ssr(yearList_tmp))
legend(continent_list, 'Location', 'Best')

% clear *tmp
