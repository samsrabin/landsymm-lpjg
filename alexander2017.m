%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get info from Alexander et al. 2017 uncertainty paper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import

table_in = readtable('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/AlexanderEtAl2017_allNonGriddedData.csv') ;

% Restrict to global runs
table_in(~strcmp(table_in.location,'Global'),:) = [] ;

% Only care about cropland and pasture
table_in(~(strcmp(table_in.land_type,'c') | strcmp(table_in.land_type,'p')),:) = [] ;

% Only care about part of the time period
yearList = 2010:10:2100 ;
Nyears = length(yearList) ;
table_in(table_in.year<min(yearList) | table_in.year>max(yearList),:) = [] ;


%% Restrict to comparables

lineWidth = 3 ;
fontSize = 18 ;
spacing = 0.07 ;
legloc = 'Southeast' ;

combos = {'SSP1', 'RCP 4.5' ;
          'SSP3', 'RCP 6.0' ;
          'SSP4', 'RCP 6.0' ;
          'SSP5', 'RCP 8.5' ;
         } ;

for c = 1:length(combos)
    
    % Get matches
    thisEcon = combos{c,1} ;
    thisClim = combos{c,2} ;
    thisTitle_tmp = sprintf('(%s-%s)',thisEcon, thisClim([5 7])) ;
    isMatch = strcmp(table_in.economicScenario,thisEcon) & strcmp(table_in.rcp,thisClim) ;
    column_runs = strcat(strcat(table_in.model(isMatch),'_'), strrep(table_in.scenarioId(isMatch),' ','')) ;
    matching_runs = unique(column_runs) ;
    matching_runs_4legend = strrep(matching_runs,'_','\_') ;
    Nmatch = length(matching_runs) ;
    if Nmatch==0
        continue
    end
    
    % Get data
    data_vyc = nan(2,Nyears,Nmatch) ;
    for i = 1:Nmatch
        thisRun = matching_runs{i} ;
        tmp = strsplit(thisRun, '_') ;
        thisModel = tmp{1} ;
        thisScenID = tmp{2} ;
        isThisRun = strcmp(table_in.model,thisModel) & strcmp(table_in.scenarioId,thisScenID) ;
        data_x = table_in.value_rebased(isThisRun) ;
        years_x = table_in.year(isThisRun) ;
        lt_x = table_in.land_type(isThisRun) ;
        thisRun_isCrop = strcmp(lt_x,'c') ;
        for y = 1:Nyears
            thisYear = yearList(y) ;
            isThisYear = years_x==thisYear ;
            if any(isThisYear)
                data_vyc(1,y,i) = data_x(isThisYear & thisRun_isCrop) ;
                data_vyc(2,y,i) = data_x(isThisYear & ~thisRun_isCrop) ;
            end
        end
    end
    data_vyc = data_vyc * 0.01 ;   % Mha to Mkm2
    data_ycv = permute(data_vyc,[2 3 1]) ;
    
    % Plot
    % Cropland
    figure('Position',figurePos,'Color','w') ;
    subplot_tight(1,2,1,spacing)
    plot(yearList, data_ycv(:,:,1), 'LineWidth', lineWidth)
    legend(matching_runs_4legend,'Location',legloc)
    set(gca,'FontSize',fontSize) ;
    thisTitle = sprintf('Cropland %s',thisTitle_tmp) ;
    title(thisTitle)
    xlabel('Year')
    ylabel('Area (Mkm^2)')
    % Pasture
    subplot_tight(1,2,2,spacing)
    plot(yearList, data_ycv(:,:,2), 'LineWidth', lineWidth)
    legend(matching_runs_4legend,'Location',legloc)
    set(gca,'FontSize',fontSize) ;
    thisTitle = sprintf('Pasture %s',thisTitle_tmp) ;
    title(thisTitle)
    xlabel('Year')
    ylabel('Area (Mkm^2)')
    
    
end