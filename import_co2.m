%%%%%%%%%%%%%%%%%%
%%% Import CO2 %%%
%%%%%%%%%%%%%%%%%%


disp('Importing CO2...')

% Baseline
[status, result] = unix(sprintf( ...
    'grep -he "^param \\"file_co2" %s/main_plum_lpjg.ins', ...
    baselineDir)) ;
if status~=0
    error('Error in unix() call')
end
tmp = strsplit(result,'"') ;
co2file = tmp{4} ;
co2file = strrep(co2file,'/home/fh1-project-lpjgpi/lr8247','/Users/Shared/lpj-guess') ;
if ~exist(co2file, 'file')
    error('co2file (%s) not found!', co2file)
end
thisTable_bl = readtable(co2file) ;
thisTable_bl(:,3) = [] ;
thisTable_bl.Properties.VariableNames = {'Year','co2'} ;
% Add extra years to beginning, if needed
while min(thisTable_bl.Year) > min(yearList_baseline)
    thisTable_bl = cat(1, {max(thisTable_bl.Year)-1 thisTable_bl.co2(1)}, thisTable_bl) ;
end
% Add extra years to end, if needed
while max(thisTable_bl.Year) < max(yearList_baseline)
    thisTable_bl = cat(1, thisTable_bl, {max(thisTable_bl.Year)+1 thisTable_bl.co2(end)}) ;
end
thisTable_bl = thisTable_bl(thisTable_bl.Year>=min(yearList_baseline) & thisTable_bl.Year<=max(yearList_baseline),:) ;
ts_co2_bl = thisTable_bl.co2 ;
clear thisTable_bl

% Future
ts_co2_yr = nan(Nyears_fu, Nruns) ;
for r = 1:Nruns
    thisDir = runDirs{r} ;
    [status, result] = unix(sprintf('grep -he "^param \\"file_co2" %s/main_plum_lpjg.ins', thisDir)) ;
    if status~=0
        error('Error in unix() call')
    end
    tmp = strsplit(result,'"') ;
    co2file = tmp{4} ;
    co2file = strrep(co2file,'/home/fh1-project-lpjgpi/lr8247','/Users/Shared/lpj-guess') ;
    if ~exist(co2file, 'file')
        error('co2file (%s) not found!', co2file)
    end
    thisTable = readtable(co2file) ;
    thisTable(:,3) = [] ;
    thisTable.Properties.VariableNames = {'Year','co2'} ;
    
    % Get future CO2
    % Add extra years to beginning, if needed
    while min(thisTable.Year) > min(yearList_future)
        thisTable = cat(1, {max(thisTable.Year)-1 thisTable.co2(1)}, thisTable) ;
    end
    % Add extra years to end, if needed
    while max(thisTable.Year) < max(yearList_future)
        thisTable = cat(1, thisTable, {max(thisTable.Year)+1 thisTable.co2(end)}) ;
    end
    % Save
    thisTable = thisTable(thisTable.Year>=min(yearList_future) & thisTable.Year<=max(yearList_future),:) ;
    ts_co2_yr(:,r) = thisTable.co2 ;

    clear thisTable co2file thisDir tmp
end; clear r
