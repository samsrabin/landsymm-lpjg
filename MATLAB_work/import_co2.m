%%%%%%%%%%%%%%%%%%
%%% Import CO2 %%%
%%%%%%%%%%%%%%%%%%


disp('Importing CO2...')

insfile = 'main_plum_lpjg.ins' ;

% Baseline
thisInsfile = sprintf('%s/%s', ...
    baselineDir, insfile) ;
if ~exist(thisInsfile, 'file')
    insfile = 'main.ins' ;
    thisInsfile = sprintf('%s/%s', ...
        baselineDir, insfile) ;
    if ~exist(thisInsfile, 'file')
        error('What ins-file should I look for file_co2 in?')
    end
end
[status, result] = unix(sprintf( ...
    'grep -he "^param \\"file_co2" %s', ...
    thisInsfile)) ;
if status~=0
    error('Error in unix() call')
end
tmp = strsplit(result,'"') ;
co2file = tmp{4} ;
co2file = get_local_co2file(co2file) ;
if ~exist(co2file, 'file')
    error('co2file (%s) not found!', co2file)
end
thisTable_bl = read_and_process_table(co2file);

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
    thisInsfile = sprintf('%s/%s', ...
        thisDir, insfile) ;
    if ~exist(thisInsfile, 'file')
        error('What ins-file should I look for file_co2 in?')
    end
    [status, result] = unix(sprintf('grep -he "^param \\"file_co2" %s', thisInsfile)) ;
    if status~=0
        error('Error in unix() call')
    end
    tmp = strsplit(result,'"') ;
    co2file = tmp{4} ;
    co2file = get_local_co2file(co2file) ;
    if ~exist(co2file, 'file')
        error('co2file (%s) not found!', co2file)
    end
    thisTable = read_and_process_table(co2file);
    
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


function co2file = get_local_co2file(co2file)

co2file = strrep(co2file,'/home/fh1-project-lpjgpi/lr8247','/Users/Shared/lpj-guess') ;
co2file = strrep(co2file,'/home/kit/imk-ifu/lr8247/input/lpj-guess','/Users/Shared/lpj-guess/input') ;

end


function table_out = read_and_process_table(co2file)

table_in = readtable(co2file, 'Delimiter', ' ', ...
        'ConsecutiveDelimitersRule', 'join', ...
        'LeadingDelimitersRule', 'ignore') ;
array_in = table2array(table_in) ;

array_out = array_in ;

if any(any(isnan(array_in)))
    nrows = size(array_in, 1) ;
    ncols = size(array_in, 2) ;
    
    bad_row = sum(isnan(array_in), 2) == ncols ;
    bad_col = sum(isnan(array_in), 1) == nrows ;
    
    array_out(bad_row,:) = [] ;
    array_out(:,bad_col) = [] ;
end

table_out = array2table(array_out, 'VariableNames', {'Year','co2'}) ;

end
