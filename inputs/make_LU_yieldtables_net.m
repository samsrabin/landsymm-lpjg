%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate land use files for yield table runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starts as NATURAL in transitions{1,1}-1. Then produces rows to ensure
% non-interpolated transitions to subsequent land uses in specified
% years.
transitions = { ...
    1970, 'CROPLAND' ;
    1980, {'CROPLAND','FOREST'} ;
    } ;

dir_out = '/Volumes/Reacher/LandSyMM/inputs/LU' ;


%% Setup

Ntrans = size(transitions,1) ;

LU_list = {'NATURAL'} ;
for t = 1:Ntrans
    tmp = transitions{t,2} ;
    if iscell(tmp)
        LU_list = [LU_list tmp] ; %#ok<AGROW>
    else
        LU_list = [LU_list {tmp}] ; %#ok<AGROW>
    end
end
LU_list = unique(LU_list, 'stable') ;
Nlu = length(LU_list) ;

yearList = ((transitions{1,1}-1):transitions{end,1})' ;
Nyears = length(yearList) ;


%% Generate LU array

LU_yv = zeros(Nyears,Nlu) ;
Nlu_now_max = -Inf ;
for t = 0:Ntrans
    if t==0
        LU_yv(1, strcmp(LU_list, 'NATURAL')) = 1 ;
    else
        LU_list_now = transitions{t,2} ;
        [~, Ilu] = intersect(LU_list, LU_list_now) ;
        if ~iscell(LU_list_now)
            LU_list_now = {LU_list_now} ;
        end
        Nlu_now = length(LU_list_now) ;
        Nlu_now_max = max(Nlu_now, Nlu_now_max) ;
        if length(Ilu) ~= Nlu_now
            error('Error finding land uses: Found %d, expected %d', ...
                length(Ilu), Nlu_now)
        end
        if t==Ntrans
            IA = Nyears ;
        else
            [~,IA] = intersect(yearList,transitions{t,1}:(transitions{t+1,1})-1) ;
        end
        LU_yv(IA, Ilu) = 1/Nlu_now ;
    end
end

data_YV = cat(2, yearList, LU_yv) ;
% array2table(data_YV, 'VariableNames', [{'Year'} LU_list'])


%% Save LU file

% Get output filename
file_out = sprintf('%s/landuse-totals', dir_out) ;
sep = '.' ;
for t = 1:Ntrans
    thisLU = transitions{t,2} ;
    if iscell(thisLU)
        thisLetter = '' ;
        for x = 1:length(thisLU)
            tmp = thisLU{x} ;
            thisLetter = [thisLetter lower(tmp(1))] ;
        end
    else
        thisLetter = lower(thisLU(1)) ;
    end
    file_out = sprintf('%s%s%s%d', ...
        file_out, sep, thisLetter, transitions{t,1}) ;
    sep = '-' ;
end
file_out = sprintf('%s.txt', ...
    file_out) ;

% Get output header and data-row format
header_out = 'Year' ;
for L = 1:Nlu
    header_out = [header_out '\t' LU_list{L}] ; %#ok<AGROW>
end
header_out = [header_out '\n'] ;
if Nlu_now_max==1
    tmp = '%d' ;
elseif Nlu_now_max==2
    tmp = '%0.1f' ;
else
    tmp = '%0.6f' ;
end
format_out = ['%d' repmat(['\t' tmp], [1 Nlu]) '\n'] ;

% Open file and save header
fid = fopen(file_out,'w') ;
fprintf(fid, header_out) ;

% Write data
fprintf(fid, ...
    format_out, ...
    data_YV' ...
    ) ;
fclose(fid) ;
disp('Done saving')


