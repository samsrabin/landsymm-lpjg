%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate land use files for yield table runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starts as NATURAL in transitions{1,1}-1. Then produces rows to ensure
% non-interpolated 100% transitions to subsequent land uses in specified
% years.
transitions = { ...
    1970, 'CROPLAND' ;
    1980, 'FOREST' ;
    } ;

dir_out = '/Volumes/Reacher/LandSyMM/inputs/LU' ;


%% Setup

Ntrans = size(transitions,1) ;

LU_list = unique([{'NATURAL'} ; transitions(:,2)], 'stable') ;
Nlu = length(LU_list) ;

yearList = ((transitions{1,1}-1):transitions{end,1})' ;
Nyears = length(yearList) ;


%% Generate LU array

LU_yv = zeros(Nyears,Nlu) ;
for t = 0:Ntrans
    if t==0
        LU_yv(1, strcmp(LU_list, 'NATURAL')) = 1 ;
    elseif t==Ntrans
        LU_yv(end, strcmp(LU_list, transitions{t,2})) = 1 ;
    else
        [~,IA] = intersect(yearList,transitions{t,1}:(transitions{t+1,1})-1) ;
        LU_yv(IA, strcmp(LU_list, transitions{t,2})) = 1 ;
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
    file_out = sprintf('%s%s%s%d', ...
        file_out, sep, lower(thisLU(1)), transitions{t,1}) ;
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
format_out = ['%d' repmat('\t%d', [1 Nlu]) '\n'] ;

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


