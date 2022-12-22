%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate land use files for yield table runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starts as NATURAL in transitions{1,1}-1. Then produces rows to ensure
% non-interpolated transitions to subsequent land uses in specified
% years.
transitions = { ...
    1970, {'CROPLAND','PASTURE'} ;
    1980, {'CROPLAND','PASTURE','FOREST'} ;
    } ;
transitions = { ...
    1970, {'CROPLAND'} ;
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



%%

LU_yv = [] ;
Nlu_now_max = -Inf ;

% Generate all possible transitions
letters = cellfun(@getletter, LU_list) ;
grossLUC_list = cellstr(cat(1,nchoosek(letters, 2),nchoosek(flip(letters), 2))) ;
NgrossLUC = length(grossLUC_list) ;

% Generate transition matrix
LU_list_prev = {'NATURAL'} ;
LU_fracs = zeros(size(LU_list)) ;
LU_fracs(strcmp(LU_list, 'NATURAL')) = 1 ;
transition_matrix = zeros(Ntrans, NgrossLUC) ;
for t = 1:Ntrans
    fprintf('transition %d\n', t)
    Nlu_prev = length(LU_list_prev) ;
    LU_list_now = transitions{t,2} ;
    Nlu_now = length(LU_list_now) ;
    for l1 = 1:Nlu_prev
        % Donor land use
        LU1 = LU_list_prev{l1} ;
        letter1 = getletter(LU1) ;
        I1 = find(strcmp(LU_list, LU1)) ;
        LU1_frac_orig = LU_fracs(I1) ;
        
%         % Skip if nothing to give
%         if LU_fracs(I1) == 0
%             continue
%         end
        
        for l2 = 1:Nlu_now
            % Recipient land use
            LU2 = LU_list_now{l2} ;
            letter2 = getletter(LU2) ;
            I2 = find(strcmp(LU_list, LU2)) ;
            
            % Skip if it's the same land use
            if strcmp(LU1, LU2)
                continue
            end
            
            LU_fracs
            
            % Get fraction of donor that goes to this recipient
            frac = LU1_frac_orig / Nlu_now ;
            
            % Get corresponding column
            thisLUC = [letter1 letter2] ;
            I = find(strcmp(grossLUC_list, thisLUC)) ;
            if length(I) ~= 1
                error('Found %d columns matching %s (expected 1)', ...
                    length(I), thisLUC)
            end
            
            % Save to transition matrix
            transition_matrix(t,I) = frac ;
            
            % Track current state
            LU_fracs(I1) = LU_fracs(I1) - frac ;
            LU_fracs(I2) = LU_fracs(I2) + frac ;
            if abs(sum(LU_fracs) - 1) > 1e-12
                warning('LU_fracs does not sum to 1')
            elseif any(LU_fracs < 0)
                warning('Negative(s) in LU_fracs')
            elseif any(LU_fracs > 1)
                warning('>1 in LU_fracs')
            end
        end
    end
    LU_fracs
%     stop
    LU_list_prev = LU_list_now ;
end
clear LU_list_prev LU_list_now
LU_fracs
transition_matrix
error('You need to expand transition matrix out to include every year, not just transition years.')

%     [~, Ilu] = intersect(LU_list, LU_list_now) ;
%     if ~iscell(LU_list_now)
%         LU_list_now = {LU_list_now} ;
%     end
%     Nlu_now = length(LU_list_now) ;
%     Nlu_now_max = max(Nlu_now, Nlu_now_max) ;
%     if length(Ilu) ~= Nlu_now
%         error('Error finding land uses: Found %d, expected %d', ...
%             length(Ilu), Nlu_now)
%     end
%     if t==Ntrans
%         IA = Nyears ;
%     else
%         [~,IA] = intersect(yearList,transitions{t,1}:(transitions{t+1,1})-1) ;
%     end
%     LU_yv(IA, Ilu) = 1/Nlu_now ;
% end

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


%% FUNCTIONS

function thisletter = getletter(name)

if strcmp(name, 'NATURAL')
    thisletter = 'v' ;
else
    thisletter = lower(name(1)) ;
end

end