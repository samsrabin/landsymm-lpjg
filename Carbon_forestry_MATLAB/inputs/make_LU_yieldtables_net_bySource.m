%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate land use files for yield table runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_lu = {'NATURAL', 'CROPLAND', 'PASTURE', 'FOREST'} ;
list_forest = {'forC', 'forT'} ;

% What year is the last in the time horizon of the project?
dest_yN = 2010 ;

% How long is the project time horizon?
dest_duration = 30 ;

% How many years does LU stay as source before transition?
% This should always be the same as dest_duration, to ensure that source
% forestry is cut immediately before transition (as long as you have
% manageforest_dec31 enabled).
source_duration = dest_duration ;

dir_out = '/Volumes/Reacher/LandSyMM/inputs/LU_bySource' ;


%% Setup

Nlu = length(list_lu) ;
Nforest = length(list_forest) ;

list_stands = [setdiff(list_lu, {'FOREST'}, 'stable') list_forest] ;
Nstands = length(list_stands) ;

if ~exist(dir_out, 'dir')
    mkdir(dir_out)
end

% Get yearList
dest_y1 = dest_yN - dest_duration ;
yearList = (dest_y1-source_duration-1):dest_y1 ;
Nyears = length(yearList) ;


%% Process all

for s = 1:Nstands
    source = list_stands{s} ;
    
    fprintf('%s...\n', source)
    
    source_is_forest = any(strcmp(list_forest, source)) ;
    
    % Set up array for file_lu
    LU_yv = zeros(Nyears, Nlu) ;

    % Fill first year (year before transition to source) with 100% NATURAL
    LU_yv(1, strcmp(list_lu, 'NATURAL')) = 1 ;
    
    % Fill source years with 100% source
    if ~source_is_forest
        LU_yv(2:end-1, strcmp(list_lu, source)) = 1 ;
    else
        LU_yv(2:end-1, strcmp(list_lu, 'FOREST')) = 1 ;
    end
    
    % Split last year among all LUs
    LU_yv(end,:) = 1 / Nlu ;
    
    % Set up array for file_luforest
    if ~source_is_forest
        forest_yv = 0.5 * ones(2+source_duration, Nforest) ;
    else
        forest_yv = zeros(2+source_duration, Nforest) ;
        forest_yv(:,strcmp(list_forest, source)) = 1 ;
        forest_yv(end,:) = 0.5 ;
    end
    
    % Sanity checks
    if any(sum(LU_yv,2) ~= 1)
        error('Some row of LU_yv does not sum to 1')
    elseif any(sum(forest_yv,2) ~= 1)
        error('Some row of forest_yv does not sum to 1')
    end
    
    % Get filename token for source LU
    if ~source_is_forest
        switch source
            case 'CROPLAND'
                thisLUtoken = 'crop' ;
            case 'PASTURE'
                thisLUtoken = 'past' ;
            case 'NATURAL'
                thisLUtoken = 'ntrl' ;
            otherwise
                error('source %s not recognized', source)
        end
    else
        thisLUtoken = source ;
    end
    
    % Save files
    file_out_mostly = sprintf('%s/source_%s_%d-%d_', ...
        dir_out, thisLUtoken, yearList(2), yearList(end-1)) ;
    save_file([file_out_mostly 'lu.txt'], LU_yv, yearList, list_lu)
    save_file([file_out_mostly 'luforest.txt'], forest_yv, yearList, list_forest)
    
end

disp('Done.')


%% FUNCTION

function save_file(file_out, fractions_yv, yearList, varNames)

data_YV = cat(2, shiftdim(yearList), fractions_yv) ;

% Get output header and data-row format
header_out = 'Year' ;
Nvar = length(varNames) ;
for v = 1:Nvar
    header_out = [header_out '\t' varNames{v}] ; %#ok<AGROW>
end
header_out = [header_out '\n'] ;
format_out = ['%d' repmat('\t%0.6f', [1 Nvar]) '\n'] ;

% Open file and save header
fid = fopen(file_out,'w') ;
fprintf(fid, header_out) ;

% Write data
fprintf(fid, ...
    format_out, ...
    data_YV' ...
    ) ;
fclose(fid) ;

end
