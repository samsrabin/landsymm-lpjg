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


%% Get gridlist

% Import previous LandSyMM gridlist and ISIMIP3 gridlist
gridlist_prev = lpjgu_matlab_read2geoArray('/Users/Shared/PLUM/input/gridlists/gridlist_62892.runAEclimOK.txt') ;
gridlist_isimip = lpjgu_matlab_read2geoArray('/Volumes/Reacher/GGCMI/inputs/phase3/ISIMIP3/landseamask-lpjg/gridlist_ggcmi_v1.1.gapfilled.lpjg.txt') ;

% Merge
gridlist.mask_YX = gridlist_prev.mask_YX | gridlist_isimip.mask_YX ;
gridlist.list2map = union(gridlist_prev.list2map, gridlist_isimip.list2map, 'stable') ;
gridlist.lonlats = union(gridlist_prev.lonlats, gridlist_isimip.lonlats, 'rows', 'stable') ;
clear gridlist_prev gridlist_isimip
Ncells = length(gridlist.list2map) ;


%% Generate LU array for a single cell

LU_list = unique([{'NATURAL'} ; transitions(:,2)], 'stable') ;
Nlu = length(LU_list) ;

yearList = nan(2*Ntrans,1) ;
Nyears = length(yearList) ;

LU_yv = zeros(Nyears,Nlu) ;
for t = 1:Ntrans
    y1 = (t-1)*2+1 ;
    y2 = (t-1)*2+2 ;
    
    if t==1
        LU_yv(y1, strcmp(LU_list, 'NATURAL')) = 1 ;
    else
        LU_yv(y1, strcmp(LU_list, transitions{t-1,2})) = 1 ;
    end
    LU_yv(y2, strcmp(LU_list, transitions{t,2})) = 1 ;
    
    yearList(y1) = transitions{t,1} - 1 ;
    yearList(y2) = transitions{t,1} ;
end

LU_yv(1, strcmp(LU_list, 'NATURAL')) = 1 ;
LU_yv(2:3, strcmp(LU_list, 'CROPLAND')) = 1 ;
LU_yv(4, strcmp(LU_list, 'FOREST')) = 1 ;

data_YV = cat(2, yearList, LU_yv) ;


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
file_out = sprintf('%s.%d.txt', ...
    file_out, Ncells) ;

% Get output header and data-row format
header_out = 'Lon\tLat\tYear' ;
for L = 1:Nlu
    header_out = [header_out '\t' LU_list{L}] ; %#ok<AGROW>
end
header_out = [header_out '\n'] ;
format_out = ['%0.2f\t%0.2f\t%d' repmat('\t%d', [1 Nlu]) '\n'] ;

% Open file and save header
fid = fopen(file_out,'w') ;
fprintf(fid, header_out) ;

% Write to file gridcell-by-gridcell
for g = 1:Ncells
    thisInd = gridlist.list2map(g) ;
    out_YV = cat(2, ...
        repmat(gridlist.lonlats(g,:), [Nyears 1]), ...
        data_YV) ;
    fprintf(fid, ...
        format_out, ...
        out_YV' ...
        ) ;
    if rem(g,10000)==0 || g==Ncells
        fclose(fid) ;
        disp(num2str(g))
        if g<Ncells
            fid = fopen(file_out,'a') ;
        end
    end
end

disp('Done')