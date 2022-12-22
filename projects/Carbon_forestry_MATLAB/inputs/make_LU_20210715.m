%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate land use files for yield table runs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use historical land use (with tiny areas of everything) until y0_past, at
% which point everything is converted to pasture (again, aside from tiny
% refugia of other LUs). After one year, the gridcell's LU is split up
% evenly among all stands.

list_lu = {'NATURAL', 'CROPLAND', 'PASTURE', 'FOREST'} ;
list_forest = {'forC', 'forT'} ;

% When is the year of conversion to 100% pasture?
y0_past = 1850 ;

% Which land use versions are we using?
remapVer = '10_g2p' ; % CerealsC3, CerealsC4, Rice, OilNfix, OilOther, 
                      % StarchyRoots, Pulses, Sugarbeet, Sugarcane, 
                      % FruitAndVeg

dir_out = '/Volumes/Reacher/LandSyMM/inputs/LU_20210715' ;


%% Setup

if ~exist(dir_out, 'dir')
    mkdir(dir_out)
end


%% Get ISIMIP3 historical land use

% Determine files to use
if contains(remapVer, '_g2p')
    LUdir = sprintf('/Volumes/Reacher/G2P/inputs/LU/remaps_v%s', ...
        remapVer) ;
else
    LUdir = sprintf('/Users/Shared/PLUM/input/remaps_v%s', ...
        remapVer) ;
end
file_lu = sprintf('%s/LU.remapv%s.txt', ...
    LUdir, remapVer);
file_cf = sprintf('%s/cropfracs.remapv%s.txt', ...
    LUdir, remapVer);

% Import and trim any unneeded years
lu_in = import_lu(y0_past, file_lu) ;
cf_in = import_lu(y0_past, file_cf) ;
Nyears = length(lu_in.yearList) ;

% Move any unneeded LU types into NATURAL
[~, I] = setdiff(lu_in.varNames, list_lu) ;
if ~isempty(I)
    lu_in.garr_xvy(:,strcmp(lu_in.varNames, 'NATURAL'),:) = ...
        lu_in.garr_xvy(:,strcmp(lu_in.varNames, 'NATURAL'),:) ...
        + sum(lu_in.garr_xvy(:,I,:), 2) ;
    lu_in.garr_xvy(:,I,:) = [] ;
    lu_in.varNames(I) = [] ;
end

% Check
if ~isequal(cf_in.lonlats, lu_in.lonlats)
    error('Need to align gridlists of land use and cropfrac files')
end

Ncells = length(lu_in.list2map) ;

% Trim any unneeded crop types
is_unneeded = nansum(cf_in.garr_xv, 1) == 0 ;
if any(is_unneeded)
    cf_in.garr_xv(:,is_unneeded) = [] ;
    cf_in.varNames(is_unneeded) = [] ;
end


%% Process

% Ensure that every cell has at least some crop
outPrec = 6 ;
minCrop = 10^-outPrec ;
lu_out = lu_in ;
cf_out = cf_in ;
noCrop_x1y = lu_out.garr_xvy(:,strcmp(lu_out.varNames, 'CROPLAND'),:) < minCrop ;
never_cropped = all(noCrop_x1y, 3) ;
cf_out.garr_xv(never_cropped,:) = 0 ;
cf_out.garr_xv(never_cropped, strcmp(cf_out.varNames, 'CerealsC3')) = 1 ;
if any(noCrop_x1y(:))
    lu_out = donate_land_to_crop(lu_out, noCrop_x1y, 'PASTURE', minCrop) ;
    noCrop_x1y = lu_out.garr_xvy(:,strcmp(lu_out.varNames, 'CROPLAND'),:) < minCrop ;
    if any(noCrop_x1y(:))
        lu_out = donate_land_to_crop(lu_out, noCrop_x1y, 'NATURAL', minCrop) ;
        noCrop_x1y = lu_out.garr_xvy(:,strcmp(lu_out.varNames, 'CROPLAND'),:) < minCrop ;
        if any(noCrop_x1y(:))
            error('There are still %d cell-years with no cropland', length(find(noCrop_x1y)))
        end
    end
end

% Add zeros for missing LUs
missing_LUs = setdiff(list_lu, lu_out.varNames) ;
if ~isempty(missing_LUs)
    lu_out.garr_xvy = cat(2, lu_out.garr_xvy, ...
        zeros(Ncells, length(missing_LUs), Nyears)) ;
    lu_out.varNames = [lu_out.varNames missing_LUs] ;
end

% Add 1 year of pasture, then all LUs
if Nyears > 0
    lu_out.yearList(end+1) = lu_out.yearList(end) + 1 ;
else
    lu_out.yearList(1) = y0_past ;
end
lu_out.garr_xvy(:,:,Nyears+1) = 0 ;
lu_out.garr_xvy(:,strcmp(lu_out.varNames,'PASTURE'),end) = 1 ;
Nlu_out = length(list_lu) ;
Nforest = length(list_forest) ;
thisDenom = Nlu_out - 1 + Nforest ;
lu_out.yearList(end+1) = lu_out.yearList(end) + 1 ;
lu_out.garr_xvy(:,:,end+1) = 1 / thisDenom ;
lu_out.garr_xvy(:,strcmp(lu_out.varNames,'FOREST'),end) = ...
    Nforest / thisDenom ; 


%% Save files

overwrite = true ;
outWidth = 1 ;
delimiter = ' ' ;
fancy = false ;

% Get output filenames
file_lu_out = sprintf('%s/remap%s_someCrop_%dpast_%dall_LU.txt', ...
    dir_out, strrep(remapVer, '_', ''), y0_past, y0_past+1) ;
file_cf_out = sprintf('%s/remap%s_someCrop_%dpast_%dall_cropfrac.txt', ...
    dir_out, strrep(remapVer, '_', ''), y0_past, y0_past+1) ;

% Get headers
lu_out_header = [{'Lon', 'Lat', 'Year'}, lu_out.varNames] ;
cf_out_header = [{'Lon', 'Lat'}, cf_out.varNames] ;

disp('Saving LU...')
lpjgu_matlab_saveTable(lu_out_header, lu_out, file_lu_out, ...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy) ;

disp('Saving cropfrac...')
lpjgu_matlab_saveTable(cf_out_header, cf_out, file_cf_out, ...
    'outPrec', outPrec, ...
    'outWidth', outWidth, ...
    'delimiter', delimiter, ...
    'overwrite', overwrite, ...
    'fancy', fancy) ;



%% FUNCTIONS

function lu = donate_land_to_crop(lu, noCrop_x1y, donor, minCrop)

% Sanity checks
if any(any(any(lu.garr_xvy < 0)))
    error('Negative value(s) in lu.garr_xvy (worst: %g)', min(min(min(lu.garr_xvy))))
elseif any(any(any(isnan(lu.garr_xvy))))
    error('NaN(s) in lu.garr_xvy')
end

I_crop = find(strcmp(lu.varNames, 'CROPLAND')) ;
if length(I_crop) ~= 1
    error('Expected 1 I_crop; found %d', length(I_crop))
end

I_donor = find(strcmp(lu.varNames, donor)) ;
if length(I_donor) ~= 1
    error('Expected 1 I_donor; found %d', length(I_donor))
end

lu_xCy = lu.garr_xvy(:,I_crop,:) ;
lu_xDy = lu.garr_xvy(:,I_donor,:) ;

okDonor_x1y = lu_xDy >= 2*minCrop ;
noCrop_okDonor_x1y = noCrop_x1y & okDonor_x1y ;

Nbad_before = length(find(noCrop_x1y)) ;

if any(noCrop_okDonor_x1y(:))

    lu_xDy(noCrop_okDonor_x1y) = ...
        lu_xDy(noCrop_okDonor_x1y) - minCrop ;
    lu_xCy(noCrop_okDonor_x1y) = ...
        lu_xCy(noCrop_okDonor_x1y) + minCrop ;
    
    Nbad_after = length(find(lu_xCy < minCrop)) ;
    
    msg = sprintf('# of no-crop cell-years: %d -> %d', ...
        Nbad_before, ...
        Nbad_after) ;
    
    if Nbad_after > Nbad_before
        error(msg) %#ok<SPERR>
    end
    
    fprintf('%s\n', msg)
    
    lu.garr_xvy(:,I_donor,:) = lu_xDy ;
    lu.garr_xvy(:,I_crop,:) = lu_xCy ;

end

% Sanity checks
if any(any(any(lu.garr_xvy < 0)))
    error('Negative value(s) in lu.garr_xvy (worst: %g)', min(min(min(lu.garr_xvy))))
elseif any(any(any(isnan(lu.garr_xvy))))
    error('NaN(s) in lu.garr_xvy')
end

end


function S = import_lu(y0_past, filename)

% Import
S = lpjgu_matlab_read2geoArray(filename) ;

% Restrict to needed years
if isfield(S, 'yearList')
    S.garr_xvy(:,:,S.yearList >= y0_past) = [] ;
    S.yearList(S.yearList >= y0_past) = [] ;
end

end