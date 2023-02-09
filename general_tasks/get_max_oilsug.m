%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Combine Oilcrops and combine Sugars %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSR 2022-12-22
% - Picks a "winner" of the constituent crops at each timestep in each
%   gridcell. For now, it does this based on which has a higher yield
%   (after applying calibration factors) at the lowest N fertilization
%   level, rainfed. One could imagine picking a winner at EACH
%   fertilizer x irrigation level instead, but I'm not sure if that
%   would make things tricky for PLUM's curve fitting. Probably worth
%   a try.
% - "Burns in" calibration factors for output Oilcrops and Sugar so that PLUM
%   doesn't have to worry about them. PLUM should thus assume calibration
%   factors of 1 for these crops.
% - The output directory, outDir_top, is what you'll be sending to PLUM.

addpath(genpath(landsymm_lpjg_path()))

% get_max_oilsug_options.m must be somewhere on your path.
% There, specify the following variables:
%    topDir: Path to "outForPLUM" directory.
%    combinations: Info about what crops to combine and how. Cell array with rows like
%                      'Combinedcrop', {'Subcrop1', 'Subcrop2'}, [calib_factor1 calib_factor2]
%                  E.g.: combinations = ...
%                         {'Sugar',    {'Sugarbeet', 'Sugarcane'}, [10.676 8.777] ;
%                          'Oilcrops', {'OilNfix', 'OilOther'},    [0.448 0.508]} ;
%    Nlevels: Text versions of the fertilization levels.
%             E.g.: Nlevels = {'0', '0060', '0200', '1000'} ;
%    irrigs: Text versions of the suffixes used to indicate irrigation type.
%            E.g.: irrigs = {'', 'i'} ;
%  OPTIONAL: OUTPUT FILES
%    delimiter: The column delimiter in output files.
%               E.g.: delimiter = ' ' ;
%    do_gzip: Whether to zip up output files or leave them as plain text files.
%             E.g.: do_gzip = true ;
%    outDir_top: Directory where the outputs will go. If not provided, uses
%                [topDir '.maxOilSug']
%    outPrec: Precision of data in output files. 
%             E.g.: outPrec = 3 ;
%    outPrec_lonlat: Precision of lon/lat values in output files.
%                    E.g., outPrec_lonlat = 2 ;
%    overwrite: Whether to overwrite existing output files.
%               E.g., overwrite = false ;


get_max_oilsug_options


%% Setup

cd(topDir)
dirList = dir('.') ;

other_files = {'done', 'runinfo_pot.tar', 'runinfo_act.tar', 'anpp.out', 'tot_runoff.out'} ;

% Process optional settings
if ~exist('outPrec', 'var')
    outPrec = 3 ;
end
if ~exist('outPrec_lonlat', 'var')
    outPrec_lonlat = 2 ;
end
if ~exist('delimiter', 'var')
    delimiter = ' ' ;
end
if ~exist('overwrite', 'var')
    overwrite = true ;
end
if ~exist('outPrec', 'var')
    outPrec = 3 ;
end
if ~exist('do_gzip', 'var')
    do_gzip = true ;
end
outWidth = 1 ;
fancy = false ;

% Get directories
topDir = strip(topDir, 'right', filesep) ;
if ~exist('outDir_top', 'var')
    outDir_top = [topDir '.maxOilSug'] ;
end
if ~exist(outDir_top, 'dir')
    mkdir(outDir_top) ;
end
fprintf('Saving to %s\n', outDir_top)


%% Process

for d = 1:length(dirList)
    thisDir = dirList(d).name ;
    if ~dirList(d).isdir || thisDir(1) == '.'
        continue
    end
    
    thisFile_yield = fullfile(topDir, thisDir, 'yield.out') ;
    thisFile_irrig = fullfile(topDir, thisDir, 'gsirrigation.out') ;

    if isempty(dir(sprintf('%s*', thisFile_yield)))
        fprintf('Skipping %s (thisFile_yield not found)\n', thisDir)
        continue
    end
    

    outDir = sprintf('%s/%s', outDir_top, thisDir) ;
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end
    outFile_yield = strrep(thisFile_yield, topDir, outDir_top) ;
    outFile_irrig = strrep(thisFile_irrig, topDir, outDir_top) ;

    if ~overwrite
        outFile_yield_test = outFile_yield ;
        outFile_irrig_test = outFile_irrig ;
        if do_gzip
            outFile_yield_test = [outFile_yield_test '.gz'] ; %#ok<AGROW>
            outFile_irrig_test = [outFile_irrig_test '.gz'] ; %#ok<AGROW>
        end
        if exist(outFile_yield_test, 'file') && exist(outFile_irrig_test, 'file')
            fprintf('    Skipping %s (max_oilsug files already exist)\n', thisDir)
            continue
        end
    end

    fprintf('%s...\n', thisDir)

    % Copy extra files
    for x = 1:length(other_files)
        thisOtherFile = fullfile(topDir, thisDir, other_files{x}) ;
        if ~exist(thisOtherFile, 'file')
            thisOtherFile2 = dir([thisOtherFile '*']) ;
            if isempty(thisOtherFile2)
                error('No file found matching %s', thisOtherFile2)
            elseif length(thisOtherFile2) > 1
                error('%s does not exist, and\n%d files found matching %s', ...
                    thisOtherFile, length(thisOtherFile2), thisOtherFile2)
            end
            thisOtherFile2 = fullfile(thisOtherFile2.folder, thisOtherFile2.name) ;
            if ~strcmp([thisOtherFile '.gz'], thisOtherFile2)
                warning('%s does not exist;\nusing %s instead', thisOtherFile, thisOtherFile2)
            end
            thisOtherFile = thisOtherFile2 ;
        end
        copyfile(thisOtherFile, outDir)
    end

    S_yield = lpjgu_matlab_read2geoArray(thisFile_yield, ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;
    Ncrops = length(S_yield.varNames) ;
    S_irrig = lpjgu_matlab_read2geoArray(thisFile_irrig, ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;

    if ~any(S_irrig.garr_xv > 0)
        error('No values > 0')
    end

    % Remove unneeded columns
    unused_columns = ~contains(S_yield.varNames, Nlevels) & ~strcmp(S_yield.varNames, 'ExtraCrop') ;
    S_yield.garr_xv(:,unused_columns) = [] ;
    S_irrig.garr_xv(:,unused_columns) = [] ;
    S_yield.varNames(unused_columns) = [] ;
    S_irrig.varNames(unused_columns) = [] ;

    if ~any(S_irrig.garr_xv > 0)
        error('No values > 0')
    end

    for x = 1:size(combinations, 1)
        mainCrop = combinations{x,1} ;
        subCrops = combinations{x,2} ;
        calib_factors = combinations{x,3} ;
        Nsub = length(subCrops) ;
        IA = contains(S_yield.varNames, subCrops) ;
        for n = 1:length(Nlevels)
            for ii = 1:length(irrigs)
                get_thisCrop = @(x) [x irrigs{ii} Nlevels{n}] ;
                thisCrop_main = get_thisCrop(mainCrop) ;
                thisCrop_sub = cellfun(get_thisCrop, subCrops, 'UniformOutput', false) ;
                [~, ~, I_sub] = intersect(thisCrop_sub, S_yield.varNames, 'stable') ;
                if length(I_sub) ~= length(thisCrop_sub)
                    tmp = sprintf('%s ', thisCrop_sub{:}) ;
                    error('Found %d matches for %d items (%s)', ...
                        length(I_sub), length(thisCrop_sub), tmp)
                end
                yield_sub_xv = S_yield.garr_xv(:,I_sub) ;
                irrig_sub_xv = S_irrig.garr_xv(:,I_sub) ;

                if ~any(S_irrig.garr_xv > 0)
                    error('No values > 0')
                end

                % Remove subcrop columns
                S_yield.garr_xv(:,I_sub) = [] ;
                S_yield.varNames(I_sub) = [] ;
                S_irrig.garr_xv(:,I_sub) = [] ;
                S_irrig.varNames(I_sub) = [] ;

                if ~any(S_irrig.garr_xv > 0)
                    error('No values > 0')
                end

                % Apply calibration factors
                for c = 1:length(subCrops)
                    thisCrop_for_cf = subCrops{c} ;
                    cf = calib_factors(c) ;
                    yield_sub_xv(:,c) = yield_sub_xv(:,c) * cf ;
                end

                % Arbitrarily: Choose between subcrops based on performance
                % at rainfed 0010
                if n==1 && ii==1
                    [~, whichMax_sub_x] = max(yield_sub_xv, [], 2) ;
                end
                yield_main_x = nan(size(yield_sub_xv, 1), 1) ;
                irrig_main_x = nan(size(yield_sub_xv, 1), 1) ;
                for s = 1:Nsub
                    thisSubIsMax = whichMax_sub_x == s ;
                    yield_main_x(thisSubIsMax) = yield_sub_xv(thisSubIsMax,s) ;
                    irrig_main_x(thisSubIsMax) = irrig_sub_xv(thisSubIsMax,s) ;
                end
                S_yield.varNames{end+1} = thisCrop_main ;
                S_yield.garr_xv(:,end+1) = yield_main_x ;
                S_irrig.varNames{end+1} = thisCrop_main ;
                S_irrig.garr_xv(:,end+1) = irrig_main_x ;

                if ~any(S_irrig.garr_xv > 0)
                    error('No values > 0')
                end
            end
        end
        
    end

    if ~any(S_irrig.garr_xv > 0)
        error('No values > 0')
    end

    lpjgu_matlab_saveTable([{'Lon', 'Lat', S_yield.varNames}], S_yield, outFile_yield,...
        'outPrec', outPrec, ...
        'outPrec_lonlat', outPrec_lonlat, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'verbose', false, ...
        'gzip', do_gzip) ;
    clear S_yield
    lpjgu_matlab_saveTable([{'Lon', 'Lat', S_irrig.varNames}], S_irrig, outFile_irrig,...
        'outPrec', outPrec, ...
        'outPrec_lonlat', outPrec_lonlat, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'verbose', false, ...
        'gzip', do_gzip) ;
    clear S_irrig

end

disp('Done.')


