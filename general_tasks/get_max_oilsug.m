%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Combine Oilcrops and combine Sugars %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% topDir = '/Users/Shared/SAI-LandSyMM/output/remap10/2022-11-05.halfdeg.20221111' ;
% outDir_top = '/Users/Shared/SAI-LandSyMM/output/remap10/2022-11-05.halfdeg.maxOilSug.20221114' ;

% topDir = '/Users/Shared/SAI-LandSyMM/output/remap10/2022-11-23.halfdeg' ;
% outDir_top = '/Users/Shared/SAI-LandSyMM/output/remap10/2022-11-23.halfdeg.maxOilSug' ;

% topDir = '/Users/Shared/SAI-LandSyMM/output/remap10/outForPLUM-2022-11-30-090159.halfdeg' ;
% outDir_top = '/Users/Shared/SAI-LandSyMM/output/remap10/outForPLUM-2022-11-30-090159.halfdeg.maxOilSug' ;

topDir = '/Users/Shared/SAI-LandSyMM/output/remap10/outForPLUM-2022-12-02-221655.halfdeg' ;
outDir_top = '/Users/Shared/SAI-LandSyMM/output/remap10/outForPLUM-2022-12-02-221655.halfdeg.maxOilSug' ;


%% Setup

cd(topDir)
dirList = dir('.') ;

file_gridlist_out = '/Users/Shared/LandSyMM/inputs/gridlists/gridlist_62892.runAEclimOK.txt' ;
gridlist_out = lpjgu_matlab_read2geoArray(file_gridlist_out, ...
    'verboseIfNoMat', false, 'force_mat_save', false, 'force_mat_nosave', true) ;

combinations = {'Sugar', {'Sugarbeet', 'Sugarcane'}, [10.676 8.777] ;
                'Oilcrops', {'OilNfix', 'OilOther'}, [0.448 0.508]} ;
Nlevels = {'0', '0060', '0200', '1000'} ;
irrigs = {'', 'i'} ;

other_files = {'done', 'runinfo_pot.tar', 'runinfo_act.tar', 'anpp.out', 'tot_runoff.out'} ;

if ~exist(outDir_top, 'dir')
    mkdir(outDir_top) ;
end

%%% Options %%%%%%%%%%%%%
outPrec = 3 ;
outPrec_lonlat = 2 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = false ;
fancy = false ;
do_gzip = true ;
%%%%%%%%%%%%%%%%%%%%%%%%%


%% Process

for d = 1:length(dirList)
    thisDir = dirList(d).name ;
    if ~dirList(d).isdir || thisDir(1) == '.'
        continue
    end
    
    thisFile_yield = sprintf('%s/%s/yield.out', topDir, thisDir) ;
    thisFile_irrig = sprintf('%s/%s/gsirrigation.out', topDir, thisDir) ;

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
        thisOtherFile = sprintf('%s/%s/%s*', topDir, thisDir, other_files{x}) ;
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


