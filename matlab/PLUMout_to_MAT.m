%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert a PLUM run's outputs to MAT files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inDir = '/Volumes/Reacher/G2P/outputs_PLUM/ggcmi/gepic_emu/s1_26' ;

fileList = {'landuse', 'cropfracs', 'fert', 'irrig', 'yield'} ;


%% Setup

% Set up output directory
outDir = [inDir '-matlab'] ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

% Get years
D = dir(inDir) ;
D_names = {D.name} ;
yearList = cellfun(@str2num, D_names(~cellfun(@isempty, regexp(D_names, ...
    '\d\d\d\d')))) ;
clear D D_names
y1 = yearList(1) ;
yN = yearList(end) ;
yStep = unique(yearList(2:end) - yearList(1:end-1)) ;
if length(yStep) ~= 1
    error('Inconsistent time step')
end
if ~isequal(yearList, y1:yStep:yN)
    error('~isequal(yearList, y1:yStep:yN). Missing year(s)?')
end
Nyears = length(yearList) ;


%% Import

% Import
fprintf('Importing:') % Line break happens later
for y = 1:Nyears
    thisYear = yearList(y) ;
    thisDir = sprintf('%s/%d', inDir, thisYear) ;
    
    if y==1
        fprintf('\n    ')
    elseif rem(y-1,10)==0
        fprintf('\n    ')
    end
    fprintf('%d... ', thisYear) ;
    
    % Import file
    detailed_tmp = lpjgu_matlab_read2geoArray( ...
        sprintf('%s/LandUse.txt', thisDir), ...
        'force_mat_save', false, 'force_mat_nosave', true, 'verboseIfNoMat', false) ;
    
    % Get info as needed
    if y==1
        
        % Get LU list
        if any(strcmp(fileList, 'landuse'))
            isincl = get_incl_lu(detailed_tmp.varNames) ;
            luList_in = detailed_tmp.varNames(isincl) ;
            Nlu_in = length(luList_in) ;
            nonnat_list = {'cropland', 'pasture', 'barren', 'urban'} ;
            [~, nonnat_inds_in] = intersect(luList_in, nonnat_list) ;
            [~, natural_inds_in] = setdiff(luList_in, nonnat_list) ;
            if isempty(natural_inds_in)
                error('No natural land covers found')
            end
            luList_out = luList_in ;
            luList_out{natural_inds_in(1)} = 'natural' ;
            if length(natural_inds_in) > 1
                luList_out(natural_inds_in(2:end)) = [] ;
            end
            [~, nonnat_inds_out] = intersect(luList_out, nonnat_list) ;
            [~, natural_inds_out] = setdiff(luList_out, nonnat_list) ;
            if length(natural_inds_out) ~= 1
                error('Error getting natural_inds_out: %d found (expected 1)', ...
                    length(natural_inds_out))
            end
            Nlu_out = length(luList_out) ;
        end
        
        % Get crop list
        if ~isempty(intersect(fileList, {'cropfracs', 'fert', 'irrig'}))
            cropList = detailed_tmp.varNames( ...
                get_incl_suffix(detailed_tmp.varNames, '_A', [])) ;
            cropList = strrep(cropList, '_A', '') ;
            Ncrops = length(cropList) ;
        end
        
    end
    
    % Establish target arrays
    if y==1
        lonlats = detailed_tmp.lonlats ;
        list2map = detailed_tmp.list2map ;
        Ncells = length(list2map) ;
        if any(strcmp(fileList, 'landuse'))
            landuse = setup_struct(lonlats, list2map, luList_out, yearList) ;
            landuse.garr_xvy = nan(Ncells, Nlu_out, Nyears) ;
        end
        if any(strcmp(fileList, 'cropfracs'))
            cropfracs = setup_struct(lonlats, list2map, cropList, yearList) ;
            cropfracs.garr_xvy = nan(Ncells, Ncrops, Nyears) ;
        end
        if any(strcmp(fileList, 'irrig'))
            irrig = setup_struct(lonlats, list2map, cropList, yearList) ;
            irrig.garr_xvy = nan(Ncells, Ncrops, Nyears) ;
        end
        if any(strcmp(fileList, 'fert'))
            fert = setup_struct(lonlats, list2map, cropList, yearList) ;
            fert.garr_xvy = nan(Ncells, Ncrops, Nyears) ;
        end
        if any(strcmp(fileList, 'yield'))
            yield = setup_struct(lonlats, list2map, cropList, yearList) ;
            yield.garr_xvy = nan(Ncells, Ncrops, Nyears) ;
        end
    else
        if ~isequal(lonlats, detailed_tmp.lonlats)
            error('lonlats mismatch')
        end
    end
    
    % Save this year's land uses
    if any(strcmp(fileList, 'landuse'))
        isincl = get_incl_lu(detailed_tmp.varNames) ;
        if length(find(isincl)) ~= Nlu_in
            error('LU mismatch')
        end
        [~, ~, IB] = intersect(luList_in, ...
            detailed_tmp.varNames(isincl), ...
            'stable') ;
        if length(IB) ~= Nlu_in
            error('LU mismatch')
        end
        tmp_xv = detailed_tmp.garr_xv(:,isincl) ;
        tmp_xv = tmp_xv(:,IB) ;
        landuse.garr_xvy(:,natural_inds_out,y) = ...
            sum(tmp_xv(:,natural_inds_in),2) ;
        landuse.garr_xvy(:,nonnat_inds_out,y) = ...
            tmp_xv(:,nonnat_inds_in) ;
        clear tmp_xv
    end
    
    % Save this year's crop fractions
    if any(strcmp(fileList, 'cropfracs'))
        [isincl, IB] = get_incl_suffix(detailed_tmp.varNames, '_A', cropList) ;
        tmp_xv = detailed_tmp.garr_xv(:,isincl) ;
        cropfracs.garr_xvy(:,:,y) = tmp_xv(:,IB) ;
    end
    
    % Save this year's irrigation
    if any(strcmp(fileList, 'irrig'))
        [isincl, IB] = get_incl_suffix(detailed_tmp.varNames, '_II', cropList) ;
        tmp_xv = detailed_tmp.garr_xv(:,isincl) ;
        irrig.garr_xvy(:,:,y) = tmp_xv(:,IB) ;
    end
    
    % Save this year's fertilizer
    if any(strcmp(fileList, 'fert'))
        [isincl, IB] = get_incl_suffix(detailed_tmp.varNames, '_FQ', cropList) ;
        tmp_xv = detailed_tmp.garr_xv(:,isincl) ;
        fert.garr_xvy(:,:,y) = tmp_xv(:,IB) ;
    end
    
    % Save this year's yield
    if any(strcmp(fileList, 'yield'))
        [isincl, IB] = get_incl_suffix(detailed_tmp.varNames, '_Y', cropList) ;
        tmp_xv = detailed_tmp.garr_xv(:,isincl) ;
        yield.garr_xvy(:,:,y) = tmp_xv(:,IB) ;
    end
        
    clear detailed_tmp isincl IB
end
fprintf('\nDone importing.\n')


%% Save

% Save land uses
thisFile = 'landuse' ;
if any(strcmp(fileList, thisFile))
    fprintf('Saving %s... ', thisFile)
    out_struct = landuse ;
    outFile = sprintf('%s/%s.garr.mat', outDir, thisFile) ;
    save(outFile, 'out_struct', '-v7.3') ;
    clear out_struct
    disp('Done.')
end

% Save crop fractions
thisFile = 'cropfracs' ;
if any(strcmp(fileList, thisFile))
    fprintf('Saving %s... ', thisFile)
    out_struct = cropfracs ;
    outFile = sprintf('%s/%s.garr.mat', outDir, thisFile) ;
    save(outFile, 'out_struct', '-v7.3') ;
    clear out_struct
    disp('Done.')
end

% Save irrigation
thisFile = 'irrig' ;
if any(strcmp(fileList, thisFile))
    fprintf('Saving %s... ', thisFile)
    out_struct = irrig ;
    outFile = sprintf('%s/%s.garr.mat', outDir, thisFile) ;
    save(outFile, 'out_struct', '-v7.3') ;
    clear out_struct
    disp('Done.')
end

% Save fertilizer
thisFile = 'fert' ;
if any(strcmp(fileList, thisFile))
    fprintf('Saving %s... ', thisFile)
    out_struct = fert ;
    outFile = sprintf('%s/%s.garr.mat', outDir, thisFile) ;
    save(outFile, 'out_struct', '-v7.3') ;
    clear out_struct
    disp('Done.')
end

% Save yield
thisFile = 'yield' ;
if any(strcmp(fileList, thisFile))
    fprintf('Saving %s... ', thisFile)
    out_struct = yield ;
    outFile = sprintf('%s/%s.garr.mat', outDir, thisFile) ;
    save(outFile, 'out_struct', '-v7.3') ;
    clear out_struct
    disp('Done.')
end



%% FUNCTIONS

function isincl = get_incl_lu(varNames)

isincl = ~endsWith(varNames, ...
    {'_A', '_FI', '_FQ', '_II', '_IQ', '_OI', '_Y'}) ;
[~,IA] = intersect(varNames, ...
    {'area', 'suitable', 'protected', 'pa_fraction'}) ;
isincl(IA) = false ;

end


function [isincl, IB] = get_incl_suffix(varNames, suffix, cropList)

isincl = endsWith(varNames, suffix) ...
    & ~contains(varNames, 'pasture') ;

IB = [] ;
if ~isempty(cropList)
    [~, ~, IB] = intersect(cropList, ...
        strrep(varNames(isincl), suffix, ''), ...
        'stable') ;
    if length(IB) ~= length(cropList)
        error('cropList mismatch')
    end
end

end


function out_struct = setup_struct(lonlats, list2map, varList, yearList)

out_struct.lonlats = lonlats ;
out_struct.list2map = list2map ;
out_struct.varNames = varList ;
out_struct.yearList = yearList ;


end
    