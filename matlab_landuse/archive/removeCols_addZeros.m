function removeCols_addZeros(in_file, del_crops, Nlevels, varargin)
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
outPrec_lonlat = 2 ;
save_every_pct = 1 ;

iscellorempty = @(x) iscell(x) | isempty(x) ;
all_are_intORempty = @(x) all(isint(x)) | isempty(x) ;
is_vector = @(x) ismatrix(x) && min(size(x))==1 ;

p = inputParser ;
addRequired(p,  'in_file', @ischar) ;
addRequired(p,  'del_crops', iscellorempty) ;
addRequired(p,  'Nlevels', is_vector) ;
addParameter(p, 'outPrec', 6, @isint) ;
addParameter(p, 'outWidth', 1, @isint) ;
addParameter(p, 'delimiter', ' ', @ischar) ;
addParameter(p, 'overwrite', true, @islogical)
addParameter(p, 'fancy', false, @islogical)
addParameter(p, 'progress_step_pct', 20, @isint)

% Parse inputs
parse(p,...
    in_file, del_crops, Nlevels, ...
    varargin{:});
pFields = fieldnames(p.Results) ;
Nfields = length(pFields) ;
for f = 1:Nfields
    thisField = pFields{f} ;
    if ~exist(thisField,'var')
        eval([thisField ' = p.Results.' thisField ' ;']) ;
    end
    clear thisField
end ; clear f
clear p

if isempty(del_crops) && isempty(Nlevels)
    error('Both del_crops and Nlevels are empty.')
end

% Set up N strings
Nwidth = ceil(log10(max(Nlevels))) + 1 ;
if ~all(isint(Nlevels))
    error('Rework function to allow non-integer N levels')
end
Nformat = ['%0' num2str(Nwidth) 'd'] ;
Nlevels_token = regexprep(num2str(Nlevels), '\s*', '-') ;

% Generate output file
out_file = sprintf('%s.himelo0.%s.txt', ...
    strrep(in_file, '.txt', ''), Nlevels_token) ;

% Read original
disp('Reading...')
% if strcmp(in_file(end-2:end),'.gz')
%     in_file = strrep(in_file,'.gz','') ;
% end
in_data = lpjgu_matlab_read2geoArray(in_file, ...
    'verboseIfNoMat', false) ;

disp('Processing...')

% Delete unnecessary columns
[~, is_del_crop] = intersect(in_data.varNames, del_crops, 'stable') ;
if isfield(in_data, 'garr_xvy')
    in_data.garr_xvy(:,is_del_crop,:) = [] ;
elseif isfield(in_data, 'garr_xv')
    in_data.garr_xv(:,is_del_crop) = [] ;
else
    error('in_data has neither garr_xvy nor garr_xv')
end
in_data.varNames(is_del_crop) = [] ;

% Generate list of columns to add
add_crops = setdiff(in_data.varNames, ...
    {'Lon'; 'Lat'; 'Year'}, ...
    'stable') ;
N_addcrops = length(add_crops) ;
N_nlevels = length(Nlevels) ;
add_cols = cell(1, N_addcrops*N_nlevels) ;
for c = 1:N_addcrops
    thisCrop = add_crops{c} ;
    for n = 1:N_nlevels
        thisN = sprintf(Nformat, Nlevels(n)) ;
        add_cols{(c-1)*N_nlevels + n} = [thisCrop thisN] ;
    end
end

% Set up outputs
out_varNames = [in_data.varNames add_cols] ;


% % Delete unnecessary years at beginning
% error('untested')
% if isfield(out_data, 'garr_xvy')
%     ndel = 0 ;
%     for y = 2:length(out_data.yearList)
%         if isequaln(out_data.garr_xvy(:,:,y), out_data.garr_xvy(:,:,y-1))
%             ndel = ndel + 1 ;
%         else
%             break
%         end
%     end
% end
% if ndel > 0
%     fprintf('Removing %d years at beginning\n', ndel)
%     out_data.garr_xvy(:,:,1:ndel) = [] ;
% end

% % Delete unnecessary years at end
% error('untested')
% if isfield(out_data, 'garr_xvy')
%     ndel = 0 ;
%     for y = 0:length(out_data.yearList)-1
%         if isequaln(out_data.garr_xvy(:,:,end-y), out_data.garr_xvy(:,:,end-(y+1)))
%             ndel = ndel + 1 ;
%         else
%             break
%         end
%     end
% end
% if ndel > 0
%     out_data.garr_xvy(:,:,1:ndel) = [] ;
% end



% % Add zeros
% if isfield(in_data, 'garr_xvy')
%     out_data.yearList = in_data.yearList ;
%     Nyears = length(out_data.yearList) ;
%     out_data.garr_xvy = cat(2, ...
%         in_data.garr_xvy, ...
%         zeros(Ncells, Nvars_out - Nvars_in, Nyears)) ;
% elseif isfield(in_data, 'garr_xv')
%     out_data.garr_xv = cat(2, ...
%         in_data.garr_xv, ...
%         zeros(Ncells, Nvars_out - Nvars_in)) ;
% else
%     error('in_data has neither garr_xvy nor garr_xv')
% end
% clear in_data

% Find columns that should be zeros
[~, justZeroCols] = intersect(out_varNames, add_cols, 'stable') ;

% Get header
if isfield(in_data, 'garr_xvy')
    in_header_cell = [{'Lon', 'Lat', 'Year'} out_varNames] ;
elseif isfield(in_data, 'garr_xv')
    in_header_cell = [{'Lon', 'Lat'} out_varNames] ;
else
    error('in_data has neither garr_xvy nor garr_xv')
end
in_header_str = [] ;

% Make out_header_str
if isempty(in_header_str)
    % Get maximum length of variable names
    varName_max_width = 0 ;
    for v = 1:length(in_header_cell)
        varName_max_width = max(varName_max_width,length(in_header_cell{v})) ;
    end
    tmp_header_cell = in_header_cell ;
    for v = 1:length(in_header_cell)
        thisVar = tmp_header_cell{v} ;
        if strcmp(thisVar,'Lon') || strcmp(thisVar,'Lat') || strcmp(thisVar,'Year')
            outWidth = 5 ;
        else
            outWidth = min(outWidth,varName_max_width) ;
        end
        while length(thisVar) < outWidth
            thisVar = [thisVar ' '] ;
        end
        tmp_header_cell{v} = thisVar ;
    end
    out_header_str = strjoin(tmp_header_cell,delimiter) ;
else
    out_header_str = in_header_str ;
end

% Get lat/lon/yr columns
[i_lat,i_lon,i_year] = lpjgu_matlab_getLatLonYrCols(in_header_cell) ;

% Get output file formatSpec
out_header_cell = in_header_cell ;
out_formatSpec = '' ;
for ii = 1:length(out_header_cell)
    %     thisCol = out_header_cell{i} ;
    if ii>1
        out_formatSpec = [out_formatSpec delimiter] ;
    end
    if ii==i_lat || ii==i_lon
        out_formatSpec = [out_formatSpec '%-' num2str(outWidth) '.' num2str(outPrec_lonlat) 'f'] ;
    elseif ii==i_year || any(justZeroCols==ii)
        out_formatSpec = [out_formatSpec '%-' num2str(outWidth) 'u'] ;
    else
        out_formatSpec = [out_formatSpec '%-' num2str(outWidth) '.' num2str(outPrec) 'f'] ;
    end
end
out_formatSpec = [out_formatSpec ' \n'] ;

% Save header
fileID_out = fopen(out_file, 'w') ;
fprintf(fileID_out,'%s \n', out_header_str) ;

% Write data cell-by-cell
Ncells = length(in_data.list2map) ;
chunkSize = floor((save_every_pct/100)*Ncells) ;
Nchunks = ceil(Ncells / chunkSize) ;

% Save
disp('Saving...')

for ii = 1:Nchunks

    % Open
    if rem(ii-1, chunkSize)
        fileID_out = fopen(out_file, 'A') ;
    end
    
    % Get indices
    c1 = (ii-1)*chunkSize + 1 ;
    cN = min(ii*chunkSize, Ncells) ;
    thisChunkSize = cN - c1 + 1 ;
    
    % Get data array
    
    if isfield(in_data, 'garr_xvy')
        
        Nyears = length(in_data.yearList) ;
        
        col_lons = transpose(repmat(in_data.lonlats(c1:cN,1), [1 Nyears])) ;
        col_lons = col_lons(:) ;
        col_lats = transpose(repmat(in_data.lonlats(c1:cN,2), [1 Nyears])) ;
        col_lats = col_lats(:) ;
        col_lonlats = [col_lons col_lats] ;
        
        col_years = repmat(shiftdim(in_data.yearList), [thisChunkSize 1]) ;
        
%         out_data_tmp = cat(2, ...
%             repmat(out_data.lonlats(ii,:), [Nyears 1]), ...
%             Nyears*ones(Nyears,1)) ;
%         out_data_tmp2 = out_data.garr_xvy(ii,:,:) ;
%         out_data_tmp2 = squeeze(out_data_tmp2) ;
%         out_data_tmp2 = out_data_tmp2' ;
%         out_data_tmp = cat(2, out_data_tmp, out_data_tmp2) ;
        
        out_data_tmp2 = in_data.garr_xvy(c1:cN,:,:) ;
        out_data_tmp2 = permute(out_data_tmp2, [3 1 2]) ;
        out_data_tmp2 = reshape(out_data_tmp2, [thisChunkSize*Nyears size(out_data_tmp2,3)]) ;
        
        out_data_tmp = cat(2, ...
            col_lonlats, ...
            col_years, ...
            out_data_tmp2, ...
            zeros(thisChunkSize*Nyears, N_addcrops*N_nlevels)) ;
    elseif isfield(in_data, 'garr_xv')
        error('make this work')
    else
        error('in_data has neither garr_xvy nor garr_xv')
    end
    
    % Some array values may be "-0", which MATLAB treats as equal to "0", but
    % which can end up being printed as -0 due to different IEEE
    % representations.
    % https://www.mathworks.com/matlabcentral/answers/91430-why-does-fprintf-put-minus-signs-in-front-of-zeros-in-my-file
    out_data_tmp(out_data_tmp==0) = 0 ;
    
    fprintf(fileID_out,out_formatSpec,out_data_tmp) ;
    clear out_data_tmp
    
    fclose(fileID_out) ;
    thisPct = round(ii/Nchunks*100) ;
    fprintf('%d%%...\n', thisPct) ;
    pause(0.1)
    
end

    


end

