%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert cmass_wood_potharv_sts.out to LandSyMM-friendly form %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infile = '/Users/Shared/landsymm_forestry/1970past/1850-2014_germany_nohistluc/output-2022-04-30-020741/cmass_wood_potharv_sts.out' ;
data = lpjgu_matlab_read2geoArray(infile) ;
whos data

ok_years = find(data.yearList >= 1970) ;
[~, ok_vars] = intersect(data.varNames, {'ntrl', 'forC'}, 'stable') ;

data.garr_xvy = data.garr_xvy(:,ok_vars,ok_years) ;
data.yearList = data.yearList(ok_years) ;
data.varNames = data.varNames(ok_vars) ;

outfile = strrep(infile, 'cmass_wood_potharv_sts', 'landsymm_pharvW') ;

in_forHead = [{'Lon', 'Lat', 'Year'} data.varNames] ;
lpjgu_matlab_saveTable(in_forHead, data, outfile, 'outPrec', 6, 'overwrite', true)


function lpjgu_matlab_saveTable(in_forHead,out_array,out_file,varargin)

all_are_intORempty = @(x) all(isint(x)) | isempty(x) ;
is_nonneg_int = @(x) isint(x) & x>=0 & numel(x)==1 ;

is_ok_struct = @(x) isstruct(x) ...
    & (all(isfield(x, {'garr_xvy', 'yearList'})) | isfield(x, 'garr_xv')) ...
    & all(isfield(x, {'lonlats', 'list2map', 'varNames'})) ;

isPositiveIntegerValuedNumeric = @(x) isint(x) && x>0 ;

% Set up & parse input arguments
p = inputParser ;
addRequired(p,'in_forHead') ;
addRequired(p,'out_array',...
    @(x) is_ok_struct(x) || istable(x) || (ismatrix(x) && ~isvector(x))) ;
addRequired(p,'out_file',@ischar) ;
addParameter(p,'outPrec',3,is_nonneg_int) ;
addParameter(p,'outPrec_lonlat',2,is_nonneg_int) ;
addParameter(p,'outWidth',5,isPositiveIntegerValuedNumeric) ;
addParameter(p,'varNameAlign','L',...
    @(x) strcmpi(x,'L') || strcmpi(x,'R')) ;
addParameter(p,'dataAlign','L',...
    @(x) strcmpi(x,'L') || strcmpi(x,'R')) ;
addParameter(p,'fancy',false,@islogical) ;
addParameter(p,'progress_step_pct',5,@isnumeric) ;
addParameter(p,'save_every_n',1000,isPositiveIntegerValuedNumeric) ;
addParameter(p,'save_every_pct',5,isPositiveIntegerValuedNumeric) ;
addParameter(p,'delimiter',' ',@ischar) ;
addParameter(p,'overwrite',false,@islogical) ;
addParameter(p,'verbose',true,@islogical) ;
addParameter(p,'justZeroCols',[],all_are_intORempty) ;
addParameter(p,'gzip',false,@islogical) ;
parse(p,in_forHead,out_array,out_file,varargin{:});

if p.Results.save_every_n>uint32(inf)
    error('p.Results.save_every_n>uint32(inf)')
end
% p.Results.save_every_n = uint32(p.Results.save_every_n) ;

% Do not overwrite if file already exists, unless told to
if ~p.Results.overwrite && exist(out_file,'file')
    warning([out_file ' already exists. Specify overwrite=true to overwrite. Skipping.'])
    return
end

% Get header
if isstruct(out_array)
    if isfield(out_array, 'garr_xvy')
        in_header_cell = [{'Lon', 'Lat', 'Year'} out_array.varNames] ;
    elseif isfield(out_array, 'garr_xv')
        in_header_cell = [{'Lon', 'Lat'} out_array.varNames] ;
    else
        error('in_data has neither garr_xvy nor garr_xv')
    end
    in_header_str = [] ;
elseif ischar(in_forHead)
    fileID_in = fopen(in_forHead) ;
    in_header_str = fgetl(fileID_in) ;
    fclose(fileID_in) ;
    in_header_cell = textscan(in_header_str,'%s') ;
    in_header_cell = in_header_cell{1} ;
elseif iscell(in_forHead)
    % Check number of elements
    if length(in_forHead) ~= size(out_array,2)
        error('length(in_forHead) ~= size(out_array,2)')
    end
    
    if length(in_forHead)==2 && (iscell(in_forHead{2}))
        in_header_str = in_forHead{1} ;
        in_header_cell = in_forHead{2} ;
    else
        in_header_cell = in_forHead ;
        in_header_str = [] ;
    end
else
    error('in_forHead not valid.')
end

% NOW DOING THIS IN THE ACTUAL FUNCTIONS
% % Some array values may be "-0", which MATLAB treats as equal to "0", but 
% % which can end up being printed as -0 due to different IEEE
% % representations.
% % https://www.mathworks.com/matlabcentral/answers/91430-why-does-fprintf-put-minus-signs-in-front-of-zeros-in-my-file
% if isstruct(out_array)
%     if isfield(out_array, 'garr_xvy')
%         out_array.garr_xvy(out_array.garr_xvy==0) = 0 ;
%     elseif isfield(out_array, 'garr_xv')
%         out_array.garr_xv(out_array.garr_xv==0) = 0 ;
%     else
%         error('in_data has neither garr_xvy nor garr_xv')
%     end
% else
%     out_array(out_array==0) = 0 ;
% end

% Write data
if p.Results.fancy
    write_fancy(in_header_str, in_header_cell, p, out_file, out_array) ;
else
    write_fast(in_header_str, in_header_cell, p, out_file, out_array, p.Results.verbose) ;
% else
%     write_quiet(in_header_str, in_header_cell, p, out_file, out_array) ;
end

% Zip data
if p.Results.gzip
    disp('Zipping...')
    gzip(out_file)
    delete(out_file)
end

if p.Results.verbose
    disp('Done.')
end

end

function write_quiet(in_header_str, in_header_cell, p, out_file, out_array)

if isstruct(out_array)
    error('Rework write_quiet() to work with structure form of out_array')
end
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
            outWidth = min(p.Results.outWidth,varName_max_width) ;
        end
        while length(thisVar) < outWidth
            thisVar = [thisVar ' '] ;
        end
        tmp_header_cell{v} = thisVar ;
    end
    out_header_str = strjoin(tmp_header_cell,p.Results.delimiter) ;
else
    out_header_str = in_header_str ;
end


% Get lat/lon/yr columns
[i_lat,i_lon,i_year] = lpjgu_matlab_getLatLonYrCols(in_header_cell) ;

% Get output file formatSpec
out_header_cell = in_header_cell ;
out_header_cell = pad(out_header_cell,p.Results.outWidth) ;
out_formatSpec = '' ;
for i = 1:length(out_header_cell)
    %     thisCol = out_header_cell{i} ;
    if i>1
        out_formatSpec = [out_formatSpec p.Results.delimiter] ;
    end
    if i==i_lat || i==i_lon
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) '.' num2str(p.Results.outPrec_lonlat) 'f'] ;
    elseif i==i_year || any(p.Results.justZeroCols==i)
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) 'u'] ;
    elseif p.Results.outPrec==0
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) 'd'] ;
    else
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) '.' num2str(p.Results.outPrec) 'f'] ;
    end
end
out_formatSpec = [out_formatSpec ' \n'] ;

% Save header
fileID_out = fopen(out_file, 'w') ;
fprintf(fileID_out,'%s \n', out_header_str) ;

% Write data line-by-line
if istable(out_array)
    out_array = table2array(out_array) ;
end
Nrows = size(out_array,1) ;

for i = 1:Nrows
    
    % Some array values may be "-0", which MATLAB treats as equal to "0", but 
    % which can end up being printed as -0 due to different IEEE
    % representations.
    % https://www.mathworks.com/matlabcentral/answers/91430-why-does-fprintf-put-minus-signs-in-front-of-zeros-in-my-file
    out_line = out_array(i,:) ;
    out_line(out_line==0) = 0 ;
    
    % Write line
    fprintf(fileID_out,out_formatSpec,out_line) ;
    
end

fclose(fileID_out) ;

end

function write_fast(in_header_str, in_header_cell, p, out_file, out_array, verbose)

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
            outWidth = min(p.Results.outWidth,varName_max_width) ;
        end
        while length(thisVar) < outWidth
            thisVar = [thisVar ' '] ;
        end
        tmp_header_cell{v} = thisVar ;
    end
    out_header_str = strjoin(tmp_header_cell,p.Results.delimiter) ;
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
        out_formatSpec = [out_formatSpec p.Results.delimiter] ;
    end
    if ii==i_lat || ii==i_lon
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) '.' num2str(p.Results.outPrec_lonlat) 'f'] ;
    elseif ii==i_year || any(p.Results.justZeroCols==ii)
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) 'u'] ;
    else
        out_formatSpec = [out_formatSpec '%-' num2str(p.Results.outWidth) '.' num2str(p.Results.outPrec) 'f'] ;
    end
end
out_formatSpec = [out_formatSpec ' \n'] ;

% Save header
fileID_out = fopen(out_file, 'w') ;
fprintf(fileID_out,'%s \n', out_header_str) ;
fclose(fileID_out) ;

% Write data
if isstruct(out_array)
    
    % Write data cell-by-cell
    Ncells = length(out_array.list2map) ;
    chunkSize = floor((p.Results.save_every_pct/100)*Ncells) ;
    Nchunks = ceil(Ncells / chunkSize) ;
    
    for ii = 1:Nchunks
        % Open
        fileID_out = fopen(out_file, 'A') ;
        
        % Get indices
        c1 = (ii-1)*chunkSize + 1 ;
        cN = min(ii*chunkSize, Ncells) ;
        thisChunkSize = cN - c1 + 1 ;
        
        % Get data array
        
        if isfield(out_array, 'garr_xvy')
            
            Nyears = length(out_array.yearList) ;
            
            col_lons = transpose(repmat(out_array.lonlats(c1:cN,1), [1 Nyears])) ;
            col_lons = col_lons(:) ;
            col_lats = transpose(repmat(out_array.lonlats(c1:cN,2), [1 Nyears])) ;
            col_lats = col_lats(:) ;
            col_lonlats = [col_lons col_lats] ;
            
            col_years = repmat(shiftdim(out_array.yearList), [thisChunkSize 1]) ;
            
%             out_array_tmp = cat(2, ...
%                 repmat(out_array.lonlats(ii,:), [Nyears 1]), ...
%                 Nyears*ones(Nyears,1)) ;
%             out_array_tmp2 = out_array.garr_xvy(ii,:,:) ;
%             out_array_tmp2 = squeeze(out_array_tmp2) ;
%             out_array_tmp2 = out_array_tmp2' ;
%             out_array_tmp = cat(2, out_array_tmp, out_array_tmp2) ;

            out_array_tmp2 = out_array.garr_xvy(c1:cN,:,:) ;
            out_array_tmp2 = permute(out_array_tmp2, [3 1 2]) ;
            try
    out_array_tmp2 = reshape(out_array_tmp2, [thisChunkSize*Nyears size(out_array_tmp2,3)]) ;
catch ME
    keyboard
end

            
            out_array_tmp = cat(2, ...
                col_lonlats, ...
                col_years, ...
                out_array_tmp2) ;
        elseif isfield(out_array, 'garr_xv')
            out_array_tmp = cat(2, ...
                out_array.lonlats(c1:cN,:), ...
                out_array.garr_xv(c1:cN,:)) ;
        else
            error('in_data has neither garr_xvy nor garr_xv')
        end
        
        % Some array values may be "-0", which MATLAB treats as equal to "0", but
        % which can end up being printed as -0 due to different IEEE
        % representations.
        % https://www.mathworks.com/matlabcentral/answers/91430-why-does-fprintf-put-minus-signs-in-front-of-zeros-in-my-file
        out_array_tmp(out_array_tmp==0) = 0 ;
        
        fprintf(fileID_out,out_formatSpec,out_array_tmp') ;
        
        fclose(fileID_out) ;
        thisPct = round(ii/Nchunks*100) ;
        if verbose
            fprintf('%d%%...\n', thisPct) ;
        end
        pause(0.1)
                
    end
    
else
    
    % Write data chunk-by-chunk
    if istable(out_array)
        out_array = table2array(out_array) ;
    end
    out_array = transpose(out_array) ;
    Nrows = size(out_array,2) ;
    Nchunks = ceil(100 / p.Results.save_every_pct) ;
    chunkSize = ceil(Nrows*100/p.Results.save_every_pct) ;
    
    for ii = 1:Nchunks
        
        % Open
        fileID_out = fopen(out_file, 'A') ;
        
        % Get indices
        c1 = (ii-1)*chunkSize + 1 ;
        cN = min(ii*chunkSize, Nrows) ;
        
        % Some array values may be "-0", which MATLAB treats as equal to "0", but
        % which can end up being printed as -0 due to different IEEE
        % representations.
        % https://www.mathworks.com/matlabcentral/answers/91430-why-does-fprintf-put-minus-signs-in-front-of-zeros-in-my-file
        out_array_tmp = out_array(:,c1:cN) ;
        out_array_tmp(out_array_tmp==0) = 0 ;
        
        fprintf(fileID_out,out_formatSpec,out_array_tmp) ;
        fclose(fileID_out) ;
        
        if ii~=Nchunks
            thisPct = round(ii/Nchunks*100) ;
            if rem(ii,10)==0
                formatSpec = '%d%%...\n' ;
            else
                formatSpec = '%d%%...' ;
            end
            fprintf(formatSpec, thisPct) ;
            pause(0.1)
        end
        
    end
end

end

function write_fancy(in_header_str, in_header_cell, p, out_file, out_array)

if isstruct(out_array)
    error('Rework write_fancy() to work with structure form of out_array')
end

% Make out_header_str
Nvars = length(in_header_cell) ;
if isempty(in_header_str)
    
    outWidths = nan(size(in_header_cell)) ;
    tmp_header_cell = in_header_cell ;
    for v = 1:Nvars
        thisVar = tmp_header_cell{v} ;
        if strcmp(thisVar,'Lon') || strcmp(thisVar,'Lat') || strcmp(thisVar,'Year')
            outWidths(v) = 5 ;
        else
            outWidths(v) = max(p.Results.outWidth,length(thisVar)) ;
        end
        while length(thisVar) < outWidths(v)
            if strcmpi(p.Results.varNameAlign,'L')
                thisVar = [thisVar ' '] ;
            elseif strcmpi(p.Results.varNameAlign,'R')
                thisVar = [' ' thisVar] ;
            else
                error('How did this happen? Fix parser.')
            end
        end
        tmp_header_cell{v} = thisVar ;
    end
    out_header_str = strjoin(tmp_header_cell,p.Results.delimiter) ;
else
    out_header_str = in_header_str ;
end


% Get lat/lon/yr columns
[i_lat,i_lon,i_year] = lpjgu_matlab_getLatLonYrCols(in_header_cell) ;

% Get output file formatSpec
out_formatSpec = '' ;
for v = 1:Nvars
    out_formatSpec = [out_formatSpec '%s'] ;
    if v < Nvars
        out_formatSpec = [out_formatSpec '\t'] ;
    end
end
out_formatSpec = [out_formatSpec '\n'] ;

% Save header
fileID_out = fopen(out_file, 'w') ;
fprintf(fileID_out,'%s \n', out_header_str) ;
fclose(fileID_out) ;

% Set up contents of eval()
fprintf_command = 'fprintf(fileID_out,out_formatSpec' ;
for v = 1:Nvars
    fprintf_command = [fprintf_command ',out_cells{' num2str(v) '}'] ;
end
fprintf_command = [fprintf_command ') ;'] ;

% Write data line-by-line, saving every p.Results.save_every_n lines
if istable(out_array)
    out_array = table2array(out_array) ;
end
Nrows = size(out_array,1) ;
progress = 0 ;
isopen = false ;
if p.Results.verbose
    disp('Saving...')
end
tic ;
j = 0 ;
k = 0 ;
k_targ = ceil(Nrows * p.Results.progress_step_pct/100) ;
for i = 1:Nrows
    
    % Update j and k
    j = j + 1 ;
    k = k + 1 ;
    
    % Re-open file if just closed
    if ~isopen
        fileID_out = fopen(out_file,'a') ;
        isopen = true ;
    end
    
    % Construct string & write line
    out_cells = cell(Nvars,1) ;
    for v = 1:Nvars
        thisVar = in_header_cell{v} ;
        thisElmt = num2str(out_array(i,v)) ;
        if strcmp(thisElmt,'0')
            thisElmt = [thisElmt '.' repmat('0',[1 p.Results.outPrec])] ;
        end
        % Pad right of decimal with zeros, if needed
        if ~strcmp(thisVar,'Lon') && ~strcmp(thisVar,'Lat') && ~strcmp(thisVar,'Year')
            C = strsplit(thisElmt,'.') ;
            L = C{1} ;
            if length(C)==2
                R = C{2} ;
            else
                R = '0' ;
            end
            clear C
            while length(R) < p.Results.outPrec
                R = [R '0'] ;
            end
            thisElmt = [L '.' R] ;
        end
        % Pad with spaces, if needed
        if length(thisElmt) < outWidths(v)
            padding = repmat(' ',[1 outWidths(v)-length(thisElmt)]) ;
            if strcmpi(p.Results.varNameAlign,'L')
                thisElmt = [thisElmt padding] ;
            elseif strcmpi(p.Results.varNameAlign,'R')
                thisElmt = [padding thisElmt] ;
            else
                error('How did this happen? Fix parser.')
            end
        end
        
        out_cells{v} = thisElmt ;
    end
    eval(fprintf_command) ;
    
    % Close file every p.Results.save_every_n lines
    %             if rem(i,p.Results.save_every_n)==0
    if j==p.Results.save_every_n
        j = 0 ;
        fclose(fileID_out) ;
        isopen = false ;
    end
    
    % Update progress every whatever%
    %             if p.Results.verbose && rem(i,ceil(Nrows*p.Results.progress_step_pct/100))==0
    if p.Results.verbose && k==k_targ
        k = 0 ;
        progress = progress + p.Results.progress_step_pct ;
        disp(['   ' num2str(progress) '% complete (' toc_hms(toc) ')'])
    end
end
if isopen
    fclose(fileID_out) ;
end
if p.Results.verbose
    disp(toc_hms(toc))
end
end



















% % % % % Write data line-by-line, saving every p.Results.save_every_n lines
% % % % out_array = table2array(out_array) ;
% % % % Nchunks = ceil(size(out_array,1) / p.Results.save_every_n) ;
% % % % disp('Saving...')
% % % % tic ;
% % % % for i = 1:Nchunks
% % % %
% % % %     % Get line numbers
% % % %     i1 = (i-1)*p.Results.save_every_n + 1 ;
% % % %     iN = i*p.Results.save_every_n ;
% % % %     if iN > size(out_array,1)
% % % %         if i == Nchunks
% % % %             iN = size(out_array,1) ;
% % % %         else
% % % %             error('iN > size(out_array,1) but not on last chunk!')
% % % %         end
% % % %     end
% % % %
% % % %     % Write lines
% % % %     tmp = out_array(i1:iN,:) ;
% % % %     precision = ['%.' num2str(p.Results.outPrec) 'f'] ;
% % % %     dlmwrite(out_file,tmp,'-append',...
% % % %              'delimiter',p.Results.delimiter,...
% % % %              'precision',precision) ;
% % % %
% % % %     % Update progress
% % % %     progress = round(100*i/Nchunks,1) ;
% % % %     disp(['   Chunk ' num2str(i) ' of ' num2str(Nchunks) ' complete (' num2str(progress) '%), ' toc_hms(toc) ')'])
% % % % end
% % % % disp('Done.')