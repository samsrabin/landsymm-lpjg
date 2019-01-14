function removeCols_addZeros(in_file, del_cols, add_cols, out_file, varargin)
outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;

iscellorempty = @(x) iscell(x) | isempty(x) ;
all_are_intORempty = @(x) all(isint(x)) | isempty(x) ;

p = inputParser ;
addRequired(p,  'in_file', @ischar) ;
addRequired(p,  'del_cols', iscellorempty) ;
addRequired(p,  'add_cols', iscellorempty) ;
addRequired(p,  'out_file', @ischar) ;
addParameter(p, 'outPrec', 6, @isint) ;
addParameter(p, 'outWidth', 1, @isint) ;
addParameter(p, 'delimiter', ' ', @ischar) ;
addParameter(p, 'overwrite', true, @islogical)
addParameter(p, 'fancy', false, @islogical)
addParameter(p, 'progress_step_pct', 20, @isint)
addParameter(p, 'justZeroCols',[],all_are_intORempty) ;

% Parse inputs
parse(p,...
    in_file, del_cols, add_cols, out_file, ...
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

if isempty(del_cols) && isempty(add_cols)
    warning('Both del_cols and add_cols are empty. Exiting.')
else
    
    % Read original
    disp('Reading...')
    if strcmp(in_file(end-2:end),'.gz')
        in_file = strrep(in_file,'.gz','') ;
    end
    in_table = lpjgu_matlab_readTable(in_file) ;
    in_header_cell = in_table.Properties.VariableNames ;
    
    % Set up outputs
    disp('Processing...')
    out_array = table2array(in_table) ;
    clear in_table
    out_header_cell = in_header_cell ;
    
    % Delete columns
    if ~isempty(del_cols)
        [~,IA] = intersect(in_header_cell, del_cols) ;
        out_array(:,IA) = [] ;
        out_header_cell(IA) = [] ;
    end
    
    % Add columns (zeros)
    if ~isempty(add_cols)
        out_array = [out_array zeros(size(out_array,1),length(add_cols))] ;
        out_header_cell = [out_header_cell add_cols] ;
    end
    
    % Find columns that should be zeros
    
    disp('Saving...')
    lpjgu_matlab_saveTable(out_header_cell, out_array, out_file,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', progress_step_pct, ...
        'justZeroCols', justZeroCols) ;
    
end


end

