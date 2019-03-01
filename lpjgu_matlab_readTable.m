function out_table = lpjgu_matlab_readTable(in_file,varargin)

p = inputParser ;
addRequired(p,'in_file',@ischar) ;
addOptional(p,'dont_save_MAT',false,@islogical) ;
addOptional(p,'do_save_MAT',false,@islogical) ;
addOptional(p,'verbose',false,@islogical) ;
addOptional(p,'verboseIfNoMat',true,@islogical) ;
addOptional(p,'dispPrefix','',@ischar) ;
parse(p,in_file,varargin{:});
dont_save_MAT = p.Results.dont_save_MAT ;
verbose = p.Results.verbose ;
dispPrefix = p.Results.dispPrefix ;
if p.Results.dont_save_MAT && p.Results.do_save_MAT
    error('dont_save_MAT and do_save_MAT can''t both be true.')
end
verboseIfNoMat = p.Results.verboseIfNoMat ;
if verboseIfNoMat && verbose
%     warning('Because verboseIfNoMat==true, setting verbose=false.')
    verbose = false ;
end

% If in_file is symlink, replace it with its target
[s,w] = unix(['[[ -L ' in_file ' ]] && echo true']) ;
if s==0 && contains(w,'true') % is symlink
    disp('Symlink; pointing to target instead.')
    [~,w] = unix(['stat -f "%Y" ' in_file]) ;
    in_file = regexprep(w,'[\n\r]+','') ; % Remove extraneous newline
end

if strcmp(in_file(end-2:end),'.gz')
    in_file = in_file(1:end-3) ;
end

in_matfile = [in_file '.mat'] ;
if exist(in_matfile,'file')
    if verbose
        disp([dispPrefix '   Loading table from MAT-file...'])
    end
    load(in_matfile) ;
else
    did_unzip = gunzip_if_needed(in_file) ;
    if verbose || verboseIfNoMat
        disp([dispPrefix '   Making table...'])
    end
    import_this_to_table() ;
    if ~dont_save_MAT
        if p.Results.do_save_MAT
            save(in_matfile,'out_table','-v7.3') ;
        else
            prompt_to_save() ;
        end
    end
    if did_unzip
        if verbose || verboseIfNoMat
            disp([dispPrefix '      Deleting unzipped in_file...'])
        end
        err2 = system(['rm "' in_file '"']) ;
        if err2~=0
            error('Error in rm.')
        end
    end
end

% Standardize column names
standardize_colnames()


    function did_unzip = gunzip_if_needed(in_file)
        did_unzip = false ;
        if ~exist(in_file,'file')
            in_file_gz = [in_file '.gz'] ;
            if exist(in_file_gz,'file')
                if verbose || verboseIfNoMat
                    disp([dispPrefix '   Unzipping in_file...'])
                end
                err1 = system(['gunzip < "' in_file_gz '" > "' in_file '"']) ;
                if err1~=0
                    error('Error in gunzip.')
                end
                did_unzip = true ;
            else
                error('in_file.mat, in_file, and in_file.gz not found.')
            end
        end
    end


    function import_this_to_table()
        % Read header to get field names
        in_header = read_header() ;
        
        % Read data
        %%% Based on http://de.mathworks.com/matlabcentral/answers/89374-importdata-fails-to-import-large-files
        if verbose || verboseIfNoMat
            disp([dispPrefix '      Reading data...'])
        end
        fid = fopen(in_file) ;
        firstline = ftell(fid) ;  % Remember where we are
        numcols = length(regexp(fgetl(fid), '\s+', 'split')) ;    %read line, split it into columns, count them
        if length(in_header)==2
            if verbose || verboseIfNoMat
                warning('Not skipping first line.')
            end
            fseek(fid, firstline, 'bof') ;                    %reposition to prior line
        end
        fmt = repmat('%f', 1, numcols) ;             %maybe %d if entries are integral
        datacell = textscan( fid, fmt, 'CollectOutput', 1) ;  %read file
        fclose(fid) ;
        A = datacell{1};
        
        % Get rid of empty columns at the end
        while ~any(~isnan(A(:,numcols)))
            A(:,numcols) = [] ;
            numcols = numcols - 1 ;
        end
        
        is_gridlist = false ;
        if size(A,2) < 2
            error('Input must be 2-d matrix.')
        elseif size(A,2) == 2
            if verbose || verboseIfNoMat
                warning('Assuming columns Lon, Lat.')
            end
            colNames = {'Lon','Lat'} ;
            is_gridlist = true ;
        else
            colNames = in_header ;
        end
        
        if length(in_header) ~= size(A,2)
            error('Length of column name list does not match horizontal dimension of array.')
        end
        
        % Make table
        varNames = strrep(colNames,'"','') ;
        varNames = strrep(varNames,'-','_') ;
        out_table = array2table(A,'VariableNames',varNames) ;
        
        % Deal with doubled data in 2005 for N_fert_rcp85_6f.out
        [~,in_name,in_ext] = fileparts(in_file) ;
        if strcmp([in_name in_ext],'N_fert_rcp85_6f.out')
            bad = find(out_table.Year==2005) ;
            out_table(bad(2:2:end),:) = [] ;
        end
        
    end


    function prompt_to_save()
        ok = false ;
        while ~ok
            if exist(in_matfile,'file')
                disp([dispPrefix '      Save, overwriting existing MAT-file? Y or [N]. 10 seconds...'])
                default_save = false ;
            else
                disp([dispPrefix '      Save to MAT-file? [Y] or N. 10 seconds...'])
                default_save = true ;
            end
            dbl = getkeywait_ssr(10) ;
            if (dbl==-1 && default_save) || strcmp(char(dbl),'y') || strcmp(char(dbl),'Y')
                ok = true ;
                if verbose || verboseIfNoMat
                    disp([dispPrefix '      Saving MAT-file...'])
                end
                save(in_matfile,'out_table','-v7.3') ;
            elseif (dbl==-1 && ~default_save) || strcmp(char(dbl),'n') || strcmp(char(dbl),'N')
                ok = true ;
            elseif verbose || verboseIfNoMat
                warning(['Input (' char(dbl) ') not recognized.'])
            end
        end
    end


    function standardize_colnames()
        varNames = out_table.Properties.VariableNames ;
        for v = 1:length(varNames)
            if strcmp(varNames{v},'lon')
                out_table.Properties.VariableNames{'lon'} = 'Lon' ;
            elseif strcmp(varNames{v},'lat')
                out_table.Properties.VariableNames{'lat'} = 'Lat' ;
            elseif strcmp(varNames{v},'year')
                out_table.Properties.VariableNames{'year'} = 'Year' ;
            end
        end
    end


    function out_header = read_header()
        fid = fopen(in_file) ;
        H = fgetl(fid) ;
        fclose(fid) ;
        in_header = regexp(H, '\s+', 'split') ;
        out_header = in_header ;
        if isempty(out_header{1})
            out_header = out_header(2:end) ;
        end
        if isempty(out_header{end})
            out_header = out_header(1:end-1) ;
        end
        for f = 1:length(out_header)
            if isempty(out_header{f})
                error('Empty variable name!')
            end
        end
    end


end