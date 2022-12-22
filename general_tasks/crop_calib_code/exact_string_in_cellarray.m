function i = exact_string_in_cellarray(in_cellarray,in_string,varargin)

verbose = true ;
all = false ;
if ~isempty(varargin)
    if ~isempty(varargin{1})
        verbose = varargin{1} ;
    end
    if length(varargin)>1
        all = varargin{2} ;
    end
end

% i = find(not(cellfun('isempty', strfind(in_cellarray,in_string)))) ;
% 
% if ~all
%     if length(i) > 1
%         for c = 1:length(i)
%             if strcmp(in_cellarray(i(c)),in_string)
%                 i = i(c) ;
%                 break
%             end
%         end
%     end
% 
%     % If no EXACT match found, then empty set
%     if length(i) > 1
%         i = [] ;
%     end
% end
% 
% if isempty(i) && verbose
%     warning(['No cells match ' in_string '.'])
% end

i = find(not(cellfun('isempty', strfind(in_cellarray,in_string)))) ;

if ~all
    if length(i) > 1
        for c = 1:length(i)
            if strcmp(in_cellarray(i(c)),in_string)
                i = i(c) ;
                break
            end
        end
    end

    % If no EXACT match found, then empty set
    if length(i) > 1
        i = [] ;
    end
end

if isempty(i) && verbose
    warning(['No cells match ' in_string '.'])
end

end