function cmds_out = get_cell_forPlot(whos_in, varName, bl_or_yr, CFTnames)

theseVars = {whos_in.name}' ;
theseVars = {theseVars{strfindw(theseVars,['ts_*' varName '_*' bl_or_yr])}}' ;

if isempty(CFTnames)
    error('CFTnames is empty!')
end

for c = 1:length(CFTnames)
    thisCrop = CFTnames{c} ;
    thisVar = theseVars{getCellsWithString(theseVars,thisCrop)} ;
    if isempty(thisVar)
        error('No matching variables found!') ;
    elseif iscellstr(thisVar)
        error('More than one matching variable found!') ;
    end
    cmds_out{c,1} = ['cell_' bl_or_yr '{c} = ' thisVar ' ;'] ;
end


%     function tf_out = find_members(A,B)
%         % Finds indices of elements in A that match elements in B
%         tf_out = false(size(A)) ;
%         for i = 1:length(B)
%             is_match = find(~cellfun(@isempty,strfind(A,B{i}))) ;
%             if isempty(is_match)
%                 error([B{i} ' not found in A!'])
%             elseif length(is_match)>1
%                 error([B{i} ' repeated in A!'])
%             end
%             tf_out(is_match) = true ;
%         end
%         
%     end

end