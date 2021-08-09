function y = get_unneeded(x)

y = cellfun(@isempty, ...
    regexp(regexprep(x,'CerealsC[34]','CerealsC'),'.*\d+')) ...
    | contains(x,'G_ic') ;

end