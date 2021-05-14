function out = getN(x) 

out = regexprep(regexprep(x, 'CerealsC[34]', ''), '^[a-zA-Z_]+', '') ;

end