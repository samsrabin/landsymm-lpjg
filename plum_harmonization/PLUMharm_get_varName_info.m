function [LUnames, maxLength_LUnames, Nlu, notBare, isAgri, LUnames_agri, ...
    isAgri_isPast, isCrop] = ...
    PLUMharm_get_varName_info(varNames)

LUnames = varNames ;
maxLength_LUnames = max(cellfun(@length,LUnames)) ;
Nlu = length(LUnames) ;
notBare = ~strcmp(LUnames,'BARREN') ;
isAgri = ~strcmp(LUnames,'NATURAL') & notBare ;
LUnames_agri = LUnames(isAgri) ;
isAgri_isPast = strcmp(LUnames_agri,'PASTURE') ;
isCrop = ~strcmp(LUnames,'NATURAL') & ~strcmp(LUnames,'PASTURE') & notBare ;

end