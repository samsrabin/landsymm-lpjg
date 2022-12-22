function [array_in, varNames_in] = trim_miscanthus(array_in, varNames_in)

if any(contains(varNames_in,'Miscanthus'))
    if any(any(array_in(:,:,contains(varNames_in,'Miscanthus'),:)>0))
        error('This script assumes zero Miscanthus area!')
    end
    array_in(:,:,contains(varNames_in,'Miscanthus'),:) = [] ;
    varNames_in(contains(varNames_in,'Miscanthus')) = [] ;
end


end