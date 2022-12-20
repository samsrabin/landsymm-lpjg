function cell_out = collate_empties(cell_in)

% C = cell(size(cell_in)) ;
% C(:) = repmat({''},size(C)) ;
% cell_out = cell(length(cell_in),2) ;
% cell_out(:,1) = cell_in ;
% cell_out(:,2) = C ;
% cell_out = cell_out' ;
% cell_out = cell_out(:) ;

C = cell(size(cell_in)) ;
C(:) = repmat({''},size(C)) ;

cell_out = cell(length(cell_in),2) ;
if iscell(cell_in)
    cell_out(:,1) = cell_in ;
else
    for i = 1:length(cell_in)
        cell_out{i,1} = cell_in(i) ;
    end
end
cell_out(:,2) = C ;
cell_out = cell_out' ;
cell_out = cell_out(:) ;



end