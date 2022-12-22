function cell_out = collate_twocells(c1, c2)

C = cell(size(c1,1),2) ;
C(:,1) = c1 ;
C(:,2) = c2 ;
C = C' ;
cell_out = C(:) ;


end