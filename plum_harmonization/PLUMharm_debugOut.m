function PLUMharm_debugOut(debug_title,debug_kind,agri_YXv,landArea_YX,debugIJ,LUnames)

outPrec = 2 ;

Nlu = length(LUnames) ;
i = debugIJ(1) ;
j = debugIJ(2) ;
outPrec = num2str(outPrec) ;

if strcmp(debug_kind,'areas')
    areas = cell(Nlu+2,1) ;
    fracs = cell(Nlu+2,1) ;
else
    areas = cell(Nlu,1) ;
    fracs = cell(Nlu,1) ;
end

for v = 1:Nlu
    area_tmp = agri_YXv(i,j,v) ;
    if area_tmp == 0
        areas{v} = '0' ;
    else
        areas{v} = sprintf(['%0.' outPrec 'e'],area_tmp) ;
    end
    thisArea = landArea_YX(i,j,min(v,size(landArea_YX,3))) ;
    frac_tmp = agri_YXv(i,j,v)/thisArea ;
    frac_tmp(thisArea==0) = 0 ;
    if frac_tmp == 0
        fracs{v} = '0' ;
    else
        fracs{v} = sprintf(['%0.' outPrec 'e'],frac_tmp) ;
    end
end

if strcmp(debug_kind,'areas')
    area_tmp = sum(agri_YXv(i,j,~contains(LUnames,{'PASTURE','NATURAL','BARREN'})),3) ;
    if area_tmp == 0
        areas{end-1} = '0' ;
    else
        areas{end-1} = sprintf(['%0.' outPrec 'e'],area_tmp) ;
    end
    frac_tmp = sum(agri_YXv(i,j,~contains(LUnames,{'PASTURE','NATURAL','BARREN'})),3)/squeeze(landArea_YX(i,j,:)) ;
    if frac_tmp == 0
        fracs{end-1} = '0' ;
    else
        fracs{end-1} = sprintf(['%0.' outPrec 'e'],frac_tmp) ;
    end
    areas{end} = sprintf(['%0.' outPrec 'e'],sum(agri_YXv(i,j,~strcmp(LUnames,'BARREN')),3)) ;
    fracs{end} = sprintf(['%0.' outPrec 'f'],sum(agri_YXv(i,j,~strcmp(LUnames,'BARREN')),3)./squeeze(landArea_YX(i,j,:))) ;
    T = table([shiftdim(LUnames);{'CROPLAND';'vegd'}],areas,fracs) ;
else
    T = table(shiftdim(LUnames),areas,fracs) ;
end

T.Properties.VariableNames{1} = debug_title ;
T.Properties.VariableNames{2} = debug_kind ;
disp(T)



end