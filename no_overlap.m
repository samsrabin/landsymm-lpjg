function bool_out = no_overlap(set1, set2)

if set1(2) < set1(1) && set2(2) < set2(1)
    set1(2) = set1(2) + 12 ;
    set2(2) = set2(2) + 12 ;
end
if set1(2)<=set2(1) && set1(1)>=set2(2)
    bool_out = true ;
else
    bool_out = false ;
end


end
