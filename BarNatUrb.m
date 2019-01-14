function [barren_new, natural_new, urban_new] = ...
    BarNatUrb(barren_old, mgmt_old, natural_old, urban_old, ...
              barren_lpj, vegd_lpj, outPrec_LC)

barren_new = barren_old ;
urban_new = urban_old ;
natural_new = natural_old ;
vegd_new = natural_old + mgmt_old;

% Not enough vegetated, before URBAN
thisCond = vegd_new < vegd_lpj ;
natural_new(thisCond) = natural_new(thisCond) + urban_new(thisCond) ;
urban_new(thisCond) = 0 ;
vegd_new = natural_new + mgmt_old ;

% Move the rest of URBAN to BARREN
barren_new = barren_new + urban_new ;
urban_new = zeros(size(barren_new)) ;

% Still not enough vegetated? Move some BARREN to NATURAL.
thisCond = vegd_new < vegd_lpj ;
stillMissing = vegd_lpj - vegd_new ;
stillMissing = stillMissing(thisCond) ;
if any(stillMissing) < 0
    error('How is stillMissing negative?')
end
if any(round(barren_new(thisCond),outPrec_LC) < round(stillMissing,outPrec_LC))
    error('Not enough BARREN to make up the difference!')
end
natural_new(thisCond) = natural_new(thisCond) + stillMissing ;
barren_new(thisCond) = barren_new(thisCond) - stillMissing ;
vegd_new = natural_new + mgmt_old ;
if any(round(vegd_new,outPrec_LC)<round(vegd_lpj,outPrec_LC))
    error('How do you still not have enough vegetated?')
end

% Not enough barren? Move some NATURAL to BARREN.
thisCond = barren_new < barren_lpj  ;
stillMissing = barren_lpj - barren_new ;
if any(stillMissing) < 0
    error('How is stillMissing negative?')
elseif any(stillMissing > 0 & round(vegd_new,outPrec_LC)<=round(vegd_lpj,outPrec_LC))
%     error('any(stillMissing > 0 & vegd_new<=vegd_lpj)')
end
toMoveNtrl2bare = stillMissing ;
excess_vegd = vegd_new - vegd_lpj ;
excess_vegd = max(0,excess_vegd) ;
toMoveNtrl2bare = min(excess_vegd,toMoveNtrl2bare) ;
toMoveNtrl2bare = min(natural_new,toMoveNtrl2bare) ;
if round(min(toMoveNtrl2bare))<0
    error('round(min(toMoveNtrl2bare))<0')
end
natural_new(thisCond) = natural_new(thisCond) - toMoveNtrl2bare(thisCond) ;
barren_new(thisCond) = barren_new(thisCond) + toMoveNtrl2bare(thisCond) ;
vegd_new = natural_new + mgmt_old ;
if any(round(vegd_new,outPrec_LC) < round(vegd_lpj,outPrec_LC))
    error('Now not enough vegd!')
end
if any(barren_lpj - barren_new > 10^-outPrec_LC)
    warning(['You have some cell(s) with less barren now than in baseline! max=' ...
        num2str(max(barren_lpj - barren_new))])
end
if any(abs(vegd_new+barren_new+urban_new-1)>10^-outPrec_LC)
    error(['Something doesn''t sum to 1! max=' num2str(max(abs(vegd_new+barren_new+urban_new))) ', min=' num2str(min(abs(vegd_new+barren_new+urban_new)))])
end
if any(round(natural_new,outPrec_LC)<0 | round(natural_new,outPrec_LC)>1)
    error('Some natural_new<0 | natural_new>1')
end
if any(round(barren_new,outPrec_LC)<0 | round(barren_new,outPrec_LC)>1)
    error('Some barren_new<0 | barren_new>1')
end
if any(round(urban_new,outPrec_LC)<0 | round(urban_new,outPrec_LC)>1)
    error('Some urban_new<0 | urban_new>1')
end


end