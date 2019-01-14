function [barren_new, natural_new, urban_new] = ...
    BarNatUrb(barren_old, mgmt_old, natural_old, urban_old, ...
              barren_lpj, vegd_lpj, outPrec_LC)

%%thisCell = 45545 ;
%%thisCell = 22054 ;
%thisCell = 47outPrec_LC+27outPrec_LC+2 ;
%disp(sprintf('min(barren_old+mgmt_old+natural_old+urban_old) = %0.9g',min(barren_old+mgmt_old+natural_old+urban_old)))
%disp(sprintf('max(barren_old+mgmt_old+natural_old+urban_old) = %0.9g',max(barren_old+mgmt_old+natural_old+urban_old)))
%disp(sprintf('barren_old(thisCell) = %0.9g',barren_old(thisCell)))
%disp(sprintf('mgmt_old(thisCell) = %0.9g',mgmt_old(thisCell)))
%disp(sprintf('natural_old(thisCell) = %0.9g',natural_old(thisCell)))
%disp(sprintf('urban_old(thisCell) = %0.9g',urban_old(thisCell)))
%disp(sprintf('barren_lpj(thisCell) = %0.9g',barren_lpj(thisCell)))
%disp(sprintf('vegd_lpj(thisCell) = %0.9g',vegd_lpj(thisCell)))
%disp(' ')

barren_new = barren_old ;
urban_new = urban_old ;
natural_new = natural_old ;
vegd_new = natural_old + mgmt_old;
%barren_new = double(barren_old) ;
%urban_new = double(urban_old) ;
%natural_new = double(natural_old) ;
%vegd_new = double(natural_old) + double(mgmt_old);

% Not enough vegetated, before URBAN
%disp('(1, before)')
%disp(['natural_new(thisCell) = ' sprintf('%0.9g',natural_new(thisCell))])
%disp(['urban_new(thisCell) = ' sprintf('%0.9g',urban_new(thisCell))])
%disp(['vegd_new(thisCell) = ' sprintf('%0.9g',vegd_new(thisCell))])
thisCond = vegd_new < vegd_lpj ;
natural_new(thisCond) = natural_new(thisCond) + urban_new(thisCond) ;
urban_new(thisCond) = 0 ;
vegd_new = natural_new + mgmt_old ;
%disp('(1, after)')
%disp(['natural_new(thisCell) = ' sprintf('%0.9g',natural_new(thisCell))])
%disp(['urban_new(thisCell) = ' sprintf('%0.9g',urban_new(thisCell))])
%disp(['vegd_new(thisCell) = ' sprintf('%0.9g',vegd_new(thisCell))])
%disp(' ')

% Move the rest of URBAN to BARREN
%disp('(2, before)')
%disp(['barren_new(thisCell) = ' sprintf('%0.9g',barren_new(thisCell))])
%disp(['urban_new(thisCell) = ' sprintf('%0.9g',urban_new(thisCell))])
barren_new = barren_new + urban_new ;
urban_new = zeros(size(barren_new)) ;
%disp('(2, after)')
%disp(['barren_new(thisCell) = ' sprintf('%0.9g',barren_new(thisCell))])
%disp(['urban_new(thisCell) = ' sprintf('%0.9g',urban_new(thisCell))])
%disp(' ')

% Still not enough vegetated
thisCond = vegd_new < vegd_lpj ;
stillMissing = vegd_lpj - vegd_new ;
%disp('(3, before)')
%disp(['stillMissing(thisCell) = ' sprintf('%0.9g',stillMissing(thisCell))])
if any(stillMissing(thisCond)) < 0
    error('How is stillMissing negative?')
end
if any(round(barren_new(thisCond),outPrec_LC) < round(stillMissing(thisCond),outPrec_LC))
    error('Not enough BARREN to make up the difference!')
end
natural_new(thisCond) = natural_new(thisCond) + stillMissing(thisCond) ;
barren_new(thisCond) = barren_new(thisCond) - stillMissing(thisCond) ;
vegd_new = natural_new + mgmt_old ;
if any(round(vegd_new,outPrec_LC)<round(vegd_lpj,outPrec_LC))
    error(['How do you still not have enough vegetated? (max diff = ' sprintf('%0.9g',max(vegd_lpj-vegd_new)) ')'])
end

% Now not enough barren
thisCond = barren_new < barren_lpj  ;
stillMissing = barren_lpj - barren_new ;
if any(stillMissing) < 0
    error('How is stillMissing negative?')
elseif any(stillMissing > 0 & round(vegd_new,6)<=round(vegd_lpj,6))
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
if any(round(vegd_new,6) < round(vegd_lpj,6))
    error('Now not enough vegd!')
end
if any(barren_lpj - barren_new > 1e-6)
    warning(['You have some cell(s) with less barren now than in baseline! max=' ...
        num2str(max(barren_lpj - barren_new))])
end
if any(abs(vegd_new+barren_new+urban_new-1)>1e-6)
    error(['Something doesn''t sum to 1! max=' num2str(max(abs(vegd_new+barren_new+urban_new))) ', min=' num2str(min(abs(vegd_new+barren_new+urban_new)))])
end
if any(round(natural_new,6)<0 | round(natural_new,6)>1)
    error('Some natural_new<0 | natural_new>1')
end
if any(round(barren_new,6)<0 | round(barren_new,6)>1)
    error('Some barren_new<0 | barren_new>1')
end
if any(round(urban_new,6)<0 | round(urban_new,6)>1)
    error('Some urban_new<0 | urban_new>1')
end


end
