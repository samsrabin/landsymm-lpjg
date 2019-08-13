function [does_overlap, hlegend] = test_for_overlap(hf)

position_to_polyshape = @(p) polyshape( ...
    [p(1) p(1) p(1)+p(3) p(1)+p(3)], ...
    [p(2)+p(4) p(2) p(2) p(2)+p(4)]) ;

hlegend = hf.Children(1) ;
if ~strcmp(hlegend.Type,'legend')
    error('hf.Children(1) does not appear to be a legend; instead %s', hlegend.Type)
end
hlegend.Units = 'normalized' ;
ps_legend = position_to_polyshape(hlegend.Position) ;

haxes_kids = hf.Children(2).Children ;

does_overlap = false ;
for c = 1:length(haxes_kids)
    thisKid = haxes_kids(c) ;
    if ~isprop(thisKid,'String') || ~any(~isnan(thisKid.Extent))
        continue
    end
    thisKid.Units = 'normalized' ;
    ps_text = position_to_polyshape(thisKid.Extent) ;
    if overlaps(ps_legend, ps_text)
        does_overlap = true ;
        break
    end
end


end