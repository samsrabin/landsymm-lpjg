function PLUMharm_debugOut_deltas(debug_title,debug_kind,y0_YXv,y1_YXv,debugIJ,LUnames)

outPrec = 2 ;

Nlu = length(LUnames) ;
i = debugIJ(1) ;
j = debugIJ(2) ;
outPrec = num2str(outPrec) ;

areas_y0_dbl = nan(Nlu+2,1) ;
% fracs_y0 = nan(Nlu+2,1) ;
areas_y1_dbl = nan(Nlu+2,1) ;
% fracs_y1 = nan(Nlu+2,1) ;

% y0
areas_y0_dbl(1:Nlu) = squeeze(y0_YXv(i,j,:)) ;
areas_y0_dbl(end-1) = sum(y0_YXv(i,j,~contains(LUnames,{'PASTURE','NATURAL','BARREN'})),3) ;
areas_y0_dbl(end)   = sum(y0_YXv(i,j,~strcmp(LUnames,'BARREN')),3) ;
% fracs_y0 = areas_y0_dbl ./ repmat(landArea_YX(i,j),size(areas_y0_dbl)) ;

% y1
areas_y1_dbl(1:Nlu) = squeeze(y1_YXv(i,j,:)) ;
areas_y1_dbl(end-1) = sum(y1_YXv(i,j,~contains(LUnames,{'PASTURE','NATURAL','BARREN'})),3) ;
areas_y1_dbl(end)   = sum(y1_YXv(i,j,~strcmp(LUnames,'BARREN')),3) ;
% fracs_y1 = areas_y1_dbl ./ repmat(landArea_YX(i,j),size(areas_y1_dbl)) ;

% Diffs
diffs_dbl = areas_y1_dbl - areas_y0_dbl ;
diffs_pct_dbl = diffs_dbl ./ areas_y0_dbl * 100 ;

% Turn into strings
areas_y0 = cell(size(areas_y0_dbl)) ;
for v = 1:length(areas_y0_dbl)
    if areas_y0_dbl(v)~=0
        areas_y0{v} = sprintf(['%0.' outPrec 'e'],areas_y0_dbl(v)) ;
    else
        areas_y0{v} = '0' ;
    end
end

areas_y1 = cell(size(areas_y1_dbl)) ;
for v = 1:length(areas_y1_dbl)
    if areas_y1_dbl(v)~=0
        areas_y1{v} = sprintf(['%0.' outPrec 'e'],areas_y1_dbl(v)) ;
    else
        areas_y1{v} = '0' ;
    end
end

diffs = cell(size(diffs_dbl)) ;
for v = 1:length(diffs_dbl)
    if diffs_dbl(v)~=0
        diffs{v} = sprintf(['%.' outPrec 'e'],diffs_dbl(v)) ;
    else
        diffs{v} = '' ;
    end
end
diffs(diffs_dbl>0) = strcat('+',diffs(diffs_dbl>0)) ;

diffs_pct = cell(size(diffs_pct_dbl)) ;
for v = 1:length(diffs_pct_dbl)
    if diffs_pct_dbl(v)~=0 && ~isnan(diffs_pct_dbl(v))
        diffs_pct{v} = sprintf(['%0.' outPrec 'f'],diffs_pct_dbl(v)) ;
    else
        diffs_pct{v} = '' ;
    end
end
diffs_pct(diffs_dbl>0) = strcat('+',diffs_pct(diffs_dbl>0)) ;
diffs_pct(diffs_dbl~=0) = strcat(diffs_pct(diffs_dbl~=0),' %') ;

T = table([shiftdim(LUnames);{'CROPLAND';'vegd'}],areas_y0,areas_y1,diffs,diffs_pct) ;
T.Properties.VariableNames{1} = debug_title ;
T.Properties.VariableNames{2} = [debug_kind '_y0'] ;
T.Properties.VariableNames{3} = [debug_kind '_y1'] ;
disp(T)



end