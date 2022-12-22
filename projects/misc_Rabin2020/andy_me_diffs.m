% diffs = {
%     '380 to 415' ;
%     '393 to 459' ;
%     '1575 to 1578' ;
%     '1585 to 1587' ;
%     '5.7 to 1.5' ;
%     '4.6 to 0.3' ;
%     '1961 to 1995' ;
%     '1983 to 2047' ;
%     '0.250 to 0.240' ;
%     '0.249 to 0.240' ;
%     '0.182 to 0.179' ;
%     '0.182 to 0.179' ;
%     '58.6 to 57.9' ;
%     '58.9 to 58.8' ;
%     '52.5 to 55.1' ;
%     '52.2 to 54.3' ;
%     '17.9 to 18.9' ;
%     '17.9 to 18.8' ;
%     '28.9 to 35.9' ;
%     '27.5 to 45.2' ;
%     '60.3 to 109.7' ;
%     '73.3 to 119.0' ;
%     '477 to 419' ;
%     '503 to 495' ;
%     '40.7 to 38.9' ;
%     '41.9 to 40.5' ;
%     } ;

diffs = {
    '369 to 415' ;
    '1408 to 1396' ;
    '2.7 to 4.9' ;
    '1780 to 1817' ;
    '0.283 to 0.260' ;
    '0.181 to 0.173' ;
    '61.3 to 58.2' ;
    '36.4 to 41.0' ;
    '4.3 to 4.6' ;
    '9.8 to 8.0' ;
    '72.6 to 68.0' ;
    '496 to 452' ;
    '42.7 to 38.6' ;
    } ;

for i = 1:length(diffs)
    thisTxt = diffs{i} ;
    splitTxt = strsplit(thisTxt) ;
    y0 = str2double(splitTxt{1}) ;
    y1 = str2double(splitTxt{3}) ;
    thisDiff_pct = (y1-y0)/y0 * 100 ;
    thisSign = ' ' ;
    if thisDiff_pct<0
        thisSign = '-' ;
    elseif thisDiff_pct>0
        thisSign = '+' ;
    end
    fprintf('%s (%s%.1f%%)\n',thisTxt,thisSign,abs(thisDiff_pct)) ;
end





