function thisCountry_short_filename = thisCountry_short_toFilename( ...
    thisCountry_short)

thisCountry_short_filename = thisCountry_short ;
thisCountry_short_filename = strrep(thisCountry_short_filename,' and ','') ;
thisCountry_short_filename = strrep(thisCountry_short_filename,' ','') ;
thisCountry_short_filename = strrep(thisCountry_short_filename,'+','') ;
thisCountry_short_filename = strrep(thisCountry_short_filename,'&','') ;
thisCountry_short_filename = strrep(thisCountry_short_filename,'(','') ;
thisCountry_short_filename = strrep(thisCountry_short_filename,')','') ;
thisCountry_short_filename = strrep(thisCountry_short_filename,'-','') ;


end