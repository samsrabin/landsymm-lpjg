prev = load('~/testing20180904.mat') ;
prev_fields = fieldnames(prev) ;

for f = 1:length(prev_fields)
    thisField = prev_fields{f} ;
    if exist(thisField,'var')
        eval(['thisVar = ' thisField ' ;']) ;
        eval(['thisVar_prev = prev.' thisField ' ;']) ;
        if isstruct(thisVar) && isfield(thisVar,'maps_YXv') && ~isequaln(thisVar.maps_YXv,thisVar_prev.maps_YXv)
            fprintf('%s.maps_YXv differs!\n',thisField)
        elseif (isnumeric(thisVar) || iscell(thisVar) || ischar(thisVar) || islogical(thisVar)) ...
            && ~isequaln(thisVar,thisVar_prev)
            fprintf('%s differs!\n',thisField)
        elseif ~(isstruct(thisVar) || isnumeric(thisVar) || iscell(thisVar) || ischar(thisVar) || strcmp(thisField,'ans') || islogical(thisVar))
            stop
        end
    end
end
disp('Done')