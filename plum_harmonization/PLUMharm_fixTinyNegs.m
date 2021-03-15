function out_YXv = PLUMharm_fixTinyNegs(in_YXv, landArea_YXv, debugIJ_2deg)
% Rounding errors can result in small negative values. Fix.

Nlu = size(in_YXv, 3) ;
do_debug = ~isempty(debugIJ_2deg) ;

tolerance = 1 ;   % m2
if min(in_YXv(:)) < -tolerance
    error('"Large" negative values of area! Min %0.1e', min(in_YXv(:)))
end

out_YXv = in_YXv ;

isBad_YX = any(out_YXv<0,3) ;

if do_debug
    dbI = debugIJ_2deg(1) ;
    dbJ = debugIJ_2deg(2) ;
    isdb_YX = false(size(in_YXv, 1:2)) ;
    isdb_YX(dbI,dbJ) = true ;
    
    fprintf('PLUMharm_fixTinyNegs() for (%d,%d):\n', dbI, dbJ)
end

if any(any(isBad_YX))
    Nbad = length(find(isBad_YX)) ;
    isBad_YXv = repmat(isBad_YX,[1 1 Nlu]) ;
    
    % Get garrays
    tmp = out_YXv(isBad_YXv) ;
    if do_debug
        isdb = isdb_YX(isBad_YX) ;
        I_db = find(isdb) ;
    end
    tmp_landArea = landArea_YXv(isBad_YXv) ;
    tmp = tmp ./ tmp_landArea ;
    if any(isnan(tmp))
        error('How do you have NaN in tmp?')
    end
    tmp_xv = reshape(tmp,[Nbad Nlu]) ;
    
    % Report initial areas in debug cell
    if do_debug
        tmp_db = tmp_xv(I_db,:) ;
        fprintf(['   Before: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], tmp_db) ;
    end
    
    % Set negatives to zero
    tmp_xv(tmp_xv<0) = 0 ;
    tmp_xSum = sum(tmp_xv,2) ;
    
    % Report updated areas in debug cell
    if do_debug
        tmp_db = tmp_xv(I_db,:) ;
        fprintf(['    After: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], tmp_db) ;
    end
    
    j = 0 ;
    while max(abs(tmp_xSum-1)) > 3*eps
        j = j+1 ;
        if j>50
            error('Possible infinite loop in fixing tiny negative areas!')
        end
        tmp_xv = tmp_xv ./ repmat(tmp_xSum,[1 Nlu]) ;
        tmp_xSum = sum(tmp_xv,2) ;
        
        if do_debug
            tmp_db = tmp_xv(I_db,:) ;
            fprintf(['         %d: ' repmat('%0.3e ',[1 length(tmp_db)]) '\n'], j, tmp_db) ;
        end
    end
    out_YXv(isBad_YXv) = tmp_xv(:) .* tmp_landArea ;
end


end