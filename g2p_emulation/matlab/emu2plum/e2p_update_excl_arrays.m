function excl_vecs_after = ...
    e2p_update_excl_arrays( ...
    varNames, cropListI_before, cropListI_after, ...
    i_thisM, i_theseABetc, excl_vecs_before)


% Set up for expanding vectors
thisDest_cropI = getbasenamei(varNames{i_thisM}) ;
theseSources_cropI = getbasenamei(varNames(i_theseABetc)) ;
[~, Ibefore, Iafter] = intersect(cropListI_before, cropListI_after) ;
Iafter_dest = find(strcmp(cropListI_after, thisDest_cropI)) ;
if length(Iafter_dest) ~= 1
    error('Expected to find 1 Idest; found %d', ...
        length(Iafter_dest))
end
% % % [~, ~, Ibefore_sources] = intersect(cropListI_before, theseSources_cropI) ;
[~, Ibefore_sources] = intersect(cropListI_before, theseSources_cropI) ;
Nsources = length(theseSources_cropI) ;
if length(Ibefore_sources) ~= Nsources
    error('Expected to find %d Ibefore_sources; found %d', ...
        Nsources, length(Ibefore_sources))
end

% Expand vectors
excl_vecs_after = cell(size(excl_vecs_before)) ;
for x = 1:length(excl_vecs_after)
    excl_vecs_after{x} = combine_excl_vecs(excl_vecs_before{x}, ...
        Ibefore, Iafter, Iafter_dest, Ibefore_sources, ...
        cropListI_after) ;
end


end


function out_xc = combine_excl_vecs(in_xc, Ibefore, Iafter, ...
    Iafter_dest, Ibefore_sources, cropListI_after)

Nsources = length(Ibefore_sources) ;

if isvector(in_xc)
    size_out_xc = size(in_xc) ;
    size_out_xc(size_out_xc > 1) = 1 + size_out_xc(size_out_xc > 1) ;
    out_xc = nan(size_out_xc) ;
    out_xc(Iafter) = in_xc(Ibefore) ;
    out_xc(Iafter_dest) = ...
        sum(in_xc(Ibefore_sources)) == Nsources ;
else
    Ncells = size(in_xc, 1) ;
    out_xc = nan(Ncells, length(cropListI_after)) ;
    out_xc(:,Iafter) = in_xc(:,Ibefore) ;
    out_xc(:,Iafter_dest) = ...
        sum(in_xc(:,Ibefore_sources), 2) == Nsources ;
end

if any(any(isnan(out_xc)))
    error('NaN in output from combine_excl_vecs()')
end

end