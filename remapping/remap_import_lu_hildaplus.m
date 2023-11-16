function [out_lu, carea_YX] = remap_import_lu_hildaplus( ...
    geodata_dir, yearList_out, gridlist)

Nyears_out = length(yearList_out) ;
Ncells = length(gridlist.list2map) ;

% Land use input file and years it contains
yearList_lu_states = 1901:2020 ;
file_lu_states = fullfile(geodata_dir, 'HILDA+', 'hildaplus_global_05deg_lpjguess', ...
    'data', 'hildaplus_netfrac_1901_2020.txt' ...
    ) ;
if min(yearList_out) < min(yearList_lu_states) || max(yearList_out) > max(yearList_lu_states)
    error('yearList_out (%d-%d) not entirely contained within yearList_lu_states (%d-%d)', ...
        min(yearList_out), max(yearList_out), min(yearList_lu_states), max(yearList_lu_states))
end

% Read
disp('Importing land uses...')
fprintf('  file_lu_states: %s\n', file_lu_states) ;
in_lu = lpjgu_matlab_read2geoArray(file_lu_states, ...
    'target', gridlist, ...
    'verboseIfNoMat', false, ...
    'force_mat_save', true, 'force_mat_nosave', false) ;
Nlu_in = length(in_lu.varNames) ;
[C, ~, IB] = intersect(yearList_out, in_lu.yearList, 'stable') ;
if ~isequal(shiftdim(C), shiftdim(yearList_out))
    error('yearList mismatch')
end
in_lu.garr_xvy = in_lu.garr_xvy(:,:,IB) ;
in_lu.yearList = in_lu.yearList(IB) ;

% Get land use names from input; match to names in output
list_LU_in = in_lu.varNames ;
list_LU_out = {'NATURAL','CROPLAND','PASTURE','BARREN'} ;
Nlu_out = length(list_LU_out) ;
map_LU_in2out = cell(size(list_LU_in)) ;
map_LU_in2out(contains(list_LU_in,...
    {'CROPLAND'})) = {'CROPLAND'} ;
map_LU_in2out(contains(list_LU_in,...
    {'PASTURE'})) = {'PASTURE'} ;
map_LU_in2out(contains(list_LU_in,...
    {'NATURAL', 'FOREST'})) = {'NATURAL'} ;
map_LU_in2out(contains(list_LU_in,...
    {'URBAN', 'BARREN'})) = {'BARREN'} ;
if any(isempty(map_LU_in2out))
    error('Remaining empty element(s) in map_LU_in2out!')
elseif ~isempty(setdiff(map_LU_in2out,list_LU_out))
    error('Some member of map_LU_in2out is not present in list_LU_out!')
end
warning('Should work out specification of mapping with check for duplicates on LHS.')

% Combine input LU types into output types
if any(any(any(isnan(in_lu.garr_xvy))))
    error('Handle NaNs in in_lu.garr_xvy')
end
out_lu = rmfield(in_lu, {'garr_xvy', 'varNames'}) ;
out_lu.varNames = list_LU_out ;
out_lu.garr_xvy = zeros(Ncells, Nlu_out, Nyears_out) ;
for v = 1:Nlu_in
    thisLU_in = list_LU_in{v} ;
    thisLU_out = map_LU_in2out{v} ;
    fprintf('    %s (%d of %d) to %s...\n',thisLU_in,v,Nlu_in,thisLU_out)
    i = strcmp(list_LU_out,thisLU_out) ;
    out_lu.garr_xvy(:,i,:) = out_lu.garr_xvy(:,i,:) + in_lu.garr_xvy(:,v,:) ;
end
if max(max(abs(sum(out_lu.garr_xvy, 2) - sum(in_lu.garr_xvy, 2)))) > 1e-9
    error('Combining input LU types into output types: Changed sum of LU fractions')
end

% Import cell area (km2); aggregate to half-degree
file_lu_etc = fullfile(landsymm_lpjg_path(), 'data', 'geodata', 'staticData_quarterdeg.nc') ;
fprintf('file_lu_etc: %s\n', file_lu_etc) ;
carea_XY = ncread(file_lu_etc,'carea') ;
carea_YX = flipud(transpose(carea_XY)) ;
carea_YX = aggregate(carea_YX,0.25,0.5) ;
carea = rmfield(gridlist, 'mask_YX') ;
carea.garr_x = carea_YX(gridlist.mask_YX) ;
carea_hd_XY = aggregate(carea_XY,0.25,0.5) ;

disp('Finishing...')

% Ensure that LU fractions sum to 1, which is required for next step as coded
max_lufrac_sum_diff_from_1 = remap_get_xvy_sumLU_max_diff(out_lu.garr_xvy, 1) ;
if max_lufrac_sum_diff_from_1 >= 1e-6
    error('LU fractions don''t appear to sum to 1. How should ice/water fraction be handled?')
end

% Add water fraction to BARREN
icwtr_YX = flipud(transpose(ncread(file_lu_etc,'icwtr'))) ;
icwtr_YX(icwtr_YX==1) = 0 ;
icwtr_hd_YX = aggregate(icwtr_YX.*flipud(carea_XY'),0.25,0.5)./flipud(carea_hd_XY') ;
icwtr_hd_x = icwtr_hd_YX(gridlist.list2map) ;
v = strcmp(list_LU_out,'BARREN') ;
Nlu = length(out_lu.varNames) ;
icwtr_hd_x1y = repmat(icwtr_hd_x, [1, 1, Nyears_out]) ;
icwtr_hd_xvy = repmat(icwtr_hd_x1y, [1, Nlu, 1]) ;
tmp_xvy = out_lu.garr_xvy ;
out_lu.garr_xvy = out_lu.garr_xvy ./ icwtr_hd_xvy ;
out_lu.garr_xvy(icwtr_hd_xvy==0) = tmp_xvy(icwtr_hd_xvy==0) ;
clear tmp_xvy
out_lu.garr_xvy(:,v,:) = out_lu.garr_xvy(:,v,:) ...
    + icwtr_hd_x1y ;

% Check again
max_lufrac_sum_diff_from_1 = remap_get_xvy_sumLU_max_diff(out_lu.garr_xvy, 1) ;
if max_lufrac_sum_diff_from_1 >= 1e-6
    error('Error in adding ice/water fraction to BARREN: Fractions do not sum to 1')
end

disp('Done.')

end
