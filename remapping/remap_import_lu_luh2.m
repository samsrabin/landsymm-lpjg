function [out_lu, carea_YX] = remap_import_lu_luh2( ...
    geodata_dir, yearList_lu_states, yearList_out, gridlist)

Nyears_out = length(yearList_out) ;
Ncells = length(gridlist.list2map) ;

% LUH2 input files and years they contain
file_luh2_states = fullfile(geodata_dir, 'LUH2', 'v2h', 'states.1850-2015.nc') ;
if min(yearList_out) < min(yearList_lu_states) || max(yearList_out) > max(yearList_lu_states)
    error('yearList_out must be entirely contained within yearList_luh2_states!')
end
file_luh2_mgmts = fullfile(geodata_dir, 'LUH2', 'v2h', 'management.1850-2015.nc') ;
yearList_luh2_mgmts = 1850:2015 ;
if min(yearList_out) < min(yearList_luh2_mgmts) || max(yearList_out) > max(yearList_luh2_mgmts)
    error('yearList_out must be entirely contained within yearList_luh2_mgmts!')
end

% Get info for reading input files
[~,yearIs_luh2_states,~] = intersect(yearList_lu_states,yearList_out) ;
starts_luh2_states = [1 1 min(yearIs_luh2_states)] ;
counts_luh2_states = [Inf Inf Nyears_out] ;
[~,yearIs_luh2_mgmts,~] = intersect(yearList_luh2_mgmts,yearList_out) ;
starts_luh2_mgmts = [1 1 min(yearIs_luh2_mgmts)] ;
counts_luh2_mgmts = [Inf Inf Nyears_out] ;

% Get land use names from input; match to names in output
finfo = ncinfo(file_luh2_states) ;
list_LU_in = setdiff({finfo.Variables.Name},...
                     {'time','lon','lat','lat_bounds','lon_bounds','secma','secmb'}) ;
list_LU_out = {'NATURAL','CROPLAND','PASTURE','BARREN'} ;
Nlu_in = length(list_LU_in) ;
Nlu_out = length(list_LU_out) ;
map_LU_in2out = cell(size(list_LU_in)) ;
map_LU_in2out(contains(list_LU_in,...
    {'c3ann','c3nfx','c3per','c4ann','c4per'})) = {'CROPLAND'} ;
map_LU_in2out(contains(list_LU_in,...
    {'pastr','range'})) = {'PASTURE'} ;
map_LU_in2out(contains(list_LU_in,...
    {'primf','primn','secdf','secdn','secma','secmb'})) = {'NATURAL'} ;
map_LU_in2out(contains(list_LU_in,...
    {'urban'})) = {'BARREN'} ;
if any(isempty(map_LU_in2out))
    error('Remaining empty element(s) in map_LU_in2out!')
elseif ~isempty(setdiff(map_LU_in2out,list_LU_out))
    error('Some member of map_LU_in2out is not present in list_LU_out!')
end
warning('Should work out specification of mapping with check for duplicates on LHS.')

disp('Importing land uses...')
fprintf('  file_luh2_states: %s\n', file_luh2_states) ;

% Import cell area (km2); aggregate to half-degree
file_luh2_etc = fullfile(landsymm_lpjg_path(), 'data', 'geodata', 'staticData_quarterdeg.nc') ;
fprintf('file_luh2_etc: %s\n', file_luh2_etc) ;
carea_XY = ncread(file_luh2_etc,'carea') ;
carea_YX = flipud(transpose(carea_XY)) ;
carea_YX = aggregate(carea_YX,0.25,0.5) ;
carea = rmfield(gridlist, 'mask_YX') ;
carea.garr_x = carea_YX(gridlist.mask_YX) ;
carea_XYy = repmat(carea_XY,[1 1 Nyears_out]) ;
carea_hd_XY = aggregate(carea_XY,0.25,0.5) ;
carea_hd_XYy = repmat(carea_hd_XY,[1 1 Nyears_out]) ;

% Land use
out_lu = rmfield(gridlist, 'mask_YX') ;
out_lu.varNames = list_LU_out ;
out_lu.yearList = yearList_out ;
for v = 1:Nlu_in
    thisLU_in = list_LU_in{v} ;
    thisLU_out = map_LU_in2out{v} ;
    fprintf('    %s (%d of %d) to %s...\n',thisLU_in,v,Nlu_in,thisLU_out)
    i = strcmp(list_LU_out,thisLU_out) ;
    lu_in_XYy = ncread(file_luh2_states,thisLU_in,starts_luh2_states,counts_luh2_states) ;
    lu_in_XYy(isnan(lu_in_XYy)) = 0 ;
    lu_out_YXy = flipud(permute( ...
        aggregate(lu_in_XYy.*carea_XYy,0.25,0.5)./carea_hd_XYy, ...
        [2 1 3])) ;
    if any(any(any(lu_out_YXy-1 > 1e-6)))
        error('Some element(s) of lu_out_YXy > 1!')
    end
    clear lu_in_XYy
    % Convert to _xy
    lu_out_xy = lpjgu_YXz_to_xz(lu_out_YXy, [Ncells, size(lu_out_YXy, 3)], gridlist.list2map) ;
    clear lu_out_YXy
    
    if v==1
        out_lu.garr_xvy = zeros([Ncells length(out_lu.varNames) size(lu_out_xy, 2)]) ;
    end
    out_lu.garr_xvy(:,i,:) = out_lu.garr_xvy(:,i,:) + permute(lu_out_xy, [1 3 2]) ;
    if any(any(any(any(out_lu.garr_xvy-1 > 1e-6))))
        error('Some element(s) of out_lu.garr_xvy > 1!')
    end
    clear lu_out_xy
end
clear carea*_YXy carea_hd_XYy

disp('Finishing...')

% Add water fraction to BARREN
icwtr_YX = flipud(transpose(ncread(file_luh2_etc,'icwtr'))) ;
icwtr_YX(icwtr_YX==1) = 0 ;
icwtr_hd_YX = aggregate(icwtr_YX.*flipud(carea_XY'),0.25,0.5)./flipud(carea_hd_XY') ;
icwtr_hd_x = icwtr_hd_YX(gridlist.list2map) ;
v = strcmp(list_LU_out,'BARREN') ;
out_lu.garr_xvy(:,v,:) = out_lu.garr_xvy(:,v,:) ...
    + repmat(icwtr_hd_x,[1 1 Nyears_out]) ;

% Mask cells with no vegetated land according to LU dataset
bad_x = out_lu.garr_xvy(:,4,1)==1 ;
if any(bad_x)
    error('Decide how you want to handle masking of cells with no vegetated land.')
    fprintf('Removing %d cells with no vegetated land...\n', length(find(bad_x))) ;
    out_lu.garr_xvy(bad_x,:,:) = [] ;
    out_lu.list2map(bad_x) = [] ;
    out_lu.lonlats(bad_x,:) = [] ;
    gridlist.mask_YX(gridlist.list2map(bad_x)) = false ;
    gridlist.list2map(bad_x) = [] ;
    gridlist.lonlats(bad_x,:) = [] ;
    Ncells = length(gridlist.list2map) ;
end

% Force all land cells to sum to 1
lu_out_x1y = sum(out_lu.garr_xvy,2) ;
j = 0 ;
while any(any(abs(lu_out_x1y-1)>1e-6))
    j = j + 1;
    if j > 50
        error('Possible infinite loop in "Force all land cells to sum to 1".')
    end
    out_lu.garr_xvy = out_lu.garr_xvy ./ repmat(lu_out_x1y,[1 Nlu_out 1]) ;
    lu_out_x1y = sum(out_lu.garr_xvy,2) ;
end
clear lu_out_x1y

disp('Done.')

end