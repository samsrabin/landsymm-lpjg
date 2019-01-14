function [...
    lu_out, lu_area_YXvy, ...
    cf_out, cf_areaT_YXvy, cf_areaI_YXvy, ...
    nf_out, nfert_tot_YXvy] = ...
    ...
    compare_hist_LUs_import( ...
    lu_file, cf_file, nf_file, ...
    LUlist, yearList, landArea_YX, cropList)

Nlu = length(LUlist) ;
Nyears = length(yearList) ;
Ncrops = length(cropList) ;

% Import land uses
lu_out = lpjgu_matlab_readTable_then2map(lu_file) ;
lu_out = trim_years(lu_out, yearList) ;

% Trim land uses
if ~isequal(lu_out.varNames,LUlist)
    [~,IA] = setdiff(lu_out.varNames,LUlist) ;
    if any(lu_out.maps_YXvy(:,:,IA,:)>0)
        error('LUs excluded from lu_out have area >0!')
    end
    clear IA
    if ~isempty(setdiff(LUlist,lu_out.varNames))
        error('lu_out missing something from LUlist!')
    end
    [~,~,IB] = intersect(LUlist,lu_out.varNames,'stable') ;
    lu_out.maps_YXvy = lu_out.maps_YXvy(:,:,IB,:) ;
    clear IB
    lu_out.varNames = LUlist ;
end

% Get land use areas
lu_area_YXvy = lu_out.maps_YXvy .* repmat(landArea_YX,[1 1 Nlu Nyears]) ;

% Import crop fractions
cf_out = import_croppy(cf_file, yearList) ;

% Combine crop fractions, if needed
if sum(contains(cf_out.varNames,{'CC3G','CC4G'}))==2
    tmpFrac_3 = cf_out.maps_YXvy(:,:,strcmp(cf_out.varNames,'CC3G'),:) ;
    tmpFrac_4 = cf_out.maps_YXvy(:,:,strcmp(cf_out.varNames,'CC4G'),:) ;
    tmpWeight_3 = tmpFrac_3 ./ (tmpFrac_3 + tmpFrac_4) ;
    tmpWeight_4 = tmpFrac_4 ./ (tmpFrac_3 + tmpFrac_4) ;
    tmpWeight_3(tmpFrac_3+tmpFrac_4 == 0) = 0 ;
    tmpWeight_4(tmpFrac_3+tmpFrac_4 == 0) = 0 ;
    cf_out.maps_YXvy(:,:,strcmp(cf_out.varNames,'CC3G'),:) = sum(cf_out.maps_YXvy(:,:,contains(cf_out.varNames,{'CC3G','CC4G'}),:),3) ;
    cf_out.varNames{strcmp(cf_out.varNames,'CC3G')} = 'ExtraCrop' ;
    cf_out.maps_YXvy(:,:,strcmp(cf_out.varNames,'CC4G'),:) = [] ;
    cf_out.varNames(strcmp(cf_out.varNames,'CC4G')) = [] ;
elseif any(contains(cf_out.varNames,{'CC3G','CC4G'}))
	error('???')
end

% Get crop areas
cf_areaT_YXvy = nan([size(landArea_YX) Ncrops Nyears]) ;
cf_areaI_YXvy = nan([size(landArea_YX) Ncrops Nyears]) ;
for c = 1:Ncrops
    thisCrop = cropList{c} ;
    thisCropI = [thisCrop 'i'] ;
    if any(strcmp(cf_out.varNames,thisCrop))
        cf_areaT_YXvy(:,:,c,:) = sum(cf_out.maps_YXvy(:,:,contains(cf_out.varNames,thisCrop),:),3) ;
    else
        warning('%s missing.',thisCrop)
        cf_areaT_YXvy(:,:,c,:) = 0 ;
    end
    if any(strcmp(cf_out.varNames,thisCropI))
        cf_areaI_YXvy(:,:,c,:) = cf_out.maps_YXvy(:,:,strcmp(cf_out.varNames,thisCropI),:) ;
    elseif any(strcmp(cf_out.varNames,thisCrop))
        warning('%s missing.',thisCropI)
        cf_areaI_YXvy(:,:,c,:) = 0*cf_out.maps_YXvy(:,:,strcmp(cf_out.varNames,thisCrop),:) ;
    else
        warning('%s missing.',thisCropI)
        cf_areaI_YXvy(:,:,c,:) = 0 ;
    end
end
cf_out.maps_YXvy = cf_areaT_YXvy ;
cf_out.varNames = cropList ;
croparea_tmp_YXvy = repmat(lu_area_YXvy(:,:,strcmp(LUlist,'CROPLAND'),:),[1 1 Ncrops 1]) ;
cf_areaT_YXvy = cf_out.maps_YXvy .* croparea_tmp_YXvy ;
cf_areaI_YXvy = cf_areaI_YXvy .* croparea_tmp_YXvy ;

% Import nfert
nf_out = import_croppy(nf_file, yearList) ;
if sum(contains(nf_out.varNames,{'CC3G','CC4G'}))==2
    if any(any(any(nf_out.maps_YXvy(:,:,contains(nf_out.varNames,{'CC3G','CC4G'}),:)>0)))
        warning('Fertilizer applied to CC3G and/or CC4G.')
        tmpNfert_3 = nf_out.maps_YXvy(:,:,contains(nf_out.varNames,{'CC3G'}),:) ;
        tmpNfert_4 = nf_out.maps_YXvy(:,:,contains(nf_out.varNames,{'CC4G'}),:) ;
        nf_out.maps_YXvy(:,:,contains(nf_out.varNames,{'CC3G'}),:) = ...
            tmpWeight_3.*tmpNfert_3 + tmpWeight_4.*tmpNfert_4 ;
        clear tmp*
    end
    nf_out.varNames{strcmp(nf_out.varNames,'CC3G')} = 'ExtraCrop' ;
    nf_out.maps_YXvy(:,:,strcmp(nf_out.varNames,'CC4G'),:) = [] ;
    nf_out.varNames(strcmp(nf_out.varNames,'CC4G')) = [] ;
elseif any(contains(nf_out.varNames,{'CC3G','CC4G'}))
	error('???')
end
for c = 1:Ncrops
    thisCrop = cropList{c} ;
    thisCropI = [thisCrop 'i'] ;
    if strcmp(thisCrop,'ExtraCrop')
        if any(any(any(nf_out.maps_YXvy(:,:,strcmp(nf_out.varNames,thisCrop),:)>0)))
            error('This code assumes 0 fertilization of CC3G/CC4G/ExtraCrop.')
        end
    elseif any(strcmp(nf_out.varNames,thisCropI))
        if ~isequaln( ...
                nf_out.maps_YXvy(:,:,strcmp(nf_out.varNames,thisCrop),:), ...
                nf_out.maps_YXvy(:,:,strcmp(nf_out.varNames,thisCropI),:))
            error('This code assumes equal fertilization of rainfed and irrigated versions of a crop.')
        end
    end
end
if isequal(setdiff(cropList,nf_out.varNames),{'ExtraCrop'})
    nf_out.maps_YXvy = cat(3,nf_out.maps_YXvy,zeros([size(landArea_YX) 1 Nyears])) ;
    nf_out.varNames = [nf_out.varNames {'ExtraCrop'}] ;
end
[~,~,IB] = intersect(cropList,nf_out.varNames,'stable') ;
if length(IB) ~= Ncrops
    error('length(IB) ~= Ncrops')
end
nf_out.maps_YXvy = nf_out.maps_YXvy(:,:,IB,:) ;
clear IB
nf_out.varNames = cropList ;

% Get Nfert totals (kg)
nfert_tot_YXvy = nf_out.maps_YXvy .* (1e6*cf_areaT_YXvy) ;


end


function struct_out = import_croppy(file_in, yearList)

struct_out = lpjgu_matlab_readTable_then2map(file_in) ;
if ~isfield(struct_out,'maps_YXvy')
    [struct_out.maps_YXv, struct_out.varNames] = trim_miscanthus(struct_out.maps_YXv, struct_out.varNames) ;
    tmp_YXv = struct_out.maps_YXv ;
    struct_out = rmfield(struct_out,'maps_YXv') ;
    struct_out.yearList = yearList ;
    struct_out.maps_YXvy = repmat(tmp_YXv,[1 1 1 length(struct_out.yearList)]) ;
    clear tmp
else
    struct_out = trim_years(struct_out, yearList) ;
    [struct_out.maps_YXvy, struct_out.varNames] = trim_miscanthus(struct_out.maps_YXvy, struct_out.varNames) ;
end

end


function struct_in = trim_years(struct_in, yearList)

Nyears = length(yearList) ;
[~,IA] = intersect(struct_in.yearList,yearList) ;
if length(IA) ~= Nyears
    error('length(IA) ~= Nyears')
end
struct_in.maps_YXvy = struct_in.maps_YXvy(:,:,:,IA) ;
struct_in.yearList = yearList ;

end


function [array_in, varNames_in] = trim_miscanthus(array_in, varNames_in)

if any(contains(varNames_in,'Miscanthus'))
    if any(any(array_in(:,:,contains(varNames_in,'Miscanthus'),:)>0))
        error('This script assumes zero Miscanthus area!')
    end
    array_in(:,:,contains(varNames_in,'Miscanthus'),:) = [] ;
    varNames_in(contains(varNames_in,'Miscanthus')) = [] ;
end


end