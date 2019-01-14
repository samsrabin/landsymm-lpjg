function [S, inds_cropTypes_orig, inds_cropTypesI_orig] = ...
    CrOp_and_CrOpi(S, structname, ...
    cropTypes, merge_or_replace, ...
    varargin)
    
% intersect() checked 2018-05-11

% is1x1logical = @(x) numel(x)==1 & (x==0 | x==1) ;
isindices = @(x) all(isint(x)) & ismatrix(x) & min(size(x))==1 ;

p = inputParser ;
addRequired(p,'S',@isstruct) ;
addRequired(p,'structname',@ischar) ;
addRequired(p,'cropTypes',@iscellstr) ;
addRequired(p,'merge_or_replace',@ischar) ;
addOptional(p,'cropfile','',@ischar) ;
addOptional(p,'cropfracs_orig',struct([]),@isstruct) ;
addOptional(p,'inds_cropTypes_cropFracsOrig',[],isindices) ;
addOptional(p,'inds_cropTypesI_cropFracsOrig',[],isindices) ;
parse(p,...
    S, structname, ...
    cropTypes, merge_or_replace, ...
    varargin{:});

pFields = fieldnames(p.Results) ;
Nfields = length(pFields) ;
for f = 1:Nfields
    thisField = pFields{f} ;
    if ~exist(thisField,'var')
        eval([thisField ' = p.Results.' thisField ' ;']) ;
    end
    clear thisField
end ; clear f
clear p

if (strcmp(structname,'cropfracs') || strcmp(structname,'cropfracs_plum7')) && isempty(cropfile)
%     error('When S = cropfracs or cropfracs_plum7, you must provide cropfile!') ;
    warning('When S = cropfracs or cropfracs_plum7, you must provide cropfile! But not really, since I''m going to error() anyway.') ;
elseif (strcmp(structname,'gsirrig') || strcmp(structname,'yield'))
    if isempty(cropfracs_orig)
        error('When S = gsirrig or yield, you must provide cropfracs_orig!') ;
    elseif isempty(inds_cropTypes_cropFracsOrig)
        error('When S = gsirrig or yield, you must provide inds_cropTypes_cropFracsOrig!') ;
    elseif isempty(inds_cropTypesI_cropFracsOrig)
        error('When S = gsirrig or yield, you must provide inds_cropTypesI_cropFracsOrig!') ;
    end
end


if strcmp(merge_or_replace,'replace')
    warning('Replacing CrOp with CrOpi; removing CrOpi.')
    [~,inds_cropTypes] = intersect(S.varNames,cropTypes) ; % Using intersect() fixes any mis-ordering of CFTs in S relative to cropTypes{}
    if length(cropTypes) ~= length(inds_cropTypes)
        error('Each CrOp in cropTypes not found exactly once in S.varNames!')
    end
    cropTypesI = strcat(cropTypes,'i') ;
    [~,inds_cropTypesI] = intersect(S.varNames,cropTypesI) ;
    if length(cropTypesI) ~= length(inds_cropTypesI)
        error('Each CrOpi in cropTypesI not found exactly once in S.varNames!')
    end
    % Delete CrOps from maps
    S.maps_YXvy(:,:,inds_cropTypes,:) = [] ;
    S.varNames(inds_cropTypes) = [] ;
    % Rename CrOpis to CrOps
    for c = 1:length(cropTypes)
        S.varNames{strcmp(S.varNames,cropTypesI{c})} = cropTypes{c} ;
    end
    % Sanity check
    [~,inds_cropTypes_new] = intersect(S.varNames,cropTypes) ;
    if nansum(nansum(nansum(nansum(S.maps_YXvy(:,:,inds_cropTypes_new,:)))))<=0
        error(['Something went wrong in ' structname ': replace_CrOp_with_CrOpi!'])
    end
elseif strcmp(merge_or_replace,'merge')
    warning('Merging CrOpi into CrOp; removing CrOpi.')
    [~,inds_cropTypes] = intersect(S.varNames,cropTypes) ;
    if length(cropTypes) ~= length(inds_cropTypes)
        error('Each CrOp in cropTypes not found exactly once in S.varNames!')
    end
    cropTypesI = strcat(cropTypes,'i') ;
    [~,inds_cropTypesI] = intersect(S.varNames,cropTypesI) ;
    inds_cropTypes_orig = inds_cropTypes ;
    inds_cropTypesI_orig = inds_cropTypesI ;
    if length(cropTypesI) ~= length(inds_cropTypesI)
        error('Each CrOpi in cropTypesI not found exactly once in S.varNames!')
    end
    % Merge
    if strcmp(structname,'fpc') || strcmp(structname,'cropfracs') || strcmp(structname,'cropfracs_plum7')
        % Merge
        S.maps_YXvy(:,:,inds_cropTypes,:) = ...
            S.maps_YXvy(:,:,inds_cropTypes,:) ...
            + S.maps_YXvy(:,:,inds_cropTypesI,:) ;
        if any(any(any(any(-1+S.maps_YXvy(:,:,inds_cropTypesI,:) ...
                > 1e-6))))
            error('Think some more about how to merge CrOpi into CrOp!')
        end
        % Delete CrOpis
        S.maps_YXvy(:,:,inds_cropTypesI,:) = [] ;
        S.varNames(inds_cropTypesI) = [] ;
        % Sanity check
        [~,inds_cropTypes_new] = intersect(S.varNames,cropTypes) ;
        if nansum(nansum(nansum(nansum(S.maps_YXvy(:,:,inds_cropTypes_new,:)))))<=0
            error(['Something went wrong in ' structname ': merge_CrOpi_into_CrOp!'])
        end
    elseif strcmp(structname,'gsirrig') || strcmp(structname,'yield')
        % Merge (weight by area)
        tmp_frac_theseCrops_YXcy = cropfracs_orig.maps_YXvy(:,:,inds_cropTypes_cropFracsOrig,:) ;
        tmp_frac_theseCropsI_YXcy = cropfracs_orig.maps_YXvy(:,:,inds_cropTypesI_cropFracsOrig,:) ;
        tmp_frac_theseBoth_YXcy = tmp_frac_theseCrops_YXcy + tmp_frac_theseCropsI_YXcy ;
        S.maps_YXvy(:,:,inds_cropTypes,:) = ...
            S.maps_YXvy(:,:,inds_cropTypes,:) .* tmp_frac_theseCrops_YXcy./tmp_frac_theseBoth_YXcy...
          + S.maps_YXvy(:,:,inds_cropTypesI,:) .* tmp_frac_theseCropsI_YXcy./tmp_frac_theseBoth_YXcy ;
        % Delete CrOpis
        S.maps_YXvy(:,:,inds_cropTypesI,:) = [] ;
        S.varNames(inds_cropTypesI) = [] ;
    else
        error(['How am I supposed to merge with this structname? (' structname ')']) ;
    end
else
    error(['merge_or_replace not recognized: ' merge_or_replace])
end

if ~exist('inds_cropTypes_orig','var')
    inds_cropTypes_orig = inds_cropTypes ;
end
if ~exist('inds_cropTypesI_orig','var')
    inds_cropTypesI_orig = inds_cropTypes ;
end


end