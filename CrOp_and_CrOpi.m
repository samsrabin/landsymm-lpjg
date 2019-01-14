function [S, ...
    inds_cropTypes_orig, inds_cropTypesI_orig, inds_cropTypes_nonCGs_orig ...
    ] = ...
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
addOptional(p,'inds_cropTypesNonCGs_cropFracsOrig',[],isindices) ;
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
elseif strcmp(structname,'gsirrig') || strcmp(structname,'yield')|| strcmp(structname,'Nfert')
    if isempty(cropfracs_orig)
        error('When S = gsirrig, yield, or Nfert, you must provide cropfracs_orig!') ;
    elseif isempty(inds_cropTypes_cropFracsOrig)
        error('When S = gsirrig, yield, or Nfert, you must provide inds_cropTypes_cropFracsOrig!') ;
    elseif isempty(inds_cropTypesI_cropFracsOrig)
        error('When S = gsirrig, yield, or Nfert, you must provide inds_cropTypesI_cropFracsOrig!') ;
    elseif isempty(inds_cropTypesNonCGs_cropFracsOrig)
        error('When S = gsirrig, yield, or Nfert, you must provide inds_cropTypesNonCGs_cropFracsOrig!') ;
    end
end

[cropTypes_found,inds_cropTypes] = intersect(S.varNames,cropTypes) ;
if length(cropTypes) ~= length(inds_cropTypes)
    error('Each CrOp in cropTypes not found exactly once in S.varNames!')
end
cropTypesI = strcat(cropTypes,'i') ;
[cropTypesI_found,inds_cropTypesI] = intersect(S.varNames,cropTypesI) ;
inds_cropTypes_orig = inds_cropTypes ;
inds_cropTypesI_orig = inds_cropTypesI ;
cropgrass_onlyRF = false ;
if length(cropTypesI) ~= length(inds_cropTypesI)
    % endsWith() instead of contains() correctly returns FALSE for CC*G_ic
    if length(find(endsWith(S.varNames,{'CC3G','CC4G'})))==2 && ~any(contains(S.varNames,{'CC3Gi','CC4Gi'}))
        warning('CC3Gi and CC4Gi not present. This SHOULD be okay.')
        cropgrass_onlyRF = true ;
    elseif length(find(strcmp(S.varNames,'ExtraCrop')))==1 && ~any(contains(S.varNames,{'ExtraCropi'}))
        warning('ExtraCropi not present. This SHOULD be okay.')
        cropgrass_onlyRF = true ;
    else
        error('Each CrOpi in cropTypesI not found exactly once in S.varNames!')
    end
end

if strcmp(merge_or_replace,'replace')
    warning('%s: Replacing CrOp with CrOpi; removing CrOpi.',structname)
    if cropgrass_onlyRF
        S.varNames(strcmp(S.varNames,'CC3G')) = {'CC3Gi'} ;
        S.varNames(strcmp(S.varNames,'CC4G')) = {'CC4Gi'} ;
        S.varNames(strcmp(S.varNames,'ExtraCrop')) = {'ExtraCropi'} ;
        [~,inds_cropTypes] = intersect(S.varNames,cropTypes) ;
    end
    % Delete CrOps from maps
    S.maps_YXvy(:,:,inds_cropTypes,:) = [] ;
    S.varNames(inds_cropTypes) = [] ;
    % Rename CrOpis to CrOps
    S.varNames = strip(S.varNames,'right','i') ;
    % Sanity check
    [~,inds_cropTypes_new] = intersect(S.varNames,cropTypes) ;
    if nansum(nansum(nansum(nansum(S.maps_YXvy(:,:,inds_cropTypes_new,:)))))<=0
        error(['Something went wrong in ' structname ': replace_CrOp_with_CrOpi!'])
    end
elseif strcmp(merge_or_replace,'merge')
    warning('%s: Merging CrOpi into CrOp; removing CrOpi.',structname)
    % Merge
    if strcmp(structname,'fpc') || strcmp(structname,'cropfracs') || strcmp(structname,'cropfracs_plum7')
        % Merge (straight sum)
        if cropgrass_onlyRF
            cropTypes_found_noCGs = strip(cropTypesI_found,'right','i') ;
            [~,inds_cropTypes_nonCGs] = intersect(cropTypes_found,cropTypes_found_noCGs) ;
            [~,inds_cropTypes_nonCGs_orig] = intersect(S.varNames,cropTypes_found_noCGs) ;
            S.maps_YXvy(:,:,inds_cropTypes_nonCGs,:) = ...
                S.maps_YXvy(:,:,inds_cropTypes_nonCGs,:) ...
                + S.maps_YXvy(:,:,inds_cropTypesI,:) ;
        else
            S.maps_YXvy(:,:,inds_cropTypes,:) = ...
                S.maps_YXvy(:,:,inds_cropTypes,:) ...
                + S.maps_YXvy(:,:,inds_cropTypesI,:) ;
        end
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
    elseif strcmp(structname,'gsirrig') || strcmp(structname,'yield') || strcmp(structname,'Nfert')
        % Merge (weight by area)
        tmp_frac_theseCrops_YXcy = cropfracs_orig.maps_YXvy(:,:,inds_cropTypesNonCGs_cropFracsOrig,:) ;
        tmp_frac_theseCropsI_YXcy = cropfracs_orig.maps_YXvy(:,:,inds_cropTypesI_cropFracsOrig,:) ;
        tmp_frac_theseBoth_YXcy = tmp_frac_theseCrops_YXcy + tmp_frac_theseCropsI_YXcy ;
        
        tmp_frac_theseBoth_YXcy(tmp_frac_theseBoth_YXcy==0 | isnan(tmp_frac_theseBoth_YXcy)) = 1 ; % Ensure no division by zero or NaN
        
        S.maps_YXvy(:,:,inds_cropTypesNonCGs_cropFracsOrig,:) = ...
            S.maps_YXvy(:,:,inds_cropTypesNonCGs_cropFracsOrig,:) .* tmp_frac_theseCrops_YXcy./tmp_frac_theseBoth_YXcy...
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
if ~exist('inds_cropTypes_nonCGs_orig','var')
    inds_cropTypes_nonCGs_orig = inds_cropTypes_orig ;
end


end