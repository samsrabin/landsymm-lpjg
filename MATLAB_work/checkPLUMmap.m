function checkPLUMmap(PLUMtoLPJG,PLUMcrops,varargin)
% Check that every PLUM crop is mapped

fake_fruitveg_sugar = false ;
if ~isempty(varargin)
    fake_fruitveg_sugar = varargin{1} ;
    if length(varargin) > 1
        error('checkPLUMmap takes at most 1 optional argument (fake_fruitveg_sugar)')
    end
end

wheatfound2x = false ;
setasidefound2x = false ;

Ncrops_lpjg = length(PLUMtoLPJG) ;
Ncfts_plum = length(PLUMcrops) ;

for cP = 1:Ncfts_plum
    thisCrop_plum = PLUMcrops{cP} ;
    found1x = false ;
    for cL = 1:Ncrops_lpjg
        Nfound = length(find(strcmp(PLUMtoLPJG{cL},thisCrop_plum))) ;
        if Nfound==1
            if ~found1x
                found1x = true ;
            elseif strcmp(thisCrop_plum,'wheat')
                if wheatfound2x
                    error('wheat found in more than two elements of PLUMtoLPJG!') ;
                else
%                     warning('wheat found in two elements of PLUMtoLPJG so far...')
                    wheatfound2x = true ;
                end
            elseif strcmp(thisCrop_plum,'setaside')
                if setasidefound2x
                    error('setaside found in more than two elements of PLUMtoLPJG!') ;
                else
%                     warning('setaside found in two elements of PLUMtoLPJG so far...')
                    setasidefound2x = true ;
                end
            else
                error([thisCrop_plum ' found in more than one element of PLUMtoLPJG!']) ;
            end
        elseif Nfound > 1
            error([thisCrop_plum ' found more than once in PLUMtoLPJG{' num2str(cL) '}!'])
        end
    end ; clear cL
    if ~found1x && ~(fake_fruitveg_sugar && (strcmp(thisCrop_plum, 'fruitveg') || strcmp(thisCrop_plum, 'sugar')))
        error([thisCrop_plum ' not found in PLUMtoLPJG!']) ;
    end
    clear thisCrop_plum found
end ; clear cP
if wheatfound2x && setasidefound2x
    disp('Every PLUM crop is mapped once (except for wheat and setaside each being mapped twice, correctly).')
elseif wheatfound2x
    disp('Every PLUM crop is mapped once (except for wheat being mapped twice, correctly).')
elseif setasidefound2x
    disp('Every PLUM crop is mapped once (except for setaside being mapped twice, correctly).')
else
    disp('Every PLUM crop is mapped once.')
end

% Check that each LPJGcrop has a corresponding element of PLUMtoLPJG
if length(PLUMtoLPJG) ~= Ncrops_lpjg
    error('length(PLUMtoLPJG) ~= Ncrops_lpjg') ;
end
disp('Each LPJG crop has a corresponding member of PLUMtoLPJG.')