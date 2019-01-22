%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read MAT-files from harmonization; write as LPJG inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_gzip = false ;

% Directories for harmonized PLUM outputs
%dirList = {...
%              'SSP1.v10.s1.harm' ;
%              'SSP3.v10.s1.harm' ;
%              'SSP4.v10.s1.harm' ;
%              'SSP5.v10.s1.harm' ;
%              } ;
dirList = {...
%    'ssp11/SSP1/s1' ;
%    'ssp11/SSP1/s2' ;
%    'ssp11/SSP1/s3' ;
%    'ssp11/SSP1/s4' ;
%    'ssp11/SSP1/s5' ;
%    'ssp11/SSP1/s6' ;
%    'ssp11/SSP1/s7' ;
%    'ssp11/SSP1/s8' ;
%    'ssp11/SSP1/s9' ;
%    'ssp11/SSP1/s10' ;
%    'ssp11/SSP1/s11' ;
%    'ssp11/SSP1/s12' ;
%    'ssp11/SSP1/s13' ;
%    'ssp11/SSP1/s14' ;
%    'ssp11/SSP1/s15' ;
%    'ssp11/SSP1/s16' ;
%    'ssp11/SSP1/s17' ;
%    'ssp11/SSP1/s18' ;
%    'ssp11/SSP1/s19' ;
%    'ssp11/SSP1/s20' ;
%    'ssp11/SSP1/s21' ;
%    'ssp11/SSP1/s22' ;
%    'ssp11/SSP1/s23' ;
%    'ssp11/SSP1/s24' ;
%    'ssp11/SSP1/s25' ;
%    'ssp11/SSP1/s26' ;
%    'ssp11/SSP1/s27' ;
%    'ssp11/SSP1/s27' ;
%    'ssp11/SSP1/s28' ;
%    'ssp11/SSP1/s29' ;
%    'ssp11/SSP1/s30' ;
%    'ssp11/SSP1/s31' ;
%    'ssp11/SSP1/s32' ;
%    'ssp11/SSP1/s33' ;
%    'ssp11/SSP1/s34' ;
%    'ssp11/SSP1/s35' ;
%    'ssp11/SSP1/s36' ;
%    'ssp11/SSP1/s37' ;
%    'ssp11/SSP1/s38' ;
%    'ssp11/SSP1/s39' ;
%    'ssp11/SSP1/s40' ;
%    'ssp11/SSP1/s41' ;
%    'ssp11/SSP1/s42' ;
%    'ssp11/SSP1/s43' ;
%    'ssp11/SSP1/s44' ;
%    'ssp11/SSP1/s45' ;
%    'ssp11/SSP1/s46' ;
%    'ssp11/SSP1/s47' ;
%    'ssp11/SSP1/s48' ;
%    'ssp11/SSP1/s49' ;
%    'ssp11/SSP1/s50' ;
    'ssp11/SSP3/s1' ;
    'ssp11/SSP3/s2' ;
    'ssp11/SSP3/s3' ;
    'ssp11/SSP3/s4' ;
    'ssp11/SSP3/s5' ;
%    'ssp11/SSP3/s6' ;
%    'ssp11/SSP3/s7' ;
%    'ssp11/SSP3/s8' ;
%    'ssp11/SSP3/s9' ;
%    'ssp11/SSP3/s10' ;
%    'ssp11/SSP3/s11' ;
%    'ssp11/SSP3/s12' ;
%    'ssp11/SSP3/s13' ;
%    'ssp11/SSP3/s14' ;
%    'ssp11/SSP3/s15' ;
%    'ssp11/SSP3/s16' ;
%    'ssp11/SSP3/s17' ;
%    'ssp11/SSP3/s18' ;
%    'ssp11/SSP3/s19' ;
%    'ssp11/SSP3/s20' ;
%    'ssp11/SSP3/s21' ;
%    'ssp11/SSP3/s22' ;
%    'ssp11/SSP3/s23' ;
%    'ssp11/SSP3/s24' ;
%    'ssp11/SSP3/s25' ;
%    'ssp11/SSP3/s26' ;
%    'ssp11/SSP3/s27' ;
%    'ssp11/SSP3/s28' ;
%    'ssp11/SSP3/s29' ;
%    'ssp11/SSP3/s30' ;
%    'ssp11/SSP3/s31' ;
%    'ssp11/SSP3/s32' ;
%    'ssp11/SSP3/s33' ;
%    'ssp11/SSP3/s34' ;
%    'ssp11/SSP3/s35' ;
%    'ssp11/SSP3/s36' ;
%    'ssp11/SSP3/s37' ;
%    'ssp11/SSP3/s38' ;
%    'ssp11/SSP3/s39' ;
%    'ssp11/SSP3/s40' ;
%    'ssp11/SSP3/s41' ;
%    'ssp11/SSP3/s42' ;
%    'ssp11/SSP3/s43' ;
%    'ssp11/SSP3/s44' ;
%    'ssp11/SSP3/s45' ;
%    'ssp11/SSP3/s46' ;
%    'ssp11/SSP3/s47' ;
%    'ssp11/SSP3/s48' ;
%    'ssp11/SSP3/s49' ;
%    'ssp11/SSP3/s50' ;
    'ssp11/SSP4/s1' ;
    'ssp11/SSP4/s2' ;
    'ssp11/SSP4/s3' ;
    'ssp11/SSP4/s4' ;
    'ssp11/SSP4/s5' ;
%    'ssp11/SSP4/s6' ;
%    'ssp11/SSP4/s7' ;
%    'ssp11/SSP4/s8' ;
%    'ssp11/SSP4/s9' ;
%    'ssp11/SSP4/s10' ;
%    'ssp11/SSP4/s11' ;
%    'ssp11/SSP4/s12' ;
%    'ssp11/SSP4/s13' ;
%    'ssp11/SSP4/s14' ;
%    'ssp11/SSP4/s15' ;
%    'ssp11/SSP4/s16' ;
%    'ssp11/SSP4/s17' ;
%    'ssp11/SSP4/s18' ;
%    'ssp11/SSP4/s19' ;
%    'ssp11/SSP4/s20' ;
%    'ssp11/SSP4/s21' ;
%    'ssp11/SSP4/s22' ;
%    'ssp11/SSP4/s23' ;
%    'ssp11/SSP4/s24' ;
%    'ssp11/SSP4/s25' ;
%    'ssp11/SSP4/s26' ;
%    'ssp11/SSP4/s27' ;
%    'ssp11/SSP4/s28' ;
%    'ssp11/SSP4/s29' ;
%    'ssp11/SSP4/s30' ;
%    'ssp11/SSP4/s31' ;
%    'ssp11/SSP4/s32' ;
%    'ssp11/SSP4/s33' ;
%    'ssp11/SSP4/s34' ;
%    'ssp11/SSP4/s35' ;
%    'ssp11/SSP4/s36' ;
%    'ssp11/SSP4/s37' ;
%    'ssp11/SSP4/s38' ;
%    'ssp11/SSP4/s39' ;
%    'ssp11/SSP4/s40' ;
%    'ssp11/SSP4/s41' ;
%    'ssp11/SSP4/s42' ;
%    'ssp11/SSP4/s43' ;
%    'ssp11/SSP4/s44' ;
%    'ssp11/SSP4/s45' ;
%    'ssp11/SSP4/s46' ;
%    'ssp11/SSP4/s47' ;
%    'ssp11/SSP4/s48' ;
%    'ssp11/SSP4/s49' ;
%    'ssp11/SSP4/s50' ;
    'ssp11/SSP5/s1' ;
    'ssp11/SSP5/s2' ;
    'ssp11/SSP5/s3' ;
    'ssp11/SSP5/s4' ;
    'ssp11/SSP5/s5' ;
%    'ssp11/SSP5/s6' ;
%    'ssp11/SSP5/s7' ;
%    'ssp11/SSP5/s8' ;
%    'ssp11/SSP5/s9' ;
%    'ssp11/SSP5/s10' ;
%    'ssp11/SSP5/s11' ;
%    'ssp11/SSP5/s12' ;
%    'ssp11/SSP5/s13' ;
%    'ssp11/SSP5/s14' ;
%    'ssp11/SSP5/s15' ;
%    'ssp11/SSP5/s16' ;
%    'ssp11/SSP5/s17' ;
%    'ssp11/SSP5/s18' ;
%    'ssp11/SSP5/s19' ;
%    'ssp11/SSP5/s20' ;
%    'ssp11/SSP5/s21' ;
%    'ssp11/SSP5/s22' ;
%    'ssp11/SSP5/s23' ;
%    'ssp11/SSP5/s24' ;
%    'ssp11/SSP5/s25' ;
%    'ssp11/SSP5/s26' ;
%    'ssp11/SSP5/s27' ;
%    'ssp11/SSP5/s28' ;
%    'ssp11/SSP5/s29' ;
%    'ssp11/SSP5/s30' ;
%    'ssp11/SSP5/s31' ;
%    'ssp11/SSP5/s32' ;
%    'ssp11/SSP5/s33' ;
%    'ssp11/SSP5/s34' ;
%    'ssp11/SSP5/s35' ;
%    'ssp11/SSP5/s36' ;
%    'ssp11/SSP5/s37' ;
%    'ssp11/SSP5/s38' ;
%    'ssp11/SSP5/s39' ;
%    'ssp11/SSP5/s40' ;
%    'ssp11/SSP5/s41' ;
%    'ssp11/SSP5/s42' ;
%    'ssp11/SSP5/s43' ;
%    'ssp11/SSP5/s44' ;
%    'ssp11/SSP5/s45' ;
%    'ssp11/SSP5/s46' ;
%    'ssp11/SSP5/s47' ;
%    'ssp11/SSP5/s48' ;
%    'ssp11/SSP5/s49' ;
%    'ssp11/SSP5/s50' ;
    } ;
          
base_year = 2010 ;
y1 = 2011 ;
yN = 2100 ;
yStep = 1 ;
                        
% Trying to avoid new crop spinup time
y1_pre = 2006 ;    % Will repeat first PLUMout year for y1_pre:(y1-1)

% Make it so that each gridcell always has at least some tiny amount of
% every crop?
outPrec = 6 ;
mincropfrac = 10^-outPrec ;
someofall = true ; 
donation_order = {'PASTURE','NATURAL','BARREN'} ;
              

%% Setup

disp('Setting up...')

% Determine which system you're on
tmp = pwd ;
if strcmp(tmp(1:5),'/User')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp

cf_kgNha_kgNm2 = 1e-4 ;

outWidth = 1 ;
delimiter = ' ' ;
overwrite = true ;
fancy = false ;
save_every_pct = 1 ;
verbose_write = true ;

if onMac
    addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')
else
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/paper02-matlab-work/')) ;
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/matlab-general/')) ;
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/matlab-general-fromshared/')) ;
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/plum_harmonization/')) ;
    addpath(genpath('/pfs/data1/home/kit/imk-ifu/lr8247/lpj-guess-crop-calibration/')) ;
end

% Get year info
yearList = shiftdim(y1:yStep:yN) ;
Nyears = length(yearList) ;
if y1_pre<y1
    yearList_xtra = y1_pre:(y1-1) ;
else
    yearList_xtra = [] ;
end
Nyears_xtra = length(yearList_xtra) ;

% Get cells present in previous LPJ-GUESS output
lpjg_in = lpjgu_matlab_readTable('/Users/Shared/PLUM/trunk_runs/LPJGPLUM_expt1.1_2006-2100_PLUM6xtra_20180412171654/rcp26/2011-2015/yield.out.gz','dont_save_MAT',true,'verboseIfNoMat',false) ;
lons_lpjg = lpjg_in.Lon ;
lats_lpjg = lpjg_in.Lat ;
clear lpjg_in

% Get LUH2 land area (m2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = 1e6*transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
landArea_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%%% Convert to half-degree
tmp = landArea_YXqd(:,1:2:1440) + landArea_YXqd(:,2:2:1440) ;
landArea_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
clear *_YXqd



%% Do it

for d = 1:length(dirList)
    
    % Get directories
    inDir = find_PLUM2LPJG_inputs(dirList{d}) ;
    inDir = removeslashifneeded(inDir) ;
    disp(inDir)
    outDir = addslashifneeded([removeslashifneeded(inDir) '.forLPJG']) ;
    unix(['mkdir -p ' outDir]) ;
        
    
    %%%%%%%%%%%%%%
    %%% Import %%%
    %%%%%%%%%%%%%%
    
    disp('Importing harmonization outputs...')
    for y = 1:Nyears
        thisYear = yearList(y) ;
        disp(['   Reading ' num2str(thisYear) '...'])
        
        % Import this year's harmonized outputs
        load(sprintf('%s/%d/LandCoverFract.base%d.mat', inDir, thisYear, base_year)) ;
        S_lu = out_y1 ; clear out_y1
        load(sprintf('%s/%d/CropFract.base%d.mat', inDir, thisYear, base_year)) ;
        S_cf = out_y1 ; clear out_y1
        load(sprintf('%s/%d/Fert.base%d.mat', inDir, thisYear, base_year)) ;
        S_nf = out_y1 ; clear out_y1
        load(sprintf('%s/%d/Irrig.base%d.mat', inDir, thisYear, base_year)) ;
        S_ir = out_y1 ; clear out_y1
        
        % Set up empty arrays and do other initialization
        if y==1
            list2map = find(~isnan(S_lu.maps_YXv(:,:,1))>0) ;
            Ncells = length(list2map) ;
            lu_in.varNames = S_lu.varNames ;
            lu_in.yearList = yearList ;
            lu_in.maps_YXvy = nan([size(S_lu.maps_YXv) Nyears]) ;
            cropList_rf = S_cf.varNames ;
            cropList_ir = strcat(cropList_rf,'i') ;
            cf_in.varNames = cropList_ir ;
            nf_in.varNames = cropList_ir ;
            ir_in.varNames = cropList_ir ;
            cf_in.yearList = yearList ;
            nf_in.yearList = yearList ;
            ir_in.yearList = yearList ;
            cf_in.maps_YXvy = nan([size(S_cf.maps_YXv) Nyears]) ;
            nf_in.maps_YXvy = nan([size(S_cf.maps_YXv) Nyears]) ;
            ir_in.maps_YXvy = nan([size(S_cf.maps_YXv) Nyears]) ;
        end
                
        % Put this year's harmonized outputs into the structure
        lu_in.maps_YXvy(:,:,:,y) = S_lu.maps_YXv ;
        cf_in.maps_YXvy(:,:,:,y) = S_cf.maps_YXv ;
        nf_in.maps_YXvy(:,:,:,y) = S_nf.maps_YXv ;
        ir_in.maps_YXvy(:,:,:,y) = S_ir.maps_YXv ;
        
    end
    
    % Get arrays
    disp('Getting arrays...')
    [lu_out, lu_header_cell] = lpjgu_matlab_maps2table(lu_in,list2map) ;
    [cf_out, cf_header_cell] = lpjgu_matlab_maps2table(cf_in,list2map) ;
    [nf_out, nf_header_cell] = lpjgu_matlab_maps2table(nf_in,list2map) ;
    [ir_out, ir_header_cell] = lpjgu_matlab_maps2table(ir_in,list2map) ;
        
    % Check for NaNs
    disp('Checking arrays...')
    if any(any(isnan(lu_out)))
        error('Some value of lu_out is NaN.')
    elseif any(any(isnan(cf_out)))
        error('Some value of cf_out is NaN.')
    elseif any(any(isnan(nf_out)))
        error('Some value of nf_out is NaN.')
    elseif any(any(isnan(ir_out)))
        error('Some value of ir_out is NaN.')
    end
    
    % Check for nonsensical values
    [~,IA,~] = intersect(nf_header_cell,cropList_ir,'stable') ;
    if min(min(lu_out(:,4:end)))<0
        error('min(lu_out)<0')
    elseif max(max(lu_out(:,4:end)))>=1+10^(-outPrec)
        error('max(lu_out)>1')
    elseif max(sum(lu_out(:,4:end),2))>=1+10^(-outPrec)
        error('max(sum_lu)>1')
    elseif min(min(cf_out(:,4:end)))<0
        error('min(cf_out)<0')
    elseif max(max(cf_out(:,4:end)))>=1+10^(-outPrec)
        error('max(cf_out)>1')
    elseif max(sum(cf_out(:,IA),2))>=1+10^(-outPrec)
        error('max(sum_cf)>1')
    elseif min(min(nf_out(:,IA)))<0
        error('min(nf_out)<0')
    elseif min(min(ir_out(:,IA)))<0
        error('min(ir_out)<0')
    elseif max(max(ir_out(:,IA)))>=1+10^(-outPrec)
        error('max(ir_out)>1')
    end
    
    disp('Additional processing...')
    
    % If doing so, make sure that every grid cell has at least some of each
    % type of cropland
    if mincropfrac>0
        
        % Some cropland
        iCrop = find(strcmp(lu_in.varNames,'CROPLAND')) ;
        lu_tmp = lu_out(:,4:end) ;
        lu_tmp = round(lu_tmp,outPrec) ;
        no_cropland = lu_tmp(:,iCrop)<mincropfrac ;
        i = 0 ;
        while(any(no_cropland))
            i = i+1 ;
            if i > length(donation_order)
                error('GET CROPLAND FROM SOMEWHERE')
            end
            this_donor = donation_order{i} ;
            iThis = find(strcmp(lu_in.varNames,this_donor)) ;
            involved = lu_tmp(:,iThis)>=mincropfrac & no_cropland ;
            if any(involved)
                warning('Giving some from %s to CROPLAND (%d cells).', this_donor, length(find(involved)))
                transfer_amt = mincropfrac-lu_tmp(involved,iCrop) ;
                lu_tmp(involved,iThis) = lu_tmp(involved,iThis) - transfer_amt ;
                lu_tmp(involved,iCrop) = lu_tmp(involved,iCrop) + transfer_amt ;
                no_cropland = lu_tmp(:,iCrop)==0 ;
            end
        end
        lu_out(:,4:end) = lu_tmp ;
        
        % Some of each type: Where there was no cropland
        cf_tmp = cf_out(:,IA) ;
        Ncrops = length(IA) ;
        cf_tmp(repmat(sum(cf_tmp,2),[1 Ncrops])==0) = 1/Ncrops ;
        if max(sum(cf_tmp,2))>=1+2*10^(-outPrec)
            error('max(round(sum_cf))>1')
        end
        
        % Some of each type: Everywhere else
        for c = 1:Ncrops
            % Find rows with 0 for this crop
            thiscrop = cf_tmp(:,c) ;
            iszerothiscrop = thiscrop<mincropfrac ;
            % Find the crop that currently has the greatest area
            maxcropfrac_Xv = repmat(max(cf_tmp,[],2),[1 Ncrops]) ;
            ismaxcropfrac_Xv = cf_tmp==maxcropfrac_Xv ;
            % Make sure there's only one ismaxcropfrac in each row
            for i = fliplr(2:Ncrops)
                tmp = ismaxcropfrac_Xv(:,i) ;
                sumtoleft = sum(ismaxcropfrac_Xv(:,1:(i-1)),2) ;
                tmp(sumtoleft>0) = false ;
                ismaxcropfrac_Xv(:,i) = tmp ;
            end
            cf_tmp(ismaxcropfrac_Xv & iszerothiscrop) = cf_tmp(ismaxcropfrac_Xv & iszerothiscrop) - mincropfrac ;
            cf_tmp(iszerothiscrop,c) = cf_tmp(iszerothiscrop,c) + mincropfrac ;
        end
        
        % Save
        cf_out(:,IA) = cf_tmp ;
%         clear cf_tmp Ncrops
        
    end
    
    
    % Check for nonsensical values
    [~,IA,~] = intersect(nf_header_cell,cropList_ir,'stable') ;
    if min(min(lu_out(:,4:end)))<0
        error('min(lu_out)<0')
    elseif max(max(lu_out(:,4:end)))>=1+10^(-outPrec)
        error('max(lu_out)>1')
    elseif max(sum(lu_out(:,4:end),2))>=1+2*10^(-outPrec)
        error('max(sum_lu)>1')
    elseif min(min(cf_out(:,4:end)))<mincropfrac
        error('min(cf_out)<mincropfrac')
    elseif max(max(cf_out(:,4:end)))>=1+10^(-outPrec)
        error('max(cf_out)>1')
    elseif max(sum(cf_out(:,4:end),2))>=1+1*10^(-outPrec)
        error('max(sum_cf)>1')
    elseif min(min(nf_out(:,IA)))<0
        error('min(nf_out)<0')
    elseif min(min(ir_out(:,IA)))<0
        error('min(ir_out)<0')
    elseif max(max(ir_out(:,IA)))>=1+10^(-outPrec)
        error('max(ir_out)>1')
    end
    
    % Convert ktN/ha to kgN/m2
    [~,IA,~] = intersect(nf_header_cell,cropList_ir,'stable') ;
    nf_out(:,IA) = cf_kgNha_kgNm2 * nf_out(:,IA) ;
    
    % Add rainfed crops (zeros)
    cf_out = [cf_out zeros(size(cf_out(:,IA)))] ;
    nf_out = [nf_out zeros(size(nf_out(:,IA)))] ;
    ir_out = [ir_out zeros(size(ir_out(:,IA)))] ;
    cf_header_cell = [cf_header_cell cropList_rf] ;
    nf_header_cell = [nf_header_cell cropList_rf] ;
    ir_header_cell = [ir_header_cell cropList_rf] ;
    
    % Remove ExtraCropi, because it receives no management inputs, and
    % did not exist in historical remap_v4.
    % It did have some area. Put that in ExtraCrop.
    cf_out(:,strcmp(cf_header_cell,'ExtraCrop')) = cf_out(:,strcmp(cf_header_cell,'ExtraCropi')) ;
    cf_out(:,strcmp(cf_header_cell,'ExtraCropi')) = [] ;
    nf_out(:,strcmp(nf_header_cell,'ExtraCropi')) = [] ;
    ir_out(:,strcmp(ir_header_cell,'ExtraCropi')) = [] ;
    cf_header_cell(strcmp(cf_header_cell,'ExtraCropi')) = [] ;
    nf_header_cell(strcmp(nf_header_cell,'ExtraCropi')) = [] ;
    ir_header_cell(strcmp(ir_header_cell,'ExtraCropi')) = [] ;
    cropList_ir(strcmp(cropList_ir,'ExtraCropi')) = [] ;
    
    % Add extra years, if doing so
    if Nyears_xtra>0
        
        lons = lu_out(1:Nyears:end,1) ;
        lons_tmp = repmat(lons,[1 Nyears+Nyears_xtra])' ;
        lons_out = lons_tmp(:) ;
        lats = lu_out(1:Nyears:end,2) ;
        lats_tmp = repmat(lats,[1 Nyears+Nyears_xtra])' ;
        lats_out = lats_tmp(:) ;
        
        years_out = repmat([yearList_xtra';yearList],[Ncells 1]) ;
        
        Nvar_lu = length(lu_in.varNames) ;
        tmp1_yxv = reshape(lu_out(:,4:end),[Nyears Ncells Nvar_lu]) ;
        tmp1_yxv = cat(1, repmat(tmp1_yxv(1,:,:),[Nyears_xtra 1 1]), tmp1_yxv) ;
        tmp1_out = reshape(tmp1_yxv,[Ncells*(Nyears+Nyears_xtra) Nvar_lu]) ;
        lu_out = [lons_out lats_out years_out tmp1_out] ;
        
        Nvar_cf = size(cf_out,2) - 3 ;
        tmp1_yxv = reshape(cf_out(:,4:end),[Nyears Ncells Nvar_cf]) ;
        tmp1_yxv = cat(1, repmat(tmp1_yxv(1,:,:),[Nyears_xtra 1 1]), tmp1_yxv) ;
        tmp1_out = reshape(tmp1_yxv,[Ncells*(Nyears+Nyears_xtra) Nvar_cf]) ;
        cf_out = [lons_out lats_out years_out tmp1_out] ;
        tmp1_yxv = reshape(nf_out(:,4:end),[Nyears Ncells Nvar_cf]) ;
        tmp1_yxv = cat(1, repmat(tmp1_yxv(1,:,:),[Nyears_xtra 1 1]), tmp1_yxv) ;
        tmp1_out = reshape(tmp1_yxv,[Ncells*(Nyears+Nyears_xtra) Nvar_cf]) ;
        nf_out = [lons_out lats_out years_out tmp1_out] ;
        tmp1_yxv = reshape(ir_out(:,4:end),[Nyears Ncells Nvar_cf]) ;
        tmp1_yxv = cat(1, repmat(tmp1_yxv(1,:,:),[Nyears_xtra 1 1]), tmp1_yxv) ;
        tmp1_out = reshape(tmp1_yxv,[Ncells*(Nyears+Nyears_xtra) Nvar_cf]) ;
        ir_out = [lons_out lats_out years_out tmp1_out] ;
        
    end
            
    % Save land cover
    disp('Saving land cover...')
    file_out = [outDir 'landcover.txt'] ;    
    lpjgu_matlab_saveTable(lu_header_cell, single(lu_out), file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'save_every_pct', save_every_pct, ...
        'verbose', verbose_write) ;
    if do_gzip
        disp('gzipping...')
        unix(['gzip ' file_out]) ;
    end
   
    % Save crop fractions
    disp('Saving crop fractions...')
    file_out = [outDir 'cropfractions.txt'] ;
    lpjgu_matlab_saveTable(cf_header_cell, single(cf_out), file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'save_every_pct', save_every_pct, ...
        'verbose', verbose_write) ;
    if do_gzip
        disp('gzipping...')
        unix(['gzip ' file_out]) ;
    end
    
    
    % Save fertilization
    disp('Saving fertilization...')
    file_out = [outDir 'nfert.txt'] ;
    lpjgu_matlab_saveTable(nf_header_cell, single(nf_out), file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'save_every_pct', save_every_pct, ...
        'verbose', verbose_write) ;
    if do_gzip
        disp('gzipping...')
        unix(['gzip ' file_out]) ;
    end
    
    % Save irrigation
    disp('Saving irrigation...')
    file_out = [outDir 'irrig.txt'] ;
    lpjgu_matlab_saveTable(ir_header_cell, single(ir_out), file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'save_every_pct', save_every_pct, ...
        'verbose', verbose_write) ;
    if do_gzip
        disp('gzipping...')
        unix(['gzip ' file_out]) ;
    end
    
end

disp('All done!')


