%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read MAT-files from harmonization; write as LPJG inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directories for harmonized PLUM outputs
dirList = {...
              'SSP1.v10.s1.harm' ;
              'SSP3.v10.s1.harm' ;
              'SSP4.v10.s1.harm' ;
              'SSP5.v10.s1.harm' ;
              } ;
          
base_year = 2010 ;
y1 = 2011 ;
yN = 2015 ;
yStep = 1 ;
                        
% Trying to avoid new crop spinup time
y1_pre = 2006 ;    % Will repeat first PLUMout year for y1_pre:(y1-1)
someofall = true ; % Make it so that each gridcell always has at least some tiny amount of every crop
              

%% Setup

disp('Setting up...')

cf_kgNha_kgNm2 = 1e-4 ;

outPrec = 6 ;
outWidth = 1 ;
delimiter = ' ' ;
overwrite = false ;
fancy = false ;

addpath('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work/')

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
    if any(isnan(lu_out(:)))
        error('Some value of lu_out is NaN.')
    elseif any(isnan(cf_out(:)))
        error('Some value of cf_out is NaN.')
    elseif any(isnan(nf_out(:)))
        error('Some value of nf_out is NaN.')
    elseif any(isnan(ir_out(:)))
        error('Some value of ir_out is NaN.')
    end
    
    % Check for nonsensical values
    [~,IA,~] = intersect(nf_header_cell,cropList_ir,'stable') ;
    if min(sum(cf_out(:,IA),2))<0
        warning('min(sum_cf)<0')
    elseif max(sum(cf_out(:,IA),2))>=1+10^(-outPrec)
        warning('min(sum_cf)>1')
    elseif min(min(nf_out(:,IA)))<0
        warning('min(nf_out)<0')
    elseif min(min(ir_out(:,IA)))<0
        warning('min(ir_out)<0')
    elseif max(max(ir_out(:,IA)))>=1+10^(-outPrec)
        warning('max(ir_out)>1')
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
    
    % Remove ExtraCropi, because it receives no management inputs.
    cf_out(:,strcmp(cf_header_cell,'ExtraCropi')) = [] ;
    nf_out(:,strcmp(nf_header_cell,'ExtraCropi')) = [] ;
    ir_out(:,strcmp(ir_header_cell,'ExtraCropi')) = [] ;
    cf_header_cell(strcmp(cf_header_cell,'ExtraCropi')) = [] ;
    nf_header_cell(strcmp(nf_header_cell,'ExtraCropi')) = [] ;
    ir_header_cell(strcmp(ir_header_cell,'ExtraCropi')) = [] ;
    cropList_ir(strcmp(cropList_ir,'ExtraCropi')) = [] ;
    
    % Add extra years, if doing so
    if Nyears_xtra>0
        Ncrops = length([cropList_rf cropList_ir]) ;
        tmp = repmat(yearList_xtra,[Ncells 1]) ;
        xtra_out = [repmat(cf_out(1:Ncells,1:2),[Nyears_xtra 1]) ...
                    tmp(:) ...
                    zeros(Nyears_xtra*Ncells,Ncrops)] ;
        cf_out = cat(1, cf_out, xtra_out) ;
        nf_out = cat(1, nf_out, xtra_out) ;
        ir_out = cat(1, ir_out, xtra_out) ;
    end
            
    % Save land cover
    disp('Saving land cover...')
    file_out = [outDir 'landcover.txt'] ;
    lpjgu_matlab_saveTable(lu_header_cell, lu_out, file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', 20, ...
        'verbose', false) ;
    
    % Save crop fractions
    disp('Saving crop fractions...')
    file_out = [outDir 'cropfractions.txt'] ;
    lpjgu_matlab_saveTable(cf_header_cell, cf_out, file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', 20, ...
        'verbose', false) ;
    
    % Save fertilization
    disp('Saving fertilization...')
    file_out = [outDir 'nfert.txt'] ;
    lpjgu_matlab_saveTable(nf_header_cell, nf_out, file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', 20, ...
        'verbose', false) ;
    
    % Save irrigation
    disp('Saving irrigation...')
    file_out = [outDir 'irrig.txt'] ;
    lpjgu_matlab_saveTable(ir_header_cell, ir_out, file_out,...
        'outPrec', outPrec, ...
        'outWidth', outWidth, ...
        'delimiter', delimiter, ...
        'overwrite', overwrite, ...
        'fancy', fancy, ...
        'progress_step_pct', 20, ...
        'verbose', false) ;
        
end

disp('All done!')


