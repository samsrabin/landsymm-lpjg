%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_cropfracs_20170108 = false ;

do_save.LU              = true ;
do_save.crops           = true ;
do_save.crops_plum7     = false ;
do_save.irrig           = true ;
do_save.water           = false ;
do_save.carbon          = false ;
do_save.mrunoff         = false ;
do_save.albedo          = false ;
do_save.bvocs           = false ;
do_save.Nflux           = true ;
do_save.yield           = true ;
do_save.yield_exp       = false ;
do_save.yield_map       = true ;
do_save.yield_exp_map   = false ;
do_save.fpc             = false ;

% inDir_list = {...
% % %     'PLUM2LPJG_26_s1_1850-2010/output-2017-12-02-203502';
% %     'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347';
% % %     'PLUM2LPJG_SSP4_RCP60_mdnPLUM/output-2017-12-07-143853';
% % %     'PLUM2LPJG_26_s1_2011-2100/output-2017-12-05-132537';
% % %     'PLUM2LPJG_26constClimCO2_s1_2011-2100/output-2017-12-05-143605';
% % %     'PLUM2LPJG_45_s1_2011-2100/output-2017-12-05-162344';
% % %     'PLUM2LPJG_60_s1_2011-2100/output-2017-12-05-161832';
% % %     'PLUM2LPJG_85_s1_2011-2100/output-2017-12-06-115808';
% % %     'PLUM2LPJG_SSP1_RCP45_mdnPLUM/output-2017-12-05-093016';
% % %     'PLUM2LPJG_SSP2_RCP45_mdnPLUM/output-2017-12-05-093222',;
% % %     'PLUM2LPJG_SSP3_RCP60_mdnPLUM/output-2017-12-05-113149',;
% % %     'PLUM2LPJG_SSP4_RCP60_mdnPLUM/output-2017-12-05-121202',;
% % %     'PLUM2LPJG_SSP5_RCP85_mdnPLUM/output-2017-12-05-122918';
% % %     'PLUM2LPJGblIrr_26_s1_2011-2100/output-2017-12-07-093537';
% % %     'PLUM2LPJGblIrr_26constClimCO2_s1_2011-2100/output-2017-12-07-090538';
% % %     'PLUM2LPJGblIrr_45_s1_2011-2100_constLU/output-2017-12-08-153814';
% % %     'PLUM2LPJGblIrr_45_s1_2011-2100/output-2017-12-07-094546';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constLU/output-2017-12-08-151957';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100/output-2017-12-07-093346';
% % %     'PLUM2LPJGblIrr_85_s1_2011-2100_constLU/output-2017-12-08-150159';
% % %     'PLUM2LPJGblIrr_SSP1_RCP45_mdnPLUM/output-2017-12-07-095646';
% % %     'PLUM2LPJGblIrr_SSP2_RCP45_mdnPLUM/output-2017-12-07-094823';
% % %     'PLUM2LPJGblIrr_SSP3_RCP60_mdnPLUM/output-2017-12-07-084408';
% % %     'PLUM2LPJGblIrr_SSP4_RCP60_mdnPLUM/output-2017-12-07-144106';
% % %     'PLUM2LPJGblIrr_SSP5_RCP85_mdnPLUM/output-2017-12-07-090531';
% % %     'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU_constCO2/output-2017-12-11-115548';
% % %     'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU/output-2017-12-10-232848';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP2LU_constCO2/output-2017-12-11-180702';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP2LU/output-2017-12-11-011415';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU_constCO2/output-2017-12-11-180709';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU/output-2017-12-11-020406';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU_constCO2/output-2017-12-11-203339';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU/output-2017-12-11-024202';
% % %     'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU_constCO2/output-2017-12-11-205323';
% % %     'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU/output-2017-12-11-102405';
% % %     'PLUM2LPJGblIrr_26constClimCO2_s1_2011-2100_constLU/output-2017-12-09-001746'...
% % %     'PLUM2LPJGblIrr_SSP2_RCP60_mdnPLUM/output-2017-12-09-141738'...
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP2LU/output-2017-12-15-015508';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU/output-2017-12-14-1250c23';
% % %     'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU/output-2017-12-15-015520';
% % %     'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU/output-2017-12-15-011713';
% %     'PLUM2LPJGblIrr_SSP1v2_RCP45/output-2017-12-25-140114';
% %     'PLUM2LPJGblIrr_SSP2v2_RCP60/output-2017-12-25-134337';
% %     'PLUM2LPJGblIrr_SSP3v2_RCP60/output-2017-12-25-133608';
% %     'PLUM2LPJGblIrr_SSP4v2_RCP60/output-2017-12-25-134011';
% %     'PLUM2LPJGblIrr_SSP5v2_RCP85/output-2017-12-25-135831';
% %     'PLUM2LPJGblIrrNoFGnoRye_26_s1_1850-2010_moreProcs/matlab_merge_20180108120140' ;
% %     'PLUM2LPJGblIrrNoFG_26_s1_1850-2010_moreProcs/matlab_merge_20180108114616' ;
% %     'PLUM2LPJGblIrr_SSP1v2_RCP45.20180109/output-2018-01-11-081646' ;
% %     'PLUM2LPJGblIrr_SSP3v2_RCP60.20180109/output-2018-01-11-023306' ;
% %     'PLUM2LPJGblIrr_SSP4v2_RCP60.20180109/output-2018-01-11-023748' ;
% %     'PLUM2LPJGblIrr_SSP5v2_RCP85.20180109/output-2018-01-11-021841' ;
% %     'PLUM2LPJGblIrr_SSP1v2_RCP45.20180109/output-2018-01-13-230513' ;
%     'PLUM2LPJGblIrr_SSP3v2_RCP60.20180109/output-2018-01-13-231353' ;
% %     'PLUM2LPJGblIrr_SSP4v2_RCP60.20180109/output-2018-01-13-231930' ;
% %     'PLUM2LPJGblIrr_SSP5v2_RCP85.20180109/output-2018-01-13-225751' ;
%     } ;
% Ncrops = 4 ;

% inDir_list = {...
% %     'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180124084436' ;
%     'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180422140716' ;
%     'PLUM2LPJG.PLUM7.SSP1v2.RCP45/output-2018-02-02-054242' ;
%     'PLUM2LPJG.PLUM7.SSP3v2.RCP60/output-2018-01-25-090934' ;
%     'PLUM2LPJG.PLUM7.SSP4v2.RCP60/output-2018-01-25-090916' ;
%     'PLUM2LPJG.PLUM7.SSP5v2.RCP85/output-2018-01-24-182830' ;
% } ;
% Ncrops = 8 ;

inDir_list = {...
% %     'PLUM2LPJGblIrr_SSP1_RCP45_mdnPLUM/output-2017-12-07-095646' ;
%     'PLUM2LPJGblIrr_SSP3_RCP60_mdnPLUM/output-2017-12-07-084408' ;
%     'PLUM2LPJGblIrr_SSP4_RCP60_mdnPLUM/output-2017-12-07-144106' ;
%     'PLUM2LPJGblIrr_SSP5_RCP85_mdnPLUM/output-2017-12-07-090531' ;
%     'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347'
    'PLUM2LPJGblIrr_45_s1_2011-2100_constSSP1LU/output-2017-12-10-232848' ;
    'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP3LU/output-2017-12-14-125023' ;
    'PLUM2LPJGblIrr_60_s1_2011-2100_constSSP4LU/output-2017-12-15-015520' ;
    'PLUM2LPJGblIrr_85_s1_2011-2100_constSSP5LU/output-2017-12-15-011713' ;

} ;
Ncrops = 4 ;

% inDir_list = {...
% %     'PLUM2LPJG_SSP1_RCP45_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-23-101319' ;
% %     'PLUM2LPJG_SSP3_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-23-100011' ;
% %     'PLUM2LPJG_SSP4_RCP60_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-23-093751' ;
% %     'PLUM2LPJG_SSP5_RCP85_v3s1_constLUmgmt_asPLUMout2011-2015/output-2018-04-23-092619' ;
%     } ;
% Ncrops = 7 ;


%% Setup

addpath(genpath(landsymm_lpjg_path()))

is_baseline_list = false(length(inDir_list),1) ;
for d = 1:length(inDir_list)
    if strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_26_s1_1850-2010/output-2017-12-02-203502') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJGblIrrNoFGnoRye_26_s1_1850-2010_moreProcs/matlab_merge_20180108120140') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJGblIrrNoFG_26_s1_1850-2010_moreProcs/matlab_merge_20180108114616') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180124084436') ...
    || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180422140716')
        is_baseline_list(d) = true ;
    end
end

pftList_noCrops = {'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE',...
            'TrBR','C3G','C4G','PC3G','PC4G',...
            'CC3G_ic','CC4G_ic','Total','Crop_sum','Pasture_sum',...
            'Natural_sum','Barren_sum'} ;
LUlist = {'CROPLAND','PASTURE','NATURAL','BARREN'} ;

% Anonymous function for checking whether any version of a file exists
anyfileexist = @(in_file) ...
    exist(in_file,'file') ...
    || exist([in_file '.maps.mat'],'file') ...
    || exist([in_file '.mat'],'file') ...
    || exist([in_file '.gz'],'file') ;

% Get CFTs and PFTs
if Ncrops==4
    cropTypes = {'TeWW','TeSW','TeCo','TrRi'} ;
    PLUMcrops_2_LPJG = readtable([thisDir 'PLUMcrops_2_LPJG_20171229.csv']) ;
    cropTypes_conv_plum7 = PLUMcrops_2_LPJG.PLUMname ;
    cropTypes_conv_lpjg = PLUMcrops_2_LPJG.LPJGname ;
    calibFactors = PLUMcrops_2_LPJG.calibFactor ;
    clear PLUMcrops_2_LPJG
elseif Ncrops==8
    cropTypes = {'CerealsC3w','CerealsC3s','Oilcrops',...
                 'StarchyRoots','Pulses','CerealsC4',...
                 'Miscanthus','Rice'} ;
    LPJGcrops_2_PLUM = readtable([thisDir 'LPJGcrops_2_PLUM_20180124.csv']) ;
    cropTypes_conv_plum7 = LPJGcrops_2_PLUM.PLUMname ;
    cropTypes_conv_lpjg = LPJGcrops_2_PLUM.LPJGname ;
    calibFactors = LPJGcrops_2_PLUM.calibFactor ;
%     clear LPJGcrops_2_PLUM
elseif Ncrops==7
    cropTypes = {'CerealsC3','Oilcrops',...
                 'StarchyRoots','Pulses','CerealsC4',...
                 'Miscanthus','Rice'} ;
    LPJGcrops_2_PLUM = readtable([thisDir 'PLUM6xtra_calib_20180423.csv']) ;
    cropTypes_conv_plum7 = LPJGcrops_2_PLUM.PLUMname ;
    cropTypes_conv_lpjg = LPJGcrops_2_PLUM.LPJGname ;
    calibFactors = LPJGcrops_2_PLUM.calibFactor ;
%     clear LPJGcrops_2_PLUM
else
    error(['Rework CFT/PFT code for Ncrops==' num2str(Ncrops)]) ;
end
pftList = [pftList_noCrops cropTypes] ;
plum7_hasMoreCrops = length(unique(cropTypes_conv_plum7)) > length(cropTypes) ;
if ~plum7_hasMoreCrops
    cropTypes_plum7 = unique(cropTypes_conv_plum7) ;
elseif length(unique(cropTypes_conv_plum7))==length(cropTypes_conv_plum7)
    cropTypes_plum7 = cropTypes_conv_plum7 ;
else
    error('PLUM7 has more crops, but not all its values are unique??')
end
Ncrops_plum7 = length(cropTypes_plum7) ;

% Import bare-soil albedo
if do_save.albedo
    baresoil_albedo_file = [thisDir 'soilmap.txt'] ;
    baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true,'verbose',false) ;
    baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
    clear baresoil_albedo
end

% Read PLUMout_gridlist
if any(is_baseline_list)
    gridlist_file = '/Users/Shared/lpj-guess/gridlists/PLUMout_gridlist.txt' ;
    PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;
end

% Import land area (km2)
landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX = tmp(1:2:720,:) + tmp(2:2:720,:) ;
if any(is_baseline_list)
    land_area_YX(~PLUMout_gridlist.mask_YX) = NaN ;
end
land_area_YX_m2 = land_area_YX*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Conversion factors
cf_cpool = 1e-3*1e-9 ;
cf_water = 1e-3*1e-3*1e-3 ;
cf_nflux = 1e-9 ;
cf_bvoc = 1e-15 ;
cf_lu = 1e-6 ;   % km2 to Mkm2



%% Loop through inDir_list

for d = 1:length(inDir_list)
    
    %%%%%%%%%%%%%
    %%% Setup %%%
    %%%%%%%%%%%%%
    
    % Get full inDir
    inDir = find_PLUM2LPJG_run(inDir_list{d}) ;
    is_baseline = is_baseline_list(d) ;
    if is_baseline
        disp('is_baseline')
    end
    
    timeseries_out = [inDir 'timeseries.mat'] ;
    timeseries_PLUMexp_out = [inDir 'timeseries_PLUMexp.mat'] ;
    lastdecade_out = [inDir 'last_decade.mat'] ;
    
    disp(inDir)
    
    % Get baseline/not info
    if is_baseline
        yearList = 1850:2010 ;
        replace_CrOp_with_CrOpi = false ;
        merge_CrOpi_into_CrOp = true ;
        if ~exist('PLUMout_mask_YX1y','var')
            PLUMout_mask_YX1y = repmat(PLUMout_gridlist.mask_YX,[1 1 1 length(yearList)]) ;
        end
    else
        yearList = 2011:2100 ;
        replace_CrOp_with_CrOpi = true ;
        merge_CrOpi_into_CrOp = false ;
    end
    yearList = transpose(yearList) ;
    Nyears = length(yearList) ;
        
    % Get land use and crop fractions input files
    clear LUfile cropfile
    if strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_26_s1_1850-2010/output-2017-11-24-175840/') ...
            || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_26_s1_1850-2010/output-2017-11-29-163216/') ...
            || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_26_s1_1850-2010/output-2017-12-02-203502/')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
        cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_1901_2014_20160906.MP.noIrr.plus_Nfert0-200-1000xIrr_factorial_empties.inclZeroFert.out' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180124084436') ...
            || strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJG_PLUM7_1850-2010/matlab_merge_20180422140716')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
        cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM7_1850-2010_20180105a.assignWWorSW.out' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJGblIrrNoFG_26_s1_1850-2010_moreProcs/matlab_merge_20180108114616')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
        cropfile = '/Users/Shared/remapping_MIRCA/crop_fractions_fromMIRCA_PLUM_1yr_20180105a.MP.out' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJGblIrrNoFGnoRye_26_s1_1850-2010_moreProcs/matlab_merge_20180108120140')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
        cropfile = '/Users/Shared/remapping_MIRCA/crop_fractions_fromMIRCA_PLUM_1yr_20180105b.MP.out' ;
    elseif strcmp_ignoreTrailSlash(inDir_list{d},'PLUM2LPJGblIrr_26_s1_1850-2010/output-2017-12-06-001347/')
        LUfile = '/project/fh1-project-lpjgpi/lr8247/input/PLUM_input/LU/lu_1850_2015_luh2_aggregate_sum2x2_midpoint_nourban_orig_v21.txt' ;
        if test_cropfracs_20170108
            cropfile = '/Users/Shared/remapping_MIRCA/crop_fractions_fromMIRCA_PLUM_1yr_20180105.MP.out' ;
        else
            cropfile = '/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_1901_2014_20160906.MP.out' ;
        end
    else
%         [x,LUcropDir_tmp] = unix([thisDir 'get_lu_dir.sh '...
%             ...'/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper2/MATLAB_work']) ;
%             inDir]) ;
%         LUcropDir_tmp = LUcropDir_tmp(1:end-1) ;
%         LUcropDir = find_PLUM2LPJG_inputs(LUcropDir_tmp) ;
        [x,LUfile] = unix([thisDir 'get_lu_file.sh ' inDir]) ;
        if x~=0
            error(['get_lu_file.sh failed with error ' num2str(x)])
        end
        LUfile_tmp = regexprep(LUfile,'[\n\r]+','') ; % Remove extraneous newline
        LUfile = find_PLUM2LPJG_input_file(LUfile_tmp) ;
        [x,cropfile] = unix([thisDir 'get_cropfrac_file.sh ' inDir]) ;
        if x~=0
            error(['get_cropfrac_file.sh failed with error ' num2str(x)])
        end
        cropfile_tmp = regexprep(cropfile,'[\n\r]+','') ; % Remove extraneous newline
        cropfile = find_PLUM2LPJG_input_file(cropfile_tmp) ;
    end
    if ~exist('LUfile','var')
        LUfile = [LUcropDir 'landcover.txt'] ;
        cropfile = [LUcropDir 'cropfractions.txt'] ;
    end
    
    % Are we also going to fudge into plum7?
    fudge2plum7 = false ;
    if Ncrops==8
        fudge2plum7 = true ;
        cropfile_plum7 = '' ;
    elseif Ncrops~=4
        error(['Am I supposed to fudge2plum7 with Ncrops = ' num2str(Ncrops) '?'])
    elseif exist([inDir 'cropfile_plum7.txt'],'file')
        fudge2plum7 = true ;
        cropfile_plum7 = [inDir 'cropfile_plum7.txt'] ;
    elseif exist('LUcropDir_tmp','var') && any(strfind(LUcropDir_tmp,'20171201'))
        LUcropDir_tmp = strrep(LUcropDir_tmp,'20171201','20171228') ;
        LUcropDir_plum7 = find_PLUM2LPJG_inputs(LUcropDir_tmp,false) ;
        if ~isempty(LUcropDir_plum7)
            fudge2plum7 = true ;
            cropfile_plum7 = [LUcropDir_plum7 'cropfractions.txt'] ;
        end
    elseif exist('LUcropDir_tmp','var') ...
        && (any(strfind(LUcropDir_tmp,'20180109')) || any(strfind(LUcropDir_tmp,'AsPLUM_20180118')))
        LUcropDir_tmp = strrep(strrep(LUcropDir_tmp,'20180109','20171228'),'AsPLUM_20180118','20171228') ;
        LUcropDir_plum7 = find_PLUM2LPJG_inputs(LUcropDir_tmp,false) ;
        if ~isempty(LUcropDir_plum7)
            fudge2plum7 = true ;
            cropfile_plum7 = [LUcropDir_plum7 'cropfractions.txt'] ;
        end
    end
    if exist('LUcropDir_tmp','var') && exist('cropfile_plum7','var') && isempty(cropfile_plum7) && any(strfind(LUcropDir_tmp,'AsPLUM_20180118'))
%         LUcropDir_tmp = strrep(LUcropDir_tmp,'AsPLUM_20180118','20171228') ;
        LUcropDir_plum7 = find_PLUM2LPJG_inputs(LUcropDir_tmp,true) ;
        cropfile_plum7 = [LUcropDir_plum7 'cropfractions.txt'] ;
        expyieldfile_plum7 = [LUcropDir_plum7 'yield.txt'] ;
    end
    
    
    % Do expected yields exist? If not, make sure do_save.yield_exp* are
    % false. Reset back to true at end of this for loop (before going to
    % next directory)
    flip_do_save_yield_exp = false ;
    if do_save.yield_exp && ~exist('LUcropDir_plum7','var')
        warning('No expected yields exist; skipping do_save.yield_exp.')
        flip_do_save_yield_exp = true ;
        do_save.yield_exp = ~do_save.yield_exp ;
    end
    flip_do_save_yield_exp_map = false ;
    if do_save.yield_exp_map && ~exist('LUcropDir_plum7','var')
        warning('No expected yields exist; skipping do_save.yield_exp_map.')
        flip_do_save_yield_exp_map = true ;
        do_save.yield_exp_map = ~do_save.yield_exp_map ;
    end
    
    % Setup if processing expected yields
    if do_save.yield_exp
        if ~exist('expyieldfile_plum7','var')
            expyieldfile_plum7 = [LUcropDir_plum7 'yield.txt'] ;
        end
        if ~exist(expyieldfile_plum7,'file')
            warning('Expected yield file not found! Setting do_save.yield_exp to false.')
            do_save.yield_exp = false ;
        end
    end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get important info %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yield = lpjgu_matlab_readTable_then2map([inDir 'yield.out'],'force_mat_save',true) ;
    list_to_map = yield.list_to_map ;
    if min(yield.yearList) == min(yearList) - 5 && max(yield.yearList)==max(yearList)
        warning('Adjusting yearList to account for 5 years'' padding at beginning.')
        yearList = yield.yearList ;
        Nyears = length(yearList) ;
    elseif ~isequal(yield.yearList,yearList)
        error('Rework so yearLists match!')
    end
    if ~(do_save.yield || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map)
        clear yield
    end
    
    
    %%%%%%%%%%%%%%%%%%
    %%% Import FPC %%%
    %%%%%%%%%%%%%%%%%%
    
    if do_save.fpc || do_save.albedo
        fpc = lpjgu_matlab_readTable_then2map([inDir 'fpc.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('Processing FPC...')
        [fpc, ~, ~] = CrOp_and_CrOpi(fpc, 'fpc', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp) ;
        % Normalize to 1
        fpc_total_YX1y = fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Total'),:) ;
        isbad_YX1y = fpc_total_YX1y > 1 ;
        i = 0 ;
        while any(isbad_YX1y(:))
            i = i + 1 ;
            if i > 5
                error('Too many iterations!')
            end
            isbad_YXvy = repmat(isbad_YX1y,[1 1 size(fpc.maps_YXvy,3) 1]) ;
            fpc_total_YXvy = repmat(fpc_total_YX1y,[1 1 size(fpc.maps_YXvy,3) 1]) ;
            fpc.maps_YXvy(isbad_YXvy) = fpc.maps_YXvy(isbad_YXvy) ./ fpc_total_YXvy(isbad_YXvy) ;
            clear fpc_total_YXvy
            fpc_total_YX1y = fpc.maps_YXvy(:,:,strcmp(fpc.varNames,'Total'),:) ;
            isbad_YX1y = fpc_total_YX1y > 1 ;
        end
        clear fpc_total_YX1y isbad* i
        if is_baseline
            list_to_map = PLUMout_gridlist.list_to_map ;
            fpc.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(fpc.varNames) 1])) = NaN ;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Import land use %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    LU = lpjgu_matlab_readTable_then2map(LUfile,'force_mat_save',true) ;
    if ~isfield(LU,'yearList') && size(LU.maps_YXv,4)==1
        % One-year LU file
        LU.maps_YXvy = repmat(LU.maps_YXv,[1 1 1 Nyears]) ;
        LU.yearList = yearList ;
    end
    if ~isequal(LU.yearList,yearList)
        if min(LU.yearList) <= min(yearList) && max(LU.yearList) >= max(yearList)
            LU.maps_YXvy = LU.maps_YXvy(:,:,:,LU.yearList>=min(yearList) & LU.yearList<=max(yearList)) ;
            LU.yearList = transpose(LU.yearList(1):max(yearList)) ;
        end
        if min(LU.yearList) == min(yearList) + 5 && max(LU.yearList)==max(yearList)
            warning('Adjusting LU to account for 5 years'' padding at beginning.')
            LU.yearList = yearList ;
            LU.maps_YXvy = cat(4,repmat(LU.maps_YXvy(:,:,:,1),[1 1 1 5]),LU.maps_YXvy) ;
        elseif ~isequal(LU.yearList,yearList)
            error('Rework LU so yearLists match!')
        end
    end
    % Check, removing extraneous if necessary
    if any(strcmp(LU.varNames,'URBAN'))
        thisMax = max(max(max(LU.maps_YXvy(:,:,strcmp(LU.varNames,'URBAN'),:)))) ;
        if thisMax > 0
            warning(['Removing LU: URBAN even though up to ' num2str(thisMax) ' is present.'])
        end
        LU.maps_YXvy(:,:,strcmp(LU.varNames,'URBAN'),:) = [] ;
        LU.varNames(strcmp(LU.varNames,'URBAN')) = [] ;
    end
    if any(strcmp(LU.varNames,'PEATLAND'))
        thisMax = max(max(max(LU.maps_YXvy(:,:,strcmp(LU.varNames,'PEATLAND'),:)))) ;
        if thisMax > 0
            warning(['Removing LU: PEATLAND even though up to ' num2str(thisMax) ' is present.'])
        end
        LU.maps_YXvy(:,:,strcmp(LU.varNames,'PEATLAND'),:) = [] ;
        LU.varNames(strcmp(LU.varNames,'PEATLAND')) = [] ;
    end
    if any(~(strcmp(LU.varNames,'PASTURE') | strcmp(LU.varNames,'CROPLAND') | strcmp(LU.varNames,'NATURAL') | strcmp(LU.varNames,'BARREN')))
        error('Something is wrong with LU.varList.')
    end
    if is_baseline
        LU.list_to_map = PLUMout_gridlist.list_to_map ;
        LU.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(LU.varNames) 1])) = NaN ;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import crop fractions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(cropfile,'/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_1901_2014_20160906.MP.noIrr.plus_Nfert0-200-1000xIrr_factorial_empties.inclZeroFert.out')
        load('/project/fh1-project-lpjgpi/lr8247/PLUM/input/LU/crop_fractions_fromMIRCA_PLUM_1901_2014_20160906.MP.noIrr.mat') ;
    else
        cropfracs = lpjgu_matlab_readTable_then2map(cropfile,'force_mat_save',true) ;
    end
    if ~isfield(cropfracs,'yearList') && size(cropfracs.maps_YXv,4)==1
        % One-year cropfracs file
        cropfracs.maps_YXvy = repmat(cropfracs.maps_YXv,[1 1 1 Nyears]) ;
        cropfracs.yearList = yearList ;
    end
    if ~isequal(cropfracs.yearList,yearList)
        if min(cropfracs.yearList) <= min(yearList) && max(cropfracs.yearList) >= max(yearList)
            cropfracs.maps_YXvy = cropfracs.maps_YXvy(:,:,:,cropfracs.yearList>=min(yearList) & cropfracs.yearList>= max(yearList)) ;
            cropfracs.yearList = transpose(cropfracs.yearList(1):max(yearList)) ;
        elseif min(cropfracs.yearList) > min(yearList) && max(cropfracs.yearList) == max(yearList) && isequaln(cropfracs.maps_YXvy(:,:,:,1),cropfracs.maps_YXvy(:,:,:,end))
            Nmissing = length(yearList(1):cropfracs.yearList(1)-1) ;
            cropfracs.maps_YXvy = cat(4,repmat(cropfracs.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),cropfracs.maps_YXvy(:,:,:,cropfracs.yearList>=min(yearList) & cropfracs.yearList>= max(yearList))) ;
            cropfracs.yearList = yearList ;
        elseif min(cropfracs.yearList) > min(yearList) && max(cropfracs.yearList) > max(yearList) && isequaln(cropfracs.maps_YXvy(:,:,:,1),cropfracs.maps_YXvy(:,:,:,end))
            Nmissing = length(yearList(1):cropfracs.yearList(1)-1) ;
            cropfracs.maps_YXvy = cat(4,...
                repmat(cropfracs.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),...
                cropfracs.maps_YXvy(:,:,:,...
                cropfracs.yearList>=min(yearList) ...
                & cropfracs.yearList<=max(yearList)...
                )) ;
            %         Nextra = max(cropfracs.yearList) - max(yearList) ;
            %         cropfracs.maps_YXvy = cropfracs.maps_YXvy(:,:,:,1:find(cropfracs.yearList==max(yearList))) ;
            cropfracs.yearList = yearList ;
        end
        if min(cropfracs.yearList) == min(yearList) + 5 && max(cropfracs.yearList)==max(yearList)
            warning('Adjusting cropfracs to account for 5 years'' padding at beginning.')
            cropfracs.yearList = yearList ;
            cropfracs.maps_YXvy = cat(4,repmat(cropfracs.maps_YXvy(:,:,:,1),[1 1 1 5]),cropfracs.maps_YXvy) ;
        elseif ~isequal(cropfracs.yearList,yearList)
            error('Rework cropfracs so yearLists match!')
        end
    end
    cropfracs_orig = cropfracs ;
    [cropfracs, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig] = ...
        CrOp_and_CrOpi(cropfracs, 'cropfracs', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, cropfile) ;
    if is_baseline
        cropfracs.list_to_map = PLUMout_gridlist.list_to_map ;
        cropfracs.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs.varNames) 1])) = NaN ;
        if exist('cropfracs_orig','var')
            cropfracs_orig.list_to_map = PLUMout_gridlist.list_to_map ;
            cropfracs_orig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs_orig.varNames) 1])) = NaN ;
        end
    end
    % Get area of each crop (km2)
    cropareas = cropfracs ;
%     cropareas.maps_YXvy = cropfracs.maps_YXvy .* repmat(LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:),[1 1 length(cropfracs.varNames) 1]) ;
    cropareas.maps_YXvy = ...
        cropfracs.maps_YXvy ...
        .* repmat(LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:),[1 1 length(cropfracs.varNames) 1]) ...
        .* repmat(land_area_YX,[1 1 length(cropfracs.varNames) Nyears]) ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import crop fractions (plum7) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if fudge2plum7
        if ~isempty(cropfile_plum7)
            % Import crop fractions
            cropfracs_plum7 = lpjgu_matlab_readTable_then2map(cropfile_plum7,'force_mat_save',true) ;
            if ~isfield(cropfracs_plum7,'yearList') && size(cropfracs_plum7.maps_YXv,4)==1
                % One-year cropfracs_plum7 file
                cropfracs_plum7.maps_YXvy = repmat(cropfracs_plum7.maps_YXv,[1 1 1 Nyears]) ;
                cropfracs_plum7.yearList = yearList ;
            end
            if ~isequal(cropfracs_plum7.yearList,yearList)
                if min(cropfracs_plum7.yearList) <= min(yearList) && max(cropfracs_plum7.yearList) >= max(yearList)
                    cropfracs_plum7.maps_YXvy = cropfracs_plum7.maps_YXvy(:,:,:,cropfracs_plum7.yearList>=min(yearList) & cropfracs_plum7.yearList>= max(yearList)) ;
                    cropfracs_plum7.yearList = transpose(cropfracs_plum7.yearList(1):max(yearList)) ;
                elseif min(cropfracs_plum7.yearList) > min(yearList) && max(cropfracs_plum7.yearList) == max(yearList) && isequaln(cropfracs_plum7.maps_YXvy(:,:,:,1),cropfracs_plum7.maps_YXvy(:,:,:,end))
                    Nmissing = length(yearList(1):cropfracs_plum7.yearList(1)-1) ;
                    cropfracs_plum7.maps_YXvy = cat(4,repmat(cropfracs_plum7.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),cropfracs_plum7.maps_YXvy(:,:,:,cropfracs_plum7.yearList>=min(yearList) & cropfracs_plum7.yearList>= max(yearList))) ;
                    cropfracs_plum7.yearList = yearList ;
                elseif min(cropfracs_plum7.yearList) > min(yearList) && max(cropfracs_plum7.yearList) > max(yearList) && isequaln(cropfracs_plum7.maps_YXvy(:,:,:,1),cropfracs_plum7.maps_YXvy(:,:,:,end))
                    Nmissing = length(yearList(1):cropfracs_plum7.yearList(1)-1) ;
                    cropfracs_plum7.maps_YXvy = cat(4,...
                        repmat(cropfracs_plum7.maps_YXvy(:,:,:,1),[1 1 1 Nmissing]),...
                        cropfracs_plum7.maps_YXvy(:,:,:,...
                        cropfracs_plum7.yearList>=min(yearList) ...
                        & cropfracs_plum7.yearList<=max(yearList)...
                        )) ;
                    cropfracs_plum7.yearList = yearList ;
                end
                if min(cropfracs_plum7.yearList) == min(yearList) + 5 && max(cropfracs_plum7.yearList)==max(yearList)
                    warning('Adjusting cropfracs_plum7 to account for 5 years'' padding at beginning.')
                    cropfracs_plum7.yearList = yearList ;
                    cropfracs_plum7.maps_YXvy = cat(4,repmat(cropfracs_plum7.maps_YXvy(:,:,:,1),[1 1 1 5]),cropfracs_plum7.maps_YXvy) ;
                elseif ~isequal(cropfracs_plum7.yearList,yearList)
                    error('Rework cropfracs_plum7 so yearLists match!')
                end
            end
            cropfracs_plum7_orig = cropfracs_plum7 ;
            [cropfracs_plum7, inds_cropTypes_plum7_cropFracsOrig, ...
                inds_cropTypesI_plum7_cropFracsOrig] = ...
                CrOp_and_CrOpi(cropfracs_plum7, 'cropfracs_plum7', cropTypes_plum7, ...
                replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, cropfile_plum7) ;
            if is_baseline
                cropfracs_plum7.list_to_map = PLUMout_gridlist.list_to_map ;
                cropfracs_plum7.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs_plum7.varNames) 1])) = NaN ;
                if exist('cropfracs_plum7_orig','var')
                    cropfracs_plum7_orig.list_to_map = PLUMout_gridlist.list_to_map ;
                    cropfracs_plum7_orig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cropfracs_plum7_orig.varNames) 1])) = NaN ;
                end
            end
        else
            cropfracs_plum7 = plum7ify_struct(...
                cropfracs, cropTypes_conv_plum7, cropTypes_conv_lpjg, ...
                cropfracs.varNames, cropTypes_plum7, 'sum', []) ;
        end
        % Get area of each crop (km2)
        cropareas_plum7 = cropfracs_plum7 ;
    %     cropareas_plum7.maps_YXvy = cropfracs_plum7.maps_YXvy .* repmat(LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:),[1 1 length(cropfracs_plum7.varNames) 1]) ;
        cropareas_plum7.maps_YXvy = ...
            cropfracs_plum7.maps_YXvy ...
            .* repmat(LU.maps_YXvy(:,:,strcmp(LU.varNames,'CROPLAND'),:),[1 1 length(cropfracs_plum7.varNames) 1]) ...
            .* repmat(land_area_YX,[1 1 length(cropfracs_plum7.varNames) Nyears]) ;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save yield (kgDM/m2 --> kgDM) (i.e., actually PRODUCTION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_save.yield || do_save.yield_exp || do_save.yield_map || do_save.yield_exp_map
        [yield, ~, ~] = CrOp_and_CrOpi(yield, 'yield', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, ...
            '', cropfracs_orig, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig) ;
        yield.maps_YXvy(:,:,strcmp(yield.varNames,'CC3G_ic')|strcmp(yield.varNames,'CC4G_ic'),:) = [] ;
        tmp = yield.varNames ;
        tmp(strcmp(tmp,'CC3G_ic')) = [] ;
        tmp(strcmp(tmp,'CC4G_ic')) = [] ;
        yield.varNames = tmp ;
        
        if do_save.yield || do_save.yield_map
            disp('   Saving yield/cropprod...')
            if is_baseline
                yield.list_to_map = PLUMout_gridlist.list_to_map ;
                yield.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(yield.varNames) 1])) = NaN ;
            end
            if do_save.yield
%                 disp('Regular:')
                for c = 1:length(cropTypes)
                    thisCrop = cropTypes{c} ;
% % %                     eval(['cropprod_ts_' thisCrop ' = getTS(yield,thisCrop,1e6*cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
                    eval(['cropprod_ts_' thisCrop '_lpjg = getTS(yield,thisCrop,1e6*cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) ;']) ;
%                     eval(['disp([''   '' thisCrop ''_LPJG = '' num2str(mean(cropprod_ts_' thisCrop '_lpjg))]) ;'])
                    save(timeseries_out,['cropprod_ts_' thisCrop '_lpjg'],v73_or_append(timeseries_out)) ;
                    clear thisCrop
                end; clear c
                clear cropprod_ts*
            end
            
        end

        if fudge2plum7
            disp('   Saving yield/cropprod (PLUM7)...')
            if do_save.yield || do_save.yield_map
                % First, get new structure with new varNames and maps_YXvy
                if is_baseline
                    [yield_plum7, outFrom] = plum7ify_struct(...
                        yield, cropTypes_conv_plum7, cropTypes_conv_lpjg, ...
                        cropareas.varNames, cropareas_plum7.varNames,'max',[]) ;
                else
                    [yield_plum7, outFrom] = plum7ify_struct(...
                        yield, cropTypes_conv_plum7, cropTypes_conv_lpjg, ...
                        cropareas.varNames, cropareas_plum7.varNames,'wtd_av',[]) ;
                end
                % Then calculate and save
                if do_save.yield
%                     disp('plum7:')
                    for c = 1:Ncrops_plum7
                        thisCrop_plum7 = cropTypes_plum7{c} ;
% % %                         thisCalibFactor = calibFactors(c) ;
                        thisCalibFactor = unique(calibFactors(strcmp(cropTypes_conv_plum7,thisCrop_plum7))) ;
                        eval(['cropprod_ts_' thisCrop_plum7 '_plum = thisCalibFactor*getTS(yield_plum7,thisCrop_plum7,1e6*cropareas_plum7.maps_YXvy(:,:,strcmp(cropareas_plum7.varNames,thisCrop_plum7),:)) ;']) ;
%                         eval(['disp([''   '' thisCrop_plum7 ''_plum = '' num2str(mean(cropprod_ts_' thisCrop_plum7 '_plum))]) ;'])
%                         eval(['disp([''   '' thisCrop_plum7 ''_plum = '' num2str((1/thisCalibFactor)*mean(cropprod_ts_' thisCrop_plum7 '_plum))]) ;'])
                        save(timeseries_out,['cropprod_ts_' thisCrop_plum7 '_plum'],v73_or_append(timeseries_out)) ;
                        clear thisCrop_plum7 thisCalibFactor
                    end; clear c
                    clear cropprod_ts*
                end
                if do_save.yield_map
                    yield_plum7.maps_YXvy = yield_plum7.maps_YXvy(:,:,:,end-9:end) ;
                    yield_plum7.yearList = yield_plum7.yearList(end-9:end) ;
                    save(lastdecade_out,'yield_plum7',v73_or_append(lastdecade_out)) ;
                end
                clear yield_plum7
            end
            
            % EXPECTED yield
            if do_save.yield_exp || do_save.yield_exp_map
                disp('   Saving yield/cropprod (PLUM7 EXPECTED)...')
                expyield_plum7 = lpjgu_matlab_readTable_then2map(expyieldfile_plum7,'force_mat_save',true) ;
                if min(expyield_plum7.yearList) == min(yearList) + 5 && max(expyield_plum7.yearList)==max(yearList)
                    warning('Adjusting expyield_plum7 to account for 5 years'' padding at beginning.')
                    expyield_plum7.yearList = yearList ;
                    expyield_plum7.maps_YXvy = cat(4,repmat(expyield_plum7.maps_YXvy(:,:,:,1),[1 1 1 5]),expyield_plum7.maps_YXvy) ;
                elseif ~isequal(expyield_plum7.yearList,yearList)
                    error('Rework expyield_plum7 so yearLists match!')
                end
                [expyield_plum7, inds_cropTypes_plum7_cropFracsOrig, inds_cropTypesI_plum7_cropFracsOrig] = ...
                    CrOp_and_CrOpi(expyield_plum7, 'yield', cropTypes_plum7, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, ...
                    cropfile_plum7, cropfracs_plum7_orig, inds_cropTypes_plum7_cropFracsOrig, inds_cropTypesI_plum7_cropFracsOrig) ;
                
                % Then calculate and save (CALIBRATION FACTORS ALREADY INCLUDED)
                if do_save.yield_exp
                    for c = 1:Ncrops_plum7
                        thisCrop_plum7 = cropTypes_plum7{c} ;
                        eval(['cropprod_ts_' thisCrop_plum7 '_plum = getTS(expyield_plum7,thisCrop_plum7,1e6*cropareas_plum7.maps_YXvy(:,:,strcmp(cropareas_plum7.varNames,thisCrop_plum7),:)) ;']) ;
                        save(timeseries_PLUMexp_out,['cropprod_ts_' thisCrop_plum7 '_plum'],v73_or_append(timeseries_PLUMexp_out)) ;
                        clear thisCrop_plum7
                    end; clear c
                    clear cropprod_ts*
                end
                if do_save.yield_exp_map
                    expyield_plum7.maps_YXvy = expyield_plum7.maps_YXvy(:,:,:,end-9:end) ;
                    expyield_plum7.yearList = expyield_plum7.yearList(end-9:end) ;
                    save(lastdecade_out,'expyield_plum7',v73_or_append(lastdecade_out)) ;
                end
                clear expyield_plum7
            end
        end
        
        if do_save.yield_map
            yield.maps_YXvy = yield.maps_YXvy(:,:,:,end-9:end) ;
            yield.yearList = yield.yearList(end-9:end) ;
            save(lastdecade_out,'yield',v73_or_append(lastdecade_out)) ;
        end
        
    end
    clear yield
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save LU and crop time series %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.LU
        LUarea_ts_ntrl = getTS(LU,'NATURAL',land_area_YX)*cf_lu ;
        LUarea_ts_bare = getTS(LU,'BARREN',land_area_YX)*cf_lu ;
        LUarea_ts_crop = getTS(LU,'CROPLAND',land_area_YX)*cf_lu ;
        LUarea_ts_past = getTS(LU,'PASTURE',land_area_YX)*cf_lu ;
        save(timeseries_out,'LUarea_ts_ntrl','LUarea_ts_bare','LUarea_ts_crop','LUarea_ts_past',v73_or_append(timeseries_out)) ;
        clear *_ts_*
    end
    if do_save.crops
        for c = 1:length(cropTypes)
            thisCrop = cropTypes{c} ;
            eval(['croparea_ts_' thisCrop '_lpjg = getTS(cropareas,thisCrop,ones(size(land_area_YX)))*cf_lu ;']) ;
            save(timeseries_out,['croparea_ts_' thisCrop '_lpjg'],v73_or_append(timeseries_out)) ;
            clear thisCrop
        end; clear c
        clear *_ts_*
    end
    if do_save.crops_plum7 && fudge2plum7
        for c = 1:length(cropTypes_plum7)
            thisCrop = cropTypes_plum7{c} ;
            eval(['croparea_ts_' thisCrop '_plum = getTS(cropareas_plum7,thisCrop,ones(size(land_area_YX)))*cf_lu ;']) ;
            save(timeseries_out,['croparea_ts_' thisCrop '_plum'],v73_or_append(timeseries_out)) ;
            clear thisCrop
        end; clear c
        clear *_ts_*
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save irrigation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.irrig
        irrig = lpjgu_matlab_readTable_then2map([inDir 'irrigation.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
%         gsirrig = lpjgu_matlab_readTable_then2map([inDir 'gsirrigation.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving irrigation...')
%         [gsirrig, ~, ~] = CrOp_and_CrOpi(gsirrig, 'gsirrig', cropTypes, replace_CrOp_with_CrOpi, merge_CrOpi_into_CrOp, ...
%             '', cropfracs_orig, inds_cropTypes_cropFracsOrig, inds_cropTypesI_cropFracsOrig) ;
%         gsirrig.maps_YXvy(:,:,strcmp(gsirrig.varNames,'CC3G_ic'),:) = [] ;
%         gsirrig.varNames(strcmp(gsirrig.varNames,'CC3G_ic')) = [] ;
%         gsirrig.maps_YXvy(:,:,strcmp(gsirrig.varNames,'CC4G_ic'),:) = [] ;
%         gsirrig.varNames(strcmp(gsirrig.varNames,'CC4G_ic')) = [] ;
        if is_baseline
            irrig.list_to_map = PLUMout_gridlist.list_to_map ;
            irrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(irrig.varNames) 1])) = NaN ;
%             gsirrig.list_to_map = PLUMout_gridlist.list_to_map ;
%             gsirrig.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(gsirrig.varNames) 1])) = NaN ;
        end
        irrig_ts = getTS(irrig,'Total',land_area_YX) * cf_water ;
        save(timeseries_out,'irrig_ts',v73_or_append(timeseries_out)) ;
%         for c = 1:length(cropTypes)
%             thisCrop = cropTypes{c} ;
%             eval(['gsirrig_ts_' thisCrop '_lpjg = getTS(gsirrig,thisCrop,cropareas.maps_YXvy(:,:,strcmp(cropareas.varNames,thisCrop),:)) * cf_water ;']) ;
%             save(timeseries_out,['gsirrig_ts_' thisCrop '_lpjg'],v73_or_append(timeseries_out)) ;
%             clear thisCrop gsirrig_ts_*
%         end; clear c
%         if fudge2plum7
%             % First, get new structure with new varNames and maps_YXvy
%             if is_baseline
%                 gsirrig_plum7 = plum7ify_struct(...
%                     gsirrig, cropTypes_conv_plum7, cropTypes_conv_lpjg, ...
%                     cropareas.varNames, cropareas_plum7.varNames,'specified',outFrom) ;
%             else
%                 gsirrig_plum7 = plum7ify_struct(...
%                     gsirrig, cropTypes_conv_plum7, cropTypes_conv_lpjg, ...
%                     cropareas.varNames, cropareas_plum7.varNames,'wtd_av',[]) ;
%             end
%             
%             % Then calculate and save
%             for c = 1:Ncrops_plum7
%                 thisCrop_plum7 = cropTypes_plum7{c} ;
%                 eval(['gsirrig_ts_' thisCrop_plum7 '_plum = getTS(gsirrig_plum7,thisCrop_plum7,cropareas_plum7.maps_YXvy(:,:,strcmp(cropareas_plum7.varNames,thisCrop_plum7),:)) * cf_water ;']) ;
%                 save(timeseries_out,['gsirrig_ts_' thisCrop_plum7 '_plum'],v73_or_append(timeseries_out)) ;
%                 clear thisCrop_plum7 gsirrig_ts_*
%             end; clear c
%             clear gsirrig_plum7
%         end
        irrig.maps_YXvy = irrig.maps_YXvy(:,:,:,end-9:end) ;
        irrig.yearList = irrig.yearList(end-9:end) ;
        if ~exist(lastdecade_out,'file')
            save(lastdecade_out,'irrig') ;
        else
            save(lastdecade_out,'irrig',v73_or_append(lastdecade_out)) ;
        end
        clear irrig
%         gsirrig.maps_YXvy = gsirrig.maps_YXvy(:,:,:,end-9:end) ;
%         gsirrig.yearList = gsirrig.yearList(end-9:end) ;
%         save(lastdecade_out,'gsirrig',v73_or_append(lastdecade_out)) ;
%         clear gsirrig
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save other annual water %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.water && anyfileexist([inDir 'awater.out'])
        awater = lpjgu_matlab_readTable_then2map([inDir 'awater.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving evaporation, transpiration, and runoff...')
        aevap_ts = getTS(awater,'Evap',land_area_YX) * cf_water ;
        save(timeseries_out,'aevap_ts',v73_or_append(timeseries_out)) ;
        aaet_ts = getTS(awater,'Transp',land_area_YX) * cf_water ;
        save(timeseries_out,'aaet_ts',v73_or_append(timeseries_out)) ;
        tot_runoff_ts = getTS(awater,'Runoff',land_area_YX) * cf_water ;
        save(timeseries_out,'tot_runoff_ts',v73_or_append(timeseries_out)) ;
        clear *_ts
        awater.maps_YXvy = awater.maps_YXvy(:,:,:,end-9:end) ;
        awater.yearList = awater.yearList(end-9:end) ;
        save(lastdecade_out,'awater',v73_or_append(lastdecade_out)) ;
        clear awater
    elseif do_save.water
        % Evaporation
        if anyfileexist([inDir 'aevap.out'])
            aevap = lpjgu_matlab_readTable_then2map([inDir 'aevap.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
            disp('   Saving evaporation...')
            if is_baseline
                aevap.list_to_map = PLUMout_gridlist.list_to_map ;
                aevap.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aevap.varNames) 1])) = NaN ;
            end
            aevap_ts = getTS(aevap,'Total',land_area_YX) * cf_water ;
            save(timeseries_out,'aevap_ts',v73_or_append(timeseries_out)) ;
            aevap.maps_YXvy = aevap.maps_YXvy(:,:,:,end-9:end) ;
            aevap.yearList = aevap.yearList(end-9:end) ;
            save(lastdecade_out,'aevap',v73_or_append(lastdecade_out)) ;
            clear aevap
        else
            warning('Evaporation missing! Skipping.')
        end
        
        % Evapotranspiration (mm/yr)
        aaet = lpjgu_matlab_readTable_then2map([inDir 'aaet.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving aaet...')
        if is_baseline
            aaet.list_to_map = PLUMout_gridlist.list_to_map ;
            aaet.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aaet.varNames) 1])) = NaN ;
        end
        aaet_ts = getTS(aaet,'Total',land_area_YX) * cf_water ;
        save(timeseries_out,'aaet_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        aaet.maps_YXvy = aaet.maps_YXvy(:,:,:,end-9:end) ;
        aaet.yearList = aaet.yearList(end-9:end) ;
        save(lastdecade_out,'aaet',v73_or_append(lastdecade_out)) ;
        clear aaet
        
        % Runoff (mm/yr)
        tot_runoff = lpjgu_matlab_readTable_then2map([inDir 'tot_runoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving runoff...')
        if is_baseline
            tot_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            tot_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(tot_runoff.varNames) 1])) = NaN ;
        end
        tot_runoff_ts = getTS(tot_runoff,'Total',land_area_YX) * cf_water ;
        save(timeseries_out,'tot_runoff_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        tot_runoff.maps_YXvy = tot_runoff.maps_YXvy(:,:,:,end-9:end) ;
        tot_runoff.yearList = tot_runoff.yearList(end-9:end) ;
        save(lastdecade_out,'tot_runoff',v73_or_append(lastdecade_out))
        clear tot_runoff
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save carbon (kgC/m2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.carbon
        cpool = lpjgu_matlab_readTable_then2map([inDir 'cpool.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving cpool...')
        if is_baseline
            cpool.list_to_map = PLUMout_gridlist.list_to_map ;
            cpool.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(cpool.varNames) 1])) = NaN ;
        end
        % cpool_table = lpjgu_matlab_readTable([inDir 'cpool.out'],'do_save_mat',true) ;
        cpool_ts_VegC = getTS(cpool,'VegC',land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_VegC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_LitterSoilC = getTS(cpool,{'LitterC','SoilC'},land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_LitterSoilC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_HarvSlowC = getTS(cpool,'HarvSlowC',land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_HarvSlowC',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool_ts_Total = getTS(cpool,'Total',land_area_YX_m2) * cf_cpool ;
        save(timeseries_out,'cpool_ts_Total',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        cpool.maps_YXvy = cpool.maps_YXvy(:,:,:,end-9:end) ;
        cpool.yearList = cpool.yearList(end-9:end) ;
        save(lastdecade_out,'cpool',v73_or_append(lastdecade_out)) ;
        clear cpool
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save monthly runoff (mm/month) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_save.mrunoff
        mon_runoff = lpjgu_matlab_readTable_then2map([inDir 'mrunoff.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        if is_baseline
            mon_runoff.list_to_map = PLUMout_gridlist.list_to_map ;
            mon_runoff.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(mon_runoff.varNames) 1])) = NaN ;
        end
        mon_runoff.maps_YXvy = mon_runoff.maps_YXvy(:,:,:,end-9:end) ;
        mon_runoff.yearList = mon_runoff.yearList(end-9:end) ;
        save(lastdecade_out,'mon_runoff',v73_or_append(lastdecade_out))
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save albedo %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.albedo
        snowdepth = lpjgu_matlab_readTable_then2map([inDir 'msnowdepth.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Getting and saving albedo...') 
        if is_baseline
            snowdepth.list_to_map = PLUMout_gridlist.list_to_map ;
            snowdepth.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(snowdepth.varNames) 1])) = NaN ;
        end
        [albedo_jan_YXy,albedo_jul_YXy] = ...
            get_albedo(fpc, snowdepth, LU, cropfracs, baresoil_albedo_YX, land_area_YX, cropTypes, pftList) ;
% % % % %         albedo_ts.jan = squeeze(nansum(nansum(albedo_jan_YXy.*(repmat(land_area_YX,[1 1 size(albedo_jan_YXy,3)]) / nansum(nansum(land_area_YX))), 2), 1)) ;
% % % % %         albedo_ts.jul = squeeze(nansum(nansum(albedo_jul_YXy.*(repmat(land_area_YX,[1 1 size(albedo_jan_YXy,3)]) / nansum(nansum(land_area_YX))), 2), 1)) ;
        albedo1_ts = squeeze(nansum(nansum(albedo_jan_YXy.*(repmat(land_area_YX,[1 1 size(albedo_jan_YXy,3)]) / nansum(nansum(land_area_YX))), 2), 1)) ;
        albedo7_ts = squeeze(nansum(nansum(albedo_jul_YXy.*(repmat(land_area_YX,[1 1 size(albedo_jan_YXy,3)]) / nansum(nansum(land_area_YX))), 2), 1)) ;

    %     disp(['Jan.: ' num2str(mean(albedo_ts.jan(end-9:end)))])
    %     disp(['Jul.: ' num2str(mean(albedo_ts.jul(end-9:end)))])
% % % % %         save(timeseries_out,'albedo_ts',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'albedo1_ts',v73_or_append(timeseries_out)) ;
        save(timeseries_out,'albedo7_ts',v73_or_append(timeseries_out)) ;
        clear albedo*_ts
        albedo.list_to_map = snowdepth.list_to_map ;
        albedo.varNames = {'January','July'} ;
        tmp_YXyv = nan([size(baresoil_albedo_YX) Nyears 2],'single') ;
        tmp_YXyv(:,:,:,1) = albedo_jan_YXy ;
        tmp_YXyv(:,:,:,2) = albedo_jul_YXy ;
        tmp_YXyv = tmp_YXyv(:,:,end-9:end,:) ;
        albedo.maps_YXvy = permute(tmp_YXyv,[1 2 4 3]) ;
        clear tmp_YXyv
        albedo.yearList = snowdepth.yearList(end-9:end) ;
        save(lastdecade_out,'albedo',v73_or_append(lastdecade_out))
        clear snowdepth albedo
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save BVOCs (isoprene, monoterpenes) (mgC/m2) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.bvocs
        if anyfileexist([inDir 'aiso_smry.out'])
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            aiso = lpjgu_matlab_readTable_then2map([inDir 'aiso.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        disp('   Saving aiso...')
        if is_baseline
            aiso.list_to_map = PLUMout_gridlist.list_to_map ;
            aiso.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(aiso.varNames) 1])) = NaN ;
        end
        aiso_ts = getTS(aiso,'Total',land_area_YX_m2) * cf_bvoc ;
        save(timeseries_out,'aiso_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        aiso.maps_YXvy = aiso.maps_YXvy(:,:,:,end-9:end) ;
        aiso.yearList = aiso.yearList(end-9:end) ;
        save(lastdecade_out,'aiso',v73_or_append(lastdecade_out))
        clear aiso
        if anyfileexist([inDir 'amon_smry.out'])
            amon = lpjgu_matlab_readTable_then2map([inDir 'amon_smry.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        else
            amon = lpjgu_matlab_readTable_then2map([inDir 'amon.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        end
        disp('   Saving amon...')
        if is_baseline
            amon.list_to_map = PLUMout_gridlist.list_to_map ;
            amon.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(amon.varNames) 1])) = NaN ;
        end
        amon_ts = getTS(amon,'Total',land_area_YX_m2) * cf_bvoc ;
        save(timeseries_out,'amon_ts',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        amon.maps_YXvy = amon.maps_YXvy(:,:,:,end-9:end) ;
        amon.yearList = amon.yearList(end-9:end) ;
        save(lastdecade_out,'amon',v73_or_append(lastdecade_out))
        clear amon
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import and save N flux (kgN/ha) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.Nflux
        nflux = lpjgu_matlab_readTable_then2map([inDir 'nflux.out'],'force_mat_save',true,'list_to_map_in',list_to_map) ;
        disp('   Saving nflux...')
        if is_baseline
            nflux.list_to_map = PLUMout_gridlist.list_to_map ;
            nflux.maps_YXvy(~repmat(PLUMout_mask_YX1y,[1 1 length(nflux.varNames) 1])) = NaN ;
        end
        nflux_ts_fert = getTS(nflux,'fert',land_area_YX_m2*1e-4) * cf_nflux ;
        save(timeseries_out,'nflux_ts_fert',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        nflux_ts_flux = getTS(nflux,'flux',land_area_YX_m2*1e-4) * cf_nflux ;
        save(timeseries_out,'nflux_ts_flux',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        nflux_ts_leach = getTS(nflux,'leach',land_area_YX_m2*1e-4) * cf_nflux ;
        save(timeseries_out,'nflux_ts_leach',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        nflux_ts_harvest = getTS(nflux,'harvest',land_area_YX_m2*1e-4) * cf_nflux ;
        save(timeseries_out,'nflux_ts_harvest',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        nflux_ts_LU_ch = getTS(nflux,'LU_ch',land_area_YX_m2*1e-4) * cf_nflux ;
        save(timeseries_out,'nflux_ts_LU_ch',v73_or_append(timeseries_out)) ;
        clear *_ts_*
        nflux.maps_YXvy = nflux.maps_YXvy(:,:,:,end-9:end) ;
        nflux.yearList = nflux.yearList(end-9:end) ;
        save(lastdecade_out,'nflux',v73_or_append(lastdecade_out))
        clear nflux
    end

    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Save FPC stuff %%%
    %%%%%%%%%%%%%%%%%%%%%%

    if do_save.fpc
        disp('   Saving FPC...')
        fpc.maps_YXvy = fpc.maps_YXvy(:,:,:,end-9:end) ;
        yearList = yearList(end-9:end) ;
        save(lastdecade_out,'fpc',v73_or_append(lastdecade_out)) ;
        clear fpc
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save land use and crop maps %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_save.crops
        disp('   Saving crop maps...')
        cropfracs.maps_YXvy = cropfracs.maps_YXvy(:,:,:,end-9:end) ;
        cropfracs.yearList = cropfracs.yearList(end-9:end) ;
        save(lastdecade_out,'cropfracs',v73_or_append(lastdecade_out)) ;
        
    end
    if do_save.crops_plum7
        disp('   Saving crop_plum7 maps...')
        cropfracs_plum7.maps_YXvy = cropfracs_plum7.maps_YXvy(:,:,:,end-9:end) ;
        cropfracs_plum7.yearList = cropfracs_plum7.yearList(end-9:end) ;
        save(lastdecade_out,'cropfracs_plum7',v73_or_append(lastdecade_out)) ;
    end
    if do_save.LU
        disp('   Saving land use maps...')
        LU.maps_YXvy = LU.maps_YXvy(:,:,:,end-9:end) ;
        LUorder = nan(1,length(LUlist),'single') ;
        for L = 1:length(LU.varNames)
            thisLU = LUlist{L} ;
            LUorder(L) = find(strcmp(LU.varNames,thisLU)) ;
        end
        LU.maps_YXvy = LU.maps_YXvy(:,:,LUorder,:) ;
        LU.varNames = LUlist ;
        LU.yearList = LU.yearList(end-9:end) ;
        save(lastdecade_out,'LU',v73_or_append(lastdecade_out)) ;
    end
    clear LU cropfracs* cropareas*
    
    
    %%%%%%%%%%%%%%%%%%%%
    %%% Housekeeping %%%
    %%%%%%%%%%%%%%%%%%%%
    
    if flip_do_save_yield_exp
        do_save.yield_exp = ~do_save.yield_exp ;
    end
    if flip_do_save_yield_exp_map
        do_save.yield_exp_map = ~do_save.yield_exp_map ;
    end
    
    disp('Done.')
end
