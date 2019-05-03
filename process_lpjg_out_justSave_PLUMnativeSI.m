%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with LPJ-GUESS outputs for impacts paper: %%%
%%% PLUM-style native %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calib_file = '/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/MATLAB_work/remap5e_v18.csv' ;

% %%% All TRUE except really unnecessary ones
% do_save.LU              = true ;
% do_save.crops           = true ;
% do_save.yield           = true ;
% do_save.yield_exp       = false ;
% do_save.yield_map       = true ;
% do_save.yield_exp_map   = false ;
% do_save.irrig           = true ;
% do_save.water           = true ;
% do_save.carbon          = true ;
% do_save.mrunoff         = true ;
% do_save.albedo          = true ;
% do_save.bvocs           = true ;
% do_save.Nflux           = true ;
% do_save.Nfert           = true ;
% do_save.fpc             = false ;

%%% All FALSE except selected
do_save.LU              = false ;
do_save.crops           = false ;
do_save.yield           = false ;
do_save.yield_exp       = false ;
do_save.yield_map       = false ;
do_save.yield_exp_map   = false ;
do_save.irrig           = false ;
do_save.water           = true ;
do_save.carbon          = false ;
do_save.mrunoff         = true ;
do_save.albedo          = false ;
do_save.bvocs           = false ;
do_save.Nflux           = true ;
do_save.Nfert           = false ;
do_save.fpc             = false ;

inDir_list = {...
%     'LPJGPLUM_1850-2010_remap6p7/output-2019-02-18-120851'
    'LPJGPLUM_2011-2100_harm3_SSP1_RCP45/output-2019-02-27-103914';
%     'LPJGPLUM_2011-2100_harm3_SSP3_RCP60/output-2019-02-27-093027';
%     'LPJGPLUM_2011-2100_harm3_SSP4_RCP60/output-2019-02-27-093259';
%     'LPJGPLUM_2011-2100_harm3_SSP5_RCP85/output-2019-02-27-104120';
%     'LPJGPLUM_2011-2100_harm3_constLU_RCP85/output-2019-03-07-164546';
    } ;



%% Setup

test_cropfracs_20170108 = false ;

tmp = pwd ;
if strcmp(tmp(1:5),'/User')
    onMac = true ;
elseif strcmp(tmp(1:5),'/pfs/')
    onMac = false ;
else
    error('What system are you on?')
end
clear tmp
if onMac
    thisDir = addslashifneeded('/Users/sam/Documents/Dropbox/LPJ-GUESS-PLUM/LPJGP_paper02_Sam/MATLAB_work') ;
    gridlist_file = '/Users/Shared/lpj-guess/gridlists/PLUMout_gridlist.txt' ;
    landarea_file = '/Users/Shared/PLUM/crop_calib_data/other/staticData_quarterdeg.nc' ;
    biomes_dir = addslashifneeded('/Users/sam/Geodata/General/WWF terrestrial ecosystems') ;
else
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/matlab-general/')) ;
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/matlab-general-fromshared/')) ;
    addpath(genpath('/home/fh1-project-lpjgpi/lr8247/lpj-guess-crop-calibration/')) ;
    thisDir = addslashifneeded('/home/fh1-project-lpjgpi/lr8247/paper02-matlab-work') ;
    gridlist_file = '/home/fh1-project-lpjgpi/lr8247/paper02-matlab-work/PLUMout_gridlist.txt' ;
    landarea_file = '/home/fh1-project-lpjgpi/lr8247/PLUM/input/LUH2/supporting/staticData_quarterdeg.nc' ;
    biomes_dir = addslashifneeded('/home/fh1-project-lpjgpi/lr8247/PLUM/input/biomes') ;
end
biomes_map_file = [biomes_dir 'wwf_terr_ecos_UnpackClip.halfDeg.tif'] ;
biomes_key_file = [biomes_dir 'wwf_terr_ecos.codes.csv'] ;

if ~exist(thisDir,'dir')
    error('thisDir does not exist')
end
cd(thisDir) ;
addpath(genpath(pwd))

is_baseline_list = false(length(inDir_list),1) ;
do_PLUMout_gridlist_adjust_list = false(length(inDir_list),1) ;
for d = 1:length(inDir_list)
    if contains(inDir_list{d},'1850-')
        is_baseline_list(d) = true ;
        do_PLUMout_gridlist_adjust_list(d) = true ;
    elseif contains(inDir_list{d},'constLU')
        do_PLUMout_gridlist_adjust_list(d) = true ;
    end
end

pftList_noCrops = {'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE',...
            'TrBR','C3G','C4G','PC3G','PC4G',...
            'CC3G_ic','CC4G_ic','ExtraCrop','Total','Crop_sum','Pasture_sum',...
            'Natural_sum','Barren_sum'} ;
LUlist = {'CROPLAND','PASTURE','NATURAL','BARREN'} ;

% Anonymous function for checking whether any version of a file exists
anyfileexist = @(in_file) ...
    exist(in_file,'file') ...
    || exist([in_file '.maps.mat'],'file') ...
    || exist([in_file '.mat'],'file') ...
    || exist([in_file '.gz'],'file') ;

% Get calibration factors
LPJGcrops_2_PLUM = readtable(calib_file) ;
cropTypes_conv = LPJGcrops_2_PLUM.Crop ;
calibFactors = LPJGcrops_2_PLUM.calibFactor ;
clear LPJGcrops_2_PLUM

% Import bare-soil albedo
if do_save.albedo
    baresoil_albedo_file = [thisDir 'soilmap.txt'] ;
    baresoil_albedo = lpjgu_matlab_readTable_then2map(baresoil_albedo_file,'force_mat_save',true,'verbose',false) ;
    baresoil_albedo_YX = baresoil_albedo.maps_YXv ;
    clear baresoil_albedo
end

% Read PLUMout_gridlist, if needed
do_PLUMout_gridlist_adjust = any(do_PLUMout_gridlist_adjust_list) ;
if do_PLUMout_gridlist_adjust
    PLUMout_gridlist = lpjgu_matlab_readTable_then2map(gridlist_file,'verbose',false,'force_mat_save',true) ;
end

% Import land area (km2 to m2)
gcel_area_YXqd = transpose(ncread(landarea_file,'carea')) ;
land_frac_YXqd = 1 - flipud(transpose(ncread(landarea_file,'icwtr'))) ;
land_area_YXqd = gcel_area_YXqd .* land_frac_YXqd ;
%%%% Convert to half-degree
tmp = land_area_YXqd(:,1:2:1440) + land_area_YXqd(:,2:2:1440) ;
land_area_YX_orig = tmp(1:2:720,:) + tmp(2:2:720,:) ;
tmp = gcel_area_YXqd(:,1:2:1440) + gcel_area_YXqd(:,2:2:1440) ;
gcel_area_YX_orig = tmp(1:2:720,:) + tmp(2:2:720,:) ;
% Convert to m2
land_area_YX_orig = land_area_YX_orig*1e6 ;
gcel_area_YX_orig = gcel_area_YX_orig*1e6 ;
clear tmp gcel_area_YXqd land_frac_YXqd land_area_YXqd

% Import biomes
biomes_YX = flipud(imread(biomes_map_file)) ;
biomes_YX(biomes_YX<0) = NaN ;
biomes_key = readtable(biomes_key_file) ;

% Conversion factors
%%% All masses in kg
%%% All areas in m2
%%% All volumes in m3
cf_cpool = 1 ;   % LPJ-GUESS outputs kgC/m2
cf_water = 1e-3 ; % LPJ-GUESS outputs mm (will be multiplied by land_area_m2 --> m3)
cf_nflux = 1e-4 ;   % LPJ-GUESS outputs kgN/ha
cf_bvoc = 1e-6 ;   % LPJ-GUESS outputs mgC/m2

% Set up for "last 30 years of century" stats
last30_statList = { ...
    'last30_mean' ; ...
    'last30_median' ; ...
    'last30_std' ; ...
    'last30_prctile05' ; ...
%     'last30_prctile25' ; ...
%     'last30_prctile75' ; ...
    'last30_prctile95' ; ...
    'last30_min' ; ...
    'last30_max' ; ...
    } ;
Nstats = length(last30_statList) ;
last30_mean = @(array_YXvy) mean(array_YXvy,4) ;
last30_median = @(array_YXvy) median(array_YXvy,4) ;
last30_std = @(array_YXvy) std(array_YXvy,1,4) ;
last30_prctile05 = @(array_YXvy) prctile(array_YXvy,5,4) ;
last30_prctile25 = @(array_YXvy) prctile(array_YXvy,25,4) ;
last30_prctile75 = @(array_YXvy) prctile(array_YXvy,75,4) ;
last30_prctile95 = @(array_YXvy) prctile(array_YXvy,95,4) ;
last30_min = @(array_YXvy) min(array_YXvy,[],4) ;
last30_max = @(array_YXvy) max(array_YXvy,[],4) ;
last30_statHandles = cell(Nstats,1) ;
for s = 1:Nstats
    last30_statHandles{s} = eval(last30_statList{s}) ;
end; clear s


%% Loop through inDir_list

loop_through_indir_list
