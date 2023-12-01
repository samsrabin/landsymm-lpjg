%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make stand and pft list for ins-files %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version for crop mappings
% thisVer = 'WithFruitVegSugar_b' ; remapVer = '7b' ;
thisVer = 'WithFruitVeg_sepSugar_sepOil' ; remapVer = '10_g2p' ;
% thisVer = 'WithFruitVeg_sepSugar_sepOil_sepC3' ; remapVer = '11_g2p' ;
% thisVer = 'ggcmi5' ; remapVer = '12_g2p' ;
% thisVer = 'ggcmi5_preBNF' ; remapVer = '13_g2p' ;

include_cropphencol = true ;


%% Set up

addpath(genpath(landsymm_lpjg_path()))

% out_dir = sprintf('/Users/Shared/G2P/inputs/LU/remaps_v%s/',remapVer) ;
out_dir = '/Users/sam/Downloads/create-ins/ins-all' ;
if ~exist(out_dir,'dir')
    mkdir(out_dir) ;
end


%% Do it

get_st_pft_lists(thisVer, remapVer, include_cropphencol, out_dir)