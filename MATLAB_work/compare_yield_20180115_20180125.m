%% Import

yield1 = lpjgu_matlab_readTable_then2map('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/PLUM2LPJGblIrrNoFG_26_s1_1850-2005_moreProcs/output-2018-01-06-071106/yield.out.gz','force_mat_save',true) ;
yield2 = lpjgu_matlab_readTable_then2map('/Volumes/WDMPP_Storage/Shared/PLUM/trunk_runs/PLUM2LPJG_PLUM7_1850-2005/output-2018-01-23-195509/yield.out','force_mat_save',true) ;


%% Map

tmp1 = mean(yield1.maps_YXvy(:,:,strcmp(yield1.varNames,'TrRi'),end-9:end),4) ;
tmp2 = mean(yield2.maps_YXvy(:,:,strcmp(yield2.varNames,'Rice'),end-9:end),4) ;
% tmp1 = mean(yield1.maps_YXvy(:,:,strcmp(yield1.varNames,'TrRii'),end-9:end),4) ;
% tmp2 = mean(yield2.maps_YXvy(:,:,strcmp(yield2.varNames,'Ricei'),end-9:end),4) ;

tmp3 = (tmp2-tmp1)./tmp1 ;
% tmp3 = (tmp2-tmp1) ;
tmp3(~(tmp1>0 & tmp2>0)) = NaN ;
tmp3(tmp3==0) = NaN ;

shademap(tmp3);

%%
tmp4 = single(tmp1==0 & tmp2>0);
tmp4(isnan(tmp3)) = NaN ;
shademap(tmp4);