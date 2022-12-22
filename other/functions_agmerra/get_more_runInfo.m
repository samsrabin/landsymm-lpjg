function runInfo = get_more_runInfo(runInfo)


runInfo.thisVer = '1.1.1';
runInfo.comment = sprintf('Crop model output for ISIMIP%s (GGCMI phase 3). GDD threshold = %0.2f.', runInfo.phase, runInfo.gdd_thresh) ;
runInfo.modelname='lpj-guess';
runInfo.institution='Karlsruhe Institute of Technology, Institute of Meteorology and Climate Research / Atmospheric Environmental Research, Garmisch-Partenkirchen, Germany';
runInfo.contact='Sam Rabin, sam.rabin@kit.edu';
if strcmp(runInfo.phase, '2a') || strcmp(runInfo.phase, '3a') || strcmp(runInfo.phase, 'XX')
   runInfo.baseyear = 1901 ;
elseif strcmp(runInfo.phase, '2b')
   runInfo.baseyear = 1661 ;
elseif strcmp(runInfo.phase, '3b')
   runInfo.baseyear = 1601 ;
elseif strcmp(runInfo.phase, 'EV')
	runInfo.baseyear = 1980 ;
else
   error('Undefined baseyear for runInfo.phase: %s', runInfo.phase)
end



end