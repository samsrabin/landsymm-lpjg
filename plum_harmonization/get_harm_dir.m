function harmDir = get_harm_dir(inDir, fruitveg_sugar_2oil, combineCrops)
% Append option strings to directory name as needed.

harmDir = inDir ;

if fruitveg_sugar_2oil
    harmDir = [harmDir '.fvs2oil'] ;
end
if combineCrops
    harmDir = [harmDir '.combineCrops'] ;
end

end