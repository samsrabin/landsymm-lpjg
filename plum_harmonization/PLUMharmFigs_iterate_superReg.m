function [table_out, table_out_relY1] = PLUMharmFigs_iterate_superReg( ...
    area_x1y, s, orig_or_harm, ...
    list_superRegs, list_regions, harm_focus_regions, ...
    biomeID_x, countries_x, countries_key, lats, lons, ...
    list2map)

thisSuperReg = list_superRegs{s} ;
max_superReg_nameLength = max(cellfun(@length,list_superRegs)) ;
max_reg_nameLength = max(cellfun(@length,list_regions)) ;
max_nameLength = max(max_superReg_nameLength, max_reg_nameLength) ;
list_superRegs_pad = pad(list_superRegs, max_nameLength+4, 'right', '-') ;
list_regions_pad = pad(list_regions, max_nameLength+4, 'left', ' ') ;
thisSuperReg_pad = list_superRegs_pad{s} ;
region_in_superReg = strcmp(harm_focus_regions(:,2), thisSuperReg) ;
subset = harm_focus_regions(region_in_superReg,:) ;
incl_regions = cell2mat(harm_focus_regions(region_in_superReg,1)) ;
incl_biomes = unique(subset(:,3), 'stable') ;

list_gain = NaN ;
list_loss = NaN ;
list_net = NaN ;
list_gain_relY1 = NaN ;
list_loss_relY1 = NaN ;
list_net_relY1 = NaN ;
list_subreg = {thisSuperReg} ;
list_biome = {'TOTAL'} ;
Ntotals = length(subset(:,1))>1 ;
superReg_gain = 0 ;
superReg_loss = 0 ;
superReg_net  = 0 ;
superReg_orig = 0 ;

for s2 = 1:length(subset(:,1))
    is_cell_in_reg = false(size(biomeID_x)) ;
    tmp_incl_countries = subset{s2,4} ;
    tmp_incl_regions = subset{s2,1} ;
    thisReg_index = find(strcmp(harm_focus_regions(:,5),subset{s2,5})) ;
    thisReg_name = list_regions{thisReg_index} ;
    thisBiome = subset{s2,3} ;
    do_biome_total = length(find(strcmp(subset(:,3), thisBiome))) > 1 ;
    if do_biome_total && s2==min(find(strcmp(subset(:,3), thisBiome)))
        Ntotals = Ntotals + 1 ;
        biome_total_index = length(list_gain) + 1 ;
        list_gain(biome_total_index,1) = 0 ;
        list_loss(biome_total_index,1) = 0 ;
        list_net(biome_total_index,1) = 0 ;
        list_gain_relY1(biome_total_index,1) = 0 ;
        list_loss_relY1(biome_total_index,1) = 0 ;
        list_net_relY1(biome_total_index,1) = 0 ;
        list_biome{biome_total_index,1} = thisBiome ;
        list_subreg{biome_total_index,1} = 'TOTAL' ;
        list_superReg_out{biome_total_index,1} = thisSuperReg ;
        biome_orig = 0 ;
    end
    if isempty(tmp_incl_countries)
        is_cell_in_reg = is_cell_in_reg | ...
            ismember(biomeID_x, tmp_incl_regions) ;
    elseif iscell(tmp_incl_countries)
        incl_countryCodes = countries_key.numCode( ...
            contains(countries_key.Country, tmp_incl_countries)) ;
        is_cell_in_reg = is_cell_in_reg | ...
            (ismember(biomeID_x, tmp_incl_regions) ...
            & ismember(countries_x, incl_countryCodes)) ;
    elseif contains(tmp_incl_countries,'&')
        % Parse
        tmp = strsplit(tmp_incl_countries, '&') ;
        bound = str2double(tmp{1}) ;
        news = tmp{2} ;
        % Restrict
        if strcmpi(news, 'n')
            is_cell_in_reg = is_cell_in_reg | ...
                (ismember(biomeID_x, tmp_incl_regions) ...
                & lats >= bound) ;
        elseif strcmpi(news, 's')
            is_cell_in_reg = is_cell_in_reg | ...
                (ismember(biomeID_x, tmp_incl_regions) ...
                & lats <= bound) ;
        elseif strcmpi(news, 'e')
            is_cell_in_reg = is_cell_in_reg | ...
                (ismember(biomeID_x, tmp_incl_regions) ...
                & lons >= bound) ;
        elseif strcmpi(news, 'w')
            is_cell_in_reg = is_cell_in_reg | ...
                (ismember(biomeID_x, tmp_incl_regions) ...
                & lons <= bound) ;
        else
            error('???')
        end
    else
        error('Unsure how to parse tmp_incl_countries')
    end
    clear tmp*
    [gain, loss, orig] = get_region_numbers(area_x1y(is_cell_in_reg,:,:)) ;
    net = gain + loss ;
    superReg_gain = superReg_gain + gain ;
    superReg_loss = superReg_loss + loss ;
    superReg_net  = superReg_net  + net ;
    superReg_orig = superReg_orig + orig ;
    if length(subset(:,1))==1
        list_gain = gain ;
        list_loss = loss ;
        list_net = net ;
    else
        list_gain(end+1,1) = gain ; %#ok<*AGROW>
        list_loss(end+1,1) = loss ;
        list_net(end+1,1) = net ;
        list_gain_relY1(end+1,1) = gain/orig ;
        list_loss_relY1(end+1,1) = loss/orig ;
        list_net_relY1(end+1,1) = net/orig ;
        list_subreg{end+1,1} = thisReg_name ;
        list_biome{end+1,1} = subset{s2,3} ;
    end
    if do_biome_total
        list_gain(biome_total_index) = list_gain(biome_total_index) + gain ;
        list_loss(biome_total_index) = list_loss(biome_total_index) + loss ;
        list_net(biome_total_index)  = list_net(biome_total_index)  + net ;
        biome_orig = biome_orig + orig ;
        if s2==max(find(strcmp(subset(:,3), thisBiome)))
            list_gain_relY1(biome_total_index) = list_gain(biome_total_index) / biome_orig ;
            list_loss_relY1(biome_total_index) = list_loss(biome_total_index) / biome_orig ;
            list_net_relY1(biome_total_index)  = list_net(biome_total_index)  / biome_orig ;
        end
    end
end

% Print to console
if false
% plussign = '' ;
% if net>=0
%     plussign = '+' ;
% end
% minussign = '' ;
% if loss==0
%     minussign = '-' ;
% end
% fprintf('%s %s:  +%0.1e  %s%0.1e = %s%0.1e\n',  ...
%     thisSuperReg_pad, orig_or_harm, gain, minussign, loss, plussign, net)
% if length(incl_regions) > 1
%     for s2 = 1:length(subset(:,1))
%         % Get region info
%         thisReg = subset{s2,1} ;
%         thisReg_index = find(strcmp(harm_focus_regions(:,5),subset{s2,5})) ;
%         thisReg_name = list_regions{thisReg_index} ;
%         thisReg_name_pad = list_regions_pad{thisReg_index} ;
%         is_cell_in_reg = ismember(biomeID_x, thisReg) ;
%         if ~any(is_cell_in_reg)
%             error('No matching cells found')
%         end
%         % Get region values
%         [gain, loss] = ...
%             get_region_numbers(area_x1y(is_cell_in_reg,:,:)) ;
%         net = gain + loss ;
%         list_gain(end+1) = gain ;
%         list_loss(end+1) = loss ;
%         list_net(end+1)  = net ;
%         list_subreg{end+1,1} = thisReg_name ;
%         list_biome{end+1,1} = subset{s2,3} ;
%         % Display result
%         plussign = '' ;
%         if net>=0
%             plussign = '+' ;
%         end
%         minussign = '' ;
%         if loss==0
%             minussign = '-' ;
%         end
%         fprintf('     %s:  +%0.1e  %s%0.1e = %s%0.1e\n',  ...
%             thisReg_name_pad, gain, minussign, loss, plussign, net)
%     end
% end
end

if length(subset(:,1))==1
    list_superReg_out = string(thisSuperReg) ;
    list_subreg = string(subset{1,2}) ;
    list_biome = string(subset{1,3}) ;
else
    list_superReg_out = cellstr(repmat(thisSuperReg, [length(subset(:,1))+Ntotals 1])) ;
    list_subreg{1} = 'TOTAL' ;
    list_gain(1) = nansum(list_gain) ;
    list_loss(1) = nansum(list_loss) ;
    list_net(1) = nansum(list_net) ;
end
list_gain_relY1(1) = superReg_gain / superReg_orig ;
list_loss_relY1(1) = superReg_loss / superReg_orig ;
list_net_relY1(1)  = superReg_net  / superReg_orig ;

table_out = table( ...
    list_superReg_out, ...
    list_biome, ...
    list_subreg, ...
    shiftdim(list_gain), ...
    shiftdim(list_loss), ...
    shiftdim(list_net)) ;
table_out_relY1 = table( ...
    list_superReg_out, ...
    list_biome, ...
    list_subreg, ...
    shiftdim(list_gain_relY1), ...
    shiftdim(list_loss_relY1), ...
    shiftdim(list_net_relY1)) ;

table_out.Properties.VariableNames = ...
    {'Super-region', 'Biome', 'Sub-region', 'Gain', 'Loss', 'Net'} ;
table_out_relY1.Properties.VariableNames = ...
    {'Super-region', 'Biome', 'Sub-region', 'Gain', 'Loss', 'Net'} ;

end



function [gain, loss, orig] = get_region_numbers(area_x1y)

orig_x = area_x1y(:,:,1) ;
diff_x = ...
    area_x1y(:,:,end) ...
    - orig_x ;
gain = sum(diff_x(diff_x>0)) ;
loss = sum(diff_x(diff_x<0)) ;
orig = sum(orig_x) ;

end


