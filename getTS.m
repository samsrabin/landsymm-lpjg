function ts_out = getTS(in_struct, in_vars, land_area)

Nyears = length(in_struct.yearList) ;
if Nyears ~= size(in_struct.maps_YXvy,4)
    error('Nyears ~= size(in_struct.maps_YXvy,4)')
end

[~,IA] = intersect(in_struct.varNames,in_vars) ;

if size(land_area,4)==1
    land_area_YX1y = repmat(land_area,[1 1 1 Nyears]) ;
elseif ~isempty(land_area)
    land_area_YX1y = land_area ;
end

in_map_YX1y = sum(in_struct.maps_YXvy(:,:,IA,:),3) ;
if ~isempty(land_area)
    if ~isequal(size(in_map_YX1y), size(land_area_YX1y))
        error('~isequal(size(in_map_YX1y), size(land_area_YX1y')
    end
    in_map_YX1y = in_map_YX1y .* land_area_YX1y ;
end
ts_out = squeeze(nansum(nansum(in_map_YX1y,2),1)) ;

end


% function ts_out = getTS(in_struct, in_vars, land_area)
% 
% if isstruct(in_struct)
%     
%     Nyears = size(in_struct.maps_YXvy,4) ;
%     
%     if Nyears ~= length(in_struct.yearList)
%         error('Nyears ~= size(in_struct.maps_YXvy,4)')
%     end
%     
%     if ~isempty(in_vars)
%         [~,IA] = intersect(in_struct.varNames,in_vars) ;
%     else
%         IA = 1:size(in_struct.maps_YXvy,3) ;
%     end
%     
%     if size(land_area,4)==1
%         land_area_YX1y = repmat(land_area,[1 1 1 Nyears]) ;
%     elseif ~isempty(land_area)
%         land_area_YX1y = land_area ;
%     end
%     
%     in_map_YX1y = sum(in_struct.maps_YXvy(:,:,IA,:),3) ;
%     if ~isempty(land_area)
%         if ~isequal(size(in_map_YX1y), size(land_area_YX1y))
%             error('~isequal(size(in_map_YX1y), size(land_area_YX1y')
%         end
%         in_map_YX1y = in_map_YX1y .* land_area_YX1y ;
%     end
%     ts_out = squeeze(nansum(nansum(in_map_YX1y,2),1)) ;
%     
%     
% else
%     
%     Nyears = size(in_struct,4) ;
%     
%     if size(land_area,4)==1
%         land_area_YX1y = repmat(land_area,[1 1 1 Nyears]) ;
%     elseif ~isempty(land_area)
%         land_area_YX1y = land_area ;
%     end
%     
%     in_map_YX1y = sum(in_struct,3) ;
%     if ~isempty(land_area)
%         if ~isequal(size(in_map_YX1y), size(land_area_YX1y))
%             error('~isequal(size(in_map_YX1y), size(land_area_YX1y')
%         end
%         in_map_YX1y = in_map_YX1y .* land_area_YX1y ;
%     end
%     ts_out = squeeze(nansum(nansum(in_map_YX1y,2),1)) ;
%     
% 
% end

