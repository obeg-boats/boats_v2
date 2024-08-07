% function_map_2_vec(array_original)
%-----------------------------------------------------------------------------------------
% convert a 2-dimensional map of sites to a vector of sites
% remove land, high latitude, and open-ocean sites (not in an LME)
%-----------------------------------------------------------------------------------------

function [array_new indlat indlon] = function_map_2_vec(array_original,mask)

 [illon illat] = meshgrid([1:size(mask,2)],[1:size(mask,1)]);
 iuse = find(mask==0); 
 array_new = array_original(iuse);
 indlat = illat(iuse);
 indlon = illon(iuse);

end % function

%----------------------------------------------------------------------------------------
% END OF SCRIPT
