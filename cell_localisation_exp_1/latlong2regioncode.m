function [regionCode] = latlong2regioncode(lat, long)
% custom 1 degree region codes
% this saves resorting everything via a 4D array!
% (360*180) + 180 % latitude stripes at each longitude
regionCode = (floor(180+long)*180) + floor(90+lat);

end

