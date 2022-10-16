
% loop through WGS84 lat/long pairs and then find the OSTN value for each.
% This should be a lot quicker than manually calculating it again each
% time!

% wont work, takes too much space

latres = 0.0001;
latrange = [49 59];
longres = 0.0001;
longrange = [-8 2];

    
%%

% todo: offset the negative in the array index

% offset the negative value as index values cant be negative
% only applies to long, this isn't australia
longnegativeoffset = 0;
if(longrange(1) < 0)
    longnegativeoffset = (0 - longrange(1)) * (1/longres);
end 

% space efficiency offset, we don't want ~50,000 zero rows!
% only applies to lat
latefficiencyoffset = (-latrange(1)*(1/latres));
longefficiencyoffset = 0;

% initialise array
WGS2OSTN = zeros( ( (longrange(2)-longrange(1)) * (1/longres)), ...
                  ( (latrange(2)-latrange(1))   *(1/latres) ), ...
                  2);

% [0 0] if out of range
for lat= latrange(1) : latres : latrange(2)
    for long = longrange(1) : longres : longrange(2)
        % find the index values for lat and long in the array
        %latidx = (latefficiencyoffset+ (lat*(1/latres)) ) + 1;
        latidx = int32( (lat*(1/latres)) - (latrange(1)*(1/latres)) ) +1;
        %longidx = ((long*(1/longres)) + longnegativeoffset) +1;
        longidx = int32((20001+ (long*(1/longres)) )); %- (longrange(1)*(1/longres)) ) +1;
        try
            WGS2OSTN(latidx, longidx, 1:2)  = OSTN15_Matlab('gps-to-grid',lat, long );
        catch
        end
        %WGS2OSTN(longidx, latidx, 1:2) = [latidx longidx];
        [lat long]
    end
end

WGS2OSTN(1,1,1:2)


%%