

% Author: Thomas Pile
% 10th March 2021


%%
% determine the number of requests to get the desired map area

% todo
clc
clear
%%
% mass image download method
% usage
% enter the start and finish lat and long locations. lat should be a
% multiple of 0.005 and long a multiple of 0.05
% file name is longitude_latitude.jpg
% 10th March 2021

EW = 1;

% latitude
start_lat = 51;
stop_lat = 51.1;

% longitude
% going west of meridian
start_lon = 0.00; % 0.00 format for filenames!
stop_lon = 0.20;

% estimate size of request first, tested :)
qty_lon = (1+(stop_lon-start_lon)/0.01);
qty_lat = (1+(stop_lat-start_lat)/0.005);
qty_requests = qty_lat*qty_lon;
qty_requests

% holding image to join a sequence in
imheight = [273,1011] ;
imwidth = [314,1246]; % 932 long is 932 pixels wide
%imsequence = ones(imheight,imwidth,3);

i = 0;
for lat = start_lat : 0.005 : stop_lat
    % we request lon to lon+0.005 but we get back lon to lon+0.01 !
    for lon = start_lon : 0.01 : stop_lon
        % generate filename
        strlon = string(lon*2); % since requesting 0.005 gives 0.01!
        strlat = string(lat);
        % append 00 to round longitudes
        if(floor(lon)==lon)
            strlon = strlon+".00";
        end
        if(floor(lat)==lat)
            strlat = strlat+".000";
        end
        % final file name
        filename = "map/"+strlon+"_"+strlat+".jpg";
        % grab the image
        if(EW < 0)
            lms=[(EW*lon+0.001) (EW*lon+0.005+0.001) (lat) (lat+0.005)]; 
        else
            lms=[(lon+0.001) (lon+0.005+0.001) (lat) (lat+0.005)]; 
        end
        axis(lms);
        [Glon1,Glat1,Gimg1]=plot_google_map('maptype','satellite','refresh',0,'autoaxis',0,'APIKey','******************************');
        clf;
        %imshow(Gimg1(imheight(1):imheight(2),imwidth(1):imwidth(2),:))

        % save
        imwrite((Gimg1(imheight(1):imheight(2),imwidth(1):imwidth(2),:)),filename,'Quality',90);

        % add to sequence
        %imsequence(:, (i*932)  1+(i*932) ,:) = Gimg1(1:738,1:932,:);
        i=i+1
    end
end

% east image
%1 = -0.0234 
%1280 = -0.0096

%%
% import image set into map structure
% 11th March 2021

% list all files in a directory
files = dir('map');
% get filename. disregard any named . or .. usually the first 2
%files(3).name

% create a struct to store the map data
mapidx = 1;
map = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)

% loop through files
% skip 1, 2
for i= 3 : sel(size(files),1,1)
    % extract name
    fname = files(i).name;
    % find the fdelimiter location
    fdelim = strfind(fname,'_');
    % extract the longitude value. +1 to skip delim, -4 to skip filetype
    flong = fname(1:fdelim-1);
    % extract the latitude value
    flat = fname(fdelim+1:sel(size(fname),1,2)-4);
    % form path to image
    path = 'map/';
    path(5:4+sel(size(fname),1,2)) = fname(1,:);
     % import respective image
    fimg = imread(path);
    % add data to map structure
    map(mapidx).lat = flat;
    map(mapidx).long = flong;
    map(mapidx).nlat = str2num(flat);
    map(mapidx).nlong = str2num(flong);
    map(mapidx).img = fimg;
    mapidx = mapidx + 1;
end

% % loop through all map structure example
% for i = 1 : sel(size(map),1,2)
%     map(i).featuremap = 0;
% imshow(map(1).img)
% end


% find features and add them to the map structure (SURF)
for i = 1 : sel(size(map),1,2)
    imgg = rgb2gray(map(i).img);
    points = detectSURFFeatures(imgg);
    %imshow(imgg); 
    %hold on;
    %plot(points.selectStrongest(250));
    map(i).featuremap = points.selectStrongest(250);
end

% plot all points on the feature map
%scatter(map(1).featuremap.Location(:,1),map(1).featuremap.Location(:,2))

%%
% sort the struct so it can be printed by going through all indices
% 11th March 2021

% loop through all map structure example
maxlat = 0;
maxlong = 0;
minlong = 999;
minlat = 999;
for i = 1 : sel(size(map),1,2)
    if(map(i).nlat > maxlat)
        maxlat =  map(i).nlat;
    end
    if(map(i).nlat < minlat)
        minlat =  map(i).nlat;
    end
    if(map(i).nlong > maxlong)
        maxlong =  map(i).nlong;
    end
    if(map(i).nlong < minlong)
        minlong =  map(i).nlong;
    end
end
maxlong
minlong

% sort by latitude
T = struct2table(map);
%sortedT = sortrows(T, 'nlat'); 
sortedT = sortrows(T, [4 5]); 
mapsorted = table2struct(sortedT);

%  find how many longitude indices there are
long_index_max = 0;
for i = 1 : 1 : sel(size(mapsorted),1,1)
    if(mapsorted(i).nlat==mapsorted(1).nlat)
        long_index_max = long_index_max + 1;
    end
end
long_index_max;

%  find how many latitude indices there are
lat_index_max = 0;
for i = 1 : 1 : sel(size(mapsorted),1,1)
    if(mapsorted(i).nlong==mapsorted(1).nlong)
        lat_index_max = lat_index_max + 1;
    end
end
lat_index_max;



%%
% simple big map display
% 11th March 2021

% square is 700m wide, 550m high

bigimg = zeros(lat_index_max*739,long_index_max*933,3);
bigimg = im2uint8(bigimg);

%1D
% display as one big map (lel)
% k = 1;
% for i=0:long_index_max
%     bigimg(1:739,(1+(i*933)):((i*933)+933),1:3) = map(k).img(:,:,:);
%     k = k+1;
% end
% imshow(bigimg)

%2D
clf;
k = 1;
% for all lats
for i=lat_index_max-1:-1:0
    % print all longs
    for j=0:long_index_max-1
        %bigimg(1:739,(1+(j*933)):((j*933)+933),1:3) = mapsorted(k).img(:,:,:);
        %bigimg( (1+(i*739)):((i*739)+739) ,(1+(j*933)):((j*933)+933),1:3) = mapsorted(k).img(:,:,:);
        bigimg( (1+(i*739)):((i*739)+739) ,(1+(j*933)):((j*933)+933),1:3) = mapsorted(k).img(:,:,:);
        k = k+1;
        
    end
end
imshow(bigimg)

%%
% big map plot of features
% 12th March 2021

%clf;
k = 1;
m = 0;
% for all lats
for i=lat_index_max-1:-1:0
    % print all longs
    for j=0:long_index_max-1
        %bigimg(1:739,(1+(j*933)):((j*933)+933),1:3) = mapsorted(k).img(:,:,:);
        %bigimg( (1+(i*739)):((i*739)+739) ,(1+(j*933)):((j*933)+933),1:3) = mapsorted(k).img(:,:,:);
        %bigimg( (1+(i*739)):((i*739)+739) ,(1+(j*933)):((j*933)+933),1:3) = mapsorted(k).img(:,:,:);
        %
        locs = mapsorted(k).featuremap.Location;
        locs = round(locs);
        % convert from 
        % loop through all points so we can mark them
        % dont bother with scatter, wont work
        for(m = 1 : 1: sel(size(locs),1,1))
            
            % draw + sign
            marker = [255,0,0];
            bigimg( (1+(i*739)+locs(m,2)), (2+(j*933)+locs(m,1))-1, : ) = marker;
            bigimg( (1+(i*739)+locs(m,2)), (1+(j*933)+locs(m,1))-1, : ) = marker;
            bigimg( (1+(i*739)+locs(m,2)), (1+(j*933)+locs(m,1)), : ) = marker;
            bigimg( (1+(i*739)+locs(m,2)), (1+(j*933)+locs(m,1))+1, : ) = marker;
            bigimg( (1+(i*739)+locs(m,2)), (2+(j*933)+locs(m,1))+1, : ) = marker;
            bigimg( (1+(i*739)+locs(m,2))-2, (1+(j*933)+locs(m,1)), : ) = marker;
            bigimg( (1+(i*739)+locs(m,2))-1, (1+(j*933)+locs(m,1)), : ) = marker;
            bigimg( (1+(i*739)+locs(m,2)), (1+(j*933)+locs(m,1)), : ) = marker;
            bigimg( (1+(i*739)+locs(m,2))+1, (1+(j*933)+locs(m,1)), : ) = marker;
            bigimg( (1+(i*739)+locs(m,2))+2, (1+(j*933)+locs(m,1)), : ) = marker;
        end
        k = k+1;
    end
end
imshow(bigimg)
imwrite(bigimg,'map4_bigiimg.jpg','Quality',100);

%%
% plot all features on one grid using different colours
% allows comparison of how unique they are
% 12th March 2021

% loop through all map structure example
locimg = ones(739,933,3);
locimg = im2uint8(bigimg);


% generate colours
marker = 255*[rand(sel(size(map),1,2),1), rand(sel(size(map),1,2),1), rand(sel(size(map),1,2),1)];
marker = round(marker);

% plot
for k = 1 : sel(size(map),1,2)
    locs = mapsorted(k).featuremap.Location;
    locs = round(locs);
    
    for(m = 1 : 1: sel(size(locs),1,1))
        locimg( (locs(m,2)), locs(m,1), : ) = marker(k,:);
    end
    
    %imshow(map(1).img)
end
imshow(locimg)
imwrite(locimg,'map4_all_locs.jpg','Quality',100);

%%
% mean square error between all sets
% 12th March 2021

% sweep map

msemat = zeros(1, sel(size(mapsorted),1,1));

% loop through and compare first matrix with the others
for k = 2 : sel(size(mapsorted),1,1)
    msemat(k) = mse(mapsorted(1).featuremap.Location, mapsorted(k).featuremap.Location);
end
plot( 1:sel(size(mapsorted),1,1), msemat(:) )
ylabel('Mean Square Error')
xlabel('location(1) vs location(k)')
% there is very little correlation between feature locations of each 

%%
% try using matchfeatures
% 12th March 2021

mfeats = zeros(1, sel(size(mapsorted),1,1));

[f1 vpts1] = extractFeatures(rgb2gray(mapsorted(1).img), mapsorted(1).featuremap);
for k = 1 : sel(size(mapsorted),1,1)
    [f2 vpts2] = extractFeatures(rgb2gray(mapsorted(k).img), mapsorted(k).featuremap);
    % metric is sum of square differences
    [idx,met] = matchFeatures(f1, f2);
    mfeats(k) = sum(met(:))/sel(size(met),1,1);
end

plot( 1:sel(size(mapsorted),1,1), mfeats(:) )
ylabel('Average Feature Match Metric')
xlabel('location(1) vs location(k)')
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "0 is perfect match")


%%
% 3D surf of the correlaton of map grid features with each other
% 12th March 2021

scale = 7;
mfeats3d = zeros(sel(size(mapsorted),1,1)/scale, sel(size(mapsorted),1,1)/scale);

% find something that divides evenly. way too slow to run on every result
sel(size(mapsorted),1,1)/7

for i = 1 : 1 : sel(size(mapsorted),1,1)/scale
    for j = 1 : 1 : sel(size(mapsorted),1,1)/scale
        [f1 vpts1] = extractFeatures(rgb2gray(mapsorted(i).img), mapsorted(i).featuremap);
        [f2 vpts2] = extractFeatures(rgb2gray(mapsorted(j).img), mapsorted(j).featuremap);
        % metric is sum of square differences
        [idx,met] = matchFeatures(f1, f2);
        mfeats3d(i,j) = sum(met(:))/sel(size(met),1,1);
        if(mfeats3d(i,j) == NaN)
            mfeats3d(i,j) = -1;
        end
        if(mfeats3d(i,j) == inf)
            mfeats3d(i,j) = -1;
        end
    end 
end

size(mfeats3d)

surf( 1:sel(size(mfeats3d),1,1), 1:sel(size(mfeats3d),1,2), mfeats3d)
zlabel('Average Feature Match Metric')
xlabel('location(i)')
ylabel('location(j)')
annotation('textbox', [0.75, 0.2, 0.1, 0.1], 'String', "0 is perfect match")

%%
% test how closeness of matching works
% - sweep left to right along the map comparing to the 3rd image, see how the
% error metric varies as it gets closer and further away
% - consider rotating the array of points and then doing the above sweep
% - make a virtual map of Britain using 10 random feature points per km or
% so and then see how long different sized searches take
% - scale by altitude

% rotate test image
imgtestrgb = imread('maps/0.00_51.1.jpg');
% scale and rotate for test
imscale = 1;
imgtest = imresize(imrotate(imgtestrgb,-20),imscale);

% put the features into a big array and step through? how is best?

% step across the pixels in increments
for i = 1 : 10 : (933*10)
    %
end

%%
%%

% generate example map
% takes a fe minutes
mapidx = 1;
map = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)

for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        pts = rand(10,2);
        pts(:,1) = pts(:,1)*933;
        pts(:,2) = pts(:,2)*739;
        % map struct
        map(mapidx).index = mapidx;
        map(mapidx).nlat = n;
        map(mapidx).nlong = w;
        map(mapidx).featuremap = pts;
        mapidx = mapidx + 1
    end
end

%%

mapidx2 = 0;

for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        
    end
end


testsamples = 100;
results = zeros(testsamples,2);
for mapidx2 = 1:testsamples
    [idx, met] = matchFeatures(  map(50).featuremap,  map(mapidx2).featuremap   );

    results(mapidx2,2) = sum(sqrt(met.^2));
end

plot( 1:1:mapidx2, results(:,2))





%% 
% try comparing the first feature image with itself, but with some points
% missing, and the whole array shifted up or across



%%
%todo
% - plot within defined limits. i.e a subset of the map for bigger maps





%%

bigimg = zeros(739,5*933,3);
bigimg = im2uint8(bigimg);

% conversion from lat/long to pixels
lat2pix = (739/0.005)
long2pix = (933/0.01)
0.01* long2pix
