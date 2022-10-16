
% Test to show how unique each number of feature points is.
% + uses random values not real map data for speed
% + 3604202 in each data set
% + 100k used from each set for computing norms etc
% + Takes a few minutes to run these things
% !!! You must have 32GB of RAM. It uses ~12B, after OS etc that is 21GB
% 14/3/2021


clc
clear
%%

% generate example map

mapidx = 1;

% mapidx size = 3604202

% 3 feature points per km, approx
mapidx = 1;
map3 = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)
for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        pts = rand(3,2);
        pts(:,1) = pts(:,1)*933;
        pts(:,2) = pts(:,2)*739;
        % map struct
        map3(mapidx).index = mapidx;
        map3(mapidx).nlat = n;
        map3(mapidx).nlong = w;
        map3(mapidx).featuremap = pts;
        mapidx = mapidx + 1;
    end
end

% 6 feature points per km, approx
mapidx = 1;
map6 = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)
for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        pts = rand(6,2);
        pts(:,1) = pts(:,1)*933;
        pts(:,2) = pts(:,2)*739;
        % map struct
        map6(mapidx).index = mapidx;
        map6(mapidx).nlat = n;
        map6(mapidx).nlong = w;
        map6(mapidx).featuremap = pts;
        mapidx = mapidx + 1;
    end
end

% 15 feature points per km, approx
mapidx = 1;
map15 = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)
for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        pts = rand(15,2);
        pts(:,1) = pts(:,1)*933;
        pts(:,2) = pts(:,2)*739;
        % map struct
        map15(mapidx).index = mapidx;
        map15(mapidx).nlat = n;
        map15(mapidx).nlong = w;
        map15(mapidx).featuremap = pts;
        mapidx = mapidx + 1;
    end
end

% 30 feature points per km, approx
mapidx = 1;
map30 = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)

for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        pts = rand(15,2);
        pts(:,1) = pts(:,1)*933;
        pts(:,2) = pts(:,2)*739;
        % map struct
        map30(mapidx).index = mapidx;
        map30(mapidx).nlat = n;
        map30(mapidx).nlong = w;
        map30(mapidx).featuremap = pts;
        mapidx = mapidx + 1;
    end
end


% 60 feature points per km, approx
mapidx = 1;
map60 = struct('index',0,'lat',0,'long',0,'nlat',0,'nlong',0,'img',0,'featuremap',0)

for n=0 : 0.005 : (300*0.05) % lat
    for w = 0 : 0.01 : (600*0.02) % long
        % feature points
        pts = rand(60,2);
        pts(:,1) = pts(:,1)*933;
        pts(:,2) = pts(:,2)*739;
        % map struct
        map60(mapidx).index = mapidx;
        map60(mapidx).nlat = n;
        map60(mapidx).nlong = w;
        map60(mapidx).featuremap = pts;
        mapidx = mapidx + 1;
    end
end


%%

% plot 2 norm of the feature error
testsamples = 100000;
referencesample = testsamples/2;
results3 = zeros(testsamples,2);
results6 = zeros(testsamples,2);
results15 = zeros(testsamples,2);
results30 = zeros(testsamples,2);
results60 = zeros(testsamples,2);
for mapidx2 = 1:testsamples
    
    [idx, met] = matchFeatures(  map3(referencesample).featuremap,  map3(mapidx2).featuremap   );
    results3(mapidx2,2) = norm(met(:),2);
    results3(mapidx2,1) = (sum(met.^2));
    [idx, met] = matchFeatures(  map6(referencesample).featuremap,  map6(mapidx2).featuremap   );
    results6(mapidx2,2) = norm(met(:),2);
    results6(mapidx2,1) = (sum(met.^2));
    
    [idx, met] = matchFeatures(  map15(referencesample).featuremap,  map15(mapidx2).featuremap   );
    results15(mapidx2,2) = norm(met(:),2);
    results15(mapidx2,1) = (sum(met.^2));
    [idx, met] = matchFeatures(  map30(referencesample).featuremap,  map30(mapidx2).featuremap   );
    results30(mapidx2,2) = norm(met(:),2);
    results30(mapidx2,1) = (sum(met.^2));
    [idx, met] = matchFeatures(  map60(referencesample).featuremap,  map60(mapidx2).featuremap   );
    results60(mapidx2,2) = norm(met(:),2);
    results60(mapidx2,1) = (sum(met.^2));
    %results(mapidx2,2) = sqrt(sum(met.^2)); % verify
end

% full plot sum of squares
% larger indicates uniqueness of results

% 15, 30 and 60 features
% more smaller results indicates locations too similar to the true location
plot( 1:1:mapidx2, results15(:,1),'green', 1:1:mapidx2, results30(:,1),'blue', 1:1:mapidx2, results60(:,1),'red')
set(gca,'color','white');
savefig('graph_generation/experiment_15_feature_point_uniqueness.fig');
xlabel('Location (50k) vs Location(1:100k)')
ylabel('Sum of square of errors')
legend('15 Features','30 Features','60 Features')

% full plot 2-norm
plot( 1:1:mapidx2, results15(:,2),'green', 1:1:mapidx2, results30(:,2),'blue', 1:1:mapidx2, results60(:,2),'red')
set(gca,'color','white');
title('Error Plot for 100K Samples - Full')
savefig('graph_generation/pts_uniqueness_line_features_15_30_60_no_zoom.fig');
xlabel('Location (50k) vs Location(1:100k)')
ylabel('2-Norm of errors')
legend('15 Features','30 Features','60 Features')

% only the lowest values, 2-norm
plot( 1:1:mapidx2, results15(:,2),'green', 1:1:mapidx2, results30(:,2),'blue', 1:1:mapidx2, results60(:,2),'red')
set(gca,'color','white');
savefig('graph_generation/pts_uniqueness_line_features_15_30_60_with_zoom.fig');
xlabel('Location (50k) vs Location(1:100k)')
ylabel('2-Norm of errors')
ylim([0 0.001])
legend('15 Features','30 Features','60 Features')


% 3, 6 and 15 features

% full plot 2-norm
plot( 1:1:mapidx2, results3(:,2),'green', 1:1:mapidx2, results6(:,2),'blue', 1:1:mapidx2, results15(:,2),'red')
set(gca,'color','white');
title('Error Plot for 100K Samples - Full')
savefig('graph_generation/pts_uniqueness_line_features_3_6_15_no_zoom.fig');
xlabel('Location (50k) vs Location(1:100k)')
ylabel('2-Norm of errors')
legend('3 Features','6 Features','15 Features')

% only the lowest values, 2-norm
plot( 1:1:mapidx2, results3(:,2),'green', 1:1:mapidx2, results6(:,2),'blue', 1:1:mapidx2, results15(:,2),'red')
set(gca,'color','white');
savefig('graph_generation/pts_uniqueness_line_features_3_6_15_with_zoom.fig');
xlabel('Location (50k) vs Location(1:100k)')
ylabel('2-Norm of errors')
ylim([0 0.001])
legend('3 Features','6 Features','15 Features')


%%
% full historgram
h = histogram(results15(:,2),1000)
h.BinLimits = [0 0.001];
histogram(results3(:,2),100)
histogram(results60(:,2),100)

clf;

% 3 features 
hf1 = histogram(results3(:,2),10000)
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 3 Features - No Zoom')
savefig('graph_generation/pts_uniqueness_histo_features_3_100k_samples_no_zoom.fig');

% 30 features 
hf4 = histogram(results30(:,2),10000)
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 30 Features - No Zoom')
savefig('graph_generation/pts_uniqueness_histo_features_30_100k_samples_no_zoom.fig');

% 60 features -  1 match, nothing else even close
hf5 = histogram(results60(:,2),10000)
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 60 Features - No Zoom')
savefig('graph_generation/pts_uniqueness_histo_features_60_100k_samples_no_zoom.fig')


% focus on the closest matches
% these are the most useful

% 3 = 22130 matches with 0 error
%subplot(1,2,1)
h1 = histogram(results3(:,2),10000)
h1.BinLimits = [0 0.001];
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 3 Features')
savefig('graph_generation/pts_uniqueness_histo_features_3_100k_samples.fig');

% 6 = 1601 matches with 0 error
h2 = histogram(results6(:,2),10000)
h2.BinLimits = [0 0.001];
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 6 Features')
savefig('graph_generation/pts_uniqueness_histo_features_6_100k_samples.fig');

% 15 features - 1 match, fair few close
h3 = histogram(results15(:,2),10000)
h3.BinLimits = [0 0.001];
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 15 Features')
savefig('graph_generation/pts_uniqueness_histo_features_15_100k_samples.fig');

% 30 features = 1 match
h4 = histogram(results30(:,2),10000)
h4.BinLimits = [0 0.001];
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 30 Features')
savefig('graph_generation/pts_uniqueness_histo_features_30_100k_samples.fig');

% 60 features -  1 match, nothing else even close
h5 = histogram(results60(:,2),10000)
h5.BinLimits = [0 0.001];
ylabel('Matches')
xlabel('2-Norm of Error Metric')
title('Location Uniqueness: 100K Samples, 60 Features')
savefig('graph_generation/pts_uniqueness_histo_features_60_100k_samples.fig');





%%
% full map

% plot 2 norm of the feature error
testsamples = 1000000; %sel(size(map60),1,2);
referencesample = round(testsamples/2); % 1802100;
results15 = zeros(testsamples,2);
results30 = zeros(testsamples,2);
results60 = zeros(testsamples,2);
for mapidx2 = 1:testsamples
    [idx, met] = matchFeatures(  map15(referencesample).featuremap,  map15(mapidx2).featuremap   );
    results15(mapidx2,2) = norm(met(:),2);
    results15(mapidx2,1) = (sum(met.^2));
    [idx, met] = matchFeatures(  map30(referencesample).featuremap,  map30(mapidx2).featuremap   );
    results30(mapidx2,2) = norm(met(:),2);
    results30(mapidx2,1) = (sum(met.^2));
    [idx, met] = matchFeatures(  map60(referencesample).featuremap,  map60(mapidx2).featuremap   );
    results60(mapidx2,2) = norm(met(:),2);
    results60(mapidx2,1) = (sum(met.^2));
end

% only the lowest values, 2-norm
plot( 1:1:mapidx2, results15(:,2),'green', 1:1:mapidx2, results30(:,2),'blue', 1:1:mapidx2, results60(:,2),'red')
set(gca,'color','white');
savefig('graph_generation/experiment_15_feature_point_uniqueness.fig');
xlabel('Location (500) vs Location(1:1000)')
ylabel('2-norm of errors')
ylim([0 0.001])
legend('15 Features','30 Features','60 Features')

% historgrams

% todo: 3 and 5 features 

% 15 features or less could be used to remove very-unlikely matches
subplot(1,2,1)
title('15 Features')
h1 = histogram(results15(:,2),10000)
h1.BinLimits = [0 0.001];

% 60 features might be enough to be unique
subplot(1,2,2)
title('60 Features')
h2 = histogram(results60(:,2),10000)
h2.BinLimits = [0 0.001];

%%
% sample data

for ii=1:5
    map15(ii).featuremap(:)
end

% map3 - 3 features
%   760.1372 % 1
%   845.1039
%   118.4787 % 2
%   674.9848
%   467.3135 % 3
%    72.0824
% 
%   259.8388
%   510.2405
%   893.3539
%   713.0526
%   116.4761
%   717.2681
% 
%   893.0368
%   452.8555
%   746.6617
%   104.8540
%   311.6816
%   676.7286
% 
%   739.1294
%   895.2064
%   611.8061
%    26.3909
%   627.5066
%   690.2210
% 
%   633.2599
%   706.9715
%   693.3426
%   289.8558
%   484.3982
%   126.5070
% 
%   % map15 - 15 features
% 
%   400.3932
%   597.4340
%   506.2681
%    76.9153
%   704.1392
%   620.9369
%   859.9423
%   381.9821
%   307.7249
%   425.9180
%   468.8137
%   218.5893
%   215.8464
%   804.1619
%   436.0129
%   738.2680
%   351.7018
%    54.3799
%   406.7368
%   726.4266
%   135.3884
%    18.1652
%   131.9972
%   580.4517
%   187.7421
%   145.5334
%   303.3449
%   153.8348
%   202.0969
%   176.0883
% 
%   779.3937
%   320.1401
%   816.5107
%   547.8675
%   782.2185
%   447.7576
%   344.9757
%   246.2231
%   623.0876
%   166.3668
%   641.5133
%   276.0672
%   707.5678
%   763.6094
%   616.0721
%   608.6774
%   413.3711
%   374.4537
%   504.5772
%   383.6884
%    62.7665
%   500.4776
%    63.1627
%   588.2661
%   520.8583
%   483.7970
%    41.7948
%   292.9061
%     5.3690
%   382.0480
% 
%   855.6338
%   388.5600
%   522.9609
%   859.9596
%   388.8830
%   564.3608
%   244.3401
%   580.2516
%   282.4291
%   790.8489
%   291.6185
%   474.2216
%   101.8159
%   263.6568
%   343.8329
%   527.1184
%   696.3051
%   687.9775
%   113.9324
%   126.0009
%   349.1395
%   125.8947
%   270.1112
%   441.4465
%   121.7906
%   126.1764
%   691.3984
%   714.0867
%   183.8273
%   558.4684
%   
%   608.8478
%   803.8363
%    20.4493
%   545.1815
%   368.8285
%   141.4242
%   642.6936
%   504.7116
%   294.0922
%   586.9015
%   145.9477
%   770.3386
%   455.7438
%   453.4900
%   197.1910
%     6.4918
%   300.4697
%   420.7376
%    49.8327
%   684.8714
%   437.0950
%   536.7314
%   532.1641
%    29.2115
%   557.5888
%   192.0754
%   592.9264
%   207.5610
%   622.9568
%   337.1425
% 
%   348.1383
%   917.2748
%   474.5515
%   344.1191
%   454.7431
%    83.4693
%   102.2893
%   503.1890
%   293.5525
%     9.2366
%   759.8418
%   815.2374
%   282.6793
%   851.2556
%   591.3767
%   328.1514
%   147.1738
%   222.9589
%    98.5438
%   374.1707
%   735.0475
%   682.9917
%   615.0008
%   505.9672
%   423.2519
%   421.8996
%   695.9929
%    17.6135
%   722.8557
%   344.1671
% 



