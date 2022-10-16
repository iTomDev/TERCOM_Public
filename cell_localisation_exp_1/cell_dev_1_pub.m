% concept:
% + midpoint between the towers. ignore signal strength as it can be very
% misleading. use 5 strongest max. good for course location for which area
% of town

% functions required
% getcellsinradiusofpoint: for test




cell_path = 'C:\Users\Tom\Documents\Work\_projects\2021_tercom_system\big_data\cell_towers\MLS-full-cell-export-2022-09-17T000000.CSV';
opts = detectImportOptions(cell_path);
cellm = readtable(cell_path);

% About the data
% It's the whole world.

% col 1
% be more efficient with cell station type
% GSM = 0, UTMS = 1, LTE = 2

% col 2: country code
% use country code to extract national data
% https://en.wikipedia.org/wiki/Mobile_country_code
% UK: 234, 235

% lat 8

% long 7

% cellid 5


% Extract the GB cell towers and save them to a CSV
% Extract 1 degree grids 10x10 from the GB cell and save them to CSVs

% 1: cellid
% 2: lat
% 3: long
% 4: station type
% 5: region code
cellgb = zeros(250000, 5);
exportpath = 'C:\Users\Tom\Documents\Work\_projects\2021_tercom_system\big_data\cell_towers\tower_by_region\twr';

% extract UK data
j = 1;
for i=1:size(cellm,1)
    % select country - by code
    if( cellm.mcc(i) == 234 || cellm.mcc(i) == 235)
        % select country - by lat/long
        if(cellm.lat(i) > 49 && cellm.lat(i) < 60 && cellm.lon(i) > -9 && cellm.lon(i) < 2)
            cellgb(j,1) = cellm.cell(i);
            cellgb(j,2) = cellm.lat(i);
            cellgb(j,3) = cellm.lon(i);
            
            % more efficient storage for cell radio type
            cellgb(j,4) = 0;
            if( strcmp(cellm.radio(i),'UTMS') )
                cellgb(j,4) = 1;
            end
            if( strcmp(cellm.radio(i),'lte') )
                cellgb(j,4) = 2;
            end
            
            % marker for which subregion in the UK it should go in
            cellgb(j,5) = latlong2regioncode(cellgb(j,2), cellgb(j,3));
            
            j = j + 1;
        end
    end
end
writematrix(cellgb, [exportpath 'gb.csv'] );


% % precompute an array of UK regtion codes. less typing :)
k = 1;
ukregioncodes = 0;
for i=-8:1:1
    for j=49:58
        ukregioncodes(1,k) = latlong2regioncode(j, i);
        ukregioncodes(2,k) = i; % long for verification
        ukregioncodes(3,k) = j; % lat for verification
        k = k + 1;
    end
end
ukregioncodes(1,:)'

% % east
% ukregioncodes = 0;
% k = 1;
% for i=0:1
%     for j=49:58
%         ukregioncodes(1,k) = latlong2regioncode(j, i);
%         ukregioncodes(2,k) = i; % long for verification
%         ukregioncodes(3,k) = j; % lat for verification
%         k = k + 1;
%     end
% end

% ukregioncodes = ...
% [95899; 95900; 95901; 95902; 95903; 95904; 95905; 95906; 95907; 95908; 96079; 96080; 96081; 96082; 96083; 96084; 96085;...
%  96086; 96087; 96088; 96259; 96260; 96261; 96262; 96263; 96264; 96265; 96266; 96267; 96268; 96439; 96440; 96441; 96442;...
%  96443; 96444; 96445; 96446; 96447; 96448; 96619; 96620; 96621; 96622; 96623; 96624; 96625; 96626; 96627; 96628; 96799;...
%  96800; 96801; 96802; 96803; 96804; 96805; 96806; 96807; 96808; 96979; 96980; 96981; 96982; 96983; 96984; 96985; 96986;...
%  96987; 96988; 97159; 97160; 97161; 97162; 97163; 97164; 97165; 97166; 97167; 97168; 32539; 32540; 32541; 32542; 32543;...
%  32544; 32545; 32546; 32547; 32548; 32719; 32720; 32721; 32722; 32723; 32724; 32725; 32726; 32727; 32728];

ukregioncodes = ...
       [31099;
       31100;
       31101;
       31102;
       31103;
       31104;
       31105;
       31106;
       31107;
       31108;
       31279;
       31280;
       31281;
       31282;
       31283;
       31284;
       31285;
       31286;
       31287;
       31288;
       31459;
       31460;
       31461;
       31462;
       31463;
       31464;
       31465;
       31466;
       31467;
       31468;
       31639;
       31640;
       31641;
       31642;
       31643;
       31644;
       31645;
       31646;
       31647;
       31648;
       31819;
       31820;
       31821;
       31822;
       31823;
       31824;
       31825;
       31826;
       31827;
       31828;
       31999;
       32000;
       32001;
       32002;
       32003;
       32004;
       32005;
       32006;
       32007;
       32008;
       32179;
       32180;
       32181;
       32182;
       32183;
       32184;
       32185;
       32186;
       32187;
       32188;
       32359;
       32360;
       32361;
       32362;
       32363;
       32364;
       32365;
       32366;
       32367;
       32368;
       32539;
       32540;
       32541;
       32542;
       32543;
       32544;
       32545;
       32546;
       32547;
       32548;
       32719;
       32720;
       32721;
       32722;
       32723;
       32724;
       32725;
       32726;
       32727;
       32728];

% split up by region
cellregion = zeros(1,4);

% loop through region codes
for regioni = 1:size(ukregioncodes,1)
    cellregion = zeros(1,4);
    regioni
    % loop through the UK tower list
    j = 1;
    for i=1:size(cellgb,1)
        % extract each degree to a file
        % lopop through the relevant degrees
        if(cellgb(i,5) == ukregioncodes(regioni))
            cellregion(j,1) = cellgb(i,1); % cell id
            cellregion(j,2) = cellgb(i,2); % lat
            cellregion(j,3) = cellgb(i,3); % long
            cellregion(j,4) = cellgb(i,5); % region code
            j = j + 1;     
        end
    end
    % exporT
    writematrix(cellregion, [exportpath num2str(ukregioncodes(regioni)) '.csv'] );
end

% plot one region cell towers
scatter( cellregion(:,3), cellregion(:,2) );

% plot all UK towers
scatter( cellgb(:,3), cellgb(:,2) );


% perform a search

% load file
importpath = 'C:\Users\Tom\Documents\Work\_projects\2021_tercom_system\big_data\cell_towers\tower_by_region\twr';
importfilepath = [importpath num2str(ukregioncodes(10)) '.csv'];
cellregion2 = readmatrix(importfilepath);

% find the last record to represent the worst case
lastrecord = cellregion2(size(cellregion2,1), :)

% search the worst case record
tic
for(i=1:size(cellregion2,1) )
    if(cellregion2(i,1)==lastrecord)
        disp('tower found')
    end
end
toc
% Elapsed time is 0.010227 seconds.


% search time - all GB

lastrecord = cellgb(size(cellgb,1), :)

% search the worst case record
tic
for(i=1:size(cellgb,1) )
    if(cellgb(i,1)==lastrecord)
        disp('tower found')
    end
end
toc
% Elapsed time is 0.057044 seconds


