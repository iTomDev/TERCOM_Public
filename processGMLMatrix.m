function [outputArg1,outputArg2] = processGMLMatrix(inputArg1,inputArg2)
%PROCESSGMLMATRIX Summary of this function goes here
%   Detailed explanation goes here
% It's more complicated than I was intending.
% Thomas Pile, Feb 2022

% - - - - - - - - - - - - -  config flags

% create a RoadLinkGeom matrix with all geometry in WGS84, not just the
% summary (first, mid, last) in the RoadLink matrix. probably slow!
CONFIG_DO_FULL_GEOM = 0;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


WGS_INT64_SCALE = 10000;


% read GML data from folder
%[RNRaw, RLRaw] = parseGMLFolder(dirpath)

% load files
tmpload = load('matlab\big_data\RLRaw.mat');
RLRaw = tmpload.RLRawA;
tmpload = load('matlab\big_data\RNRaw.mat');
RNRaw = tmpload.RNRawA;
clear tmpload;

% find the real size of the RodeNode and RoadLink
szRNRaw = 0;
for i = 1:size(RNRaw,1)
    if ( strcmp(RNRaw(i,1),'') == 1 ) % 
        % found size
        break;
    else
        szRNRaw = szRNRaw + 1;
    end
end
szRLRaw = 0;
for i = 1:size(RLRaw,1)
    if ( strcmp(RLRaw(i,1),'') == 1 ) % 
        % found size
        break;
    else
        szRLRaw = szRLRaw + 1;
    end
end


% initialise number matrices  
% matrix format
% 1:3   ID-[1:3]
% 4:6   endID-[1:3]
% 7:9   startID-[1:3]
% 10    length
% 11    heading
% 12    classification
% 13:14    coordmin [x,y]
% 15:16    coordmid [x,y]
% 17:18    coordmax [x,y]
% 19        sourceFileID
% 20        road number
% 21        road number letter component e.g A, B
RoadLink = int64(zeros(size(RLRaw,1),21));
% 1:3    ID-[1:3]
% 4:5    coord [x,y]
% 6      type of node
% 7      sourceFileID
RoadNode = int64(zeros(size(RNRaw,1),7));
% 1:100     geom [x,y]
% 101:103   link ID       
RoadLinkGeom = int64(zeros(szRLRaw, 103));

% 
% conv1 = (RNRaw(1,2));
% upper = str2num(extractBetween(conv1,1,"."));
% %lower = str2num(extractBetween(conv1,".",strlength(conv1)));
% result = (upper*POSITION_SCALE)+lower;

% convert from BGN to OSTN
% eas = 328243.0;
% nor = 185698.15;
% [a,b] = OSTN15_Matlab('grid-to-gps',eas,nor)
% 
% 
% eas = 463772;
% nor = 1213537;
% [a,b] = OSTN15_Matlab('grid-to-gps',eas,nor)

% - - - - - - Node - - - - - - 

% convert matrix values from strings to numbers - Node
for i = 1: szRNRaw
    % ----- ID -----
    % 8E554699-810F-4A20-8017-E0FA8D9540E9
    temps = convertStringsToChars(RNRaw(i,1));
    % extract first section
    % 8E554699 = 2387953305
    val1 = temps([1:8]);
    RoadNode(i,1) = int64(hex2uint64( val1 )); 
    % second section
    % 810F-4A20-8017
    val1 = temps([10:13 15:18 20:23]);
    RoadNode(i,2) = int64(hex2uint64( val1 ));
    % third section
    % E0FA8D9540E9 = 247366721814761
    val1 = temps([25:36]);
    RoadNode(i,3) = int64(hex2uint64( val1 ));
    
    % ----- geometry pairs (North/East) -----
    % grab the nor, eas position and convert
    temps = convertStringsToChars(RNRaw(i,2));
    tempn = int64(zeros(1,2));
    start = 1;
    for j=1:strlength(temps)
        % if next number, delimited by ' '
        if (strcmp(' ',temps([j]) ) == 1)
            % extract the number and convert to int
            val = temps(start:j-1);
            % drop the fractional part
            val2 = extractBetween(val,1,'.');
            if( strcmp(val2{1,1}, '') )
                val2 = val;
            end
            n = sscanf(string(val2),'%lu');
            %n2 = double(n);
            % dont add the first value
            
            % they need to be split into pairs
            tempn( 1, 1) ... 
                = int64(n); 
            start = j+1;
        end
    end % position
    % grab the final value
    val = (temps(start:strlength(temps)));
    val2 = extractBetween(val,1,'.');
    if( strcmp(val2{1,1}, '') )
        val2 = val;
    end
    n = sscanf(string(val2),'%lu');
    %n2 = double(n);
    %tempn(mod(paircount,2)+1, 10) = n2;
    tempn(1,2) = int64(n); 
    % finally halve paircount as it's double what it should be!
    
    % also convert to WGS84
    [wgslat, wgslong] = OSTN15_Matlab('grid-to-gps',double(tempn(1,1)) ,double(tempn(1,2)) );
    wgslat = int64(wgslat*WGS_INT64_SCALE);     % scale up for int64
    wgslong = int64(wgslong*WGS_INT64_SCALE);
    RoadNode(i,4) = wgslat; RoadNode(i,5) = wgslong;
    
    % ----- Node Classification -----
%     1	junction 	eu-legal 	Valid
%     2	roundabout 	eu-legal 	Valid
%     3	road end 	eu-legal 	Valid
%     4	enclosed traffic area 	eu-legal 	Valid
%     5	level crossing 	eu-legal 	Valid
%     6	pseudo node 	eu-legal 	Valid
%     7	road service area 	eu-legal 	Valid
%     8	traffic square
    if( strcmp(RNRaw(i,3), '>junction')==1 )
        RoadNode(i,6) = int64(1);
    end
    if( strcmp(RNRaw(i,3), '>roundabout')==1 )
        RoadNode(i,6) = int64(2);
    end
    if( strcmp(RNRaw(i,3), '>road end')==1 )
        RoadNode(i,6) = int64(3);
    end
    if( strcmp(RNRaw(i,3), '>pseudo node')==1 )
        RoadNode(i,6) = int64(6);
    end
    
    % source data classification (region)
    % HP HT HU HY HZ    1:5
    % NA NB NC ND NF    6:10
    % NG NH NJ NK NL    11:15
    % NM NN NO NR NS    16:20
    % NT NU NW NX NY    21:25
    % NZ SD SE SH SJ    26:30
    % SK SM SN SO SP    31:35
    % SR SS ST SU SV    36:40
    % SW SX SY SZ TA    41:45
    % TF TG TL TM TQ    46:50
    % TR TV             51,52
    % Look up the string in the array. The index is the new value
    DataSourceIDLookup = ['HP', 'HT', 'HU', 'HY', 'HZ',...
                          'NA', 'NB', 'NC', 'ND', 'NF',...
                          'NG', 'NH', 'NJ', 'NK', 'NL',...
                          'NM', 'NN', 'NO', 'NR', 'NS',...
                          'NT', 'NU', 'NW', 'NX', 'NY',...
                          'NZ', 'SD', 'SE', 'SH', 'SJ',...
                          'SK', 'SM', 'SN', 'SO', 'SP',... 
                          'SR', 'SS', 'ST', 'SU', 'SV',... 
                          'SW', 'SX', 'SY', 'SZ', 'TA',...
                          'TF', 'TG', 'TL', 'TM', 'TQ',...
                          'TR', 'TV'];
%     DataSourceIDLookup = [HP, HT, HU, HY, HZ,...
%                           NA, NB, NC, ND, NF,...
%                           NG, NH, NJ, NK, NL,...
%                           NM, NN, NO, NR, NS,...
%                           NT, NU, NW, NX, NY,...
%                           NZ, SD, SE, SH, SJ,...
%                           SK, SM, SN, SO, SP,... 
%                           SR, SS, ST, SU, SV,... 
%                           SW, SX, SY, SZ, TA,...
%                           TF, TG, TL, TM, TQ,...
%                           TR, TV];                  
    % loop
    for j=1:size(DataSourceIDLookup)
        if( strcmp(RNRaw(i,3), [DataSourceIDLookup( ((j*2)-1):(j*2) )] )==1 )
            RoadNode(i,7) = int64(j);
        end
    end
end

save('matlab\big_data\RoadNode.mat','RoadNode','-v7.3')


% - - - - - - Link - - - - - - 

% convert matrix values from strings to numbers - Link
% loop through the records 
for i = 1:szRLRaw
    
    % ----- ID -----
    % 8E554699-810F-4A20-8017-E0FA8D9540E9
    temps = convertStringsToChars(RLRaw(i,1));
    % extract first section
    % 8E554699 = 2387953305
    val1 = temps([1:8]);
    val2 = int64(hex2uint64( val1 )); 
    RoadLink(i,1) = val2;
    % second section
    % 810F-4A20-8017
    val1 = temps([10:13 15:18 20:23]);
    val2 = int64(hex2uint64( val1 ));
    RoadLink(i,2) = val2;
    % third section
    % E0FA8D9540E9 = 247366721814761
    val1 = temps([25:36]);
    val2 = int64(hex2uint64( val1 ));
    RoadLink(i,3) = val2;
    
    % ----- Start ID -----
    temps = convertStringsToChars(RLRaw(i,4));
    % extract first section
    % 8E554699 = 2387953305
    val1 = temps([1:8]);
    val2 = int64(hex2uint64( val1 )); 
    RoadLink(i,4) = val2;
    % second section
    % 810F-4A20-8017
    val1 = temps([10:13 15:18 20:23]);
    val2 = int64(hex2uint64( val1 ));
    RoadLink(i,5) = val2;
    % third section
    % E0FA8D9540E9 = 247366721814761
    val1 = temps([25:36]);
    val2 = int64(hex2uint64( val1 ));
    RoadLink(i,6) = val2;
    
    % ----- End ID -----
    temps = convertStringsToChars(RLRaw(i,3));
    % extract first section
    % 8E554699 = 2387953305
    val1 = temps([1:8]);
    val2 = int64(hex2uint64( val1 )); 
    RoadLink(i,7) = val2;
    % second section
    % 810F-4A20-8017
    val1 = temps([10:13 15:18 20:23]);
    val2 = int64(hex2uint64( val1 ));
    RoadLink(i,8) = val2;
    % third section
    % E0FA8D9540E9 = 247366721814761
    val1 = temps([25:36]);
    val2 = int64(hex2uint64( val1 ));
    RoadLink(i,9) = val2;
    
    % ----- Length -----
    % handle zero length road. its possible how?!?!
    if ( strcmp(RLRaw(i,6),'') == 1 )
        RoadLink(i,10) = 0 ;
    else
        RoadLink(i,10) = int64(str2num( RLRaw(i,6) ));
    end
    
    % ----- geometry pairs (North/East) -----
    % extract the string to a temp vector
    temps = convertStringsToChars(RLRaw(i,2));
    % temp vector for the numbers
    tempn = int64(zeros(500,2));
    % extract string of lat/long points
    start = 1;
    paircount = 0;
    for j=1:strlength(temps)
        % if next number, delimited by ' '
        if (strcmp(' ',temps([j]) ) == 1)
            % extract the number and convert to int
            val = temps(start:j-1);
            % drop the fractional part
            val2 = extractBetween(val,1,'.');
            if( strcmp(val2{1,1}, '') )
                val2 = val;
            end
            n = sscanf(string(val2),'%lu');
            %n2 = double(n);
            % dont add the first value
            if (paircount > 1)
                % they need to be split into pairs
                tempn( round(((paircount)/2)-0.001), ... % new row every other time, -0.001 rounds 0.5 down instead
                    mod(round(paircount),2)+1) ... 
                    = int64(n); 
            end
            start = j+1;
            paircount = paircount + 1;
        end
    end
    % grab the final value
    val = (temps(start:strlength(temps)));
    val2 = extractBetween(val,1,'.');
    if( strcmp(val2{1,1}, '') )
        val2 = val;
    end
    n = sscanf(string(val2),'%lu');
    %n2 = double(n);
    %tempn(mod(paircount,2)+1, 10) = n2;
    tempn( round(((paircount)/2)-0.001), mod(round(paircount),2)+1) = int64(n); 
    % finally halve paircount as it's double what it should be!
    paircount = round(round((paircount-1) / 2));
    
    % process the array of WGS coords into the final array
    % why? sometimes there are 100+ values, we cant have an array that big
    % also convert to WGS84
    [wgslat, wgslong] = OSTN15_Matlab('grid-to-gps',double(tempn(1,1)) ,double(tempn(1,2)) );
    wgslat = int64(wgslat*WGS_INT64_SCALE);     % scale up for int64
    wgslong = int64(wgslong*WGS_INT64_SCALE);
    RoadLink(i,13) = wgslat; RoadLink(i,14) = wgslong;
    
    % midpoint
    [wgslat, wgslong] = OSTN15_Matlab('grid-to-gps',double(tempn(round(paircount/2),1)) ,double(tempn(round(paircount/2),2)) );
    wgslat = int64(wgslat*WGS_INT64_SCALE);     % scale up for int64
    wgslong = int64(wgslong*WGS_INT64_SCALE);
    RoadLink(i,15) = wgslat; RoadLink(i,16) = wgslong;
    
    % end
    [wgslat, wgslong] = OSTN15_Matlab('grid-to-gps',double(tempn(paircount,1)) ,double(tempn(paircount,2)) );
    wgslat = int64(wgslat*WGS_INT64_SCALE);     % scale up for int64
    wgslong = int64(wgslong*WGS_INT64_SCALE);
    RoadLink(i,17) = wgslat; RoadLink(i,18) = wgslong;
    
    % full geom conversion to WGS-84
    if (CONFIG_DO_FULL_GEOM == 1)
        % if less than 50, just put them all in
        if (paircount < 50)
            for j=1:paircount % through the nor/est cols
                [wgslat, wgslong] = OSTN15_Matlab('grid-to-gps',double(tempn(j,1)) ,double(tempn(j,2)) );
                wgslat = int64(wgslat*WGS_INT64_SCALE);     % scale up for int64
                wgslong = int64(wgslong*WGS_INT64_SCALE);
                RoadLinkGeom(i,2*j) = wgslat; RoadLink(i,(2*j)+1) = wgslong;
                % todo: add link ID
            end
        % if more than 50, sort a range
        else
            % gen 50 numbers between 1 and paircount, then make ascending
            % easiest scalable way to select the sample points
            samples = sort(randperm(paircount,50));
            % work through the sample points
            for j=1:50 
                [wgslat, wgslong] = OSTN15_Matlab('grid-to-gps',double(tempn(samples(j),1)) ,double(tempn(samples(j),2)) );
                wgslat = int64(wgslat*WGS_INT64_SCALE);     % scale up for int64
                wgslong = int64(wgslong*WGS_INT64_SCALE);
                RoadLinkGeom(i,2*j) = wgslat; RoadLink(i,(2*j)+1) = wgslong;
                % LinkID
                % NOTE: geom and link indexes are the same! You can use this 
                RoadLinkGeom(i,101) = RoadLink(i,1);
                RoadLinkGeom(i,101) = RoadLink(i,1);
                RoadLinkGeom(i,101) = RoadLink(i,1);         
            end
        end
    end
    
    % calculate heading
    % working clockwise results would be... (*57.2958)
    % atan2(y,x) = (2,2) = 45
    % atan2(y,x) = (-2,2) = -45
    % atan2(y,x) = (-2,-2) = -135
    % atan2(y,x) = (2,-2) = 135
    heading = atan2( double(tempn(paircount,2)-tempn(1,2)), double(int64(tempn(paircount,1))-int64(tempn(1,1))) ) *57.2958;
    RoadLink(i,11) = int64(heading);
    
    % Classification
%     1	Motorway	
%     2	A Road	
%     3	B Road	
%     4	Classified Unnumbered	
%     5	Unclassified	
%     6	Not Classified	
%     7	Unknown	
    if( strcmp(RLRaw(i,5), 'Motorway')==1 )
        RoadLink(i,12) = int64(1);
    end
    if( strcmp(RLRaw(i,5), 'A Road')==1 )
        RoadLink(i,12) = int64(2);
    end
    if( strcmp(RLRaw(i,5), 'B Road')==1 )
        RoadLink(i,12) = int64(3);
    end
    if( strcmp(RLRaw(i,5), 'Classified Unnumbered')==1 )
        RoadLink(i,12) = int64(4);
    end
    if( strcmp(RLRaw(i,5), 'Unclassified')==1 )
        RoadLink(i,12) = int64(5);
    end
    if( strcmp(RLRaw(i,5), 'Not Classified')==1 )
        RoadLink(i,12) = int64(6);
    end
    if( strcmp(RLRaw(i,5), 'Unknown')==1 )
        RoadLink(i,12) = int64(7);
    end
    
    % 19, source file ID
    DataSourceIDLookup = ['HP', 'HT', 'HU', 'HY', 'HZ',...
                          'NA', 'NB', 'NC', 'ND', 'NF',...
                          'NG', 'NH', 'NJ', 'NK', 'NL',...
                          'NM', 'NN', 'NO', 'NR', 'NS',...
                          'NT', 'NU', 'NW', 'NX', 'NY',...
                          'NZ', 'SD', 'SE', 'SH', 'SJ',...
                          'SK', 'SM', 'SN', 'SO', 'SP',... 
                          'SR', 'SS', 'ST', 'SU', 'SV',... 
                          'SW', 'SX', 'SY', 'SZ', 'TA',...
                          'TF', 'TG', 'TL', 'TM', 'TQ',...
                          'TR', 'TV'];           
    % loop
    for j=1:size(DataSourceIDLookup)
        if( strcmp(RLRaw(i,8), [DataSourceIDLookup( ((j*2)-1):(j*2) )] )==1 )
            RoadLink(i,19) = int64(j);
        end
    end
    
    % road number - 20/21st field
    if strcmp(RLRaw(i,7),'')
        % empty field
        RoadLink(i,20) = int64(0);
        RoadLink(i,21) = int64(0);
    else
        % first letter 
        % A4, M4, B1234, A1(M)
        roadname = char(RLRaw(i,7));
        % sort by letter
        if strcmp(roadname(1),'M')
            RoadLink(i,21) = int64(1);
            RoadLink(i,20) = str2num( roadname(2:numel(roadname)) );
        end
        if strcmp(roadname(1),'A')
            RoadLink(i,21) = int64(2);
            % if it's an A1(M) type situation
            if isempty(str2num( roadname(2:numel(roadname)) ))
                % iterate through numerical bit only
                for k=1:numel(roadname)
                    % end of numerical bit
                    if( strcmp(roadname(k),'(') )
                        % substring out the numerical part
                        RoadLink(i,20) = str2num( roadname(2:(k-1)) );
                    end
                end 
            else
                % not an A1(M) edge case
                RoadLink(i,20) = str2num( roadname(2:numel(roadname)) );
            end
        end % 'A'
        if strcmp(roadname(1),'B')
            RoadLink(i,21) = int64(3);
            RoadLink(i,20) = str2num( roadname(2:numel(roadname)) );
        end
    end
    i
end

% save
save('matlab\big_data\RoadLink.mat','RoadLink','-v7.3')
if (CONFIG_DO_FULL_GEOM == 1)
    save('matlab\big_data\RoadLinkGeom.mat','RoadLinkGeom','-v7.3')
end

%%
% strcuture of new reissue ID mnatrices i.e RoadLink2 and RoadNode2

% RoadNode2
% 1			ID
% 2:3		coord [x,y]
% 4			type of node
% 5			highest classification of connected link
RoadNode2 = int64(zeros(4096000,5));

% RoadLink2
% 1			ID
% 2			endID
% 3			startID
% 4         roadID - for joining links to form a road
% 5			length
% 6			heading
% 7			classification
% 8:9		coordmin [x,y]
% 10:11		coordmid [x,y]
% 12:13		coordmax [x,y]
RoadLink2 = int64(zeros(4096000,13)); 

% reindex nodes
newidx = 1;
% for all nodes
for i=1:size(RoadNode,1)
	% reissue ID
	RoadNode2(i,1) = newidx;
	newidx = newidx + 1;
    % copy existing data across
    RoadNode2(i,2) = RoadNode(i,4); % lat
    RoadNode2(i,3) = RoadNode(i,5); % long
    RoadNode2(i,4) = RoadNode(i,6); % type of node
end
RoadNodeIndexMax = newidx-1;

% for all links
for i=1:size(RoadLink,1)
	% reissue link ids
	RoadLink2(i,1) = newidx;
	newidx = newidx + 1;
    
    % substitute the old style start and end indexes for the new style
    % loop through RoadNodes
    for j=1:size(RoadNode2,1)
        % look up start and end node Ids and substitute in the new issued IDs
        % end ID matches
        if(RoadLink(i,4) == RoadNode(i,1) && RoadLink(i,4) == RoadNode(i,1))
            RoadLink2(i,2) = RoadNode2(j,1);
            break;
        end
        % start ID matches
        if(RoadLink(i,7) == RoadNode(i,1) && RoadLink(i,8) == RoadNode(i,1))
            RoadLink2(i,3) = RoadNode2(j,1);
            break;
        end
    end
    
    % copy over the remaining data
    RoadLink2(i,5) = RoadLink(i,10); % length
    RoadLink2(i,6) = RoadLink(i,11); % heading 
    RoadLink2(i,7) = RoadLink(i,12); % classification
    RoadLink2(i,8) = RoadLink(i,13); % coordmin [x,y]
    RoadLink2(i,9) = RoadLink(i,14);
    RoadLink2(i,10) = RoadLink(i,15); % coordmid [x,y]
    RoadLink2(i,11) = RoadLink(i,16);
    RoadLink2(i,12) = RoadLink(i,17); % coordmax [x,y]
    RoadLink2(i,13) = RoadLink(i,18);
    i
end

% find the highest link classification for each node
% for all links
for i=1:size(RoadLink2,1)
	% if start link is a node
    if(RoadLink2(i,2) < RoadNodeIndexMax)
        % if that nodes classification level is less than this one then
        % increase it. use the startID to index in RoadNode2
        if(RoadLink2(i,7) > RoadNode2(RoadLink2(i,2),5))
            RoadNode2(RoadLink2(i,2),5) = RoadLink2(i,7);
        end
    end
end




%%
% Experiment: Extract only the motorway links

% tmpload = load('matlab\big_data\RoadLink.mat');
% RoadLink = tmpload.RoadLink;
% tmpload = load('matlab\big_data\RoadNode.mat');
% RoadNode = tmpload.RoadNode;
% MotorwayNodeLUT = zeros(1000,2);
% clear tmpload;
% 
% % generate a tree structure of la/long so regions can be searched and
% % plotted
% % 1:3   ID-[1:3]
% % 4:6   endID-[1:3]
% % 7:9   startID-[1:3]
% % 10    length
% % 11    heading
% % 12    classification
% % 13:14    coordmin [x,y]
% % 15:16    coordmid [x,y]
% % 17:18    coordmax [x,y]
% % digital degrees are used here e.g -178.9999
% 
% WGS84_TO_REGION_LUT_SCALE = 100;
% 
% % loop through the RoadLink dataset
% RoadLUT0 = zeros(1200,1100,580); % RoadLUT0(lat,long,n-entries). entry=[count,road1,road2,...]): lowest level of road LUT. find roads in a geo area
% MotorwayLUT = zeros(7300,2);
% motorwaycount = 0;
% 
% latstartoffset = -49;
% longstartoffset = +9;
% 
% for i=1:szRLRaw
%     %round(60.8,0) % 61
%     %round(60.8,-1) %60
%     %round(-60.8,0) % -61
%     %round(-60.8,-1) %-60
%     %round(601234,-2) % 601200
%     %round(601255,-2) % 
%     
%     tv = RoadLink(i,15); %608267
%     lat = ( ((latstartoffset*WGS_INT64_SCALE) + round(double(tv),-2)) / WGS84_TO_REGION_LUT_SCALE) + 1;
%     tv = RoadLink(i,16); %-8086;
%     long = ( ((longstartoffset*WGS_INT64_SCALE) + round(double(tv),-2)) / WGS84_TO_REGION_LUT_SCALE ) + 1;
%     [lat long];
%     
%     % add to the motorway count
%     if(RoadLink(i,12)==1)
%         motorwaycount = motorwaycount + 1;
%         MotorwayLUT(motorwaycount,1) = i;
%     end
%     
%     % increase count
%     lutcount = RoadLUT0(lat,long,1);
%     lutcount = lutcount + 1;
%     RoadLUT0(lat,long,1) = lutcount;
%     RoadLUT0(lat,long, lutcount+1) = i;
% end
% 
% save('matlab\big_data\RoadLUT0.mat','RoadLUT0','-v7.3')
% save('matlab\big_data\MotorwayLUT.mat','MotorwayLUT','-v7.3')

%%


%%
% todo
% - combine small road components into full road lengths
% data set is to big because the roads 
% - reissue IDs for all nodes and links and then propagate through 
% - sort geographically into lat/long tree
% This can be used to look up potential components in any given geo region



outputArg1 = inputArg1;
outputArg2 = inputArg2;
end
