% Thomas Pile
% April 2022
% 2022/5/01
% - functionally works but has bodges etc
% v8 todo
% - use RoadNode to store processing info about each node so we can process
% multiple junctions
% - inter frame tracking
% - tidy and simplify
% - test more broadly

% bug: removestartpoints is called on a straight line with a single right
% angle pixel at the end. This should only be called when reaching an edge
% start point. Also it doesn't actually work as the startpoints array iusnt
% passed to the loop!
% that pixel makes it believe its at an edge for some reason :s

% also
% - remember offsets on padimg



%%
% Data Structures

% RoadNode
% Junctions and edgepoints in a road network
% (1) parent node ID
% (2) number of exits
% (3) number which have been explored
% (4) is junction==1, edge marker==0
roadnodeidx = 1;
roadnodemax = 100;
RoadNode(1,1).roadnodeidx = roadnodeidx;
RoadNode(1,roadnodeidx).parentNode = 0;
RoadNode(1,roadnodeidx).exitCount = 0;
RoadNode(1,roadnodeidx).exitsExplored = 0;
RoadNode(1,roadnodeidx).isEdge = 0;
RoadNode(1,roadnodeidx).ExitAngles = [0 0 0 0];
RoadNode(1,roadnodeidx).CentrePixel = [0 0];
RoadNode(1,roadnodeidx).ExitPixelOffsets = [0 0 0 0;
                                          0 0 0 0];     % x
RoadNode(1,roadnodeidx).EntryPixels = [0 0]; % pixels the junction was approached from (so we don't go back down there)


roadnodeidx = 0;
%roadnodeidx = roadnodeidx + 1;

% RoadVector
% Edges between nodes
% (1) startNode
% (2) endNode
% (3) pixelchainidx
% (4) explored
% (5) 
% (6) 
% (7)
roadvectoridx = 1;
roadvectormax = 100;
% use a struct 
RoadVector(1,1).roadvectoridx = roadvectoridx;
RoadVector(1,roadvectoridx).startNode = 0;
RoadVector(1,roadvectoridx).endNode = 0;
RoadVector(1,roadvectoridx).pixelchainidx = 0;
RoadVector(1,roadvectoridx).PixelChain = zeros(50,2);
RoadVector(1,roadvectoridx).explored = 0;

% RoadPixelChain

% IMPORTANT => ABOUT THE VARIABLES!
% current_node: for navigating the roadnode structure
% roadnodeidx: for adding new values to the structure



% StartPoints = zeros(1000,2)

%%



% test environment
% work through a long image to simulate push bvroom effect
imgpb = imread('road_compare_and_sequence\road_sequence_512_1.png');
%imgpb = imread('road_compare_and_sequence\road_sequence_512_2.png');
imgpb = rgb2gray(imgpb);
imgpb = imbinarize(imgpb);

% id y1 x1 y2 x2  
% RoadFinal = zeros(100,5);11
% roadfinalidx = 1;x
szframe = [512 512]; % h, w
padoffset = 6;

% loop through the pushbroom image
%for fi = size(imgpb,1):-1:512
for fi = size(imgpb,1)
    % take a subimage
    img = imgpb(fi-511:fi, :,:);
    imshow(uint8(255*img));
    
    % pad the image 
    padimg = logical(ones(size(img,1)+12, size(img,1)+12));
    %padimg(searchradius+7:size(padimg,1)-6, searchradius+7:size(padimg,2)-6) = (img(:,:));
    padimg(7:size(img,1)+6, 7:size(img,2)+6 ) = (img(:,:));
    imshow(padimg)
    
    % ------------------ method 1 -----------------------------------------
    % 1
    % work around the edges of the image and find pixels with roads in them
    % first ring
    %startpts = zeros(50,3); startptsidx = 1;
    % roadvector[roads_max, pixels_max, 3 (x/y/start edge)] 
    % start edge: 1 lhs, 2 rhs, 3 top, 4 bottom
    roads_max = 512*4; % maximum number of trackable roads in a scene
    pixels_max = 1000; % maximum number of pixels that can make up a road
    %roadvector = zeros(roads_max, pixels_max, 2);
    StartPoints = zeros(1000,2);
    startptsidx = 1;
    %roadvectoridx = ones(roads_max,1); 
    
    % get the start points from around the edge
    [startptsidx, StartPoints] = getStartPoints(img, szframe, startptsidx, StartPoints);
    for i=1:startptsidx
        StartPoints(i,1:2) = StartPoints(i,1:2) + [padoffset padoffset];
    end
    
    % loop through the start points
    for i=1:startptsidx
        left_the_edge = 0;
        end_of_line = 0;
        % loop through pixels
        %while end_of_line == 0 && pixelchainidx(i) < pixels_max 
            start_pixel = StartPoints(i,1:2); % y,x
            startptsidx;
            [roadnodeidx, RoadNode, RoadVector, roadvectoridx] = process_next_pixel(padimg, start_pixel, roadnodeidx, szframe, padoffset, RoadVector, roadvectoridx, StartPoints);
            disp ("--------")
        %end % while not end of line
    end% startpts
end% pushbroom loop

%draw_map_graph(RoadFinal, roadfinalidx)

% the state process to handle detecting potential junctions, finding out
% how many exits, turning this into nodes and edges, associating edges with
% pixel sequences, the tree structure that develops, and eliminating
% duplicates and loops! Better than doing it in loop above
% implemented as a state machine
%
% currentnode when working with data, roadnodeidx when adding a new record
function [roadnodeidx, RoadNode, RoadVector, roadvectoridx] = process_next_pixel(padimg, start_pixel, roadnodeidx, szframe, padoffset, RoadVector, roadvectoridx, StartPoints)
    % state defs
    START_NEXT_PIXEL = 0;
    EDGE_TESTING = 1;
    PROCESS_PIXEL_FAST = 2;
    
    PROCESS_NOT_EDGE_DEPARTED_NO = 5;
    PROCESS_NOT_EDGE_DEPARTED_YES = 6;
    PROCESS_IS_EDGE_DEPARTED_NO = 7; CROSSING_SCENE_IN_A_STRAIGHT_LINE = 12;
    PROCESS_IS_EDGE_DEPARTED_YES = 8;
    
    TEST_JUNCTION_EXIT_COUNT = 17;
    PROCESS_JUNCTION_EXIT = 18;  % replaced by PROCESS_NEXT_ROAD_OR_NODE
    PROCESS_NEXT_ROAD_OR_NODE = 19;
    
    ALL_DONE = 99;
    
    % state variables
    pixel_state = 0;
    left_the_edge = 0;
    end_of_line = 0;
    % tree structure tracking. We traverse the structure!
    current_node = 1;
    current_vector = 1;
    current_pixel = [1 1]; % y, x
    previous_pixel = [1 1]; % y, x
    junction_pixel = [0 0]; % location of the junction centre
    junction_entry_pixel = [0 0]; % the pixel we came from to get to the junction
    AvoidPixels = zeros(8,2);
    %
    pixelchainlengthmax = 100;
    PixelChain = zeros(pixelchainlengthmax,2);
    pixelchainidx = 1;
    pixel_processing = 1;
    loop_count = 0; % which number loop is this
    start_edge = 'N';
    junction_counter_skip = 0; % used to skip the other side of a junction. prevents repeat detections!
    
    while pixel_processing > 0
        pause(0.1)
        disp("----------------------------")
        pixel_state
        roadnodeidx
        current_node
        loop_count = loop_count + 1;
        % process 
        switch pixel_state
            %--------------------------------------------------------------
            case START_NEXT_PIXEL % 0
                % get next pixel
                if start_pixel(1) > 0
                    start_pixel;
                    current_pixel(:) = start_pixel(:);
                    current_pixel;
                    previous_pixel(:) = current_pixel(:);
                end
                %  create a start edge node for completeness
                % intentionally not incrementing the index!!!
                %roadnodeidx = RoadNode(1,1).roadnodeidx;
                roadnodeidx = roadnodeidx + 1;
                RoadNode(1,1).roadnodeidx = roadnodeidx;
                RoadNode(1,roadnodeidx).parentNode = 0;
                RoadNode(1,roadnodeidx).exitCount = 0;
                RoadNode(1,roadnodeidx).exitsExplored = 0;
                RoadNode(1,roadnodeidx).isEdge = 1;
                RoadNode(1,roadnodeidx).ExitAngles = [0 0 0 0];
                RoadNode(1,roadnodeidx).CentrePixel = [0 0];
                RoadNode(1,roadnodeidx).ExitPixelOffsets = [0 0 0 0;
                                                            0 0 0 0];    
                RoadNode(1,roadnodeidx).EntryPixels = [0 0];
                %roadnodeidx = roadnodeidx + 1;                                  
                % next state
                pixel_state = EDGE_TESTING;

            %--------------------------------------------------------------
            % while at the edge, and haven't left the edge
            case EDGE_TESTING % 1
                % have we left the edge yet
                near_edge = edgeCheck(current_pixel, szframe, padoffset);
                %pixel_state = PROCESS_PIXEL_FAST;
                if near_edge > 0
                    disp ("near edge")
                    % we moved away from an edge and are now returning
                    if left_the_edge > 0
                        % have left the edge then gone to a new edge
                        % end of road segment
                        %pixel_state = PROCESS_PIXEL_FAST;
                        pixel_state = PROCESS_IS_EDGE_DEPARTED_YES;
                    else % near edge, haven't left
                        pixel_state = PROCESS_IS_EDGE_DEPARTED_NO;
                        disp ("along edge")
                    end
                else
                    if left_the_edge > 0
                        pixel_state = PROCESS_NOT_EDGE_DEPARTED_YES;
                    else 
                        pixel_state = PROCESS_NOT_EDGE_DEPARTED_NO;
                    end
                    disp ("away edge")
                    left_the_edge = 1;
                end
                
            %--------------------------------------------------------------
%             case PROCESS_PIXEL_FAST % 
%                 % check if its an edge
%                 if edgeCheck(current_pixel, szframe, padoffset)
%                 else                
%                 end
                
            %--------------------------------------------------------------
            % Not left the edge yet, but not at the edge
            case PROCESS_NOT_EDGE_DEPARTED_NO % 5
                % by definition no longer at the edge!
                left_the_edge = 1;
                pixel_state = EDGE_TESTING;
            %--------------------------------------------------------------
            % Normal processing
            % find next pix, add to list, wor back through nodes etc
            case PROCESS_NOT_EDGE_DEPARTED_YES % 6
                 % not at end edge, continue
                [current_pixel, previous_pixel_2, possible_junction] = getNextPixel(padimg, current_pixel, previous_pixel, RoadNode, current_node);
                [current_pixel, previous_pixel_2]
                previous_pixel = previous_pixel_2;
                
                % found next pixel so not end of line
                if current_pixel(1) > 0 && current_pixel(2) > 0
                    % add to pixel list
                    PixelChain(pixelchainidx,1:2) = current_pixel(1:2);
                    %pixel_state = PROCESS_PIXEL_FAST;
                    pixel_state = EDGE_TESTING;
                    % is there a potential junction
                    if possible_junction > 0 && junction_counter_skip == 0
                        pixel_state = TEST_JUNCTION_EXIT_COUNT;
                    end
                    if (junction_counter_skip > 0)
                        junction_counter_skip = junction_counter_skip - 1;
                    end
                else
                    % end of line / dead end
                    % new end node
                    roadnodeidx = RoadNode(1,1).roadnodeidx;
                    roadnodeidx = roadnodeidx + 1;
                    RoadNode(1,1).roadnodeidx = roadnodeidx;
                    
                    RoadNode(1,roadnodeidx).parentNode = current_node;
                    RoadNode(1,roadnodeidx).exitCount = 1;
                    RoadNode(1,roadnodeidx).exitsExplored = 1;
                    RoadNode(1,roadnodeidx).isEdge = 1;
                    RoadNode(1,roadnodeidx).CentrePixel = current_pixel;
                    RoadNode(1,roadnodeidx).ExitPixelOffsets = [0 0 0 0;
                                                                0 0 0 0];    
                    RoadNode(1,roadnodeidx).EntryPixels = [0 0];
                    
                    % add nodes and pix chain as a new road                            
                    RoadVector(roadvectoridx).startNode = RoadNode(roadnodeidx-1).parentNode; %parent_node;
                    RoadVector(roadvectoridx).endNode = roadnodeidx - 1; %current_node;
                    RoadVector(roadvectoridx).pixelchainidx = pixelchainidx;
                    RoadVector(roadvectoridx).PixelChain = PixelChain;
                    RoadVector(roadvectoridx).explored = 1;
                    roadvectoridx = roadvectoridx + 1;
                    pixelchainidx = 1;
                    
                    % prevent the edge case of the road ending before
                    % repeat prevention filter is reset
                    if (junction_counter_skip > 0)
                        junction_counter_skip = 0;
                    end

                    end_of_line = 1;
                    pixel_state = PROCESS_NEXT_ROAD_OR_NODE;
                end
                    
            %--------------------------------------------------------------        
            % Moving along the edge 
            % needs to be able top resume moving along the edge if it
            % leaves for a junction but more hozizontal rd is leftr
            case PROCESS_IS_EDGE_DEPARTED_NO % 7
                % once it has moved one pixel, see which edge is closestr
                % but nottouching. log this as the source edge. Then when
                % it touches the other side...
                if loop_count == 1
                    % find which edge is closest but not next to
                    if current_pixel(2) == 1 % y
                        if current_pixel(1) < (512/2)
                            % departed top
                            start_edge = 'L';
                        else
                            % departed bottom
                            start_edge = 'R';
                        end
                    end
                    if current_pixel(1) == 1 || current_pixel(1) < 512 % x 
                        if current_pixel(2) < (512/2)
                            % departed top
                            start_edge = 'T';
                        else
                            % departed bottom
                            start_edge = 'B';
                        end
                    end
                end
                % after the start edge has been set
                if loop_count > 1
                    % has it crossed along the edge and reached the other
                    % side?
                    crossing_complete = 0;
                    switch start_edge
                        case 'L'
                            if current_pixel(2) < 2
                                crossing_complete = 1;
                            end
                        case 'R'
                            if 512-current_pixel(2) < 2
                                crossing_complete = 1;
                            end
                        case 'B'
                            if current_pixel(1) < 2
                                crossing_complete = 1;
                            end
                        case 'T'
                            if 512-current_pixel(1) < 2
                                crossing_complete = 1;
                            end
                        otherwise
                            disp('Error in: PROCESS_IS_EDGE_DEPARTED_NO')
                    end
                    % paths from that start point fully explored
                    if crossing_complete > 0
                        % add edge node
                        roadnodeidx = RoadNode(1,1).roadnodeidx;
                        roadnodeidx = roadnodeidx + 1;
                        RoadNode(1,1).roadnodeidx = roadnodeidx;
                        RoadNode(1,roadnodeidx).parentNode = current_node;
                        RoadNode(1,roadnodeidx).exitCount = 0;
                        RoadNode(1,roadnodeidx).exitsExplored = 1;
                        RoadNode(1,roadnodeidx).isEdge = 1;
                        %roadnodeidx = roadnodeidx + 1;
                        % remove pixel from startpoints list so we don't repeat it
                        [StartPoints startptsidx] = removeStartPoint(current_pixel, StartPoints, startptsidx);
                        % set road as explored (close out road vector)
                        RoadVector(current_vector).explored = 1;
                        % end for this start point
                        pixel_state = ALL_DONE;
                    end
                end
            %--------------------------------------------------------------    
            % End of road segment in the scene (but continues overleaf)
            case PROCESS_IS_EDGE_DEPARTED_YES % 8
                % close out the RoadVector but mark incomplete
                % add to list of incomplete road vectors
                
                % at an edge so can't go further
                % add edge node
                roadnodeidx = RoadNode(1,1).roadnodeidx;
                roadnodeidx = roadnodeidx + 1;
                RoadNode(1,1).roadnodeidx = roadnodeidx;
                
                RoadNode(1,roadnodeidx).parentNode = current_node;
                RoadNode(1,roadnodeidx).exitCount = 0;
                RoadNode(1,roadnodeidx).exitsExplored = 1;
                RoadNode(1,roadnodeidx).isEdge = 1;
                %roadnodeidx = roadnodeidx + 1;
                % remove pixel from startpoints list so we don't repeat it
                [StartPoints startptsidx] = removeStartPoint(current_pixel, StartPoints, startptsidx);
                % set road as explored
                RoadVector(current_vector).explored = 1;
                % next pixel or node
                pixel_state = PROCESS_NEXT_ROAD_OR_NODE;
                
            %--------------------------------------------------------------
            % junction handling
            case TEST_JUNCTION_EXIT_COUNT % 17
                % count the number of exits
                searchradius = 4;
                [Jmat, jmati] = junction_hunter_exit_counter(padimg, current_pixel, searchradius);
                % if it is definitely a junction, create a new RoadNode (junction)
                if(jmati > 0)
                    % close out existing road
                    %RoadVector(road_current).explored = 1;
                    RoadVector(current_vector).explored = 1;
                    
                    % mark junction centre
                    junction_pixel = current_pixel;
                    
                     % search Jmat for the entry pixerl and remove it
                    % the full list is still in RoadNode
%                     for i=1:jmati
%                         if (Jmat(i,1:2)==previous_pixel(1:2))
%                             % just zero it, since zero values aren't
%                             % processed anyway. Saves an array copy
%                             Jmat(i,1:2)= [0 0];
%                         end
%                     end
                    
                    junction_entry_pixel = previous_pixel;

                    % create new junction node. mark them unexplored
                    roadnodeidx = RoadNode(1,1).roadnodeidx;
                    roadnodeidx = roadnodeidx + 1;
                    RoadNode(1,1).roadnodeidx = roadnodeidx;
                    RoadNode(1,roadnodeidx).parentNode = current_node;
                    RoadNode(1,roadnodeidx).exitCount =  jmati; %  1                           
                    RoadNode(1,roadnodeidx).exitsExplored = 0;
                    RoadNode(1,roadnodeidx).isEdge = 0;
                    RoadNode(1,roadnodeidx).ExitAngles = Jmat(:,1);
                    RoadNode(1,roadnodeidx).ExitPixelOffsets = [Jmat(:,2),...    % y
                                                              Jmat(:,3)];     % x
                    RoadNode(1,roadnodeidx).CentrePixel = current_pixel(1:2);        % guess
                    RoadNode(1,roadnodeidx).EntryPixels = previous_pixel(1:2);
                    
                    %roadnodeidx = roadnodeidx + 1;
                    % set that junction as the current junction
                    current_node = RoadNode(1,1).roadnodeidx;
                    
                    junction_counter_skip = 4;
                    
                    % start processing the exits
                    pixel_state = PROCESS_NEXT_ROAD_OR_NODE;
                else
                    % not a junction after all, keep pixel hunting
                    %pixel_state = PROCESS_PIXEL_FAST;
                    pixel_state = EDGE_TESTING;
                end
            
            %--------------------------------------------------------------    
            case CROSSING_SCENE_IN_A_STRAIGHT_LINE % 12
                if pixelchainidx > 511
                    % if pixelchainidx is frame width then road has crossed
                    % probably all done
                    % exits along the way will have been processed by
                    % next_pixel
                else
                    % crossing but not crossed yet
                    % 
                end
            %--------------------------------------------------------------
            case PROCESS_NEXT_ROAD_OR_NODE % 19
                % current road/node reached a dead end. First check current
                % parent node for more unprocessed exit roads to follow, then
                % keep going back through parent nodes and repeating the
                % process. when all are processed return so the process can run
                % again on a new start point.

                % set new node or vector current
                % go to the process step
                disp(['number of exits explored: ' num2str( RoadNode(current_node).exitsExplored )])                % print
                disp(['is it an edge node? ' num2str( RoadNode(current_node).isEdge )])
                
                % if its an end node
                if (RoadNode(1,current_node).isEdge > 0)
                    %
                    current_node = RoadNode(current_node).parentNode;
                    pixel_state = PROCESS_NEXT_ROAD_OR_NODE;
                    
                    % end condition
                    if (current_node == 0)                                      % was 1
                        pixel_state = ALL_DONE;
                    else
                        % iterate again by keeping the current state
                        pixel_state = PROCESS_NEXT_ROAD_OR_NODE;
                    end
                else
                    % are there any unexplored exits (roads) from current
                    % node/junction?
                    disp(['exploring exit ' num2str(RoadNode(current_node).exitsExplored) ' of ' num2str(RoadNode(current_node).exitCount) ' total exits' ])                % print
                    
                    % if explored all the exits, or if its an edge node
                    if ( (RoadNode(current_node).exitCount > RoadNode(current_node).exitsExplored + 1) ) 

                        % ======================================================================================================================================
                        % select next road to explor
                        exitexploredidx = RoadNode(current_node).exitsExplored;
                        exitexploredidx = exitexploredidx + 1;
                        RoadNode(current_node).exitsExplored = exitexploredidx;
                        exitexploredidx
                        %disp(['RoadNode(current_node).ExitPixelOffsets: '])
                        %RoadNode(current_node).ExitPixelOffsets                                                         % print
                        previous_pixel = junction_pixel;
                        
                        % skip index of the pixel which mnatches the entry
                        % pixel
                        if (RoadNode(current_node).ExitPixelOffsets(exitexploredidx,1:2)==junction_entry_pixel(1:2))
                            exitexploredidx = exitexploredidx + 1;
                        end
                        
                        current_pixel = RoadNode(current_node).ExitPixelOffsets( exitexploredidx,...
                                                                             1:2);
                        % mask out the pixels of the road we came in on
                        % by adding previous pixel to AvoidPixels
                        % this is a duplicate of down usually so you can
                        % check and remove the other if you want
                        disp(['curr pix at state 19 exits found ' num2str(current_pixel(1:2)) ])
                        %AvoidPixels(1,1:2) = previous_pixel(1:2);
                        % add all exits explored to far to the pixels structrure
                        for j=1:exitexploredidx-1
                            AvoidPixels(j,1:2) = RoadNode(current_node).ExitPixelOffsets(j,1:2);
                        end
                        
                        disp(['next exit pixel : ' num2str(current_pixel(1)) ' ' num2str(current_pixel(2)) ])
                        
                        junction_counter_skip = 4;
                        
                        % when jumping back to a earlier node, set a new
                        % previous pixel
                        previous_pixel = RoadNode(current_node).CentrePixel; %junction_pixel;
                        
                        % back to pixel processing
                        %pixel_state = PROCESS_PIXEL_FAST;
                        pixel_state = EDGE_TESTING;
                    else
                        % no more exits at this node. Are there anymore nodes to
                        % explore?
                        % go to parent node, see if there are more exits
                        current_node = RoadNode(current_node).parentNode;
                        parent_node_at_end = RoadNode(current_node).parentNode;
                        parent_node_at_end                                                                              % print

                        % reset avoid pixels array
                        AvoidPixels(1,1:2) = [0 0];
                        junction_pixel = [0 0];
                        junction_entry_pixel = [0 0];
                        
                        
                        
                        % end condition
                        if (current_node == 1)
                        else
                            % iterate again by keeping the current state
                            pixel_state = PROCESS_NEXT_ROAD_OR_NODE;
                        end
                    end
                end
            %--------------------------------------------------------------
            case ALL_DONE % 99
                pixel_processing = 0;
            %--------------------------------------------------------------
            otherwise
                disp('Error in: Pixel state switch')
        end
        % for pixel
        % is end of line?
        % is junction?
        % how many exits
        % end road line, start new ones as many as are exits
    end
end


% find all starting locations around the edge of the frame
function [startptsidx, StartPoints] = getStartPoints(img, szframe, startptsidx, StartPoints)
    % work around the edges of the image and find pixels with roads in them
    % first ring
    
    % down rows - left
    for i=1:szframe(2) 
        if (img(i,1) < 1)
            % add to start locations list
            StartPoints(startptsidx,1:3) = [i 1 1]; % y, x, LHS start
            startptsidx = startptsidx + 1;
        end
    end 
    % right
    for i=1:szframe(2) 
        if (img(i,szframe(2)) < 1)
                % add to start locations list
                StartPoints(startptsidx,1:3) = [i szframe(2) 2];
                startptsidx = startptsidx + 1;
        end
    end 
    % across cols - top
    for i=2:szframe(1)-1
        if (img(szframe(2),i) < 1)
            % add to start locations list
            StartPoints(startptsidx,1:3) = [szframe(2) i 3];
            startptsidx = startptsidx + 1;
        end
    end
    % bottom
    for i=2:szframe(1)-1 
        if (img(1,i) < 1)
            % add to start locations list
            StartPoints(startptsidx,1:3) = [1 i 4];
            startptsidx = startptsidx + 1;
        end
    end
    if(startptsidx > 1)
        startptsidx = startptsidx - 1;
    end
end

% will a surrounding pixel be at the edge?
% starts off 
function is_at_edge = edgeCheck(current_pixel, szframe, offset)
    y = current_pixel()-offset;
    x = current_pixel()-offset;
    is_at_edge = 0;
    if x<2
        is_at_edge = 1;
    end
    if y<2
        is_at_edge = 1;
    end
    if x>szframe(2)-1
        is_at_edge = 1;
    end
    if y>szframe(1)-1
        is_at_edge = 1;
    end
end

% todo
% - needs to not go backwards
function [current_pixel, previous_pixel, possible_junction] = getNextPixel(padimg, current_pixel, previous_pixel, RoadNode, current_node)
    % use a weighted mask subtraction
    y = current_pixel(1);
    x = current_pixel(2);
    
    temp = [padimg(y-1:1:y+1, x-1), padimg(y-1:1:y+1, x), padimg(y-1:1:y+1, x+1)];
    
    mask = [2   64 4;
            128 0  256;
            8   32 16];

    % get pixel locations for current node
    %roadnodeidx = RoadNode(1,1).roadnodeidx;
    junction_pixel = RoadNode(1,current_node).CentrePixel;
    junction_entry_pixel = RoadNode(1,current_node).EntryPixels;
    AvoidPixels = RoadNode(1,current_node).ExitPixelOffsets(:,:);
    
    % if  previous_pixel(1:2) > 0 then we are moving along a line and need
    % to mask out a previous pixel. Otherwise we are masking out junction
    % exits in which case AvoidPixels should be used
    disp(['AvoidPixels(1,1): ' num2str(AvoidPixels(1,1))   ])
    disp(['Current Pixel: ' num2str(current_pixel(1)) ' ' num2str(current_pixel(2))   ])
    disp(['Previous Pixel: ' num2str(previous_pixel(1)) ' ' num2str(previous_pixel(2))   ])
    disp(['Junction Pixel: ' num2str(junction_pixel(1)) ' ' num2str(junction_pixel(2))   ])
    
    if abs(current_pixel-junction_pixel(1)) <2
        AvoidPixels(1,1)
        AvoidPixels(1,2)
        % apply mask out to zeros pixels we don't was
        for i=1:size(AvoidPixels,1)
            if (AvoidPixels(i,1) > 0)
                %mask( AvoidPixels(i,1), AvoidPixels(i,2) ) = 0;
                disp(['Pixel to mask : ' num2str((AvoidPixels(i,1)-junction_pixel(1))+2) ' ' num2str((AvoidPixels(i,2)-junction_pixel(2))+2)  ])  
                mask( ((AvoidPixels(i,1)-junction_pixel(1)))+2,...
                      ((AvoidPixels(i,2)-junction_pixel(2)))+2  ) = 0;
            else
                %i=size(AvoidPixels,1);
                break;
            end
        end
        
        % mask out the already covered exit nodes at the junction
        %RoadNode
        % better apporoach btw!  
    else
       % mask out the previous pixel so we don't go there
        % +2 to shift (0,0) to (2,2) i./e centre pixel etc
        mask( (previous_pixel(1)-current_pixel(1))+2, (previous_pixel(2)-current_pixel(2))+2 ) = 0;
    end
    
    disp(['L589 temp array : '])                  % print
    temp
    mask
    
    result = sum(sum((ones(3,3)-temp).*mask));
    previous_pixel = current_pixel;
    % no further pixels
    if (result <= 1)
        current_pixel = [0 0 0]; 
    end
    % 2 up left
    if (result > 1)
        current_pixel = [y-1 x-1 0]; 
    end
    % 4 up right
    if (result > 3)
        current_pixel = [y-1 x+1 0]; 
    end
    % 8 down left
    if (result > 7)
        current_pixel = [y+1 x-1 0]; 
    end
    % 16 down right
    if (result > 15)
        current_pixel = [y+1 x+1 0]; 
    end
    % 32 down  
    if (result > 31)
        current_pixel = [y+1 x 0]; 
    end
    % 64 up
    if (result > 63)
        current_pixel = [y-1 x 0]; 
    end
    %  128 left
    if (result > 127)
        current_pixel = [y x-1 0]; 
    end
    % 255 right. This is fine, but keeps it 8 bit which is nice
    if (result > 254)
        current_pixel = [y x+1 0]; 
    end
    
    % testing for junction
    % optionally check for junction bulge, bit quicker!
    possible_junction = 0;
    temp = [padimg(y-1:1:y+1, x-1), padimg(y-1:1:y+1, x), padimg(y-1:1:y+1, x+1)];
    regionintegral = sum(sum(temp));
    if regionintegral < 6 % 6      max 9
        possible_junction = 1;
    end
end

% generate a set of locations which make up the outside of the circle, for
% different radius size circles. Avoids repeatedly computing circles and
% fitting to the grid.
function precompute_circle_circumference_points()

    searchradius = 30; % pixels. radius to search around junction
    circleimg = logical(ones(2+(2*searchradius)));
    imshow(circleimg)
    circlecentre = 1+(size(circleimg,1)/2);
    %drawcircle('Color','k','Center',[circlecentre circlecentre],'Radius',searchradius,'InteractionsAllowed', 'none','LineWidth', 1, 'FaceAlpha',0);
    circleimg2 = insertShape(uint8(circleimg),'Circle',[circlecentre circlecentre searchradius],'Color','white','LineWidth', 1);
    imshow(circleimg2)
end


function draw_map_graph(RoadFinal, roadfinalidx)

    % if [max, min] undefined, determine
    szmap = [512 512];  % y x
    
    
    % [x0 y0 x1 y1]

%     LineMap = zeros(100,4);
%     LineMap(1, 1) = 0;
%     LineMap(1, 2) = 0;
%     LineMap(1, 3) = 204;
%     LineMap(1, 4) = 204;
%     %
%     LineMap(2, 1) = 5;
%     LineMap(2, 2) = 5;
%     LineMap(2, 3) = 80;
%     LineMap(2, 4) = 5;
%     %
%     LineMap(3, 1) = 20;
%     LineMap(3, 2) = 20;
%     LineMap(3, 3) = 30;
%     LineMap(3, 4) = 30;

    I = uint8(zeros(512, 512));
    imshow(I)
    for j=1:roadfinalidx
        drawline('Position', [RoadFinal(j, 3), RoadFinal(j, 2); RoadFinal(j, 5), RoadFinal(j, 4)],...
                 'InteractionsAllowed', 'none', ...
                 'LineWidth', 1,...
                 'color', 'red');
    end
end




% performs higher confidence circle test on a specific junction pixel
function [jmat, jmati] = junction_hunter_exit_counter(padimg, current_pixel, searchradius)
    jmat = zeros(10,5);
    jmati = 1;
    
    y = current_pixel(1);
    x = current_pixel(2);
    
    if (searchradius==3)
        % radius 3
        dispimg = [padimg(y-2:1:y+2,x-2),...
                   padimg(y-2:1:y+2,x-1),...
                   padimg(y-2:1:y+2,x),...
                   padimg(y-2:1:y+2,x+1),...
                   padimg(y-2:1:y+2,x+2)]
        nodeID = 0;
        % top
        if ( padimg(y-2,x-1)<1 || padimg(y-2,x)<1 || padimg(y-2,x+1)<1 ) 
            disp("Up")
            % detected up
            jmat(jmati,:) = [0 y-1 x y-2 x];
            jmati = jmati + 1;
        end
        % right
        if ( padimg(y-1,x-2)<1 || padimg(y,x-2)<1 || padimg(y+1,x-2)<1 ) 
            disp("Right")
            jmat(jmati,:) = [90 y x+1 y x+2];
            jmati = jmati + 1;
        end
        % bottom
        if ( padimg(y+2,x-1)<1 || padimg(y+2,x)<1 || padimg(y+2,x+1)<1 ) 
            disp("Down")
            jmat(jmati,:) = [180 y+1 x y+2 x];
            jmati = jmati + 1;
        end
        % left
        if ( padimg(y-1,x-2)<1 || padimg(y,x-2)<1 || padimg(y+1,x-2)<1 ) 
            disp("Left")
            jmat(jmati,:) = [270 y x-1 y -x2];
            jmati = jmati + 1;
        end
    end

    % padding 3
    % too smallt o effectively detect diagonals unless modified further
    if (searchradius==4)
        padimg(y-3:1:y+3,x)
        dispimg = [padimg(y-3:1:y+3,x-3),...
                   padimg(y-3:1:y+3,x-2),...
                   padimg(y-3:1:y+3,x-1),...
                   padimg(y-3:1:y+3,x),...
                   padimg(y-3:1:y+3,x+1),...
                   padimg(y-3:1:y+3,x+2),...
                   padimg(y-3:1:y+3,x+3)]
        % radius 4
        % top
        if ( padimg(y-3,x)<1 ) 
            disp("Up");
            jmat(jmati,:) = [0 y-1 x y-2 x];
            jmati = jmati + 1;
        end
        % right
        if ( padimg(y,x+3)<1 ) 
            disp("Right");
            jmat(jmati,:) = [90 y x+1 y x+2];
            jmati = jmati + 1;
        end
        % bottom
        if ( padimg(y+3,x)<1 ) 
            disp("Down");
            jmat(jmati,:) = [180 y+1 x y+2 x];
            jmati = jmati + 1;
        end
        % left
        if ( padimg(y,x-3)<1 ) 
            disp("Left");
            jmat(jmati,:) = [270 y x-1 y x-2];
            jmati = jmati + 1;
        end
        % TR
        if ( padimg(y-2,x+2)<1)
            disp("Top Right");
            jmat(jmati,:) = [45 y-1 x+1 y-2 x+2];
            jmati = jmati + 1;
        end
        % BR
        if ( padimg(y+2,x+2)<1) 
            disp("Bottom Right");
            jmat(jmati,:) = [135 y+1 x+1 y+2 x+2];
            jmati = jmati + 1;
        end
        % BL
        if ( padimg(y+2,x-2)<1) 
            disp("Bottom Left");
            jmat(jmati,:) = [225 y+1 x-1 y+2 x-2];
            jmati = jmati + 1;
        end
        % TL
        if ( padimg(y-2,x-2)<1) 
            disp("Top Left");   %
            jmat(jmati,:) = [315 y-1 x-1 y-2 x-2];
            jmati = jmati + 1;
        end
    end

    % padding 5
    % radius 6
    % work work, needs a longer pixel exit sequence
    % als the JMAT values are wrong, fix like above versions!
    if(searchradius==6)
        % top
        if ( padimg(y-5,x)<1 ) 
            disp("Up");
            jmat(jmati,:) = [0 y x];
            jmati = jmati + 1;
        end
        % bottom
        if ( padimg(y+5,x)<1 ) 
            disp("Down");
            jmat(jmati,:) = [180 y x];
            jmati = jmati + 1;
        end
        % left
        if ( padimg(y,x-5)<1 ) 
            disp("Left");
            jmat(jmati,:) = [270 y x];
            jmati = jmati + 1;
        end
        % right
        if ( padimg(y,x-5)<1 ) 
            disp("Right");
            jmat(jmati,:) = [90 y x];
            jmati = jmati + 1;
        end

        % TR-T
        if ( padimg(y-5,x+2)<1  ) 
            disp("TR-T");
            jmat(jmati,:) = [22 y x];
            jmati = jmati + 1;
        end
        % TR
        if ( padimg(y-4,x+4)<1 ) 
            disp("TR");
            jmat(jmati,:) = [45 y x];
            jmati = jmati + 1;
        end
        % TR-B
        if ( padimg(y-2,x+5)<1 ) 
            disp("TR-B");
            jmat(jmati,:) = [67 y x];
            jmati = jmati + 1;
        end

        % BR-T
        if ( padimg(y+5,x+2)<1 )
            disp("BR-T");
            jmat(jmati,:) = [112 y x];
            jmati = jmati + 1;
        end
        % BR
        if ( padimg(y+4,x+4)<1 ) 
            disp("BR");
            jmat(jmati,:) = [135 y x];
            jmati = jmati + 1;
        end
        % BR-B
        if ( padimg(y+2,x+5)<1 ) 
            disp("BR-B");
            jmat(jmati,:) = [157 y x];
            jmati = jmati + 1;
        end

        % BL-T
        if ( padimg(y+5,x-2)<1 ) 
            disp("TL-T");
            jmat(jmati,:) = [247 y x];
            jmati = jmati + 1;
        end
        % BL
        if ( padimg(y+4,x-4)<1 ) 
            disp("TL");
            jmat(jmati,:) = [225 y x];
            jmati = jmati + 1;
        end
        % BL-B
        if ( padimg(Y+2,x-5)<1 ) 
            disp("TL-T");
            jmat(jmati,:) = [202 y x];
            jmati = jmati + 1;
        end

        % TL-T
        if ( padimg(y-5,x-2)<1 )
            disp("TL-T");
            jmat(jmati,:) = [337 y x];
            jmati = jmati + 1;
        end
        % TL
        if ( padimg(y-4,x-4)<1 )
            disp("TL");
            jmat(jmati,:) = [315 y x];
            jmati = jmati + 1;
        end
        % TL-B
        if ( padimg(y-2,x-5)<1 )
            disp("TL-B");
            jmat(jmati,:) = [292 y x];
            jmati = jmati + 1;
        end
    end
    
    % ignore junctions with 2 or less exits because that's a corner :)
    jmati = jmati - 1;
    if jmati < 3
        % disregard as a junction
        jmati = 0;
    end
end

% remove a location from StartPoints
% this should remove any neighbouring pixel from the sta1rtpoints too
% remember the offset on x,y
function [StartPoints2 startptsidx2] = removeStartPoint(remove_pixel, StartPoints, startptsidx)
    y = remove_pixel(1);
    x = remove_pixel(2);
    startptsidx2 = 1;
    StartPoints2 = zeros(size(StartPoints));
    
    % loop through start points vector
    for i = 1:startptsidx
        % compare to start points
        for j=-1:1:1
            for k=-1:1:1
                % remove
                if (StartPoints(i,2)==y+i && StartPoints(i,2)==x+j)
                    % no nothing
                else
                    StartPoints2(startptsidx2,1:2) = StartPoints(i,1:2);
                    startptsidx2 = startptsidx2 + 1;
                end
            end
        end
    end
end

% returns the index of the node at that pixel location [y x], and 0 if it
% cant be found
function idx = getNodeByPixel(RoadNode, y, x)
    idx = 0;
    maxidx = RoadNode(1,1).roadnodeidx;
    for i=1:maxidx
        if RoadNode(1,i).CentrePixel == [y x];
            idx = i;
            break;
        end
    end
end

% returns 1 if a road node exists at the pixel, 0 otherwise
function exists = getNodeByPixelExists(RoadNode, y, x)
    exists = 0;
    maxidx = RoadNode(1,1).roadnodeidx;
    for i=1:maxidx
        if RoadNode(1,i).CentrePixel == [y x];
            exists = 1;
            break;
        end
    end
end




