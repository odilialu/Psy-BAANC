% Elevated plus maze data analysis after deeplabcut training 
% Written by Odilia Lu, last edited October 7, 2022


%% variables to manipulate: 
path = "D:\EPMAcute\Videos\males\sal";

videoType = 'wmv'; % avi is recorded from biobserve, wmv from logitech webcam. 
framerate = 30; %30 fps 
start_min = 0; 
end_min = 10; 

epmLength = 35; % length of an open arm in centimeters

neckX_columnNum = 5; % which column of the DLC data file has the nose X coordinates? 
bodyX_columnNum = 14; % which column of the DLC data file has the body X coordinates? 
pawBackLeftX_columnNum = 8; 
pawBackRightX_columnNum = 11; 


%% Create a string array of files of interest
strVids = strArray(path, videoType); 
strFiles = strArray(path, 'csv'); 

%% pre-allocate variables
startExp = start_min*60*framerate+1; 
lengthExp = end_min*60*framerate;

bodyTime_openCenter = NaN(length(strVids), 1);
distance = NaN(length(strFiles), 1); 

    timeOpen = NaN(length(strVids), 2); 
    visitsOpen = NaN(length(strVids), 2); 
    distanceTotal = NaN(length(strVids), 2); 
    distanceAverage = NaN(length(strVids), 2); 
    distanceMax = NaN(length(strVids), 2); 
    peeps = NaN(length(strVids), 2); 
    timePeeps = NaN(length(strVids), 2); 
    timeOpenBody = NaN(length(strVids), 2); 
    visitsOpenBody = NaN(length(strVids), 2); 

    timeClosed = NaN(length(strVids), 2); 
    visitsClosed = NaN(length(strVids), 2); 
    timeClosedBody = NaN(length(strVids), 2); 
    visitsClosedBody = NaN(length(strVids), 2); 

%% Calibration - pixels to centimeters. 
figure()
    video = VideoReader(strVids(1, 1)); % read the first video in the list. 
    frame = read(video, 1); %grab the first frame from each video
    imshow(frame)

    hold on 
    fprintf(['[Pixels to Centimeters calibration] \n' ...
        '1. Left click to define a coordinate marking one end of the open arm. Hit enter. \n' ...
        '2. Left click to define a coordinate marking the other end of the open arm. Hit enter. \n'])
    [pixelCoords(1,1), pixelCoords(1, 2)] = getpts();
    [pixelCoords(2, 1), pixelCoords(2,2)] = getpts();
    pixelDistance = pdist(pixelCoords, 'euclidean');
    pixelsToCM = epmLength/pixelDistance;

%% analysis
for i = 1:length(strVids)

    % extract coordinates of interest 
    coordinates = readmatrix(strFiles(i, 1)); 

    pawBackLeft = round(coordinates(startExp:lengthExp, pawBackLeftX_columnNum:pawBackLeftX_columnNum+1)); 
    pawBackRight = round(coordinates(startExp:lengthExp, pawBackRightX_columnNum:pawBackRightX_columnNum+1)); %
    neck = round(coordinates(startExp:lengthExp, neckX_columnNum:neckX_columnNum+1)); 
    body = round(coordinates(startExp:lengthExp, bodyX_columnNum:bodyX_columnNum+1)); 

    % distance travelled
    bodyDiff = diff(body); 
    bodyDistanceChange = NaN(length(bodyDiff), 1); 
    for k = 1:length(bodyDiff)
        bodyDistanceChange(k, 1) = sqrt((bodyDiff(k, 1)^2) + (bodyDiff(k, 2)^2));
    end
    distance(i, 1) = sum(bodyDistanceChange)*pixelsToCM; 

    %%%%%%%%%%%%%%%%%%%%%% DRAW ZONES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw zones for open arm 1 and 2 if you haven't already
    pathSplit = split(strVids(i, 1), "/"); 
    openOneString = strcat('openOne-', pathSplit(end, 1), '.mat'); 
    openTwoString = strcat('openTwo-', pathSplit(end, 1), '.mat'); 

    if isfile(openOneString)
        trackAreaOne = importdata(openOneString);
        trackAreaTwo = importdata(openTwoString); 
    else

    figure()
    video = VideoReader(strVids(i, 1)); % read the first video in the list. 
    frame = read(video, 1); %grab the first frame from each video
    imshow(frame)

    hold on
    fprintf('Please draw polygon of open arm 1 (no center) ...\n')
    trackAreaOne = drawpolygon;
    save(openOneString, "trackAreaOne")

    hold on
    fprintf('Please draw polygon of open arm 2 (no center)...\n')
    trackAreaTwo = drawpolygon; 
    save(openTwoString, "trackAreaTwo")

    end

    % define open+ center zone. 
    openCenterString = strcat('openCenter-', pathSplit(end, 1), '.mat'); 

    if isfile(openCenterString)
        trackAreaOpenCenter = importdata(openCenterString); 
    else
    
    figure()
    video = VideoReader(strVids(i, 1)); 
    frame = read(video, 1); 
    imshow(frame)

    hold on 
    fprintf('Please draw polygon of open + center zones ... \n')
    trackAreaOpenCenter = drawpolygon; 
    save(openCenterString, "trackAreaOpenCenter")

    end

    % Define closed arm one and two if you haven't already. 
    closedOneString = strcat('closedOne-', pathSplit(end, 1), '.mat');
    closedTwoString = strcat('closedTwo-', pathSplit(end, 1), '.mat');

    if isfile(closedOneString)
        trackAreaClosedOne = importdata(closedOneString);
        trackAreaClosedTwo = importdata(closedTwoString); 

    else
    
    figure()
    video = VideoReader(strVids(i, 1)); % read the first video in the list. 
    frame = read(video, 1); %grab the first frame from each video
    imshow(frame)

    hold on
    fprintf('Please draw polygon of closed arm 1...\n')
    trackAreaClosedOne = drawpolygon;
    save(closedOneString, "trackAreaClosedOne")

    hold on
    fprintf('Please draw polygon of closed arm 2...\n')
    trackAreaClosedTwo = drawpolygon; 
    save(closedTwoString, "trackAreaClosedTwo")
 
    end


    %%%%%%%%%%%%%%%%%%%% create and plot masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    video = VideoReader(strVids(i, 1));

    maskOne = poly2mask(trackAreaOne.Position(:,1), trackAreaOne.Position(:,2), video.Height, video.Width); %make mask to filter out area outside polygon
    polyShapeOne = polyshape(trackAreaOne.Position); 
    maskTwo = poly2mask(trackAreaTwo.Position(:, 1), trackAreaTwo.Position(:, 2), video.Height, video.Width); 
    polyShapeTwo = polyshape(trackAreaTwo.Position); 

    maskOpenCenter = poly2mask(trackAreaOpenCenter.Position(:, 1), trackAreaOpenCenter.Position(:, 2), video.Height, video.Width); 
    polyShapeOpenCenter = polyshape(trackAreaOpenCenter.Position); 

    maskClosedOne = poly2mask(trackAreaClosedOne.Position(:,1), trackAreaClosedOne.Position(:,2), video.Height, video.Width); %make mask to filter out area outside polygon
    polyShapeClosedOne = polyshape(trackAreaClosedOne.Position); 
    maskClosedTwo = poly2mask(trackAreaClosedTwo.Position(:,1), trackAreaClosedTwo.Position(:,2), video.Height, video.Width); %make mask to filter out area outside polygon
    polyShapeClosedTwo = polyshape(trackAreaClosedTwo.Position);

    figure()
    frame = read(video, 1); 
    imshow(frame)
    hold on
    plot(body(:, 1), body(:, 2)); 
    plot(neck(:, 1), neck(:, 2));
    plot(pawBackLeft(:, 1), pawBackLeft(:, 2));
    plot(pawBackRight(:, 1), pawBackRight(:, 2));
    hold on
    plot(polyShapeOne)
    plot(polyShapeTwo)
    plot(polyShapeOpenCenter)
    plot(polyShapeClosedOne)
    plot(polyShapeClosedTwo)

    %% main analysis 

    %%%%%% detailed analysis for each arm of the EPM %%%%%%%%%
    [timeOpen(i, 1), visitsOpen(i,1), distanceTotal(i,1), distanceAverage(i,1), distanceMax(i,1), peeps(i,1), timePeeps(i,1), timeOpenBody(i,1), visitsOpenBody(i,1)] = epm_func(maskOne, lengthExp, framerate, pixelsToCM, pawBackLeft, pawBackRight, neck, body);
    [timeOpen(i, 2), visitsOpen(i,2), distanceTotal(i,2), distanceAverage(i,2), distanceMax(i,2), peeps(i,2), timePeeps(i,2), timeOpenBody(i,2), visitsOpenBody(i,2)] = epm_func(maskTwo, lengthExp, framerate, pixelsToCM, pawBackLeft, pawBackRight, neck, body);
    [timeClosed(i, 1), visitsClosed(i,1), ~, ~, ~, ~, ~, timeClosedBody(i,1), visitsClosedBody(i,1)] = epm_func(maskClosedOne, lengthExp, framerate, pixelsToCM, pawBackLeft, pawBackRight, neck, body);
    [timeClosed(i, 2), visitsClosed(i,2), ~, ~, ~, ~, ~, timeClosedBody(i,2), visitsClosedBody(i,2)] = epm_func(maskClosedTwo, lengthExp, framerate, pixelsToCM, pawBackLeft, pawBackRight, neck, body);

    %%%%%%%%%%%%%%%% determine when the animal's body is in the Open+ Center Zones %%%%%%%%%%%%%%%%%%%%%%%
    [yCoordMaskOpenCenter, xCoordMaskOpenCenter] = find(maskOpenCenter); 
    coordMaskOpenCenter = cat(2, xCoordMaskOpenCenter, yCoordMaskOpenCenter); % x and y coordinates of the open / center zones. 
    bodyTime_openCenter(i, 1) = length(find(ismember(body, coordMaskOpenCenter, 'rows')))/framerate; 
    
end


%% summary data

    ntimeOpen = sum(timeOpen, 2); 

    nVisitsOpen = sum(visitsOpen, 2); 
   
    nDistanceOpen = sum(distanceTotal, 2); 
    nDistanceAverage = zeros(size(distanceAverage, 1), 1); 
    for k = 1:size(distanceAverage, 1) 
        if sum(distanceAverage(k, :)) == 0
            continue
        else
        nDistanceAverage(k, 1) = mean(nonzeros(distanceAverage(k, 1:2)), 'all');
        end
    end
    nDistanceMax = max(distanceMax,[], 2);

    nPeeps= sum(peeps, 2) - nVisitsOpen; 
    nTimePeeps = sum(timePeeps, 2); 
  
    nTimeOpenBody = sum(timeOpenBody, 2); 
    nVisitsOpenBody = sum(visitsOpenBody, 2); 

    nTimeOpenNorm = ntimeOpen./distance; 
    nTimeOpenBodyNorm = nTimeOpenBody./distance; 

    nTimeClosed = sum(timeClosed, 2); 
    nTimeClosedBody = sum(timeClosedBody, 2); 
    nVisitsClosed = sum(visitsClosed, 2); 
    nVisitsClosedBody = sum(visitsClosedBody, 2); 

    nTotalVisits = nVisitsOpen + nVisitsClosed;
    nVisitsOpenNorm = nVisitsOpen./nTotalVisits*100;
    nTotalVisitsBody = nVisitsOpenBody + nVisitsClosedBody;
    nVisitsOpenNormBody = nVisitsOpenBody./nTotalVisitsBody*100;

    nTimeCenter = bodyTime_openCenter - nTimeOpenBody;

summary = table(ntimeOpen, nVisitsOpen, nDistanceOpen, nDistanceAverage, nDistanceMax, nPeeps, nTimePeeps, nVisitsOpenNorm, nTimeOpenNorm, distance);
summaryBody= table(nTimeOpenBody, nVisitsOpenBody, nTimeOpenBodyNorm, nVisitsOpenNormBody, nTimeCenter, bodyTime_openCenter, distance);

%% Functions 
% 1. get filepaths of interest in a folder, based on filetype 
function strFiles = strArray(path, fileType)
sprintfDir = strcat(path, "\%s");
filesDirectory = dir(path); 

listFiles = strings(length(filesDirectory), 1); 
% only get files with correct file type
for i = 1:length(filesDirectory)
    if endsWith(filesDirectory(i).name, fileType) == 1
        listFiles(i, 1) = filesDirectory(i).name; 
    end
end
listFiles(cellfun('isempty', listFiles)) = []; 

% create a string array with paths to files of interest
strFiles = strings([length(listFiles), 1]); 
for i = 1:length(listFiles)
    file = listFiles(i, 1); 
    strFiles(i, 1) = sprintf(strrep(sprintfDir, '\', '/'), file);
end 
end

% 2. epm main analysis 
function [timeOpen, visits, distanceTotal, distanceAverage, distanceMax, peeps, timePeeps, timeOpenBody, visitsBody] = epm_func(mask, lengthExp, framerate, pixelsToCM, pawBackLeft, pawBackRight, neck, centerMass)

    [yCoordMaskOne, xCoordMaskOne] = find(mask); 
    coordMaskOne = cat(2, xCoordMaskOne, yCoordMaskOne); % x and y coordinates of zone one 

    % strict definition: animal is defined as within the zone if both of its back paws and neck are within the zone. 
    ind_pawCoordMaskOne_temp = find(ismember(pawBackLeft, coordMaskOne, 'rows')); % find the frames of pawBackLeft that are within the zone. 
    ind_pawCoordMaskOne_temp2 = find(ismember(pawBackRight, coordMaskOne, 'rows')); % find the frames of pawBackRight that are within the zone.
    ind_pawCoordMaskOne_temp3 = intersect(ind_pawCoordMaskOne_temp, ind_pawCoordMaskOne_temp2); % frames where both paws are in the zone. 
    ind_pawCoordMaskOne_neckInc = find(ismember(neck, coordMaskOne, 'rows'));% frames where neck is in the zone. 
    ind_pawCoordMaskOne = intersect(ind_pawCoordMaskOne_temp3,ind_pawCoordMaskOne_neckInc); %frames have to include when neck is also in the zone. 

    timeOpen = length(ind_pawCoordMaskOne)/framerate; 

    %%%%% separate data into exploration bouts, for visits / distance in arm quantifications %%%%% 
    % find the body coordinates of the animal when it is in the zone. 
    bodyCoordMaskOne = NaN(lengthExp, 2); 
    for p = 1:length(ind_pawCoordMaskOne)
        bodyCoordMaskOne(ind_pawCoordMaskOne(p, 1), 1:2) = centerMass(ind_pawCoordMaskOne(p, 1), 1:2); 
    end

    % Set a threshold of 1 sec for animal to be in / out of a zone. 
    for q = 2:(length(bodyCoordMaskOne)-framerate)
        if isnan(bodyCoordMaskOne(q, 1))
            if isnan(bodyCoordMaskOne(q:q+framerate, 1)) % if there are continuous NaNs for one second
                continue
            else 
                bodyCoordMaskOne(q, 1:2) = bodyCoordMaskOne(q-1, 1:2); %fillmissing with the last available coordinate. 
            end
        end
    end

    % separate data into bouts of exploration periods. 
    propsOne = regionprops(~isnan(bodyCoordMaskOne), bodyCoordMaskOne, 'PixelValues'); 
    for l = 1:length(propsOne)
        propsOne(l).PixelValues = reshape(propsOne(l).PixelValues, length(propsOne(l).PixelValues)/2, 2); 
    end
    
    if isempty(propsOne) == 0
        % determine the number of visits
        visits = length(propsOne);
    
        %determine the distance travelled in the arm, stored in disOne_tot_sum
        propsOneCell = struct2cell(propsOne); 
        disOne_tot = NaN(length(propsOneCell), 1); 
        for m = 1:length(propsOneCell)
            disOne = hypot(diff(propsOneCell{1, m}(:, 1)), diff(propsOneCell{1, m}(:, 2))); 
            disOne_tot(m, 1) =  sum(disOne); 
        end
        distanceTotal = sum(disOne_tot)*pixelsToCM; 
        distanceAverage = mean(disOne_tot)*pixelsToCM; 
        distanceMax = max(disOne_tot)*pixelsToCM;
    else
        visits = 0;
        distanceTotal = 0; 
        distanceAverage = 0; 
        distanceMax = 0; 
    end
    
    
    %%%%%%%%%%%%%%%% Determining peeps %%%%%%%%%%%%%%%%%%%%%%%%%%
    ind_neckNoPaws = setdiff(ind_pawCoordMaskOne_neckInc, ind_pawCoordMaskOne_temp3); % frames when neck but not back legs are in the zone. 

    % find the neck coordinates of the animal when it is in the zone. 
    neckCoordMaskOne = NaN(lengthExp, 2); 
    for p = 1:length(ind_neckNoPaws)
        neckCoordMaskOne(ind_neckNoPaws(p, 1), 1:2) = neck(ind_neckNoPaws(p, 1), 1:2); 
    end

    % Set a threshold of 1 sec for animal to be in / out of a zone. 
    for q = 2:(length(neckCoordMaskOne)-framerate)
        if isnan(neckCoordMaskOne(q, 1))
            if isnan(neckCoordMaskOne(q:q+framerate, 1))
                continue
            else 
                neckCoordMaskOne(q, 1:2) = neckCoordMaskOne(q-1, 1:2); 
            end
        end
    end

    % separate data into bouts of exploration periods. 
    propsOne_neck = regionprops(~isnan(neckCoordMaskOne), neckCoordMaskOne, 'PixelValues'); 
    for l = 1:length(propsOne_neck)
        propsOne_neck(l).PixelValues = reshape(propsOne_neck(l).PixelValues, length(propsOne_neck(l).PixelValues)/2, 2); 
    end

    % determine the number of peeps
    if isempty(propsOne_neck) ==0
        peeps = length(propsOne_neck);
    else
        peeps = 0; 
    end

    % determine length of peeps 
    propsOne_neckCell = struct2cell(propsOne_neck);
    timePerPeep = NaN(length(propsOne_neckCell), 1);
    for m = 1:length(propsOne_neckCell)
        timePerPeep(m, 1) = size(propsOne_neckCell{1, m}, 1);
    end
    timePeeps = sum(timePerPeep)/framerate;

    %%%%%% Analysis with loose definition of open arm entry: animal is defined as within the zone if its center mass is in the zone. 
    ind_bodyCoordMaskOne = find(ismember(centerMass, coordMaskOne, 'rows'));
    timeOpenBody = length(ind_bodyCoordMaskOne)/framerate; 

    % find the body coordinates of the animal when it is in the zone. 
    bodyCoordLoose = NaN(lengthExp, 2); 
    for p = 1:length(ind_bodyCoordMaskOne)
        bodyCoordLoose(ind_bodyCoordMaskOne(p, 1), 1:2) = centerMass(ind_bodyCoordMaskOne(p, 1), 1:2); 
    end

    % Set a threshold of 1 sec for animal to be in / out of a zone. 
    for q = 2:(length(bodyCoordLoose)-framerate)
        if isnan(bodyCoordLoose(q, 1))
            if isnan(bodyCoordLoose(q:q+framerate, 1)) % if there are continuous NaNs for one second
                continue
            else 
                bodyCoordLoose(q, 1:2) = bodyCoordLoose(q-1, 1:2); % fillmissing with the last available coordinate. 
            end
        end
    end

    % separate data into bouts of exploration periods. 
    propsLoose = regionprops(~isnan(bodyCoordLoose), bodyCoordLoose, 'PixelValues'); 
    for l = 1:length(propsLoose)
        propsLoose(l).PixelValues = reshape(propsLoose(l).PixelValues, length(propsLoose(l).PixelValues)/2, 2); 
    end
    
    if isempty(propsLoose) == 0
        visitsBody = length(propsLoose);
    else
        visitsBody= 0; 

    end


end
