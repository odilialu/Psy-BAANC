% Object exploration video analysis
% Odilia Lu, last edited 8/29/22 

%{ 
This script quantifies object exploration, defined as when the nose is within an object
zone, but not the body (which implies the animal is on top of the object). 

Instructions: 
1. Predict nose and body markers using DeepLabCut. 
2. Then, change parameters in the variables to change section below.
3. View summaryStats for relevant quantifications. Each row is an animal in
the folder (alphabetically organized). 

%}


%% variables to change
path = "D:\NOEPers\1-10"; % path to folder with all the videos and DLC analysis files

videoType = 'mp4';
framerate = 15; %30 fps 
start_min = 0; 
end_min = 10; 

% which column has the x coordinates for the following markers? 
noseX_columnNum = 2; 
bodyX_columnNum = 17;

% These parameters should be adjusted based on how large objects are. 
% Lammel lab radius: higher res = 43; lower res = 23; cms = 34, 32; 31; 
% Lammel lab radius exclude: higher res = 31; lower res = 16; cms = 29,27; 26
radiusLeft = 43;
radiusRight = 43; 
radiusLeftExclude = 31; 
radiusRightExclude = 31;

%% loop through all the CSV and AVI files in the folder. 
strFiles = strArray(path, 'csv');
strVids = strArray(path, videoType);

%% Analysis 

startExp = framerate*60*start_min+1; 
lengthExp = framerate*60*end_min+1;

% pre-allocate arrays
bodyLeftFrames = cell(1, length(strFiles)); 
bodyRightFrames = cell(1, length(strFiles)); 

noseLeftFrames = cell(1, length(strFiles)); 
noseRightFrames = cell(1, length(strFiles)); 

objectX = cell(1, length(strFiles)); 
objectY = cell(1, length(strFiles)); 

% get frames that animal has nose in object zone and frames when animal's
% body is on top of the object. 
for i = 1:length(strVids)

    video = VideoReader(strVids(i, 1)); 
    frame = read(video, 1); % extract the first frame of every video

    % read in csv data, get coordinates for body and nose
    data = readmatrix(strFiles(i,1)); 
    data = data(startExp:lengthExp, :);
    coordinatesNose(:, 1) = data(:, noseX_columnNum); % x data for nose tip
    coordinatesNose(:, 2) = data(:, noseX_columnNum+1); % y data for nose tip
    coordinatesBody(:, 1) = data(:, bodyX_columnNum); 
    coordinatesBody(:, 2) = data(:, bodyX_columnNum+1);

    % find the x and y coordinates for objectLeft and objectRight
    pathSplit = split(path, ":"); 
    objectXString = strcat('objectXCoord-', pathSplit(end, 1), '.mat'); 
    objectXString = strrep(objectXString, "\", "-"); 
    objectYString = strcat('objectYCoord-', pathSplit(end, 1), '.mat'); 
    objectYString = strrep(objectYString, "\", "-"); 

    if isfile(objectXString) % if you've previously defined objects, import locations
        objectX = importdata(objectXString);
        objectY = importdata(objectYString); 
        figure()
        imshow(frame)
        hold on
    else % if you've never defined object locations, get object locations. 
        figure()
        imshow(frame)
        hold on
        disp("1) left click on center of first object, 2) left click on center of second object, 3) hit enter key.")
        [objectX{1, i}, objectY{1, i}] = getpts; % select coordinate of the objects
    end

    objectLeftCenter = [objectX{1, i}(1, 1), objectY{1, i}(1, 1)];
    objectRightCenter = [objectX{1, i}(2, 1), objectY{1, i}(2, 1)];

    %time exploring object is defined as when the nose tip is within a
    %certain radius of the center of the object
    noseLeftDistance = (coordinatesNose(:, 1)-objectLeftCenter(1,1)).^2 + (coordinatesNose(:, 2)-objectLeftCenter(1, 2)).^2; 
    noseRightDistance = (coordinatesNose(:, 1)-objectRightCenter(1,1)).^2 + (coordinatesNose(:, 2)-objectRightCenter(1, 2)).^2; 
    bodyLeftDistance = (coordinatesBody(:, 1)-objectLeftCenter(1, 1)).^2 + (coordinatesBody(:, 2)-objectLeftCenter(1, 2)).^2;
    bodyRightDistance = (coordinatesBody(:, 1) - objectRightCenter(1, 1)).^2 + (coordinatesBody(:, 2)-objectRightCenter(1, 2)).^2; 

    noseLeftFrames{1, i} = find(noseLeftDistance<= (radiusLeft)^2 & noseLeftDistance >= (radiusLeftExclude)^2); 
    noseRightFrames{1, i} = find(noseRightDistance <= (radiusRight)^2 & noseRightDistance >= (radiusRightExclude)^2); 
    bodyLeftFrames{1, i} = find(bodyLeftDistance <= (radiusLeftExclude)^2); 
    bodyRightFrames{1, i} = find(bodyRightDistance <= (radiusRightExclude)^2); 

    % visualization
    plot(coordinatesNose(:, 1), coordinatesNose(:, 2))
    plot(objectX{1, i}(1,1), objectY{1, i}(1,1), 'd')
    plot(objectX{1, i}(2, 1), objectY{1, i}(2, 1), 'd')
    viscircles(objectLeftCenter, radiusLeft)
    viscircles(objectRightCenter, radiusRight)
    viscircles(objectLeftCenter, radiusLeftExclude)
    viscircles(objectRightCenter, radiusRightExclude)
end

    save(objectXString, 'objectX')
    save(objectYString, 'objectY')

%% final results - secondsLeftTotal, secondsRightTotal, secondsTotal.

secondsLeftTotal = NaN(length(strFiles), 1); 
secondsRightTotal = NaN(length(strFiles), 1); 

% eliminate frames where both the body and nose are in bounds.
for i = 1:length(strFiles)
    framesExclude = length(intersect(noseRightFrames{1, i}, bodyRightFrames{1, i})); 
    secondsRightTotal(i, 1) = (length(noseRightFrames{1, i}) - framesExclude)/framerate; 

    framesExclude = length(intersect(noseLeftFrames{1, i}, bodyLeftFrames{1, i})); 
    secondsLeftTotal(i, 1) = (length(noseLeftFrames{1, i}) - framesExclude)/framerate; 
end

% get additional quantifications
discrIndex = (secondsRightTotal-secondsLeftTotal)./(secondsLeftTotal+secondsRightTotal); 
secondsTotal = secondsLeftTotal + secondsRightTotal; 
percentLeft = (secondsLeftTotal)./secondsTotal*100; 
percentRight = (secondsRightTotal)./secondsTotal*100; 

summaryStats = table(secondsLeftTotal, secondsRightTotal, discrIndex, secondsTotal, percentLeft, percentRight); 

%% Functions
% get filepaths of interest in a folder, based on filetype 
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