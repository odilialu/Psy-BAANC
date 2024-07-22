%% Open field test analysis 
% Odilia Lu 
% Last edited December 5, 2022 5PM

%% variables to change 
path = "D:\CMS_Exp\Males\NOR\cropped\stress_D2"; 

videoType = 'avi';
framerate = 30; %30 fps 
start_min = 0; 
end_min = 7; 
lengthCM = 50; % length/width of open field (in cm) 
fractionCenterZone = 0.75; % what fraction of the open field should be designated the center? 
bodyX_columnNum = 17; % for OFT = 26; for NOE = 17; 

%% get paths of all videos and csv files in the folder. 
strVids = strArray(path, videoType); 
strFiles = strArray(path, 'csv');

%% used to define zones and pixel to cm conversion: 

% look to see if zones have already been defined, if not, define zones
pathSplit = split(path, ":");
rectCoordsString = strcat('rectCoords-', pathSplit(end, 1), '.mat');
rectCoordsString = strrep(rectCoordsString, "\", "-");
if isfile(rectCoordsString)
    rectCoords = importdata(rectCoordsString);
else
    video = VideoReader(strVids(1, 1));
    frame = read(video, 1);
    imshow(frame)
    disp("click and drag to define the perimeter of the open field box base")
    rectCoords = getrect; % define the area in which the mouse moves.
    save(rectCoordsString, 'rectCoords')
end

%define zones
scaler = (1-fractionCenterZone)/2;
insideCoords = [rectCoords(1)+(scaler*rectCoords(3)),  rectCoords(2)+ (scaler*rectCoords(4)), fractionCenterZone*rectCoords(3), fractionCenterZone*rectCoords(4)]; % defines the center zone 
% pixel to centimeter conversion 
lengthPixels = (rectCoords(3) + rectCoords(4))/2; 

%% Analysis 
startExp = framerate*60*start_min+1; 
lengthExp = framerate*60*end_min; 

% preallocate arrays 
distanceArray = cell(length(strFiles), 1); 
distancePerSecond = cell(length(strFiles), 1); 

distanceTotal = NaN(length(strFiles), 1); 
averageVelocity = NaN(length(strFiles), 1); 
timeMoving = NaN(length(strFiles), 1); 
timeRunning = NaN(length(strFiles), 1); 
timeWalking = NaN(length(strFiles), 1); 
timeFreezing = NaN(length(strFiles), 1); 
timeInside = NaN(length(strFiles), 1); 
velocityWhileMoving = NaN(length(strFiles), 1); 

for i = 1:length(strVids)
    % read in body coordinates
    data = readmatrix(strFiles(i,1)); 
    data = data(startExp:lengthExp, :); 

    body2 = nan(length(data), 2);
    body2(:, 1) = data(:, bodyX_columnNum);
    body2(:, 2) = data(:, bodyX_columnNum+1); 

    % visualization of track and inner zone against video
    video = VideoReader(strVids(i, 1)); 
    frame = read(video, 1); % extract the first frame of every video
    figure()
    imshow(frame)
    hold on
    plot(body2(:, 1), body2(:, 2))
    rectangle('Position', rectCoords)
    rectangle('Position', insideCoords)

    % distance travelled and average velocity 
    diffBody2 = diff(body2); 
    distanceArray{i, 1} = hypot(diffBody2(:, 1), diffBody2(:, 2))/lengthPixels*lengthCM; 
    distanceTotal(i, 1) = sum(distanceArray{i, 1}); 
    averageVelocity(i, 1) = distanceTotal(i, 1)/(end_min-start_min)/60; 

    % determine velocity over time (cm/s) (1 second resolution)
    n = 0; 
    for j = 1:framerate:(length(distanceArray{i, 1})-framerate)
        n = n+1; % seconds counter
        distancePerSecond{i, 1}(n) = sum(distanceArray{i, 1}(j:j+(framerate-1), 1)); 
    end
 
    % find time spent doing the following activities (freezing < 5 cm/s, walking = 5-20 cm/s, running > 20 cm/s) 
    timeMoving(i, 1) = sum(distancePerSecond{i, 1}> 5); 
    timeRunning(i, 1) = sum(distancePerSecond{i, 1}>20); 
    timeWalking(i, 1) = timeMoving(i, 1) - timeRunning(i, 1);  
    timeFreezing(i, 1) = sum(distancePerSecond{i, 1}<5);

    % find average velocity during time moving 
    velocityWhileMoving(i, 1) = mean(distancePerSecond{i, 1}(distancePerSecond{i, 1} > 5)); 

    % find time spent in the center. 
    xInsideBottom = find(body2(:, 1)>insideCoords(1));
    xInsideTop = find(body2(:, 1)<(insideCoords(1)+insideCoords(3))); 
    xInside = intersect(xInsideBottom, xInsideTop); 
    yInsideLeft = find(body2(:, 2)>insideCoords(2));
    yInsideRight = find(body2(:, 2)<(insideCoords(2)+insideCoords(4))); 
    yInside = intersect(yInsideLeft, yInsideRight);
    body2Inside = intersect(xInside, yInside); 

    timeInside(i, 1) = length(body2Inside)/framerate; % time in seconds. 

end

summaryStats = table(averageVelocity, distanceTotal, timeInside, timeFreezing, timeWalking, timeRunning, timeMoving, velocityWhileMoving); 

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
