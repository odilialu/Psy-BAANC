%% Social interaction test analysis
% Odilia Lu, edited 07/22/24

%{ 

This script quantifies mouse sociability in the 3-chamber social interaction test in two ways: 
[1] based on when the nose, but not body, is within a social interaction cup zone; 
[2] based on when the body is in the empty vs. social chambers. 

Instructions: 
1. Predict nose and body markers using DeepLabCut. 
2. Then, change parameters in the variables to change section below.
3. Run the script, and follow instructions that are displayed in the command prompt to define non-social and social interaction zones / chambers. 
4. View summaryStats for relevant quantifications. Each row is an animal in the folder (alphabetically organized). 

%}

%% variables to change
path = "F:\SM_SIT\230203 SIT Pers\cropped\analyze"; % path with all videos and DLC files

videoType = 'avi';
framerate = 30; %30 fps 
start_min = 0;
end_min = 10;

pCutoff = 0.95; % Define pCutoff based on DLC file. This excludes coordinates below this value that may not be in the frame. 
noseX_columnNum = 2; % which column of the DLC data file has the nose X coordinates? 
bodyX_columnNum = 5; % which column of the DLC data file has the body X coordinates? 

excludeSize = 60; % defines the size of the empty / social cup zones that should be excluded because it implies the animal is on top of the object, not exploring it. 
semiAxes = [104, 95]; % defines ellipse dimensions for masks around social and empty cup. Adjust based on visual inspection of figures. 

%% loop through all the CSV and AVI files in the folder. 
strVids = strArray(path, videoType); 
strFiles = strArray(path, 'csv'); 

%% pre-allocate variables. 
startExp = framerate*60*start_min+1; 
lengthExp = framerate*60*end_min;

zoneX = cell(1, length(strFiles)); 
zoneY = cell(1, length(strFiles)); 

zoneTimeOne = NaN(length(strFiles), 1); 
zoneTimeTwo = NaN(length(strFiles), 1); 

noseTimeMaskOne = NaN(length(strFiles), 1); 
noseTimeMaskTwo = NaN(length(strFiles), 1); 

MaskOneExcludeTime = NaN(length(strFiles), 1); 
MaskTwoExcludeTime = NaN(length(strFiles), 1); 

%% main analysis

for i = 1:length(strVids)

    figure()
    video = VideoReader(strVids(i, 1)); 
    frame = read(video, 1); % extract the first frame of every video
    ax = imshow(frame);

    hold on

    % read in the nose and body coordinates. 
    data = readmatrix(strFiles(i,1)); 
    data = data(startExp:lengthExp, :); 

    coordinatesNose = NaN(length(data), 2); 
    coordinatesBody = NaN(length(data), 2); 
    for k = 1:length(data)
        if data(k, 4) > pCutoff
            coordinatesNose(k, 1) = round(data(k, noseX_columnNum)); 
            coordinatesNose(k, 2) = round(data(k, noseX_columnNum+1)); 
        end

        if data(k, 7) > pCutoff
            coordinatesBody(k, 1) = round(data(k, bodyX_columnNum)); 
            coordinatesBody(k, 2) = round(data(k, bodyX_columnNum+1));
        end
    end

    plot(coordinatesNose(:, 1), coordinatesNose(:, 2))

    % get masks for familiar and novel cups, if not defined already.  
    strSplit = split(strVids, "/"); 
    socialOneString = strcat('socialOneShape-', strSplit(i, end), '.mat'); 
    socialTwoString = strcat('socialTwoShape-', strSplit(i, end), '.mat'); 

    if isfile(socialOneString)
        roiOne = importdata(socialOneString);
        roiOne = drawellipse('Center', roiOne.Center, 'SemiAxes', roiOne.SemiAxes, 'Color', 'r', 'RotationAngle', roiOne.RotationAngle);
        maskOne = roiOne.createMask(ax); 
        roiTwo = importdata(socialTwoString); 
        roiTwo = drawellipse('Center', roiTwo.Center, 'SemiAxes', roiTwo.SemiAxes, 'Color', 'r', 'RotationAngle', roiTwo.RotationAngle);
        maskTwo = roiTwo.createMask(ax); 
    else
    fprintf(['[Define empty cup zone] \n' ...
        '1) In the pop up figure, left click on the center of the empty cup zone. Hit enter. \n' ... 
        '2) Then, edit the ellipse to define the empty cup zone. Hit enter.'])
    [x, y] = getpts; 
    pointOne = cat(2, x, y); 
    roiOne = drawellipse('Center', pointOne, 'SemiAxes', semiAxes, 'Color', 'r'); 
    pause; 
    save(socialOneString, 'roiOne'); 
    maskOne = roiOne.createMask(ax); 

    fprintf(['[Define social cup zone] \n' ...
        '1) In the pop up figure, left click on the center of the social cup zone. Hit enter. \n' ... 
        '2) Then, edit the ellipse to define the social cup zone. Hit enter.'])
    [x2, y2] = getpts; 
    pointTwo = cat(2, x2, y2); 
    roiTwo = drawellipse('Center', pointTwo, 'SemiAxes', semiAxes, 'Color','r'); 
    pause; 
    save(socialTwoString, 'roiTwo'); 
    maskTwo = roiTwo.createMask(ax);
    end
    
    [yCoordMaskOne, xCoordMaskOne] = find(maskOne); 
    coordMaskOne = cat(2, xCoordMaskOne, yCoordMaskOne);

    [yCoordMaskTwo, xCoordMaskTwo] = find(maskTwo); 
    coordMaskTwo = cat(2, xCoordMaskTwo, yCoordMaskTwo); 

    % find when nose coordinate is within each of the masks
    noseMaskOne = (find(ismember(coordinatesNose, coordMaskOne, 'rows')));
    noseMaskTwo = (find(ismember(coordinatesNose, coordMaskTwo, 'rows')));

    noseTimeMaskOne(i, 1) = (length(noseMaskOne))/framerate; 
    noseTimeMaskTwo(i, 1) = (length(noseMaskTwo))/framerate;

    % mouse on top of object exclude
    roiOneExclude = drawellipse('Center', roiOne.Center, 'SemiAxes', (roiOne.SemiAxes-excludeSize), 'Color', 'g', 'RotationAngle', roiOne.RotationAngle);
    maskOneExclude = roiOneExclude.createMask(ax);
    [yCoordMaskOneExclude, xCoordMaskOneExclude] = find(maskOneExclude); 
    coordMaskOneExclude = cat(2, xCoordMaskOneExclude, yCoordMaskOneExclude);

    roiTwoExclude = drawellipse('Center', roiTwo.Center, 'SemiAxes', (roiTwo.SemiAxes-excludeSize), 'Color', 'g', 'RotationAngle', roiTwo.RotationAngle); 
    maskTwoExclude = roiTwoExclude.createMask(ax); 
    [yCoordMaskTwoExclude, xCoordMaskTwoExclude] = find(maskTwoExclude); 
    coordMaskTwoExclude = cat(2, xCoordMaskTwoExclude, yCoordMaskTwoExclude); 

    bodyMaskOneExclude = find(ismember(coordinatesBody, coordMaskOneExclude, 'rows')); 
    noseMaskOneExclude = find(ismember(coordinatesNose, coordMaskOneExclude, 'rows')); 
    MaskOneExclude = cat(1, bodyMaskOneExclude, noseMaskOneExclude); 
    MaskOneExclude = unique(MaskOneExclude); 
    MaskOneExcludeTime(i, 1) = length(find(ismember(MaskOneExclude, noseMaskOne)))/framerate; 

    bodyMaskTwoExclude = find(ismember(coordinatesBody, coordMaskTwoExclude, 'rows')); 
    noseMaskTwoExclude = find(ismember(coordinatesNose, coordMaskTwoExclude, 'rows')); 
    MaskTwoExclude = cat(1, bodyMaskTwoExclude, noseMaskTwoExclude); 
    MaskTwoExclude = unique(MaskTwoExclude); 
    MaskTwoExcludeTime(i, 1) = length(find(ismember(MaskTwoExclude, noseMaskTwo)))/framerate; 

    %%%%%%%%%%%%% find x coordinates marking the edge of the zones %%%%%%%%%%%%%%%% 
    pathSplit = split(path, ":"); 
    socialZoneXString = strcat('socialZoneXCoord-', pathSplit(end, 1), '.mat'); 
    socialZoneXString = strrep(socialZoneXString, "\", "-"); 
    socialZoneYString = strcat('socialZoneYCoord-', pathSplit(end, 1), '.mat'); 
    socialZoneYString = strrep(socialZoneYString, "\", "-"); 

    if isfile(socialZoneXString)
        zoneX = importdata(socialZoneXString);
        zoneY = importdata(socialZoneYString); 
    else
        fprintf(['[Define empty and social chambers] \n' ...
            '1. Left click to designate the X-coordinate that delineates the empty chamber. \n' ...
            '2. Left click to designate the X-coordinate that delineates the social chamber. \n' ...
            '3. Hit enter. '])
        [zoneX{1, i}, zoneY{1, i}] = getpts; 
    end

    % time exploring zone is defined as when body coordinate is within that
    % zone. 
    if zoneX{1, i}(1, 1) < zoneX{1, i}(2, 1)
        zoneTimeOne(i,1) = length(find(coordinatesBody(:, 1) < zoneX{1, i}(1, 1)))/framerate; 
        zoneTimeTwo(i, 1) = length(find(coordinatesBody(:, 1) > zoneX{1, i}(2, 1)))/framerate; 
    else 
        zoneTimeOne(i, 1) = length(find(coordinatesBody(:, 1) > zoneX{1, i}(1, 1)))/framerate; 
        zoneTimeTwo(i, 1) = length(find(coordinatesBody(:, 1) < zoneX{1, i}(2, 1)))/framerate; 
    end

    % visualization
    xline(zoneX{1, i}(1, 1))
    xline(zoneX{1, i}(2, 1))

end

    save(socialZoneXString, 'zoneX')
    save(socialZoneYString, 'zoneY')

%% final results - secondsLeftTotal, secondsRightTotal, secondsTotal
explorationTimeOne = noseTimeMaskOne - MaskOneExcludeTime; 
explorationTimeTwo = noseTimeMaskTwo - MaskTwoExcludeTime; 
explorationTimeIndex = (explorationTimeTwo-explorationTimeOne)./(explorationTimeOne + explorationTimeTwo); 

zoneIndex = (zoneTimeTwo-zoneTimeOne)./(zoneTimeTwo+zoneTimeOne); 

summaryStats = table(explorationTimeOne, explorationTimeTwo, explorationTimeIndex, zoneTimeOne, zoneTimeTwo, zoneIndex); 


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