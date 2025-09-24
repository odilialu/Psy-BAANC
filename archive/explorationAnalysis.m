%% General use - object exploration analysis
% Odilia Lu, edited 07/30/24

%{ 

This script quantifies mouse exploration of various objects, including
those used in object exploration tasks and social interaction tests. 

An animal is considered exploring the object if the nose, but not body, is within a defined zone; 

Instructions: 
1. Predict nose and body markers using DeepLabCut. 
2. Then, change parameters in the variables to change section below.
3. Run the script, and follow instructions that are displayed in the command prompt to define exploration and exclusion zones. 
4. View summaryStats for relevant quantifications. Each row is an animal in the folder (alphabetically organized). 

For social interaction test, roiOne should mark the empty cup, roiTwo should mark the social cup. 

For novel object recognition, roiOne should mark the familiar object, roiTwo should mark the novel object. 

For novel object exploration, roiOne should mark the left object, roiTwo should mark the right object. 

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

objectShapeOne = "circle"; % options = "ellipse", "circle", "rectangle", "polygon" 
objectShapeTwo = "circle"; 

experiment = "SIT_Pers"; % Name your experiment, used to retrieve saved ROIs. 

%% loop through all the CSV and AVI files in the folder. 
strVids = strArray(path, videoType); 
strFiles = strArray(path, 'csv'); 

%% pre-allocate variables. 
startExp = framerate*60*start_min+1; 
lengthExp = framerate*60*end_min;

noseTimeMaskOne = NaN(length(strFiles), 1); 
noseTimeMaskTwo = NaN(length(strFiles), 1); 

MaskOneExcludeTime = NaN(length(strFiles), 1); 
MaskTwoExcludeTime = NaN(length(strFiles), 1); 

%% main analysis
strSplit = split(strVids, "/"); 

for i = 1:length(strVids)

    identifier = strSplit(i, end); 

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

    % draw masks for object zones, if not done already. 
    coordMaskOne = getMaskCoordinates("-roiOne-", identifier, experiment, objectShapeOne, ax); 
    coordMaskTwo = getMaskCoordinates("-roiTwo-", identifier, experiment, objectShapeTwo, ax); 
    coordMaskOneExclude = getMaskCoordinates("-roiOneExclude-", identifier, experiment, objectShapeOne, ax); 
    coordMaskTwoExclude = getMaskCoordinates("-roiTwoExclude-", identifier, experiment, objectShapeTwo, ax); 

    % find when nose coordinate is within each of the masks
    noseMaskOne = (find(ismember(coordinatesNose, coordMaskOne, 'rows')));
    noseMaskTwo = (find(ismember(coordinatesNose, coordMaskTwo, 'rows')));

    noseTimeMaskOne(i, 1) = (length(noseMaskOne))/framerate; 
    noseTimeMaskTwo(i, 1) = (length(noseMaskTwo))/framerate;

    % mouse on top of object exclude
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

end

%% final results - secondsLeftTotal, secondsRightTotal, secondsTotal
explorationTimeOne = noseTimeMaskOne - MaskOneExcludeTime; 
explorationTimeTwo = noseTimeMaskTwo - MaskTwoExcludeTime; 
explorationTimeIndex = (explorationTimeTwo-explorationTimeOne)./(explorationTimeOne + explorationTimeTwo); 
totalExploration = explorationTimeOne + explorationTimeTwo; 

summaryStats = table(explorationTimeOne, explorationTimeTwo, explorationTimeIndex, totalExploration); 

%% Functions
%%%%%%%%% 1.  get filepaths of interest in a folder, based on filetype 
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


%%%%%%%%% 2. Create masks and get coordinates of mask location. 
function coordMask = getMaskCoordinates(roiName, identifier, experiment, objectShape, ax)
    roiString = strcat(experiment, roiName, identifier, '.mat'); 

    if isfile(roiString)
        roi = importdata(roiString);
        if objectShape == "ellipse"
            roi = drawellipse('Center', roi.Center, 'SemiAxes', roi.SemiAxes, 'RotationAngle', roiOne.RotationAngle);
        elseif objectShape == "circle" 
            roi = drawcircle('Center', roi.Center, 'Radius', roi.Radius); 
        elseif objectShape == "rectangle" 
            roi = drawrectangle('AspectRatio', roi.AspectRatio, 'Position', roi.Position);
        elseif objectShape == "polygon"
            roi = drawpolygon("Position", roi.Position); 
        end


    else
        disp("1) Draw shape to define zone for " + roiName + ". 2) Hit enter when satisfied with shape.")
        if objectShape == "ellipse" 
            roi = drawellipse; 
        elseif objectShape == "circle" 
            roi = drawcircle; 
        elseif objectShape == "rectangle"
            roi = drawrectangle; 
        elseif objectShape == "polygon"
            roi = drawpolygon;
        end
        pause; 

    end
    
    save(roiString, 'roi');  
    mask = roi.createMask(ax);
    [yCoordMask, xCoordMask] = find(mask); 
    coordMask = cat(2, xCoordMask, yCoordMask);

end
