%% fear conditioning freezing analysis 
%{ 
Odilia Lu. 
Last edited 11.25.23

1) Move all videos you want to analyze to a folder on your computer.
Specify this folder's path into the path variable below. Each video should have a unique
file name. 

2) Change the variables in the "Variables to edit" section as necessary. 

3) Run the script. The portion requiring user input is done at the
beginning. Specifically, you must 1) draw a mask around a region of the image
where you can detect the cue light going on and off, and  2) type in the 
command prompt a y-cutoff value that determines cue frames. 
Directions are also explained in the command window. You should only have
to do this once per video as relevant variables are saved to your current
working directory. 

4) Analysis should then run automatically. Please check the quality of freezing
graphs and re-run, adjusting parameters as necessary. You can also check
quality by uncommenting the imshow lines.

Main outputs of analysis are located in the 'summaryFreezing' structure.
For arrays in the structure, each row is a different, successive cue. 
Each column is a different animal (goes by alphabetical order of filenames). 
- DuringCue = freezing during the cue 
- PreCue = freezing in ITI just before cue, time is "timeWithoutCue" 
- Habituation = freezing during habituation 

%}

%% Variables to edit 
path = "Y:\fearConditioningExp\retrieval\males\JL\231102\mp4"; 

videoType = 'mp4';
nCues = 2; 
lengthCue = 30; % in seconds 
lengthHabituation = 300; % in seconds (300 sec for encoding, 180 sec for retrieval & extinction)
framerate = 30; 
framerateNew = 3; % FreezeFrame samples at ~ 3 Hz
timeWithoutCue = 10; % (in seconds, for plots, how much time before and after do you want to show?) 

encoding = 0; % Is this an encoding session? 1 (yes) or 0 (no)

baselineTime = 300; % in seconds
flickerThresh = 15; % how much does each pixel tend to change from frame to frame? (flicker level) 
thresholdMove = 2; % number of seconds that animal must be continuously immobile to be considered freezing
thresholdPixels = 10; % number of pixels that have to be different between frames to be considered 'moving' 

%% loop through all the video files in the folder to get a list of strings of vid paths. 
directoryVids = dir(path);
sprintfDirVids = strcat(path, "\%s"); 
[~, strVids] = strArray(path, sprintfDirVids, videoType); 

%% preallocate variables
prctFreezing = NaN(nCues, length(strVids));
prctFreezingPreCue = NaN(nCues, length(strVids)); 
prctFreezingBaseline = NaN(1, length(strVids)); 

%% Use houselight ON to determine when the cue happens. 

% get masks for where the houselight is if you haven't already 
for i = 1:length(strVids)
    video = VideoReader(strVids(i, 1)); 
    subjectID = split(strVids(i, 1), "/");  
    subjectID = subjectID(end); 

    houselightMaskString = strcat('houselightMask_', subjectID, '.mat');
    if isfile(houselightMaskString) == 0
        frame = read(video, 1);
        figure()
        imshow(frame)
        disp("draw a rectangle mask around a region of the image where you" + ...
            " can detect change in brightness due to cue light going on and off")
        houselightRect = drawrectangle;
        houselightMask = houselightRect.createMask; 
        clear maskIndex
        [maskIndex(:, 1), maskIndex(:, 2)] = find(houselightMask == 1, 50, 'first'); 
        save(houselightMaskString, 'maskIndex')
    end
end

% get frames of when cue is on if you haven't already
for i = 1:length(strVids) 
    subjectID = split(strVids(i, 1), "/");  
    subjectID = subjectID(end); 

    % define where houselight is 
    houselightMaskString = strcat('houselightMask_', subjectID, '.mat');
    if isfile(houselightMaskString)
        clear maskIndex
        maskIndex = importdata(houselightMaskString);
    else
        error('Please define mask zone first')
    end

    % use houselight to define when cue is on
    cueFramesString = strcat('cueFrames_', subjectID, '.mat');
    if isfile(cueFramesString) == 0
        video = VideoReader(strVids(i, 1)); 
        startFrame = (lengthHabituation-10)*framerate; % skip habituation frames, slows down analysis.
        houselightBrightness =  nan(video.NumFrames-startFrame+1, 1); % preallocate vector for speed.
        idx = 0;
        for j = startFrame:video.NumFrames
            idx = idx+1;
            houselightBrightFrame = 0;
            frame = rgb2gray(read(video, j));
            for k = 1:length(maskIndex)
                houselightBrightFrame = houselightBrightFrame + double(frame(maskIndex(k, 1), maskIndex(k, 2)));
            end
            houselightBrightness(idx, 1) = abs(houselightBrightFrame);
        end
        figure()
        plot(houselightBrightness)

        prompt = "At what y-value can you distinguish between houselight on vs. off?";
        houselightBrightnessCutoff = input(prompt);

        houselightONFrames = find(houselightBrightness > houselightBrightnessCutoff);
        houselightONFrames = houselightONFrames + startFrame;

        if length(houselightONFrames)/framerate > (nCues*lengthCue+1) && length(houselightONFrames)/framerate < (nCues*lengthCue-1)
            error('Imporoper determination of houselightONFrames');
        end
        save(cueFramesString, 'houselightONFrames')
    end
end

close all

%% main analysis
for i = 1:length(strVids)
    video = VideoReader(strVids(i, 1)); 
    time = string(datetime);
    disp("analyzing "+ strVids(i, 1) + " " + time)

    %%%%%%%%%%% Determine cue times %%%%%%%%%%%%%%%%%%%%%
    % import frames when cue comes on.
    subjectID = split(strVids(i, 1), "/");
    subjectID = subjectID(end);
    cueFramesString = strcat('cueFrames_', subjectID, '.mat');
    houselightONFrames = importdata(cueFramesString);

    % separate houselightONFrames into distinct cues
    for j = 1:length(houselightONFrames)-1
        if houselightONFrames(j+1) - houselightONFrames(j) > 30*framerate
            houselightONFrames(j+1) = nan;
        end
    end
    index=find(~isnan(houselightONFrames));
    idx=find(diff(index)~=1);
    A=[idx(1);diff(idx);numel(index)-idx(end)];
    cueTimes=mat2cell(houselightONFrames(~isnan(houselightONFrames)),A,1);

    % Remove noise.
    if length(cueTimes)~= nCues % if the correct number of cues was not detected
        for j = 1:length(cueTimes)
            if length(cueTimes{j})<(lengthCue*framerate-framerate)
                cueTimes{j}= []; % remove cue times where the length of cue is less than how long a cue should be
            end
        end
    end
    cueTimes = cueTimes(~cellfun('isempty',cueTimes));

    % Throw an error if there is still not the correct number of cues
    % detected.
    if length(cueTimes)~= nCues
        error('Improper determination of number of cues')
    end

    %%%%%%%%%%%% GET FRAMES OF INTEREST %%%%%%%%%%%%%%%%%%

    % frames of interest during cue + and - timeWithoutCue
    framesOfInterest = nan(nCues, (timeWithoutCue*2+lengthCue)*framerate+1);
    for j = 1:length(cueTimes)
        framesOfInterestStart = cueTimes{j}(1) - (framerate*timeWithoutCue);
        framesOfInterestEnd = cueTimes{j}(1) + (framerate*(lengthCue+timeWithoutCue));
        framesOfInterest (j, :) = framesOfInterestStart:framesOfInterestEnd;
    end

    % frames of interest to calculate baseline freezing levels during
    % habituation
    framesOfInterestBaseline = 1:(framerate/framerateNew):(baselineTime-timeWithoutCue)*framerate;

    %%%%%%%%% NUMBER DIFF PIXELS DURING CUE %%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the number of different pixels for frames of interest (CUE)
    framesOfInterest = framesOfInterest(:, 1:(framerate/framerateNew):length(framesOfInterest));
    nDiffPixels = nan(nCues, length(framesOfInterest)-1, 1);
    for j = 1:length(cueTimes)
        idx = 0;
        for k = framesOfInterest(j, 1:end-1)
            idx = idx+1;
            frameFirst = imresize(rgb2gray(read(video, k)), 0.25);
            frameNext = imresize(rgb2gray(read(video, k+(framerate/framerateNew))), 0.25);
            frameDiffTemp = abs(frameNext-frameFirst);
            frameDiffTemp(frameDiffTemp<flickerThresh) =0;
            nDiffPixels(j, idx) = nnz(frameDiffTemp);
%             imshow(frameDiffTemp) % uncomment if you want to see moving pixels
        end
    end

    %%%%%%%%%%%%% NUMBER DIFF PIXELS DURING BASELINE %%%%%%%%%%%%%%%%%%%%%%%%%
    nDiffPixelsBaseline = nan(1, length(framesOfInterestBaseline)-1);
    idx = 0;
    for k = framesOfInterestBaseline(1:end-1)
        idx =  idx + 1;
        frameFirst = imresize(rgb2gray(read(video, k)), 0.25);
        frameNext = imresize(rgb2gray(read(video, k+(framerate/framerateNew))), 0.25);
        frameDiffTemp = abs(frameNext-frameFirst);
        frameDiffTemp(frameDiffTemp<flickerThresh) =0;
        nDiffPixelsBaseline(1, idx) = nnz(frameDiffTemp);
        % imshow(frameDiffTemp) % uncomment if you want to see moving pixels
    end

    %%%%%%%%%% PLOTS USED FOR QUALITY CONTROL CHECK OF THRESHOLD VALUES %%%%%%%%%%%%%%%%

    if nCues == 10
        figure()
        figIdx = 0;
        for l = 1:nCues
            figIdx = figIdx +1;
            subplot(2, nCues/2, figIdx)
            plot(nDiffPixels(l, :))
            ylim([-50, 1000])
            xline(timeWithoutCue*framerateNew, 'r')
            xline((timeWithoutCue+lengthCue)*framerateNew, 'r')
        end
    else
        figure()
        figIdx = 1;
        subplot(1, nCues+1, figIdx)
        plot(nDiffPixelsBaseline)
        title('baseline')
        ylim([-50, 4000])
        
        for l = 1:nCues
            figIdx = figIdx+1;
            subplot(1, nCues+1, figIdx)
            plot(nDiffPixels(l, :))
            ylim([-50, 4000])
            xline(timeWithoutCue*framerateNew, 'r')
            if encoding == 1
                xline((timeWithoutCue+lengthCue-2)*framerateNew, 'r')
            end
            xline((timeWithoutCue+lengthCue)*framerateNew, 'r')
            title('cue_' + string(l))
        end
    end

    %%%%%%%%% APPLY THRESHOLDS AND CALCULATE FREEZING (CUES) %%%%%%%%%
    % Apply threshold for what number of pixels can be different to still
    % be considered freezing
    nDiffPixels(nDiffPixels<thresholdPixels) = 0;

    % Mouse is freezing if it has remained still for thresholdMove seconds.
    freezeArray = zeros(size(nDiffPixels, 1), size(nDiffPixels, 2));
    for j = 1:size(nDiffPixels, 1)
        for k = 1:size(nDiffPixels, 2) - framerateNew*thresholdMove
            if sum(nDiffPixels(j, k:k+(framerateNew*thresholdMove))) == 0
                freezeArray(j, k:k+(framerateNew*thresholdMove)) = 1;
            end
        end
    end

    % freezing during cues
    for j = 1:nCues
        prctFreezing(j, i) = (sum(freezeArray(j, timeWithoutCue*framerateNew+1:(timeWithoutCue+lengthCue)*framerateNew)))/framerateNew/lengthCue*100;
        prctFreezingPreCue(j, i) = (sum(freezeArray(j, 1:timeWithoutCue*framerateNew)))/framerateNew/lengthCue*100;
    end

    %%%%%%%%% APPLY THRESHOLDS AND CALCULATE FREEZING (BASELINE) %%%%%%%
    % Apply threshold for what number of pixels can be different to still
    % be considered freezing
    nDiffPixelsBaseline(nDiffPixelsBaseline<thresholdPixels) = 0;

    % Mouse is freezing if it has remained still for thresholdMove seconds.
    freezeArrayBaseline = zeros(size(nDiffPixelsBaseline, 1), size(nDiffPixelsBaseline, 2));
    for k = 1:size(nDiffPixelsBaseline, 2) - framerateNew*thresholdMove
        if sum(nDiffPixelsBaseline(1, k:k+(framerateNew*thresholdMove))) == 0
            freezeArrayBaseline(1, k:k+(framerateNew*thresholdMove)) = 1;
        end
    end

    prctFreezingBaseline(1, i) = (sum(freezeArrayBaseline, 'all'))/framerateNew/baselineTime*100;
end

summaryFreezing.DuringCue = prctFreezing; 
summaryFreezing.PreCue = prctFreezingPreCue; 
summaryFreezing.Habituation = prctFreezingBaseline; 