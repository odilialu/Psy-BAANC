%% Crop and cut videos 
% Odilia Lu 
% Last edited June 28, 2024
% Crops and cuts all videos in the folder designated in dirString
% Writes to your current directory 

% variables to change
dirString = "Y:\paperExperiments\social CPP\SkinnerBoxExperiment\052424_post"; % Insert path of folder with videos that need to be cropped.
videoType = 'mkv'; % current video format (i.e., mp4, avi, wmv, etc.)
framerate = 30; % desired fps of output video
start_min = 0; % for trimming video, set start time, in min (If no trim to beginning, set start_min = 0). 
end_min = 30; % for trimming video, how long do you want the video to be, in min? 
cropConsistent = 1; % do you want to crop all videos in the folder consistently (1 = yes, 0 = no; crop each video individually). 

%% generate a string array of the video paths 

sprintfDir = strcat(dirString, "\%s");
filesDirectory = dir(dirString); 

listVids = strings(length(filesDirectory), 1); 
% only get files with correct video type
for i = 1:length(filesDirectory)
    if endsWith(filesDirectory(i).name, videoType) == 1
        listVids(i, 1) = filesDirectory(i).name; 
    end
end
listVids(cellfun('isempty', listVids)) = []; 

% create a string array with paths to vids.
strVids = strings([length(listVids), 1]); 
for i = 1:length(listVids)
    fileVid = listVids(i, 1); 
    strVids(i, 1) = sprintf(strrep(sprintfDir, '\', '/'), fileVid);
end 

%% get filenames only, used for naming new vid 
filename = split(strVids, "/");

if length(strVids) == 1
    filename = split(filename(end, 1), ".");
    filename = filename(:, 1); 
else
    filename = split(filename(:, end), ".");
end

%% cropping and cutting 

% set cutting parameters
startExp = framerate*start_min+1;
lengthExp = framerate*end_min*60; 

% set cropping parameters
if cropConsistent == 1
    video = VideoReader(strVids(1, 1)); 
    frame = read(video, 1);
    imshow(frame)
    rectCoords = getrect;
end

% main
for i = 1:length(strVids)

    video = VideoReader(strVids(i, 1)); 

    if cropConsistent == 0
        frame = read(video, 1);
        imshow(frame)
        rectCoords = getrect;
    end

    writer = VideoWriter(filename(i, 1)); 
    writer.FrameRate = framerate;
    
    open(writer); 

    for j = startExp:lengthExp % define what part of vid to keep 
        frame = read(video, j); % read frame
        frame = imcrop(frame, rectCoords); % does the cropping
        writeVideo(writer, frame); 
    end

    close(writer); 

end
