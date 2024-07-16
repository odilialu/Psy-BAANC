%% Crop and cut videos 
% Christine Liu GUI for Odilia Lu videoCrop.m script 
% Last edited July 16,2024
% Writes to your current directory 

clear all

%% User input settings 
inputs = inputdlg({"Video Type: (i.e. mp4, avi, wmv, etc.)", "Frame Rate:", "Minutes to Analyze:", "Crop videos individually (0) or all at once (1)?"},"Settings",[1 40; 1 40; 1 40; 1 40],{'avi','10','10','0'});
videoType = inputs{1};
framerate = str2double(inputs{2});
end_min = str2double(inputs{3});
cropConsistent = str2double(inputs{4});

dirString = string(uigetdir('C:\',"Load folder containing videos"));

%% Get videos 
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

% set cropping parameters
if cropConsistent == 1
    video = VideoReader(strVids(1, 1)); 
    frame = read(video, 1);
    imshow(frame)
    title('Crop all videos in dir at once');
    rectCoords = getrect;
    answer = inputdlg({'Time start (s):', 'Add suffix to filename:'});
    start_s = answer{1};
    suffix = answer{2};
    
    f = waitbar(0,'Cropping Videos...');

    for i = 1:length(strVids)
        video = VideoReader(strVids(i, 1)); 
        writer = VideoWriter(strcat(filename(i, 1), '_',suffix)); 
        writer.FrameRate = framerate;
        
        open(writer); 
        waitbar(i/length(strVids),f,'Cropping Videos...')
        startExp = framerate*start_s+1;
        lengthExp = (framerate*end_min*60)+startExp; 
    
        for j = startExp:lengthExp % define what part of vid to keep 
            frame = read(video, j); % read frame
            frame = imcrop(frame, rectCoords); % does the cropping
            writeVideo(writer, frame); 
        end
    
        close(writer); 
    end
    close(f);
    d = msgbox("Done Cropping");
end
% main


if cropConsistent == 0
    [indx, tf] = listdlg('ListString',strVids,'SelectionMode','single','ListSize',[800,300]);

    while tf == 1
        video = VideoReader(strVids(indx, 1));

        frame = read(video, 1);
        imshow(frame)
        title(filename(indx,1), 'Interpreter','none');
        answer = inputdlg({'Enter time start (s):', 'Filename suffix:'});
        if isempty(answer)
            [indx, ~] = listdlg('ListString',strVids,'SelectionMode','single','ListSize',[800,300]);
            frame = read(video, 1);
            imshow(frame)
            title(filename(indx,1), 'Interpreter','none');
            answer = inputdlg({'Enter time start (s):', 'Filename suffix:'});
        end

        start_s = str2double(answer{1});
        suffix = answer{2};
        title(strcat(filename(indx, 1), '_',suffix), 'Interpreter','none');
        rectCoords = getrect;

        writer = VideoWriter(strcat(filename(indx, 1), '_',suffix)); 
        writer.FrameRate = framerate;

        open(writer); 
        f = waitbar(.5,'Cropping Video...');
        
        startExp = framerate*start_s+1;
        lengthExp = (framerate*end_min*60)+startExp; 
        
        for j = startExp:lengthExp % define what part of vid to keep 
            frame = read(video, j); % read frame
            frame = imcrop(frame, rectCoords); % does the cropping
            writeVideo(writer, frame); 
        end
        
        close(writer); 
        close(f);
        [indx, tf] = listdlg('ListString',strVids,'SelectionMode','single','ListSize',[800,300],'InitialValue',indx);
    end

    if tf == 0
        close all
    end
end
