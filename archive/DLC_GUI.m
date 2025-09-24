% DeepLabCut Coordinate Analysis GUI
% Christine Liu based off Odilia Lu NOE_test_240701.m script
% Last updated 9/8/24

%{ 
This script obtains all bodyparts from a DLC csv output file and prompts
the user to select relevant bodyparts for analysis. 

Zones can be defined where select bodyparts must be within the ROI, and
other bodyparts cannot be in the ROI for the script to consider the animal
"interacting" or "in" the object. Zones can also be summed together at the
end.

ROIs and other variables are saved as "shapes.mat" and a summary output
file is saved as "analyzedData.csv". 

Users have the option to create a plotted video where the trace of a select bodypart is
overlaid each frame of the video, and the marker turns red and stays on the
video if the bodypart is "included" in a select zone. There may be some
bugs in the workflow if a user is loading in previously analyzed data to
plot the videos if not all of the relevant variables were previously saved
in "shapes.mat".

To-do:
- Create zones for exclusion (i.e. the center of an object when the
inclusion zone should only be the edge of an object)

- Allow users to batch write videos or watch as they get plotted
individually. For now, users must change line 461 to be a parfor loop for
batch writing, and a for loop for individual videos.
%}

clear all

%% User input settings 
inputs = inputdlg({"Video Type: (i.e. mp4, avi, wmv, etc.)","Folder name for this analysis"},"Settings",[1 40; 1 40],{'avi','analysis1'});
videoType = inputs{1};
saveFolder = inputs{2};

path = string(uigetdir('C:\',"Load folder containing videos and csv containing DLC coordinates")); % path to folder with all the videos and DLC analysis files
mkdir(strcat(path,'\',saveFolder));

% get names of videos and csv files
strFiles = strArray(path, 'csv');
strVids = strArray(path, videoType);
vidIDs = erase(strVids,(wildcardPattern+"/"));

if isfile(fullfile(path,saveFolder,'analyzedData.csv'))
    analyzeOption = questdlg('These videos have been analyzed previously. Reanalyze or skip and use saved data?','analyzedData file detected','Yes, reanalyze','Skip','Skip');
    if contains(analyzeOption,'Skip') %load saved data
        load(strcat(path,"\",saveFolder,"\shapes.mat"));
    end
end
if ~isfile(fullfile(path,saveFolder,'analyzedData.csv')) || contains(analyzeOption,'Yes')
    % select relevant bodyparts
    data1 = readcell(strFiles(1,1)); %get first csv
    bodyIdx = find(strcmp(data1(3,:),'x'));
    bodyList = data1(2,bodyIdx);
    [select, ~] = listdlg('PromptString',{'Select all bodyparts of interest.'},'ListString',bodyList,'SelectionMode','multiple');
    include = bodyList(select);
    incXcol = bodyIdx(select);
    
    % measure distance?
    analyzeDistance = questdlg('Analyze distance traveled?','Option','Yes','No','Yes');
    if contains(analyzeDistance,'Yes')  % calibrate units if analyzing distance
        video = VideoReader(strVids(1, 1));
        frame = read(video, video.NumFrames); %load last frame
        figure;
        imshow(frame)
        title("Calibrate distance: Create a rectangle of known size")
        set(gca,'Position',[0 0 1 0.9])
        draw = drawrectangle;
        arena = customWait(draw);
        disp("Double click to confirm ROI");
        setArena = inputdlg({"height (cm)?", "width (cm)?"},"Set rectangle size to calibrate pixel/distance ratio",[1 40; 1 40], {'1','1'});
        cmX = arena.Position(4)/str2double(setArena{1}); % width in cm / x value pixels 
        cmY = arena.Position(3)/str2double(setArena{2}); % height in cm / y value pixels  
        cmPix = mean([arena.Position(4) arena.Position(3)]) / mean([str2double(setArena{1}) str2double(setArena{2})]); % average conversion for both dimensions
    
        [distBody, ~] = listdlg('PromptString',{'Select bodypart for tracking distance'},'ListString',bodyList,'SelectionMode','single');
        distBodypart = bodyList{distBody};
        distBodyXcol = bodyIdx(distBody);
    end
    
    
    % create final data table for export
    dataOut = table(vidIDs);
   
    %% Analysis options
    
    % any bodyparts to exclude from animal time spent in zone?
    [bodyCriteriaIdx,~] = listdlg('PromptString',{'Select all bodyparts that must be inside each object/zone', 'to meet criteria for inclusion (e.g. the nose and ears must be in the zone)'},'ListString',include,'SelectionMode','multiple');
    bodyCriteria = include(bodyCriteriaIdx);
    temp = erase(include,bodyCriteria(bodyCriteriaIdx));
    excludeList = temp(~cellfun('isempty',temp));
    if ~isempty(excludeList)
        [excludeIdx,~] = listdlg('PromptString',{'Select any bodyparts that must be excluded from each object/zone', 'to meet criteria for inclusion (e.g. the centerMass and tail must not be in the zone)'},'ListString',excludeList,'SelectionMode','multiple');
        bodyExclude = excludeList(excludeIdx);
    end

    [plotBodyIdx,~] = listdlg('PromptString',{'Select bodypart to be plotted'},'ListString',bodyCriteria,'SelectionMode','single','OKString','Select','CancelString','Do not create videos');
    plotBody = bodyCriteria(plotBodyIdx);
    
    %% draw ROIs
    video = VideoReader(strVids(1, 1));
    videos = struct();
    shapes = struct();
    
    %check if ROIs were previously created
    if isfile(strcat(path,"\",saveFolder,"\shapes.mat")) 
        useROIs = questdlg('ROIs were previously generated for this batch of videos. Would you like to use the same ROIs for analysis?','ROIs detected','Yes, apply previous settings','Edit previous ROIs','Create New ROIs','Create New ROIs');
    end
    
    if ~isfile(strcat(path,"\",saveFolder,"\shapes.mat")) || contains(useROIs,'Create') %create new ROIs
        %create shapes on first video
        video = VideoReader(strVids(1, 1)); 
        frame = read(video, video.NumFrames); %load last frame
        figure;
        imshow(frame)
        title("Create objects/zones")
        set(gca,'Position',[0 0 1 0.9])

        hold on
        videos(1).data = readmatrix(strFiles(1));
        bodyXcol = bodyIdx(plotBodyIdx);
        bodyYcol = bodyXcol + 1;
        xcoords = videos(1).data(:,bodyXcol);
        ycoords = videos(1).data(:,bodyYcol);
        plot(xcoords,ycoords,'Color',[0.3 0.3 0.3 0.5])
        [newshape, tf] = listdlg('PromptString',{'Select shape to draw ROI'},'ListString',{'round','rectangular'},'SelectionMode','single','OKString','Create New Shape', 'CancelString','Done');
        while tf == 1
            switch newshape
                case 1 %round
                    shapeinputs = inputdlg({"Object/Zone Name:", "Preset height (in pixels)?", "Preset width (in pixels)?"},"Settings",[1 40; 1 40; 1 40], {'zone1','50','50'});
                    disp("Double click to confirm ROI");
                    draw = drawellipse('Label', shapeinputs{1},'Center',[100,100],'SemiAxes',[str2double(shapeinputs{2}),str2double(shapeinputs{3})]);
                    roi = customWait(draw);
                    if ~isfield(shapes,'name')
                        shapes(1).name = shapeinputs{1};
                    else 
                        shapes(end+1).name = shapeinputs{1};
                    end
                    shapes(end).shapeType = 'round';
                    shapes(end).roi = roi;
                    disp(strcat(shapeinputs{1}," created! height = ",string(roi.SemiAxes(1)), " width = ", string(roi.SemiAxes(2))));
                    [newshape, tf] = listdlg('PromptString',{'Select shape to draw ROI'},'ListString',{'round','rectangular'},'SelectionMode','single','OKString','Create New Shape', 'CancelString','Done');
                case 2 %rect
                    shapeinputs = inputdlg({"Object/Zone Name:", "Preset height (in pixels)?", "Preset width (in pixels)?"},"Settings",[1 40; 1 40; 1 40], {'zone1','50','50'});
                    disp("Double click to confirm ROI");
                    draw = drawrectangle('Label',shapeinputs{1},'Position',[100,100,str2double(shapeinputs{3}),str2double(shapeinputs{2})]);
                    roi = customWait(draw);
                    if ~isfield(shapes,'name')
                        shapes(1).name = shapeinputs{1};
                    else 
                        shapes(end+1).name = shapeinputs{1};
                    end
                    shapes(end).shapeType = 'rect';
                    shapes(end).roi = roi;
                    disp(strcat(shapeinputs{1}," created! height = ",string(roi.Position(4)), " width = ", string(roi.Position(3))));
                    [newshape, tf] = listdlg('PromptString',{'Select shape to draw ROI'},'ListString',{'round','rectangular'},'SelectionMode','single','OKString','Create New Shape', 'CancelString','Done');
            end
        end
        videos(1).shapes = shapes;
        %repeat for remaining videos
        for v = 2:length(strVids)
            video = VideoReader(strVids(v,1)); %load first video
            frame = read(video, video.NumFrames); % extract the last frame of every video
            figure;
            imshow(frame)
            title(strcat("Set objects/zones for video",string(v)))
            set(gca,'Position',[0 0 1 0.9])
            hold on
            
            videos(v).data = readmatrix(strFiles(v,1));
            bodyXcol = bodyIdx(plotBodyIdx);
            bodyYcol = bodyXcol + 1;
            xcoords = videos(v).data(:,bodyXcol);
            ycoords = videos(v).data(:,bodyYcol);
            plot(xcoords,ycoords,'Color',[0.3 0.3 0.3 0.5])
            
            %copy rois from last video but allow repositioning
            videos(v).shapes = videos(v-1).shapes;
            for o = 1:length(videos(v-1).shapes)
                if contains(videos(v-1).shapes(o).shapeType,'rect')
                    draw = drawrectangle('Position',videos(v-1).shapes(o).roi.Position,'Label',videos(v-1).shapes(o).roi.Label,'Color','red');
                    disp("Double click to confirm ROI");
                    roi = customWait(draw);
                    videos(v).shapes(o).roi = roi;
                elseif contains(videos(v-1).shapes(o).shapeType,'round')
                    draw = drawellipse('Center',videos(v-1).shapes(o).roi.Center,'SemiAxes',videos(v-1).shapes(o).roi.SemiAxes,'RotationAngle',videos(v-1).shapes(o).roi.RotationAngle,'Label',videos(v-1).shapes(o).roi.Label,'Color','red');
                    disp("Double click to confirm ROI");
                    roi = customWait(draw);
                    videos(v).shapes(o).roi = roi;
                end
            end
        end
        save(strcat(path,"\",saveFolder,"\shapes.mat"),"videos","shapes","bodyCriteria","bodyIdx",'-mat');
    
    elseif contains(useROIs,'Yes') %load existing ROIs
        load(strcat(path,"\",saveFolder,"\shapes.mat"));

    elseif contains(useROIs,'Edit')
        editROIoptions = questdlg('Resize previous ROIs or redraw?','Edit ROI options','Resize','Redraw','Resize');
        if contains(editROIoptions,'Resize')
                resizeROI = inputdlg({'Enter a fraction to resize the ROI, keeping the center the same as before:'},"Options for resizing ROIs",[1 40],"1.1");
        end
        load(strcat(path,"\",saveFolder,"\shapes.mat"));
        for v = 1:length(strVids)
            video = VideoReader(strVids(v,1)); %load first video
            frame = read(video, video.NumFrames); % extract the last frame of every video
            figure;
            imshow(frame)
            title(strcat("Set objects/zones for video",string(v)))
            set(gca,'Position',[0 0 1 0.9])
            hold on
            videos(v).data = readmatrix(strFiles(v,1));
            bodyXcol = bodyIdx(plotBodyIdx);
            bodyYcol = bodyXcol + 1;
            xcoords = videos(v).data(:,bodyXcol);
            ycoords = videos(v).data(:,bodyYcol);
            plot(xcoords,ycoords,'Color',[0.3 0.3 0.3 0.5])
            
            %load prev rois
            if contains(editROIoptions,'Redraw')
                for o = 1:length(videos(v).shapes)
                    if contains(videos(v).shapes(o).shapeType,'rect')
                        draw = drawrectangle('Position',videos(v).shapes(o).roi.Position,'Label',videos(v).shapes(o).roi.Label,'Color','red');
                        disp("Double click to confirm ROI");
                        roi = customWait(draw);
                        videos(v).shapes(o).roi = roi;
                    elseif contains(videos(v).shapes(o).shapeType,'round')
                        draw = drawellipse('Center',videos(v).shapes(o).roi.Center,'SemiAxes',videos(v).shapes(o).roi.SemiAxes,'RotationAngle',videos(v).shapes(o).roi.RotationAngle,'Label',videos(v).shapes(o).roi.Label,'Color','red');
                        disp("Double click to confirm ROI");
                        roi = customWait(draw);
                        videos(v).shapes(o).roi = roi;
                    end
                end
            elseif contains(editROIoptions,'Resize')
                for o = 1:length(videos(v).shapes)
                    if contains(videos(v).shapes(o).shapeType,'rect')
                        centerX = videos(v).shapes(o).roi.Position(1) + videos(v).shapes(o).roi.Position(3)/2;
                        centerY = videos(v).shapes(o).roi.Position(2)+ videos(v).shapes(o).roi.Position(4)/2;
                        newWidth = videos(v).shapes(o).roi.Position(3)*str2double(resizeROI{1});
                        newHeight = videos(v).shapes(o).roi.Position(4)*str2double(resizeROI{1});
                        newX = centerX - newWidth/2;
                        newY = centerY - newHeight/2;
                        videos(v).shapes(o).roi.Position = [newX, newY, newWidth, newHeight];
                        draw = drawrectangle('Position',videos(v).shapes(o).roi.Position,'Label',videos(v).shapes(o).roi.Label,'Color','red');
                    elseif contains(videos(v).shapes(o).shapeType,'round')
                        videos(v).shapes(o).roi.SemiAxes = videos(v).shapes(o).roi.SemiAxes*str2double(resizeROI{1});
                        draw = drawellipse('Center',videos(v).shapes(o).roi.Center,'SemiAxes',videos(v).shapes(o).roi.SemiAxes,'RotationAngle',videos(v).shapes(o).roi.RotationAngle,'Label',videos(v).shapes(o).roi.Label,'Color','red');
                    end
                end
            end
        end
        save(strcat(path,"\",saveFolder,"\shapes.mat"),"videos","shapes","bodyCriteria","bodyIdx",'-mat');
    end
    
    

    
    %% Get binary arrays for each timepoint whether each included bodypart is inside each ROI
    
    flag = 0;
    
    for v = 1:length(strVids)
        videos(v).data = readmatrix(strFiles(v,1)); 
      
        for i = 1:length(include)
            bodyXcol = bodyIdx(i);
            bodyYcol = bodyXcol + 1;
            xcoords = videos(v).data(:,bodyXcol); 
            ycoords = videos(v).data(:,bodyYcol);
    
            % check whether the trace appears to move over 10% of the image between
            % timepoints to detect tracing errors
    
            xthresh = max(xcoords)/3;
            ythresh = max(ycoords)/3;
            if flag == 0
                for x = 2:length(xcoords)
                    if xcoords(x)-xcoords(x-1) > xthresh || ycoords(x)-ycoords(x-1) > ythresh
                        flag = 1;
                    end
                end
            end
            
            if flag == 1
                video = VideoReader(strVids(v));
                frame = read(video,1);
                fig = figure();
                h = imshow(frame);
                set(gca,'Position',[0 0 1 0.9])
                title(vidIDs(v),'Interpreter','none')
                hold on
                plot(xcoords,ycoords);
                fix = questdlg('Possible tracing errors detected: Coordinates jump over 33% of the image between timepoints. Would you like to correct the traces by smoothing?','Possible Tracing Errors Detected','Yes','No','Yes');
                if contains(fix,'Yes')
                    fixThresh = inputdlg({sprintf('Set pixel thresholds for maximum distance that a bodymarker could reasonably move between timepoints to smooth the trace.\n\nIf bodymarker coordinates for adjacent timepoints exceed this threshold, the coordinate is replaced with the preceding coordinate as if the animal has not moved.\n\nThe prefilled values below correspond to 33%% of the image in pixels.\n\n This correction will apply to the entire batch of videos.\n\n\nX threshold:'),'Y threshold'},'Set Smoothing Thresholds',[1 80; 1 80],{string(xthresh),string(ythresh)});
                    xthresh = str2double(fixThresh{1});
                    ythresh = str2double(fixThresh{2});
                    xcoordsFix = fixTrace(xcoords,xthresh);
                    ycoordsFix = fixTrace(ycoords,ythresh);
                    traceFixed = plot(xcoordsFix, ycoordsFix);
                    accept = questdlg('Accept smoothed trace?','Verify Trace','Yes, apply to all videos','No, reset threshold','Yes, apply to all videos');
                    while contains(accept,'No')
                        fixThresh = inputdlg({sprintf('Set pixel thresholds for maximum distance that a bodymarker could reasonably move between timepoints to smooth the trace.\n\nIf bodymarker coordinates for adjacent timepoints exceed this threshold, the coordinate is replaced with the preceding coordinate as if the animal has not moved.\n\nThe prefilled values below correspond to 10%% of the image in pixels.\n\n This correction will apply to the entire batch of videos.\n\n\nX threshold:'),'Y threshold'},'Set Smoothing Thresholds',[1 80; 1 80],{string(xthresh),string(ythresh)});
                        xthresh = str2double(fixThresh{1});
                        ythresh = str2double(fixThresh{2});
                        xcoordsFix = fixTrace(xcoords,xthresh);
                        ycoordsFix = fixTrace(ycoords,ythresh);
                        set(traceFixed,'XData',xcoordsFix,'YData',ycoordsFix);
                        accept = questdlg('Accept smoothed trace?','Verify Trace','Yes, apply to all videos','No, reset threshold','Yes, apply to all videos');
                    end
                    videos(v).data(:,incXcol(i)) = fixTrace(xcoords,xthresh);
                    videos(v).data(:,incXcol(i)+1) = fixTrace(ycoords,ythresh);
                end
                close(fig);
                flag = 2;
            elseif flag == 2
                    videos(v).data(:,incXcol(i)) = fixTrace(xcoords,xthresh);
                    videos(v).data(:,incXcol(i)+1) = fixTrace(ycoords,ythresh);
            end
    
            for o = 1:length(videos(v).shapes)
                if ~contains(videos(v).shapes(o).shapeType,'summed')
                    inside = inROI(videos(v).shapes(o).roi, videos(v).data(:,incXcol(i)), videos(v).data(:,incXcol(i)+1));
                    videos(v).shapes(o).(include{i}).logical = inside;
                    videos(v).shapes(o).(include{i}).seconds = sum(inside)/video.FrameRate;
        
                    if contains(analyzeDistance,'Yes') && contains(include{i},distBodypart)
                        marker(:,1) = videos(v).data(:,incXcol(i));
                        marker(:,2) = videos(v).data(:,incXcol(i)+1);
                        diffPoints = diff(marker);
                        videos(v).shapes(o).distanceArray{:,v} = hypot(diffPoints(:,1),diffPoints(:,2))/arena.Position(4)*str2double(setArena{1});
                        videos(v).shapes(o).distanceTotal{:,v} = sum(videos(v).shapes(o).distanceArray{:,v});
                        videos(v).shapes(o).averageVelocity{:,v} = videos(v).shapes(o).distanceTotal{:,v}/video.Duration; %cm/s
                    end
                end
            end
        end
        
        %sum up the bodyparts to be included, and subtract the parts that are excluded
        for o = 1:length(videos(v).shapes)
            if ~contains(videos(v).shapes(o).shapeType,'summed')
                for b = 1:length(bodyCriteria)
                    critLog(:,b) = videos(v).shapes(o).(bodyCriteria{b}).logical;
                end
                videos(v).shapes(o).criteriaMet.logical = all(critLog,2);
                videos(v).shapes(o).criteriaMet.seconds = sum(all(critLog,2))/video.FrameRate;
                if ~isempty(excludeList)
                    for e = 1:length(bodyExclude)
                        excLog(:,e) = videos(v).shapes(o).(bodyExclude{e}).logical;
                    end   
                    videos(v).shapes(o).criteriaMet.logical = all(critLog,2) & ~any(excLog,2);
                    videos(v).shapes(o).criteriaMet.seconds = sum(all(critLog,2) & ~any(excLog,2))/video.FrameRate;
                end
    
                dataOut.(strcat(shapes(o).name,'TotalSeconds'))(v) = videos(v).shapes(o).criteriaMet.seconds;
                if contains(analyzeDistance,'Yes')
                    dataOut.(strcat(shapes(o).name,'TotalDistance'))(v) = videos(v).shapes(o).distanceTotal;
                    dataOut.(strcat(shapes(o).name,'AvgVelocity'))(v) = videos(v).shapes(o).averageVelocity;
                end
            end
        end
    end
    
    %% sum up any zones
    summary = struct();
    [zoneIdx,tf] = listdlg('PromptString',{'Select any zones to sum'},'ListString',{shapes.name},'SelectionMode','multiple','OKString','Next','CancelString','Cancel');
    if tf == 1
        zoneSumName = inputdlg({"Give the summed zones a name"},"Zone Name",[1 40],{"sumZone1"});
        summary(1).name = zoneSumName{1};
        for i = 1:length(zoneIdx)
            summary(1).zones{i} = shapes(zoneIdx(i)).name;
        end
    end
    
    while tf == 1
        [zoneIdx,tf] = listdlg('PromptString',{'Select any zones to sum'},'ListString',{shapes.name},'SelectionMode','multiple','OKString','Next','CancelString','Done');
        if tf == 1
            zoneSumName = inputdlg({"Give the summed zones a name"},"Zone Name",[1 40],{"sumZone2"});
            summary(end+1).name = zoneSumName{1};
            for i = 1:length(zoneIdx)
            summary(end).zones{i} = shapes(zoneIdx(i)).name;
            end
        end
    end
    
    if isfield(summary,'name')
        for v = 1:length(strVids)
            for s = 1:length(summary)
                for r = 1:length(summary(s).zones)
                    for o = 1:length(shapes)
                        if any(contains(summary(s).zones(r),videos(v).shapes(o).name)) %if this shape is part of this summary zone
                            summary(s).sumZones(:,r) = videos(v).shapes(o).criteriaMet.logical;
    
                            if contains(analyzeDistance,'Yes')
                                summary(s).sumDist(:,r) = videos(v).shapes(o).distanceTotal{:,v};
                                summary(s).avgVelocity(:,r) = videos(v).shapes(o).averageVelocity{:,v};
                            end
    
                        end
                    end
                end
            summary(s).logical{:,v} = any(summary(s).sumZones,2);
            summary(s).seconds{:,v} = sum(any(summary(s).sumZones,2))/video.FrameRate;
    
            if contains(analyzeDistance,'Yes')
                summary(s).distanceTotal{:,v} = sum(summary(s).sumDist);
                summary(s).averageVelocity{:,v} = mean(summary(s).avgVelocity);
            end
            end
        end
        
        for s = 1:length(summary)
            dataOut.(strcat(summary(s).name,'Seconds')) = cell2mat(summary(s).seconds)';
            if contains(analyzeDistance,'Yes')
                dataOut.(strcat(summary(s).name,'TotalDistance')) = cell2mat(summary(s).distanceTotal)';
                dataOut.(strcat(summary(s).name,'AvgVelocity')) = cell2mat(summary(s).averageVelocity)';
            end
    
            shapes(end+1).name = summary(s).name;
            shapes(end).shapeType = 'summed';
            shapes(end).zones = summary(s).zones;
            shapes(end).criteriaMet = summary(s).logical;

            videos(v).shapes(end+1).name = summary(s).name;
            videos(v).shapes(end).shapeType = 'summed';
            videos(v).shapes(end).zones = summary(s).zones;
            videos(v).shapes(end).criteriaMet = summary(s).logical;
        end
    save(strcat(path,"\",saveFolder,"\shapes.mat"),"videos","shapes","bodyCriteria","bodyIdx",'-mat');
    end
    
    %save data table
    writetable(dataOut,fullfile(path,saveFolder,'analyzedData.csv'));
    disp('Data analyzed and saved!');
end

%% create video with traces and color coded interaction timepoints
createVid = questdlg('Create plotted videos?','Video Options','Yes','No','Yes');
if contains(createVid,'Yes')
    [plotBodyIdx,~] = listdlg('PromptString',{'Select bodypart to be plotted'},'ListString',bodyCriteria,'SelectionMode','single','OKString','Select','CancelString','Do not create videos');
    vidOption = inputdlg({'Enter duration(s) of bodyMarker trail'},"Video Output Settings",[1 40],{'5'});
    trailDur = str2double(vidOption{1});
    plotBody = bodyCriteria(plotBodyIdx);

    [plotROI,~] = listdlg('PromptString',{'Select a region of interest for the plotted bodypart'},'ListString',{shapes.name},'SelectionMode','single','OKString','Select','CancelString','Do not create videos');
    
    vidPath = strcat(path,'\',saveFolder,'\plotted');
    
    if ~isfolder(vidPath)
        mkdir(vidPath);
        createVid = 'new';
    elseif isfolder(vidPath)
        createVid = questdlg('Overwrite existing videos, create new files, or skip?','Existing videos detected','Yes, create new files','Yes, overwrite existing videos','Skip','Skip');
    end

    % parallelVid = questdlg('View and plot videos individually or batch plot without viewing?','Video Writing Options','View individual','Batch write','Batch write');
    % 
    % if contains(parallelVid,'View')
    %     parpool()
    % end

    parfor v = 1:length(strVids)
        vidFilename = fullfile(vidPath,strcat(extractBefore(vidIDs(v,1),'.'),'_plotted.avi'));
        if isfile(vidFilename) && contains(createVid,'Skip') %do not make video
            write = 0;
            disp(strcat('Skipping video: ', vidFilename));
        elseif contains(createVid,'new') %create new video
            addTime = string(datetime('now','Format','yyyy-MM-dd_HHmm'));
            writerObj = VideoWriter(strcat(fullfile(vidPath,extractBefore(vidIDs(v,1),'.')),'_',addTime));
            write = 1;
        elseif ~isfile(vidFilename) || contains(createVid,'overwrite') %create new or overwrite
            writerObj = VideoWriter(vidFilename);
            write = 1;
        end

        if write == 1
            %create video
            video = VideoReader(strVids(v, 1)); 
            writerObj.FrameRate = video.FrameRate;
            open(writerObj);
            disp(strcat('Creating video: ',writerObj.Filename));
            frame = read(video,1);
            fig = figure();
            h = imshow(frame);
            set(gca,'Position',[0 0 1 0.9])
            title(vidIDs(v),'Interpreter','none')
            hold on

            %data = readmatrix(strFiles(v,1)); 

            %get coordinates of bodypart to be plotted
            bodyXcol = bodyIdx(plotBodyIdx);
            bodyYcol = bodyXcol + 1;
            xcoords = videos(v).data(:,bodyXcol);
            ycoords = videos(v).data(:,bodyYcol);

            %make video with chosen bodypart for plot, changes color when criteria is met    
            %initialize frames
            frameArray = repmat(getframe(gcf),length(xcoords),1);

            for t = 1:length(xcoords)
                frame = read(video,t);
                if t == 1 %first frame of video shows whole plot
                   h = imshow(frame);
                else
                    set(h,'CData',frame);
                end

                %plot trail up to this time point
                trail = trailDur*video.FrameRate; %how long the trail should be in terms of time
                if t <= trail
                    p = plot(xcoords(1:t),ycoords(1:t),'Color','black');
                elseif t == length(xcoords) %last frame of video has entire plot
                    p = plot(xcoords(1:t),ycoords(1:t),'Color','black');
                else
                    p = plot(xcoords(t-trail:t),ycoords(t-trail:t),'Color','black');
                end

                if shapes(plotROI).criteriaMet{1,v}(t) == 0
                    m = plot(xcoords(t),ycoords(t),'o','MarkerFaceColor','blue','MarkerSize',5,'MarkerEdgeColor','white');
                    F = getframe(gcf);
                    frameArray(t) = F;
                    delete(m);
                elseif shapes(plotROI).criteriaMet{1,v}(t) == 1
                    m = plot(xcoords(t),ycoords(t),'o','MarkerFaceColor','red','MarkerSize',5,'MarkerEdgeColor','white');
                    F = getframe(gcf);
                    frameArray(t) = F;
                end
                delete(p);
            end
            writeVideo(writerObj,frameArray)
            close(writerObj)
            close(fig)
        end
    end
    disp('Finished creating videos');
    close all
end

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


%% if sequential timepoints have coordinates that exceed the preset threshold, set the timepoint to the one previous to avoid tracking errors
function trace = fixTrace(coord, thresh)
    trace = coord;
    for i = 1:length(trace)-1
        if abs(trace(i+1) - trace(i)) > thresh
            trace(i+1) = trace(i);
        end
    end
end

%% ROI functions
function ROI = customWait(ROI)
    l = addlistener(ROI,'ROIClicked',@clickCallback);
    uiwait;
    delete(l);
end

function clickCallback(~,evt)
    if strcmp(evt.SelectionType,'double')
        uiresume;
    end
end