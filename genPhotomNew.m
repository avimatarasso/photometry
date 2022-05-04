% general photometry processing script for preprocessed TDT files
% preprocessed .mat files should contain "Dts", "data1" 
%  
% please use tdtExtract_clean to extract code
% 
% please email akmat@uw.edu with any questions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% add the path to all code and data directories
codeDir = 'C:\Users\avima\OneDrive\Documents\MATLAB\Starter Code - photometry'; %CUSTOMIZE
%dataDir = 'G:\PhotomData\NEED TO ANALYZE\STIM ASSAYS\VTAterminals\20hz3s_5hz30s\Avi_autoStim-210907-110207'; %CUSTOMIZE
PhotomPack = ('F:\code\newPackage_020922');

addpath(PhotomPack); addpath(codeDir); %addpath(dataDir);

    
%% CUSTOMIZE this whole section
%do you only want to look at the averaged data of *some* of the subjects?
% if you have subjects 0,1,2,4 in the directory you're working in, 
% and you want to look at 1 and 4, type [2 4];

%paths must be in a cell array
path_to_data = {'F:\PhotomData\NEED TO ANALYZE\vCA1\scent\Avi_stim-210913-102052\2MT'};

%%%%%%%% change the savename to get each timelocked trial (usually of stim)
doYouWantAll = 1; % make 0 if you want to only choose some of the subs %and define which ones you want through whichSub
% if doing multip[le paths, want to be careful and check whichSub you choose
whichSub = [1:4]; % which subjects do you want to look at/for multiple regions include all

%%% events what event do you want to look at
eventLabel  = 1; %i usually do 2 for tail lift and 1 for stim %-2 = leaveClosedTS, 1 = loco.startTS
manualEvent = 0; %1 if the scoring was done manually, 0 if TTL

%How many trials and do you want all?
allTrials = 1; %if you don't want all trials, make allTrials = 0
if ~allTrials 
    trials   = 3; % how many of eventLabel# did you do?
    trialsOI = 1:trials; % 1:10 for example, which trials are you interested in
end

%customize these lines + defineSubject.m for your specific metadata + style
%%%% commonStr == common string between your experiments of interest
region = 'grabNE_vCA'; commonStr = region; %put in format of grabNE_CA1
%regionON = 0; %mark as 0 if region not included in name
savename = [region '_30s_2MT']; %3
a = findstr(region,'_'); region2 = region; region2(a) = '-'; %
titleName = [region2 ' 2MT for 30s'];

%%% do you want to z-score? if yes, make zQuest = 1, if no, make zQuest = 0
zQuest = 1;
% do you want first5 only?
first5 = 0; nn = 1; %can change nn to be any first nn events
normalQ = 1; %do you want to normalize to time immediately prior to event

%%% do you want each subject's timelocked and averaged traces and heat map?
% if averageQuest = 1, they will be shown when you run the code
averageQuest = 0;

%%% how much time (in s) after the event do you want to look at 
%timewindow = 30; % +/- seconds from event %usually 120
%specialCase = 2; %multiplier for timewindow to end. [-timewindow specialCase*timewindow]
% If only want to see time within +/- timewindow, keep specialCase as 1
%older options above, use timeBefore and after below instead
timeBefore = 30; %10 to 30 for short, 30 to 120
timeAfter  = 180; 
xlims = [-timeBefore timeAfter];

%CUSTOMIZE
dsFactor = 500; % how much do you want to downsample by?
baselineTime = 180; % in s
FS = 1017.25; % sampling frequency of synapse

%% Define Directory and constant variables

workdir = defineDir(path_to_data,commonStr);

%Define the working directories based on whichSub
if ~doYouWantAll
    workdir = workdir(whichSub);
end

%initialize baseline variables and timing arrays
baselineSD = zeros(length(workdir),1);
baselineMu = zeros(length(workdir),1);
subName = cell(1, length(workdir));
allTimings    = cell(length(workdir),1);
lengthOfBL = floor(abs(timeBefore)*FS);
lengthOfData = floor((timeAfter+timeBefore)*FS);
trials = []; 

% Allows us to ignore temporary files in the directory 
workdir = workdir(~startsWith({workdir.name}, '._'));
%% save current script
FileNameAndLocation = [mfilename('fullpath')];
if ~strcmp(path_to_data(end),'\')
    path_to_data = [path_to_data '\'];
end
scriptName = strsplit(FileNameAndLocation,'\'); scriptName = strjoin([path_to_data scriptName{end}]); 
currentTime = datestr(clock, 30); currentTime = [currentTime(1:8) '_' currentTime(10:13)];
newbackup=sprintf(['%s' '_dt' currentTime '.m'],scriptName);
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);


%% Set up events and align them

for subjectIdx = 1:length(workdir) 

    clearvars -except sess k itchStack workdir savename trialname ...
         tempPath path_to_data sTs tT whichSub zQuest baselineTime baselineSD baselineMu...
        subName subStr averageQuest firstQuest dsFactor trialsOI region first5 nn ...
        eventLabel allTrials timewindow specialCase normalQ lengthOfBL lengthOfData trials ...
        timeBefore timeAfter proc time sem subjectIdx manualEvent FS bl nData timingIdxs allData titleName xlims 
    
    %scroll through each photometry folder 
    cd(workdir(subjectIdx).folder)
    photoname = workdir(subjectIdx).name;
    load(photoname)
    data1 = dfF; % CUSTOMIZE if you want to use a different parameter for your analysis
    raw470 = dat470; raw405 = dat405; 
    %CUSTOMIZE the labels and subject names you may want
    [subStr, subNumb,subjectLabel, txtName] = defineSubject(photoname, region,subjectIdx);
    subName{subjectIdx} = subjectLabel;
    % txtName = [photoname '.txt'];
    rawtimings = dlmread(txtName); %CUSTOMIZE
      
    % find event times
    eventTime = rawtimings(:,1);
    eventtype = round(rawtimings(:,2));
    eventTime = eventTime(eventtype==eventLabel); 
    eventTime=eventTime(eventTime>5); %eventTime must be after 5s for some sort of baselin
    del_idx = [0; (diff(eventTime) < 20)]; %get rid of any events with less than 10 s in between
    
    % don't need to look at the timings if there is no event
    if isempty(eventTime)
        continue
    end
    eventTime = eventTime(~del_idx); % deletes events with less than 10 s in between
    
    % If you want to limit how many of the events you use change the if statement below!
    if first5
        if length(eventTime)> nn
             eventTime = eventTime(1:nn);
        end
    end
    
    % if manually scored and using synapse, you want to adjust for
    % nonconstant frame rate
    if manualEvent
        correctTime = max(Dts);
        vidDir = photoname(1:end-4); 
        cd(vidDir); 
        vidName = dir(['Avi*' vidDir '*.avi']); vD = vidName.name;
        if strcmp(vD(1:2), '._')
            vD = vidName(2).name;
        end
        vid = VideoReader(vD); vidTime = vid.Duration;
        vidRatio = correctTime/vidTime;
        cd('..')
        allTimings{subjectIdx} = eventTime*vidRatio; %CUSTOMIZE
    else
        allTimings{subjectIdx} = eventTime;
    end
    
    %find baseline for each session
    baselineMu(subjectIdx) = mean(data1(1:ceil(baselineTime*FS)));
    baselineSD(subjectIdx) = std(data1(1:ceil(baselineTime*FS)));
    
    % Z score photometry data
    if zQuest
        data1 = (data1-mean(data1))./std(data1);
        warning('You z-scored your data')
    else
        warning('NOT z-scored')
    end    
    
    %Align the data to events 
    [nData, timingIdxs, sDiff] = alignEvent(data1, FS, eventTime, timeBefore, timeAfter);

    time = linspace(1/FS,length(data1)/FS,floor(length(data1)));
    numActualTrials  = size(nData,2);
    
    if allTrials
        trials   = [trials numActualTrials]; % how many of eventLabel# did you do?
        trialsOI = 1:numActualTrials; % 1:10 for example, which trials are you interested in
    else
        trials   = [trials numActualTrials]; % how many of eventLabel# did you do?
        trialsOI = 1:numActualTrials; % 1:10 for example, which trials are you interested in
    end
    
    if numActualTrials < length(trialsOI) 
        cols = (subjectIdx-1)*length(trialsOI)+trialsOI;
        cols = cols(1:numActualTrials);
    else
        cumTrial = cumsum(trials);
        if subjectIdx > 1 
            cols = cumTrial(subjectIdx-1) + trialsOI;
        else
            cols = trialsOI;
        end
    end
    
    if averageQuest
        heatFig = 200+subjectIdx;
        avgFig = 300+subjectIdx;
        plotHeatAndIndAvg(nData, heatFig, avgFig, timeBefore, timeAfter, subStr, subjectIdx, dsFactor,FS)
    end
    
    if normalQ
        [nData] = centAndNormData(nData, FS, timeBefore); %sDiff, eventTime); %normalized to pre-event
        %NcData = cData/std(data1(1:round(180*FS))); 
    end
    
    % save your data by each subject for plotting 
    allData(:,cols) = nData(:,trialsOI); 
    raw.raw405(subjectIdx,:) = raw405;
    raw.raw470(subjectIdx,:) = raw470;
    raw.fit405(subjectIdx,:) = fit405;
    raw.dfF(subjectIdx,:)    = dfF;
    raw.subName{subjectIdx} = subStr;
    
    %plot individual averages and heat map
    figureN = subjectIdx+100;
    plotIndividual(figureN, eventTime, raw, data1, subName, subjectIdx, dsFactor, time)
    
end


%% plots

dataForPlots.allData = allData;
dataForPlots.timeBefore = timeBefore;
dataForPlots.timeAfter = timeAfter;
dataForPlots.dsFactor = dsFactor;
dataForPlots.FS = FS;
dataForPlots.savename = savename;
dataForPlots.titleName = titleName;
dataForPlots.workdir = workdir;
dataForPlots.trials = trials;
dataForPlots.trialsOI = trialsOI;
dataForPlots.subName= subName;

if exist('xlims','var')
dataForPlots.xlims = xlims;
end

% create plots, might want to CUSTOMIZE
createPlots(dataForPlots)


%% save variables in folder
timevec=linspace(-timeBefore,timeAfter,length(allData));
if isfolder('varsAndFigs')
    cd('varsAndFigs')
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots','-append')
    else
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots')
    end
    cd('..')
else
    mkdir('varsAndFigs')
    cd('varsAndFigs')
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots','-append')
    else
        save(savename, 'allData','timingIdxs','timevec','path_to_data','path_to_data','dataForPlots')
    end
    cd('..')
    % close all
end
%}

%% Save All Figures
answer = questdlg('All graphs are saved. Would you like to save the individual figures too?', ...
	'Saving?', ...
	'Yes','No','No');
switch answer
    case 'Yes'
        figureSaveName = savename;
        if exist('varsAndFigs','dir')
            cd('varsAndFigs')
        else
            mkdir('varsAndFigs')
            cd('varsAndFigs')
        end
        tempdir = pwd;
        FolderName = tempdir;   % Your destination folder
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure'); FigList2 = FigList(1:end);
        savefig(FigList2, fullfile(FolderName,[figureSaveName '.fig']));
       % save([savename '.mat']);
        disp('Figures have been saved!')        
        cd('..')
    case 'No'
        disp('You may manually save figures if you want.')
end

