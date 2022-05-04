function createPlots(dataForPlots)
% You might want to customize each section of this code. It will create
% plots in the style that Avi typically likes them in

allData   = dataForPlots.allData;
timeBefore = dataForPlots.timeBefore;
timeAfter = dataForPlots.timeAfter;
dsFactor  = dataForPlots.dsFactor;
FS        = dataForPlots.FS;
titleName = dataForPlots.titleName;
workdir   = dataForPlots.workdir;
trials    = dataForPlots.trials;
savename  = dataForPlots.savename;
subName  = dataForPlots.subName;
trialsOI    = dataForPlots.trialsOI;
xlims = [-timeBefore timeAfter];

%% event triggered average
figure(61)
hold on

% make an x axis 
xlong = linspace(-timeBefore,timeAfter,length(allData));
x = decimate(linspace(-timeBefore,timeAfter,length(allData)),dsFactor);

% get rid of any nan points
if any(any(isnan(allData)))
    cols2getridof = isnan(allData(1,:));
    allData(:,cols2getridof) = [];
    warning('some columns were deleted!')
end
y = mean(allData,2); centerY = mean(mean(allData(1:round(timeBefore*FS),:)));
y = (y - centerY); %center the data
[peakY,ebIdx] = findpeaks(y'); [maxY,maxIdx] = max(peakY); ebIdx = ebIdx(maxIdx);
xMax = xlong(ebIdx); 
if y(ebIdx) ~= maxY
    warning('you may have different max value') 
end
y = decimate(y,dsFactor);
%y = smoothdata(y,'SmoothingFactor',0.1);
sem = std(allData,0,2)./sqrt(size(allData,2)); % sem = std/sqrt(n)
semPeak = sem(ebIdx);
eb = decimate(sem',dsFactor);
lineProps.col{1} = [0 0.5 0];
mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 timeAfter])
end
ylim([-1 5])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score')
title([titleName ', n = ' num2str(length(workdir))])
hold off

if exist('varsAndFigs','dir')
    cd('varsAndFigs')
else
    mkdir('varsAndFigs')
    cd('varsAndFigs')
end
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_avged.svg'],'svg'), saveas(gcf,[savename '_avged.fig'],'fig'), 

save([savename '_peak_SEM.mat'],'maxY','semPeak','xMax')
cd('..')


%% plot only some events, such as the last stims

%{
% event triggered average
figure(63)
allData = ad(:,[2 3 5 6 8 9]);
hold on
xlong = linspace(-timeBefore,timeAfter,length(nData));
x = decimate(linspace(-timeBefore,timeAfter,length(nData)),dsFactor);
if any(any(isnan(allData)))
    cols2getridof = isnan(allData(1,:));
    allData(:,cols2getridof) = [];
    warning('some columns were deleted!')
end
y = mean(allData,2); centerY = mean(mean(allData(1:round(timeBefore*FS),:)));
y = (y - centerY); %center the data
[peakY,ebIdx] = findpeaks(y'); [maxY,maxIdx] = max(peakY); ebIdx = ebIdx(maxIdx);
xMax = xlong(ebIdx); 
if y(ebIdx) ~= maxY
    warning('you may have different max value') 
end
y = decimate(y,dsFactor);
%y = smoothdata(y,'SmoothingFactor',0.1);
sem = std(allData,0,2)./sqrt(size(allData,2)); % sem = std/sqrt(n)
semPeak = sem(ebIdx);
eb = decimate(sem',dsFactor);
lineProps.col{1} = [0 0.5 0];
mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 timeAfter])
end
ylim([-1 5])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score')
title([titleName ', n = ' num2str(length(workdir))])
hold off

if exist('varsAndFigs','dir')
    cd('varsAndFigs')
else
    mkdir('varsAndFigs')
    cd('varsAndFigs')
end
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_avged_lastStims.svg'],'svg'), saveas(gcf,[savename '_avged_lastStims.fig'],'fig'), 

save([savename '_peak_lastStims_SEM.mat'],'maxY','semPeak','xMax')
cd('..')
%}

%% WORK IN PROGRESS - early (all mice) to late (all mice) heatmap
figure(60)

% Until hold on, CUSTOMIZE for your dataset 
% i.e.: trials = 7    15     7     9     7    13     8    12     5
cst = cumsum(trials); dcst = diff([1 cst]);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);

cstNew = [1 1+cst]; newIdxs = cstNew(1:end-1);
for ii = 1:max(dcst)-1
    newIdxs = [newIdxs cstNew+ii];
end
newIdxs = unique(newIdxs,'stable');
newIdxs(newIdxs>size(allData,2)) = [];
allChrono = allData(:,newIdxs);

hold on
% make a time vector for the heat map and plot the heat map
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allChrono,2),allChrono')
L = line([0 0],[0 size(allChrono,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1)
xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ' chronologically, n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
colormap('parula')
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allChrono,2)+0.5])
%set(gca, 'yticklabel', subName, 'ytick', cst);%trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
hold off


set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_Chrono_heatmap.svg'],'svg'),saveas(gcf,[savename '_heatmap.fig'],'fig'),

%% ascend heatmap
figure(62)

%uncomment these if you prefer to get 
%[~,idx1] = sort(mean(LickTrig1(:,90000:150000),2),'ascend');
%LickTrig1 = LickTrig1(idx1,:);

cst = cumsum(trials);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);

%creating a matrix that has ascending peak response
[~,idx] = sort(mean(allData,1),'ascend');
allData2=allData(:,idx);

hold on
% make a time vector for the heat map and plot the heat map
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allData2,2),allData2')
L = line([0 0],[0 size(allData2,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1)
xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ', n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
colormap('parula')
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allData2,2)+0.5])
%set(gca, 'yticklabel', subName, 'ytick', cst);%trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
hold off

set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_heatmap.svg'],'svg'),saveas(gcf,[savename '_heatmap.fig'],'fig'),


%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(65)

cst = cumsum(trials);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);
% uncomment if you want your data to be descending
%{
for ii = 1:length(workdir)

if ii == 1
    [~,idx] = sort(mean(allData(:,1:cst(ii)),1),'descend');
    allData2(:,1:cst(ii))=allData(:,idx);
else
    [~,idx] = sort(mean(allData(:,cst(ii-1)+1:cst(ii)),1),'descend');
    allData2(:,cst(ii-1)+1:cst(ii))=allData(:,cst(ii-1)+idx);
end
end
%}
allData2 =allData;
hold on
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allData2,2),allData2')
L = line([0 0],[0 size(allData2,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1.5)

% plot horizontal lines
for ii = 1:length(workdir)
L2 = line([-150 size(allData2,1)],repmat(cst'+0.5, [1 2]));%([trials*ii + 0.5]));
set(L2,'Color','white')
set(L2,'LineWidth',1.5)
end

xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ', n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allData2,2)+0.5])
tickPts = [0 cst] + diff([1 cst 0])/2; tickPts= tickPts(1:end-1); 
yticks(trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
set(gca, 'yticklabel', subName, 'ytick', tickPts); 
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
hold off

cd('varsAndFigs')
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_SubjectHeatmap.svg'],'svg'),saveas(gcf,[savename '_heatmap.fig'],'fig'),

cd('..')
end