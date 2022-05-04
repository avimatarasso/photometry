function plotHeatAndIndAvg(data, heatFig, avgFig, timeBefore, timeAfter, subName, sess, dsFactor,FS)
    sem = std(data,0,2)./sqrt(size(data,2)); % sem = std/sqrt(n)
    data = data';    
%Figure 20#
    figure(avgFig)

    hold on
    %eb1 = errorbar(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,sem,'Color',[0.2,0.2,0.5]);
    x = decimate(linspace(-timeBefore,timeAfter,length(data)),dsFactor);
    if size(data,1) > 1
        y = mean(decimate(data',dsFactor),1); centerY = mean(data(1:ceil(timeBefore*FS)));
        y = y - centerY; %center the data
        %y = smoothdata(y,'SmoothingFactor',0.005);
        eb = decimate(sem',dsFactor);

        lineProps.col{1} = 'blue';
        mseb(x,y,eb,lineProps,1);

        L = line([0 0],[-3 3]);
        %axis([-timewindow timewindow -1.5 1.5])
        set(L,'Color','black')
        xlim([-20 timeAfter])
        ylim([-3 3])
        xlabel('Peri-Event Time (sec)')
        ylabel('Z score (averaged)')
        title(['Event-wise average Z-score with SEM for ' subName{sess}])
        box off
        hold off


        % Subject heat map (x dim: time, ydim: event, zdim: deltaF/F) %deltaF
        figure(heatFig)
        numbBouts = size(data,1);
        hold on
        imagesc(linspace(-timeBefore, timeAfter,length(data)),1:size(data,2),data)  %4*timewindow,length(photoPerLick)),1:size(LickTrig,1),LickTrig)
        L = line([0 0],[0 numbBouts+1]);
        set(L,'Color','black')
        L2 = line([30 30],[0 numbBouts+1]);
        set(L2,'Color','black')
        xlabel('Peri-Event Time (sec)')
        ylabel('Bout Number')
        %ytick(1:size(LickTrig,1))
        title(['Z-score around event time for each bout in ' subName{sess}])
        cb = colorbar;
        title(cb,'Z score')
        caxis([-2 2])
        colormap(flipud(brewermap([],'YlGnBu')))
        xlim([-20 timeAfter])%4*timewindow])
        ylim([0.5 numbBouts+0.5])%ylim([0 length(licks)+1])
        box off
        hold off
    else
        warning([subName{sess} ' had only one event, so no averaged heat map'])
    end
end