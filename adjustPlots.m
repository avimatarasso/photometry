% adjustPlots
% This script will allow you to adjust fonts and parts of your plot to 
%  optimize workflow in illustrator

%% avged

figure(1)
%{
if timeBefore ~= 180
    xlim([-10 30])
    ylim([-0.5 1])
else
    xlim([-30 180])
ylim([-0.5 2])
end
%}
xlim([-30 180])
ylim([-1.4 1])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off
if exist('figureSaveName')
    savename = [figureSaveName '_avged'];
else
    savename = uigetfile({'*.svg';'*.fig';'*.'}, 'File Selector')
end

saveas(gcf,[savename],'svg')

%% color map

caxis([-5 3])
xlim([-10 180])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off
savename = uigetfile({'*.svg';'*.fig';'*.'}, 'File Selector')
saveas(gcf,[savename],'svg')
