xlim([1,128]);
ylim([-2000 ,2000]);
figure(27);

set(gca, 'box', 'off');
set(gca,'LineWidth',2);
set(gca,'TickLength',[0.01 0.01]);
set(gca,'XTick', [0, 400, 800, 1200 1600]);
set(gca,'YTick', [0, 400, 800, 1200]);
set(gca,'YTickLabels', {'0.1', '1', '10', '100'});
set(gca, 'YColor', 'none');
set(gca,'FontSize',23, 'FontWeight','normal'); %'bold' or 'normal'