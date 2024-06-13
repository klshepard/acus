data = xlsread('/Users/sinan/Library/CloudStorage/Dropbox/NatElec_dualShank_2023/Source Data/Fig4_b.xlsx', '4b');

xlswrite('/Users/sinan/Library/CloudStorage/Dropbox/NatElec_dualShank_2023/Source Data/Fig3_bb.xlsx', saveData_Fig3b, '3b');

% figure(34); %Fig5b
% plot(data(3:66,2), '-.', 'LineWidth',2, 'color', '#d95f02'); hold on;
% plot(data(3:66,4), '-.', 'LineWidth',2, 'color', '#7570b3');
% plot(data(3:66,6), '-.', 'LineWidth',2, 'color', '#1b9e77');
% 
% plot(data(:,3), 'LineWidth',2, 'color', '#d95f02');
% plot(data(:,5), 'LineWidth',2, 'color', '#7570b3');
% plot(data(:,7), 'LineWidth',2, 'color', '#1b9e77');
% 
% set(gca, 'box', 'off');
% set(gca,'LineWidth',2);
% ylim([-300, 1200]);
% xlim([3, 66]);
% set(gca,'TickLength',[0.01 0.01]);
% set(gca,'YTick', [0, 400, 800, 1200]);
% set(gca,'XTick', [3, 19, 35, 50, 64]);
% set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
% set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
% legend({'Time = 0.0 sec', 'Time = 3.0 sec', 'Time = 4.5 sec', 'Time = 0.0 sec', 'Time = 3.0 sec', 'Time = 4.5 sec'}, 'location', 'northeast', 'numColumns', 2);

% figure(33); %Fig3f
% plot(data(:,2), 'LineWidth',2, 'color', '#ff4a4a'); set(gca,'LineWidth',2); hold on;
% plot(data(:,3), 'LineWidth',2, 'color', '#3f42ff'); set(gca,'LineWidth',2); hold on;
% set(gca, 'box', 'off');
% ylim([0, 15000]);
% xlim([1, 128]);
% set(gca,'XTick', [500/25-1, 1800/25-1, 3100/25-1]);
% set(gca,'XTickLabels', {'500', '1800', '3100'});
% set(gca,'FontSize', 16, 'FontWeight','normal'); %'bold' or 'normal'
% legend('Row 1', 'Row 2');

% figure(32); %Fig4b
% plot(data(:,2), '--', 'LineWidth',2, 'color', '#d95f02'); set(gca,'LineWidth',2); hold on;
% plot(data(:,4), '--', 'LineWidth',2, 'color', '#7570b3');
% plot(data(:,6), '--', 'LineWidth',2, 'color', '#1b9e77');
% 
% plot(data(:,3), 'LineWidth',2, 'color', '#d95f02');
% plot(data(:,5), 'LineWidth',2, 'color', '#7570b3');
% plot(data(:,7), 'LineWidth',2, 'color', '#1b9e77');
% 
% set(gca, 'box', 'off');
% set(gca,'LineWidth',2);
% ylim([-5000, 70000]);
% xlim([1, 64]);
% set(gca,'TickLength',[0.01 0.01]);
% set(gca,'YTick', [0, 35000, 70000]);
% set(gca,'XTick', [1, 17, 33, 48, 62]);
% set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
% set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
% legend({'No Slice', 'Slice with No Neuron', 'Slice with Single Neuron', 'Row 1', 'Row 2', '.'}, 'location', 'northwest', 'numColumns', 2);
