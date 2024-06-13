% figure(803);
% subplot(3,2,2); imagesc(physImg); set(gca,'Ydir','normal');
% set(gca,'YTick', [1,2]); yticklabels({'Row 1', 'Row 2'});
% set(gca,'XTick',[20,70,120]); xticks = 25.67*[20,70,120]; xticklabels(strsplit(num2str(xticks)));
% set(gca,'FontSize',12, 'FontWeight','normal'); %'bold' or 'normal'
% set(gca,'LineWidth',2);
% set(gca,'TickLength',[0 0]);
% 
% subplot(3,2,4); plot(physImg','LineWidth',2); xlim([1,length(physImg(1,:))]); legend('Row 1','Row 2');
% set(gca,'LineWidth',2);
% set(gca,'FontSize',12, 'FontWeight','normal'); %'bold' or 'normal'
% set(gca,'XTick', [1, 16, 31, 46, 61]); xticks = 25.67*[20,40,60]; xticklabels(strsplit(num2str(xticks))); % xticklabels([500, 1000, 1500]); 
% set(gca, 'box', 'off');

figure(805); %Fig4b-5b
subplot(2,3,1); physImg = data(:,2:3)'; imagesc(physImg); set(gca,'Ydir','normal'); caxis([0, 70000]);
set(gca,'YTick', [1,2]); yticklabels({'Row 1', 'Row 2'}); set(gca,'XTick', [1, 17, 33, 48, 62]); set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
%set(gca,'XTick',[20,70,120]); xticks = 25.67*[20,70,120]; xticklabels(strsplit(num2str(xticks)));
set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
set(gca,'LineWidth',2);
set(gca,'TickLength',[0 0]);

subplot(2,3,4); plot(physImg(1,:),'LineWidth',2,'color', '#ff4a4a'); hold on; plot(physImg(2,:),'LineWidth',2,'color', '#3f42ff');
xlim([1,length(physImg(1,:))]); ylim([-5000, 70000]); legend('Row 1','Row 2', 'location', 'northwest');
set(gca,'LineWidth',2);
set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
set(gca,'YTick', [0, 35000, 70000]); set(gca,'XTick', [1, 17, 33, 48, 62]); set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
set(gca, 'box', 'off');

subplot(2,3,2); physImg = data(:,4:5)'; imagesc(physImg); set(gca,'Ydir','normal'); caxis([0, 70000]);
set(gca,'YTick', [1,2]); yticklabels({'Row 1', 'Row 2'}); set(gca,'XTick', [1, 17, 33, 48, 62]); set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
%set(gca,'XTick',[20,70,120]); xticks = 25.67*[20,70,120]; xticklabels(strsplit(num2str(xticks)));
set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
set(gca,'LineWidth',2);
set(gca,'TickLength',[0 0]);

subplot(2,3,5); plot(physImg(1,:),'LineWidth',2,'color', '#ff4a4a'); hold on; plot(physImg(2,:),'LineWidth',2,'color', '#3f42ff');
xlim([1,length(physImg(1,:))]); ylim([-5000, 70000]); legend('Row 1','Row 2', 'location', 'northwest');
set(gca,'LineWidth',2);
set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
set(gca,'YTick', [0, 35000, 70000]); set(gca,'XTick', [1, 17, 33, 48, 62]); set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
set(gca, 'box', 'off');

subplot(2,3,3); physImg = data(:,6:7)'; imagesc(physImg); set(gca,'Ydir','normal'); caxis([0, 70000]);
set(gca,'YTick', [1,2]); yticklabels({'Row 1', 'Row 2'}); set(gca,'XTick', [3, 19, 35, 50, 64]); set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
%set(gca,'XTick',[20,70,120]); xticks = 25.67*[20,70,120]; xticklabels(strsplit(num2str(xticks)));
set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
set(gca,'LineWidth',2);
set(gca,'TickLength',[0 0]);

subplot(2,3,6); plot(physImg(1,:),'LineWidth',2,'color', '#ff4a4a'); hold on; plot(physImg(2,:),'LineWidth',2,'color', '#3f42ff');
xlim([1,length(physImg(1,:))]); ylim([-5000, 70000]); legend('Row 1','Row 2', 'location', 'northwest');
set(gca,'LineWidth',2);
set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
set(gca,'YTick', [0, 35000, 70000]); set(gca,'XTick', [1, 17, 33, 48, 62]); set(gca,'XTickLabels', {'0', '400', '800', '1200', '1600'});
set(gca, 'box', 'off');