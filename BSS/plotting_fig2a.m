meanImages_nobead = zeros(length(spadphase),128);
for i=1:91
    
    %image = Images(i,:);
    image = raw2SepRows_NASP_only(Images(i,:));
    
    figure(91);
    meanImages_nobead(i,:) = mean(image, 1);
    semilogy(mean(image, 1), 'LineWidth',2);
    ylim([1, 3e7]); xlim([1,128]); set(gca, 'box', 'off'); set(gca,'LineWidth',2); %set(gca,'FontSize',16, 'FontWeight','normal'); %'bold' or 'normal'
    
    legend(['Spadphase: ', int2str(i)]);
    saveas(gcf,['Fig2a_nobead', int2str(i), '.png'])
end