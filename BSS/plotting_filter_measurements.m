figure(908);
semilogy(lam, ECI_only_0./noFilter, 'LineWidth', 2, 'color', '#2b8cbe'); hold on;
%semilogy(lam, ECI_only_15./noFilter, 'LineWidth', 2, 'color', '#');
semilogy(lam, ECI_only_30./noFilter, 'LineWidth', 2, 'color', '#');
semilogy(lam, ECI_only_45./noFilter, 'LineWidth', 2, 'color', '#');

set(gca, 'box', 'off');
set(gca,'LineWidth',2);
set(gca,'TickLength',[0.01 0.01]);
set(gca,'FontSize', 24, 'FontWeight','normal'); %'bold' or 'normal'
legend('Normal', '30°', '45°', 'location', 'best');
ylim([0.006, 1]);
xlabel('Wavelength (nm)'); ylabel('Filter Transmission');