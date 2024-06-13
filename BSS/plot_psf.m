%% create one bead in 3D space

% L=1;
% [XYZloc,xloc,x_int] = create_random_beads_3D(X,Y,Z,L,X(1),Y(1),Z(1));
% 
% XINTG = reshape(x_int,length(X),length(Y),length(Z));
% figure(28); clf;
% scatter3(XYZloc(:,1),XYZloc(:,2),XYZloc(:,3),'g', 'o', 'filled');
% axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]);
% 
% %superposition of these beads' shank response will create the final measurement.
% % using detectionField:
% linearizedShankResponse = detectionField*x_int;
% figure(20); subplot(L+1,1,1);
% plot(linearizedShankResponse); xlim([1 NPixel]); title(sprintf('Total Response'));
% 
% % if we'd like to see individual beads' effects:
% beadResponses = zeros(NPixel, L);
% for k=1:L
%     beadResponses(:,k) = detectionField(:,xloc(k));
%     figure(20); subplot(L+1,1,k+1); plot(beadResponses(:,k)); xlim([1 NPixel]); title(sprintf('Bead Location: X=%i um, Y=%i um, Z=%i um',XYZloc(k,1),XYZloc(k,2),XYZloc(k,3)));
% end

%% PSF calculations for a new bead
% xloc ======> 35807...... X 550, Y 90, Z 50
%PSF = detectionField*x_int;
%PSF = pseudoInverseDecField_tol*PSF;
%PSF_3D = reshape(PSF,length(X),length(Y),length(Z));

%[bead_x,bead_y,bead_z] = ind2sub(size(XINTG),xloc);
%XYZ_loc_indices = [bead_x,bead_y,bead_z];

%figure(16); imagesc(rot90(squeeze(PSF_3D(:,:,XYZ_loc_indices(3)))));
%xticklabels({'y = 140','y = 90','y = 40','y = -10','y = -60'});
%yticks([11 14 17]);
%yticklabels({'z = 100','z = 70','z = 40'});
%colorbar;

%% Gaussian-fit and FWHMs for each voxel in 3D
% FWHM_min_x_y_z = NaN(sizes);
% 
% for b=1:N
%     x_int = zeros(N,1);
%     x_int(b) = 1;
%     
%     PSF = detectionField*x_int;
%     PSF = pseudoInverseDecField*PSF;
%     PSF_3D = reshape(PSF,length(X),length(Y),length(Z));
%     
%     [x,y,z]=ind2sub(sizes,b);
%     
%     %gauss1 fit
%     data_x = PSF_3D(:,y,z);
%     data_y = PSF_3D(x,:,z);
%     data_z = PSF_3D(x,y,:);
% 
%     gfit_x = fit(X.', squeeze(data_x), 'gauss1');
%     gfit_y = fit(Y.', squeeze(data_y).', 'gauss1');
%     gfit_z = fit(Z.', squeeze(data_z), 'gauss1');
% 
%     fit_data_x =  gfit_x.a1.*exp(-((X-gfit_x.b1)./gfit_x.c1).^2);
%     fit_data_y =  gfit_y.a1.*exp(-((Y-gfit_y.b1)./gfit_y.c1).^2);
%     fit_data_z =  gfit_z.a1.*exp(-((Z-gfit_z.b1)./gfit_z.c1).^2);
% 
%     [peak_p_x, peak_loc_x, peak_width_x, peak_prom_x] = findpeaks(fit_data_x,X,'WidthReference','halfheight');
%     [peak_p_y, peak_loc_y, peak_width_y, peak_prom_y] = findpeaks(fit_data_y,Y,'WidthReference','halfheight');
%     [peak_p_z, peak_loc_z, peak_width_z, peak_prom_z] = findpeaks(fit_data_z,Z,'WidthReference','halfheight');
%     
%     minimum_widths = [peak_width_x, peak_width_y, peak_width_z];
%     
%     try
%         FWHM_min_x_y_z(x,y,z) = min(minimum_widths);
%     catch
%         FWHM_min_x_y_z(x,y,z) = FWHM_min_x_y_z(x-1,y,z);
%     end
%     
% end

%% Plot FWHMs
for i=1:sizes(3)
    figure(30+i);
    imagesc(rot90(squeeze(FWHM_min_x_y_z(:,:,i))));
    set(gca,'FontSize',12, 'FontWeight','normal','LineWidth',2,'TickLength',[0.01 0.01]);
    xticks([41 111 181 251 321]); xticklabels({'300','1000','1700','2400','3100'});
    yticks([5 10 15 20 25]); yticklabels({'130','80','30','-20','-70'});
    colormap(gca,jet); %colorbar;
    caxis([10, 140]);
    axis off;
    saveas(gcf, ['FWHM_axis_off_' int2str(Z(i)) '.png']);
    %figure(60+i);
    %histogram(FWHM_min_x_y_z(:,:,i), 'BinWidth',  2.5); xlim([10 140]); ylim([0 2100]);
    %saveas(gcf, ['Histogram_FWHM_' int2str(Z(i)) '.png']);
end