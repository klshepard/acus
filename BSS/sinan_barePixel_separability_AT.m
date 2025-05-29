addpath(genpath('/users/syilmaz/MATLAB/20230503_quad2020_in_vivo_bss/20220817_BSS_with_angle_sensitivity'))
addpath(genpath('/users/syilmaz/MATLAB/20230503_quad2020_in_vivo_bss/pca_ica/pca_ica'))
addpath(genpath('/users/syilmaz/MATLAB/20230503_quad2020_in_vivo//users/syilmaz/MATLAB/20230420_quad2020_in_vivo'))

%analyzed fluorescent dye (only green):
    %FluoroSpheres, tau = 6ns, Ex = 505, Em = 515, QY = 0.8, SS = 30nm
    
ncolor = 50;
greencolormapHC = ([54*linspace(0,1,ncolor)' 255*linspace(0,1,ncolor)' zeros(ncolor,1)]/255).^1; 

%% step1: make the D-matrix
deadangle = pi/2;
clusterSize = 16;
xpitch = 25.67; % xpitch = 25.67; for dualShank2018........xpitch = 24.93; for quadShank2020
ypitch = 73.22; % ypitch = 73.22; for dualShank2018........ypitch = 97.00; for quadShank2020
NPixel = 128;
DProbe = create_Probe(NPixel,xpitch,ypitch,0,0);

% %for the 2D image
xres = 10;
yres = 10;
zres = 10;
X = (-100+min(DProbe(:,1))):xres:(max(DProbe(:,1))+100); % extension is 100um for both X..
Y = (-100+min(DProbe(:,2))):yres:(max(DProbe(:,2))+100); % .. and Y.
Z = 20:zres:200;

sizes = [length(X) length(Y) length(Z)];
N = prod(sizes);
fprintf('N = %i\n', N);

% for dualShank2018 angle sensitivity. Normalized.
%[angles, angleSenstivity] = loadAngleSensitivity('angle_sensitivity_data_quad2020/air_sum_data', 1, 3);
[angles, angleSenstivity] = loadAngleSensitivity('angle_sensitivity_data_dualShank2018/di_water___air_sum_data', 1, 3);

angleData = [angles; angleSenstivity]';

%in free space
detectionField = measured_volume(X,Y,Z,angleData, NPixel, DProbe);
maxOfEachPixel = max(detectionField,[],2);

detectionFieldColumnNormalized = detectionField./max(detectionField);
detectionFieldAllNormalized = detectionField./max(max(detectionField));

%% step2: create multiple fluorosphere point sources within the volume and measure the shank response.
% %load('angle_sensitivity_data_dualShank2018/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x25_y15_z30.mat')
% %load('angle_sensitivity_data_dualShank2018/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x10_y10_z10.mat')
% 
% L=3;
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
% 
%% step 3: add noises and take multiple frames/samples
% 
% %load('angle_sensitivity_data/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x25_y15_z30.mat')
% %load('angle_sensitivity_data_dualShank2018/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x10_y10_z10.mat')
% 
% %add the time-gate to the sources
% cntrate = 1e6; %non-dcr counts/sec
% %scale by 1e12 for correct usage of imnoise 
% cntscaler = cntrate/1e12;
% simSIG = 1; % simulate noises for signal data, iterates for the number of samples.
% simMAP = 1; % simulate noises for mapping data, iterates for the number of cubic units in the volume.
% 
% FPS = 400;
% time = 5; %sec
% numOfSamples = FPS*time;
% 
% if simSIG
%     temporal_data = createFramesAddNoises(linearizedShankResponse, numOfSamples);
% end
% 
% % to create a map of all cubic units in the volume
% if simMAP
%     mappedDetField = zeros(N, numOfSamples); % [numOfCubicUnits x numOfSamples]
%     for n=1:N
% 
%         % creating frames with different noise levels for detectionField.
%         noised_detField = createFramesAddNoises(detectionField(:,n), numOfSamples);
%         
%         %fprintf('Mapping continues. n = %i\n', n);
%         mappedDetField(n,:) = fastICA(noised_detField, 1, 'kurtosis', 0);
%     
%     end
% end
%     
%% step4: reconstruct BSS from SIG_blurred (simulation) or measured data N_SPAD(257:512) (empirical)
% % % % % % % % 
% % % % % % % % %load('angle_sensitivity_data_quad2020/Plots/simulated_shank/detection_field_mapped_200fps_30sec_simSIG_simMAP_resolution_x10_y10_z10.mat')
% % % % % % % % 
% % % % % % % % temporal_data = delta_mean_stim;
% % % % % % % % [Lhat, PCA_scores, Zica, Zica2, W, Wunnoised] = bss_with_pca_and_ica(temporal_data, 95);
% % % % % % % % 
% % % % % % % % XYZlocDetected = zeros(Lhat,3);
% % % % % % % % xlocDetected = zeros(Lhat,1);
% % % % % % % % for b=1:Lhat
% % % % % % % %     mappedDiff = mappedDetField - Zica2(b);
% % % % % % % %     RMSE = sqrt(sum(mappedDiff.^2,2));
% % % % % % % %     [val, ind] = min(RMSE);
% % % % % % % %     xlocDetected(b) = ind;
% % % % % % % %     
% % % % % % % %     [x, y, z] = ind2sub(sizes,ind);
% % % % % % % %     XYZlocDetected(b,1) = X(x); XYZlocDetected(b,2) = Y(y); XYZlocDetected(b,3) = Z(z);
% % % % % % % % end
% % % % % % % % 
% % % % % % % % figure(28); hold on; grid;
% % % % % % % % scatter3(XYZlocDetected(:,1),XYZlocDetected(:,2),XYZlocDetected(:,3),'r', '*');
% % % % % % % % axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]);

%% step4_v2: Custom BSS

%load('angle_sensitivity_data_quad2020/Plots/simulated_shank/detection_field_mapped_200fps_30sec_simSIG_simMAP_resolution_x10_y10_z10.mat')

%load('2bead_data/test_data.mat'); linearizedShankResponse = z50u_y400u_xONSPAD(1:128)';
%load('20230423_in_vivo_workdir/1-3_stim_delta_image.mat'); linearizedShankResponse = delta_mean_stim;

%load('20230512_GFP_workdir/savedata/20230603_all_data.mat');

load('20230909_in_vivo_GFP_workdir/savedata/day_2_animal_1_moving_shank_4_all_data.mat'); fps = 40;
%load('20230909_in_vivo_GFP_workdir/savedata/day_2_animal_1_moving_shank_4_all_data_400fps.mat'); fps = 400;

%load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20230519_10um_bead_video_workdir/savedata/10um_bead_video_day6_3.mat');

%load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20230519_10um_bead_video_workdir/20241126_revisions_four_beads/N_SPAD_four_beads.mat'); N_SPAD_all = N_SPAD; N_SPAD_all_diff = N_SPAD;

%load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_2_28-29_40fps.mat'); % N_SPAD_all_diff
%load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_2_36&38_40fps.mat'); % N_SPAD_all_diff
%load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_2_28-29_400fps_detrended_filtered.mat'); % N_SPAD_all_diff

%%%%detection

smoothWindow = 4; %4-30
minPkProm = 0.1; 
threshold = 1.5; %1-1.25-1.5
spatialLPF = 40; %N -> 40

% Reconstruct the bead video

% in vivo structural imaging
XYZlocDetected_single_bead_all = BSS_bead_video(N_SPAD_all_diff, NPixel, detectionFieldColumnNormalized, sizes, X,Y,Z, smoothWindow, minPkProm, threshold, spatialLPF); %12, 0.35, 1.5, 20);
%save('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20230519_10um_bead_video_workdir/savedata/day2_animal_1_moving_shank_4_40fps_bss_smooth20_prom04_thres1_last40.mat', '-v7.3');
% in vivo functional imaging
%save('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_1_102_400fps_detrended3_filtered_3D_detected_all_voxels.mat','N_SPAD_all_diff', 'XYZlocDetected_single_bead_all', 'all_voxels_RMSE', '-v7.3');

%% step5_to create a video from a previously analyzed array, through plotting figures
% % % fps = 40;
% % % fpsReductionRatio = 1;
% % % 
% % % v1 = VideoWriter(['in_vivo_sst_gcamp6s_animal_1_102_deltaF_video_xyz_' num2str(fps/fpsReductionRatio) 'fps_smoWin' num2str(smoothWindow) '_thres' num2str(threshold) '_spaLPF' num2str(spatialLPF) '.avi']);
% % % v1.FrameRate = fps/fpsReductionRatio;
% % % open(v1);
% % % %v2 = VideoWriter(['in_vivo_sst_gcamp6s_animal_1_102_deltaF_video_xy_' num2str(fps/fpsReductionRatio) 'fps.avi']);
% % % %v2.FrameRate = fps/fpsReductionRatio;
% % % %open(v2);
% % % 
% % % XYZlocDetected_single_bead_all = blockproc(XYZlocDetected_single_bead_all, [fpsReductionRatio 3], @(x) mean(x.data, 1));
% % % 
% % % for ff = 1: size(XYZlocDetected_single_bead_all, 1)
% % %     disp(ff);
% % %     try
% % %         figure(28);
% % %         scatter3(XYZlocDetected_single_bead_all(ff,1,1),XYZlocDetected_single_bead_all(ff,2,1),XYZlocDetected_single_bead_all(ff,3,1), 200, 'g', 'filled', 'hexagram'); grid;
% % %         axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid; % For 3D view, view([-10 15]); grid;
% % % % % if you'd like to plot the second bead coming in and out of the view from time to time, uncomment below
% % % % %         %hold on; scatter3(XYZlocDetected_single_bead_all(ff,1,2),XYZlocDetected_single_bead_all(ff,2,2),XYZlocDetected_single_bead_all(ff,3,2), 200, 'g', 'filled', 'pentagram'); grid;
% % % % %         %axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;
% % %         frame = getframe(gcf);
% % %         writeVideo(v1, frame);
% % %         
% % %         %view([0 90]); grid;
% % %         %frame = getframe(gcf);
% % %         %writeVideo(v2, frame);
% % %         clf(28,'reset');
% % %         
% % %     catch
% % %         try
% % %             clf(28,'reset');
% % %             figure(28);
% % %             scatter3(XYZlocDetected_single_bead_all(ff,1,1),XYZlocDetected_single_bead_all(ff,2,1),XYZlocDetected_single_bead_all(ff,3,1), 200, 'g', 'filled', 'hexagram'); grid;
% % %             axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;
% % %             frame = getframe(gcf);
% % %             writeVideo(v1, frame);
% % %             
% % %             %view([0 90]); grid;
% % %             %frame = getframe(gcf);
% % %             %writeVideo(v2, frame);
% % %         catch
% % %             clf(28,'reset');
% % %             figure(28);
% % %             scatter3(NaN,NaN,NaN, 200, 'g', 'filled', 'hexagram'); grid;
% % %             axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;
% % %             frame = getframe(gcf);
% % %             writeVideo(v1, frame);
% % %             
% % %             %view([0 90]); grid;
% % %             %frame = getframe(gcf);
% % %             %writeVideo(v2, frame);
% % %         end
% % %     end
% % %         
% % % end
% % % 
% % % close(v1);
% % % %close(v2);
% % % 
% % % %% to create a video from a previously analyzed array in temporal domain, through plotting various intensities in the same location
% % % NPixel = 128;
% % % fps = 400;
% % % 
% % % v1 = VideoWriter(['in_vivo_sst_gcamp6s_animal_1_102_deltaF_video_xyz_' num2str(fps) 'fps_detected_intensityMask_sized.avi']);
% % % v1.FrameRate = fps;
% % % open(v1);
% % % for i=70:10:70
% % %     for j =30:10:30
% % %         for k=150:10:150
% % %             XYZ_location = [i, j, k]; index = sub2ind(sizes,find(X==XYZ_location(1)), find(Y==XYZ_location(2)), find(Z==XYZ_location(3))); % where the neuron is detected previously.
% % %             
% % %             %pseudoInverseDecField = pinv(detectionFieldColumnNormalized);
% % %             %sourceBackpropagated_temporal = pseudoInverseDecField(index,:)*N_SPAD_all_diff;
% % %             %figure; plot((1:length(sourceBackpropagated_temporal))./fps,movmean(sourceBackpropagated_temporal,50)); title(num2str(XYZ_location));
% % % 
% % %             sourceBackpropagated_temporal = detectionFieldColumnNormalized(:,index)'*N_SPAD_all_diff;
% % %             sourceBackpropagated_temporal_normalized = normalize(sourceBackpropagated_temporal, 'range', [1e-9 1]);
% % %         end
% % %     end
% % % end
% % % 
% % % for ff = 1: length(sourceBackpropagated_temporal_normalized)
% % %     disp(ff);
% % %     figure(28);
% % %     scatter3(XYZ_location(1),XYZ_location(2),XYZ_location(3), 200*sourceBackpropagated_temporal_normalized(ff), 'g', 'filled', 'hexagram', 'MarkerFaceAlpha', sourceBackpropagated_temporal_normalized(ff));
% % %     axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]);
% % %     grid on;
% % %     frame = getframe(gcf);
% % %     writeVideo(v1, frame);
% % %     clf(28,'reset');
% % % end
% % % close(v1);