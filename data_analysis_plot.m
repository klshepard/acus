% to plot the in vivo data section that is relevant vs. time
%% to get the whole data into a single matrix

cdir = '/users/syilmaz/MATLAB/20230420_quad2020_in_vivo';
addpath(genpath(cdir));

%workdir = '20230423_in_vivo_workdir';
%workdir = '20230909_in_vivo_GFP_workdir';
workdir = '20241105_revisions_in_vivo_sparse_GCaMP_workdir';

AllData = 102;
%kstestResults = double.empty(0,3);
%videoDataList = [25, 26, 28, 29, 36, 38, 41, 48 49];

runList = AllData;
for i = 1:length(runList)
%close all;
run = runList(i); 
%expname = 'in_vivo_gcamp_insertion_1'; %expname = 'imaging_test'; % if filename == 'imaging_test_1', we delete the rawdata created. Not saved.
%expname = 'animal_1_moving_shank';
expname = 'in_vivo_sst_gcamp6s_animal_1';
rawdir = [workdir '/rawdata/'];
tsdir = [workdir '/timestamp/'];
filename = [expname '_' num2str(run)];

NPixel = 128;
%NPixel = 1024;

%plotOut = 0;

timepoints = 400*15; % FPS: 500 Hz ==> 1 sec = 10 MB.
usedframes = 80; % #frames. First ~10 are contaminated by RAM.10

framePeriods = 2500;
fspadon = 80e6;
fps = 1/(1/fspadon * framePeriods * usedframes); %fps = fps*3/4;
spadduty  = 50; %RANGE 1 to 99
spadphase = 180; numPhases = length(spadphase);

%quad2020
quadtype = 2;     %2020QP
shankNumber = 1;
numPatterns = 1; % no change to this. This resets everything in Python and then takes a whole another pattern again. Not continuous data.
usedframesConsidNPixel = usedframes / (1024/NPixel);
patternFrames = usedframesConsidNPixel * timepoints;
restPeriods = 0; soulPeriods = 0; oleddiv = 10;

sizeKB_per_TimePoint = 2*usedframesConsidNPixel;
sizeMB_total = sizeKB_per_TimePoint*timepoints/1024;

numsnaps = 1;
N_SPAD_all = zeros(NPixel, timepoints);
fullData_all = zeros(NPixel, usedframes, timepoints);
nframes_all = zeros(timepoints,1);

allSnaps = zeros(NPixel, timepoints, numsnaps);

for j=1:numsnaps
    
    ignoreframes = 2*(1024/NPixel); %ignore first 10, to compensate for RAM
   
    %Open timestamp data into matlab vector
    fileID_t = fopen([tsdir filename '_timestamp.txt']);
    TS = textscan(fileID_t,'%s','delimiter','\n'); 
    TS = str2double(TS{1});
    fclose(fileID_t);
    
   %Open intensity data into matlab 2D-matrix
    fileID = fopen([rawdir filename '_raw']);
    B = fread(fileID,'uint16');
    fclose(fileID);
    %B = dec2bin(A); slow method. Faster with numbers.
    
    % get rid of the ignoreframes here.    
    cutRows = NPixel*ignoreframes;
    B = B(cutRows+1:end,:);
    B_divided = reshape(B,[],timepoints); % 16 bits
    for pp = 1:numPhases
        for nn = 1:numPatterns 
            for kk=1:timepoints
                pattlen = size(B_divided,1)/(numPatterns*numPhases);
                [ N_SPAD, nframes, lambdam, lambdav, fullData, ~] = ASPdataValidationv5_sinan(B_divided((pattlen*numPatterns*(pp-1)+pattlen*(nn-1)+1):(pattlen*numPatterns*(pp-1)+pattlen*nn),kk),quadtype, NPixel, shankNumber);
            
                N_SPAD_all(:,kk) = N_SPAD;
                fullData_all(1:size(fullData,1), 1:size(fullData,2),kk) = fullData;
                nframes_all(kk) = nframes;
            end
        end
    end
    allSnaps(:,:,j) = N_SPAD_all;
    %hot pixels
end

%save([workdir '/savedata/' filename '_400fps.mat'], 'N_SPAD_all','fullData_all','-v7.3')

%% summed frames in time
    % plotting only the region of interest counts, for GFP in vivo tests.
    roi =1:128; % 512, 256, 116, 92 //////// 512, 256, 68, 56
    
    afterFiber = N_SPAD_all(roi,:);
    afterFiber_temporal_spike_data = sum(afterFiber, 1);
    
    figure(30); set(gcf, 'Position', [400 600 500 200]);
    plot((1:length(afterFiber_temporal_spike_data))./fps, afterFiber_temporal_spike_data);
    title(['Total Count in ROI vs. Time [at ' num2str(fps) '  fps]']); xlabel('Time (sec)'); ylabel('Total Count');
    
%% deltaF plotting only the region of interest counts, for in vivo tests.
%     roi =1:48; % 512, 256, 116, 92 //////// 512, 256, 68, 56
%     
%     afterFiber = N_SPAD_all(roi,:);
%     avgImage = mean(afterFiber, 2);
%     afterFiber_avg_subtracted = afterFiber - avgImage;
%     afterFiber_temporal_spike_data = sum(afterFiber_avg_subtracted, 1);
%
%     figure(21);
%     plot((1:length(afterFiber_temporal_spike_data))./fps, afterFiber_temporal_spike_data);
%     title(['Total Count in ROI vs. Time [at ' num2str(fps) '  fps]']); xlabel('Time (sec)'); ylabel('Total Count');
    
%% deltaF plotting only the end of the shank with fiber on, for in vivo tests.
    end_NPixel = 128; % 512, 256, 116, 92 //////// 512, 256, 68, 56
    
    afterFiber = N_SPAD_all(NPixel-end_NPixel+1:end,fps*0+1:end);
    
    avgImage = mean(afterFiber(:,fps*0+1:end), 2); % while averaging, ignore first a few seconds
    %avgImage(avgImage<max(avgImage)*0.1) = 0; % dark side of the moon
    
    afterFiber_avg_subtracted = afterFiber - avgImage;
    afterFiber_avg_subtracted_normalized = (afterFiber - avgImage) ./ avgImage;
    afterFiber_avg_subtracted_normalized(~isfinite(afterFiber_avg_subtracted_normalized)) = NaN; [afterFiber_avg_subtracted_normalized, mask] = fillmissing(afterFiber_avg_subtracted_normalized,'constant', 0);
    afterFiber_temporal_spike_data = sum(afterFiber_avg_subtracted, 1);%, 'omitnan'); % omiting nan, but infinites are there. discontinuities.
    
    lp_fpass = 25;
    afterFiber_temporal_spike_data = detrend(afterFiber_temporal_spike_data, 3); % 1 or 3
    %afterFiber_temporal_spike_data = lowpass(afterFiber_temporal_spike_data', lp_fpass, fps, 'Steepness', 0.9)';
    %[b,a] = butter(4, 25/(fps/2));
    %afterFiber_temporal_spike_data = filtfilt(b,a,afterFiber_temporal_spike_data);
    
    afterFiber_temporal_spike_data = movmean(afterFiber_temporal_spike_data',50)'; %50 for 400fps, 5 for 40 fps.
    afterFiber_temporal_spike_data = detrend(afterFiber_temporal_spike_data', 3)';
    
    %%%%%%%%%% Two-sample Kolmogorov-Smirnov test: kstest2
    %[h, p, ks2stat] = kstest2(afterFiber_temporal_spike_data, noise_a1_34, 'alpha', 0.000001, 'Tail', 'smaller');
    %kstestResults = [kstestResults; h, p, ks2stat];
    
    figure(23); set(gcf, 'Position', [750 600 700 200]);
    plot((1:length(afterFiber_temporal_spike_data))./fps, afterFiber_temporal_spike_data, 'color', 'blue');
    title(['Total Count on imager vs. Time [at ' num2str(fps) '  fps]  - Movmean & Detrended']); xlabel('Time (sec)'); ylabel('\DeltaF');

%%%%%%%%%%%%%% Separate into channels and plot them row by row.
    %NPixel = 128;
    pixelN_per_channel = 16;
    n_channels = NPixel/pixelN_per_channel;
    
    channeled_N_SPAD_all  = reshape(N_SPAD_all(:,fps*0+1:end), pixelN_per_channel, n_channels, timepoints);
    channeled_avgImage = reshape(avgImage, pixelN_per_channel, n_channels);
    channeled_subtracted = channeled_N_SPAD_all - channeled_avgImage;
    channeled_subtracted_normalized = (channeled_N_SPAD_all - channeled_avgImage) ./ channeled_avgImage;
    channeled_temporal_spike_data = squeeze(sum(channeled_subtracted, 1));
    y = channeled_temporal_spike_data;

    figure(25); set(gcf, 'Position', [500 50 700 700]);
    
    % scaling each channel based on its noise std to mimick uniform illumination scenario
    hpfiltered_y = highpass(y', 2, fps)';
    channel_scalar_list_from_std = std(hpfiltered_y(1,:)) ./ std(hpfiltered_y, [], 2);
    y_scaled_for_same_noise = y .* channel_scalar_list_from_std;
    channel_dist = max(peak2peak(y_scaled_for_same_noise,2));
    for c=1:n_channels
        plot((1:timepoints)./fps, y_scaled_for_same_noise(c,:)+channel_dist*(c-1), 'color', 'black');
        hold on;
    end
    title('Original Data'); ylabel('Channel \DeltaF'); xlabel('Time (sec)'); 
    ylim('padded'); set(gca, 'box', 'off'); set(gca, 'YColor', 'none');  %ylim([-2,30]); set(gca,'YTick', [0, 4, 8, 12, 16, 20, 24, 28]);
    %newcolors = ["#0B0" "#00F" "#50F" "#F04" "#0D6" "#C7F" "#35F" "#A9A"]; colororder(gca, newcolors);
    
%%%%%%%%%%%%%% Detrending and Filtering the temporal data
    %y_scaled_for_same_noise = detrend(y_scaled_for_same_noise', 3)';
    %y_scaled_for_same_noise = lowpass(y_scaled_for_same_noise', lp_fpass, fps, 'Steepness', 0.9)';
    %[b,a] = butter(4, 25/(fps/2));
    %y_scaled_for_same_noise = filtfilt(b,a,y_scaled_for_same_noise')';
    
    y_scaled_for_same_noise = movmean(y_scaled_for_same_noise',50)'; %50 for 400fps, 5 for 40 fps.
    y_scaled_for_same_noise = detrend(y_scaled_for_same_noise', 3)';
    
    figure(27); set(gcf, 'Position', [750 50 700 700]);
    
    channel_dist = max(peak2peak(y_scaled_for_same_noise,2));
    for c=1:n_channels
        plot((1:timepoints)./fps, y_scaled_for_same_noise(c,:)+channel_dist*(c-1), 'color', 'black');
        hold on;
    end
    title(['Data filtered with Movmean and Detrended']); ylabel('Channel \DeltaF'); xlabel('Time (sec)');
    ylim('padded'); set(gca, 'box', 'off'); set(gca, 'YColor', 'none'); %ylim([-2,30]); set(gca,'YTick', [0, 4, 8, 12, 16, 20, 24, 28]);
    %newcolors = ["#0B0" "#00F" "#50F" "#F04" "#0D6" "#C7F" "#35F" "#A9A"]; colororder(gca, newcolors);

%%%%%%%%%%%%%% one-plot to see if there is anything at all

figure(190); set(gcf, 'Position', [0 0 400 780]); t = tiledlayout(6,1); title(t, ['Run: ' num2str(run) ' [' num2str(fps) ' FPS]']);
figure(23); fig23ax = gca;
figure(190); sub1 = nexttile; %set(gca, 'YColor', 'none');
sub1b = copyobj(fig23ax.Children, sub1); title(['Sum of All Channels [Movmean & Detrended]']); ylabel('\DeltaF'); ylim('padded'); set(gca,'YTick', [0]);

figure(27); fig27ax = gca;
figure(190); sub2 = nexttile([5,1]); set(gca, 'YColor', 'none'); 
sub2b = copyobj(fig27ax.Children, sub2);
title(['Channel Array [2x' num2str(pixelN_per_channel/2) ']']);
ylim('padded'); %ylim([-2,30]); set(gca,'YTick', [0, 4, 8, 12, 16, 20, 24, 28]);

hline = findobj(190, 'type', 'line'); set(hline,'LineWidth',1);
%exportgraphics(gcf, [workdir '/savedata/' filename '_deltaF_Movmean_Detrended_' num2str(fps) '_fps.png'], 'Resolution', 600);
%save([workdir '/savedata/' filename '.mat'],  "-nocompression");

%% to display data conventionally
% % %     fig_list_to_close = [991, 4, 802, 25, 27, 190];
% % %     for f=1:length(fig_list_to_close)
% % %         close(findobj('type', 'figure', 'number', fig_list_to_close(f)))
% % %     end
% % %     
% % %     plot_timepoints = [1.9, 5]; %seconds
% % %     plot_timepoints = uint16(plot_timepoints .* fps);
% % %     
% % %     for p=1:length(plot_timepoints)
% % %         
% % %         NPixel = 128;
% % %         end_NPixel = 128; % 512, 256, 116, 92 //////// 512, 256, 68, 58
% % % 
% % %         figure(991)
% % %         set(gcf, 'Position', [50 50 700 700])
% % %         ax1 = subplot(4,1,1);
% % %         scatter(repmat((1:NPixel)',nframes_all(plot_timepoints(p)),1), reshape(fullData_all(:,1:nframes_all(plot_timepoints(p)),plot_timepoints(p)),1,[]), 10);  hold on;
% % %         axis([NPixel-end_NPixel+1 NPixel 0 31])
% % %         title('Scatterplot of All 5-bit Captures');
% % % 
% % %         ax2 = subplot(4,1,2);
% % %         scatter(repmat((1:NPixel)',nframes_all(plot_timepoints(p)),1), reshape(fullData_all(:,1:nframes_all(plot_timepoints(p)),plot_timepoints(p)),1,[]), 10);  hold on;
% % %         axis([NPixel-end_NPixel+1 NPixel 0 10])
% % %         title('Scatterplot of All 5-bit Captures');
% % % 
% % %         ax3 = subplot(4,1,3);
% % %         bar(NPixel-end_NPixel+1:NPixel, N_SPAD_all(NPixel-end_NPixel+1:end, plot_timepoints(p)), 2); hold on;
% % %         set(gca,'YScale','log')
% % %         xlim([NPixel-end_NPixel+1 NPixel])
% % %         %xline(456-384, '-.b', 'Fiber', 'LabelHorizontalAlignment', 'left')
% % %         %xline(444-384, 'black', 'Epoxy', 'LabelHorizontalAlignment', 'left')
% % %         title('Accumulation, log');
% % % 
% % %         ax4 = subplot(4,1,4);
% % %         plot(NPixel-end_NPixel+1:NPixel, N_SPAD_all(NPixel-end_NPixel+1:NPixel, plot_timepoints(p)));
% % %         hold on;
% % %         xlim([NPixel-end_NPixel+1 NPixel])
% % %         %xline(456-384, '-.b', 'Fiber', 'LabelHorizontalAlignment', 'left')
% % %         %xline(444-384, 'black', 'Epoxy', 'LabelHorizontalAlignment', 'left')
% % %         title('Accumulation, linear'); 
% % %         
% % %         figure(4); set(gcf, 'Position', [50 270 700 200]);
% % %         plot(NPixel-end_NPixel+1:NPixel, N_SPAD_all(NPixel-end_NPixel+1:NPixel, plot_timepoints(p)));
% % %         xlim([NPixel-end_NPixel+1 NPixel]); hold on;
% % %         
% % %         %[physImg] = raw2SepRows_quad2020_in_vivo(N_SPAD_all(NPixel-end_NPixel+1:NPixel, plot_timepoints(p)));
% % %     end
% % %     figure(991);
% % %     axes(ax1); legend(string(plot_timepoints./fps));
% % %     axes(ax2); legend(string(plot_timepoints./fps));
% % %     axes(ax3); legend(string(plot_timepoints./fps)); %xline(4, '-.b', 'Fiber', 'LabelHorizontalAlignment', 'left'); 
% % %     axes(ax4); legend(string(plot_timepoints./fps)); %xline(4, '-.b', 'Fiber', 'LabelHorizontalAlignment', 'left');
% % %     
% % %     figure(4);
% % %     plot(NPixel-end_NPixel+1:NPixel, avgImage);
% % %     legend([string(plot_timepoints./fps) 'Average Image']);
    
%% Prepare data for BSS
% % %     % max difference with average: 22.64
% % %     stim_timepoints = 19.385:1/fps:19.41;
% % %     stim_timepoints = uint16(stim_timepoints .* fps);
% % %     
% % %     maxStimTimePoint = 12*fps;
% % %     
% % %     all_stim = N_SPAD_all(NPixel-end_NPixel+1:NPixel, stim_timepoints);
% % %     meanStim = mean(all_stim, 2);
% % %     delta_mean_stim = meanStim - avgImage;
% % %     figure(30); plot(meanStim); hold on; plot(avgImage); plot(meanStim - avgImage); legend('Average Stim TP', 'Average Image', 'Delta'); xlim([0 128]);
% % %     figure(31); plot(delta_mean_stim); xlim([0 128]);
% % %     
% % %     highestStim = N_SPAD_all(NPixel-end_NPixel+1:NPixel, maxStimTimePoint);
% % %     delta_highestStim = highestStim - avgImage;
% % %     figure(33); plot(highestStim); hold on; plot(avgImage); plot(highestStim-avgImage); legend('Highest Stim TP', 'Average Image'); xlim([0 128]);
% % %     figure(34); plot(delta_highestStim); xlim([0 128]);
    
%% to create a video from a previously analyzed array, through plotting figures

% v1 = VideoWriter([workdir '/savedata/' filename '_deltaF_LPF' num2str(lp_fpass) '_' num2str(fps) 'fps_1D_non-uniform.avi']); % For moving_shank_in_vivo: v = VideoWriter('day_2_animal_1_moving_shank_4_differential_to_initial_position_400fps_1D.avi');
% v1.FrameRate = fps;
% open(v1);
% v2 = VideoWriter([workdir '/savedata/' filename '_deltaF_LPF' num2str(lp_fpass) '_' num2str(fps) 'fps_2D_non-uniform.avi']); % For moving_shank_in_vivo: v = VideoWriter('day_2_animal_1_moving_shank_4_differential_to_initial_position_400fps_2D.avi');
% v2.FrameRate = fps;
% open(v2);
% 
% avgImage = mean(N_SPAD_all, 2); %avgImage(avgImage<max(avgImage)*0.01) = 0; % dark side of the moon
% N_SPAD_all_diff = N_SPAD_all - avgImage;
% %N_SPAD_all_diff = (N_SPAD_all - avgImage)./avgImage; N_SPAD_all_diff(~isfinite(N_SPAD_all_diff)) = NaN; [N_SPAD_all_diff, mask] = fillmissing(N_SPAD_all_diff,'constant', 0);
% 
% %%%%%%%%%% scaling each pixel based on its noise std to mimic uniform illumination scenario
% %hpfiltered_N_SPAD_all_diff = highpass(N_SPAD_all_diff', 2, fps)';
% %pixel_scalar_list_from_std = std(hpfiltered_N_SPAD_all_diff(1,:)) ./ std(hpfiltered_N_SPAD_all_diff, [], 2);
% %N_SPAD_all_diff = N_SPAD_all_diff .* pixel_scalar_list_from_std;
% 
% %Detrended and Filtered each pixel separately in temporal domain.
% N_SPAD_all_diff = detrend(N_SPAD_all_diff', 3)';
% N_SPAD_all_diff = lowpass(N_SPAD_all_diff', lp_fpass, fps, 'Steepness', 0.9)';
% 
% lim_min = min(N_SPAD_all_diff, [], 'all');
% lim_max = max(N_SPAD_all_diff, [], 'all');
% for ff = 1: timepoints
%     % 1D
%     figure(28); set(gcf, 'Position', [800 200 350 350]);
%     plot(1:NPixel, N_SPAD_all_diff(:, ff), 'LineWidth',1); xlim([1 NPixel]); title(['Time = ' num2str(ff/fps)]); ylim([lim_min lim_max]);
%     
%     frame = getframe(gcf);
%     writeVideo(v1, frame);
%     clf(28,'reset');
%     
%     % 2D
%     raw2SepRows_quad2020_in_vivo(N_SPAD_all_diff(:, ff));
%     figure(802); subplot(2,1,1); caxis([lim_min, lim_max]); title(['Time = ' num2str(ff/fps)]);
%     figure(802); subplot(2,1,2); ylim([lim_min, lim_max]); set(gca,'LineWidth',1); set(gcf, 'Position', [200 200 350 500]);
%     set(subplot(2,1,1), 'Position', [0.1, 0.8, 0.8, 0.1]);
%     set(subplot(2,1,2), 'Position', [0.1, 0.1, 0.8, 0.6]);
%     
%     frame = getframe(gcf);
%     writeVideo(v2, frame);
%     clf(802,'reset');
% end
% close(v1);
% close(v2);

%% deltaF + movmean + detrend + backpropagation

avgImage = mean(afterFiber, 2); %avgImage(avgImage<max(avgImage)*0.01) = 0; % dark side of the moon
N_SPAD_all_diff = afterFiber - avgImage;
N_SPAD_all_diff = movmean(N_SPAD_all_diff',50)';
N_SPAD_all_diff = detrend(N_SPAD_all_diff', 3)';

all_voxels_backpropagated = ones(size(N_SPAD_all_diff,2), size(detectionFieldColumnNormalized,2)).*(NaN);
left_psinv = pinv(detectionFieldColumnNormalized);%, 1e-5);
for t=1:timepoints
    all_voxels_backpropagated(t,:) = left_psinv*N_SPAD_all_diff(:,t);
end

save([workdir '/savedata/' filename '_400fps_detrended3_filtered_backprop.mat'], 'N_SPAD_all','N_SPAD_all_diff','fullData_all','all_voxels_backpropagated','-v7.3');

%% plot backpropagated voxel vs. time
load('/users/syilmaz/MATLAB/20230503_quad2020_in_vivo_bss/20220817_BSS_with_angle_sensitivity/angle_sensitivity_data_quad2020/Plots/simulated_shank/detection_field_mapped_200fps_30sec_simSIG_simMAP_resolution_x10_y10_z10.mat')
load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_1_102_400fps_detrended3_filtered_3D_detected_all_voxels_and_backProp.mat')

XYZ_location_center = [70, 30, 150];
spatial_sum = zeros(size(all_voxels_backpropagated, 1), 1);
resolution=10;
figure;
for rb = 0:10:40
    range_border = rb;

    voxel_count=0;
    for p=-range_border:resolution:range_border
        for q=-range_border:resolution:range_border
            for r=-range_border:resolution:range_border
                XYZ_location = [XYZ_location_center(1)+p, XYZ_location_center(2)+q, XYZ_location_center(3)+r];
                index = sub2ind(sizes,find(X==XYZ_location(1)), find(Y==XYZ_location(2)), find(Z==XYZ_location(3)));
                spatial_sum = spatial_sum + all_voxels_backpropagated(:,index);
                voxel_count = voxel_count+1;
            end
        end
    end

    spatial_avg = spatial_sum/voxel_count;
    subplot(5,1,rb/10+1);
    plot(detrend(movmean(spatial_avg,100),3)); title(['-' num2str(range_border) ' < Voxel < ' num2str(range_border)]); hold on;
    plot(spatial_avg);
    ylim([-0.02, 0.02]);
end
sgtitle(['Backprop. Voxel Coordinates: ' num2str(XYZ_location_center) 'Data' i])

end