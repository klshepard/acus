%% Separate into channels and plot them row by row.
    NPixel = 128;
    
    pixelN_per_channel = 16;
    n_channels = NPixel/pixelN_per_channel;
    fps = 400;
    
    %channeled_N_SPAD_all  = reshape(N_SPAD_all, pixelN_per_channel, n_channels, timepoints);
    
    %channeled_avgImage = reshape(avgImage, pixelN_per_channel, n_channels);
    %channeled_subtracted = channeled_N_SPAD_all - channeled_avgImage;
    %channeled_subtracted_normalized = (channeled_N_SPAD_all - channeled_avgImage) ./ channeled_avgImage;
    
    %channeled_temporal_spike_data = squeeze(sum(channeled_N_SPAD_all, 1));
    
    hp_fpass = 0.1;
    
    %y = highpass(channeled_temporal_spike_data', hp_fpass, fps)';
    y = data(:,2:end);
    
    figure(25); set(gcf, 'Position', [50 50 700 700]);
    for c=1:n_channels
        plot((1:size(y,2))./fps, y(c,1:end)*power(c,1)+max(y(:,1:end)+8000, [], 'all')*(c-1), 'linewidth', 1.5); ylim([-10000 100000]);
        hold on;
    end
    title(['Data filtered with ' num2str(hp_fpass) ' Hz HPF' ]); ylabel('Photon Count'); xlabel('Time (sec)');