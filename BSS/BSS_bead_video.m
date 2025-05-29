function [XYZlocDetected_single_bead_all, all_voxels_RMSE] = BSS_bead_video(N_SPAD_all_diff, NPixel, detectionFieldColumnNormalized, sizes, X,Y,Z, smoothWindow, minPkProm, threshold, spatialLPF)

singleAllowed = 1; twoAllowed = 0;
if singleAllowed == 1
    XYZlocDetected_single_bead_all = ones(size(N_SPAD_all_diff,2), 3).*(NaN); % only single bead is allowed to be detected.
    all_voxels_RMSE = ones(size(N_SPAD_all_diff,2), size(detectionFieldColumnNormalized,2)).*(NaN); % we can save all voxel information if only single bead is allowed to be detected.
elseif twoAllowed == 1
    XYZlocDetected_single_bead_all = ones(size(N_SPAD_all,2), 3, 2).*(NaN); % two beads are even detected here.
    all_voxels_RMSE = ones(size(N_SPAD_all_diff,2), size(detectionFieldColumnNormalized,2)).*(NaN); % we can save all voxel information for the highest peaks.
end
previousLocs = -1; % this is needed if we're not doing video.
    for ff=1:1:size(N_SPAD_all_diff, 2)
        disp(ff); close all;
        linearizedShankResponse = N_SPAD_all_diff(:, ff);
        %linearizedShankResponse = rollingWindowHotPixRemoval(yes_neuron,5,3);
        %linearizedShankResponse = yes_neuron;
        %linearizedShankResponse = rollingWindowHotPixRemoval(z50u_y200u_xONSPAD',3,3);
        %linearizedShankResponse = z50u_y200u_xONSPAD';

        %figure(7); plot(linearizedShankResponse); xlim([1 NPixel]); hold on; title('Linearized Shank Response');
        %[autoCorr, w, diffHeight] = autocorrAndDiff(linearizedShankResponse, L);
        
        % check if the mean of the data is below zero. This makes the code a lot faster. Discard all measurements with mean below zero.
        if mean(linearizedShankResponse)>0
            smoothData = smooth(linearizedShankResponse, smoothWindow);
            normalizedSmoothData = smoothData./max(smoothData);
            %figure(1); plot(normalizedSmoothData); title('Normalized Smooth Data'); xlim([1 NPixel]); %hold on;
        
            [peaks, peakLocs, peakWidths, peakProminences] = findpeaks(normalizedSmoothData,'SortStr','descend','MinPeakProminence',minPkProm,'MinPeakHeight',0,'MinPeakDistance',12,'WidthReference','halfheight','Annotate','extents');
            findpeaks(normalizedSmoothData,'SortStr','descend','MinPeakProminence',minPkProm,'MinPeakHeight',0,'MinPeakDistance',12,'WidthReference','halfheight','Annotate','extents');
        else
            continue
        end
        
        % if the normalized data has a peak in the first a few pixels, 'findpeaks'  cannot find it. Add it manually.
        if normalizedSmoothData(1) == 1 || normalizedSmoothData(2) == 1 || normalizedSmoothData(3) == 1 || normalizedSmoothData(4) == 1
            peaks = [1; peaks];
            peakLocs = [1; peakLocs];
            idx = find(normalizedSmoothData<0.25);  peakWidths = [idx(1); peakWidths];
            try peakProminences = [1-peaks(2); peakProminences]; catch peakProminences = [1; peakProminences]; end
        end
        
        NPeaks = length(peaks);
        fprintf('We have %i number of neurons found.\n', NPeaks);

        % create a window centered around each peak with width to extract individual peaks
        signalsExtracted = zeros(NPeaks, NPixel);
        signalsExtractedNormalized = zeros(NPeaks, NPixel);
        XYZlocDetected = zeros(NPeaks,3);
        xlocDetected = zeros(NPeaks,1);

        thres = threshold; % adjustable. So that there is at least one peak satisfying this threshold condition if prominence criterion is met.
        for i=1:NPeaks
            window = zeros(NPixel,1);
            if peakLocs(i)-ceil(peakWidths(i)) >= 1 && peakLocs(i)+ceil(peakWidths(i)) <= NPixel
                window(peakLocs(i)-ceil(peakWidths(i)):peakLocs(i)+ceil(peakWidths(i))) = 1;
            elseif peakLocs(i)-ceil(peakWidths(i)) < 1
                try window(1:peakLocs(i)+ceil(peakWidths(i)*2)) = 1; catch window(1:peakLocs(i)+ceil(peakWidths(i)*2)+10) = 1; end %window(1:peakLocs(i)+ceil(peakWidths(i)*2)+10) = 1;
            elseif peakLocs(i)+ceil(peakWidths(i)) > NPixel
                try window(peakLocs(i)-ceil(peakWidths(i)*2):end) = 1; catch window(peakLocs(i)-ceil(peakWidths(i)-10):end) = 1; end%window(peakLocs(i)-ceil(peakWidths(i)-10):end) = 1;
            end

            signalsExtracted(i,:) = linearizedShankResponse .* window(1:NPixel);
            if mean(signalsExtracted(i,:))>0
                %figure(21); subplot(NPeaks,1,i); sgtitle('Extracted Bead Windows'); plot(signalsExtracted(i,:)); xlim([1 NPixel]);hold on;% title(sprintf('Bead Location: X=%i um, Y=%i um, Z=%i um',XYZloc(k,1),XYZloc(k,2),XYZloc(k,3)));
                signalsExtractedNormalized(i,:) = signalsExtracted(i,:) ./ max(signalsExtracted(i,:));
                diff = detectionFieldColumnNormalized' - signalsExtractedNormalized(i,:);
                RMSE = sqrt(sum(diff.^2,2)); %Euclidean distance (norm-2)
                if i==1
                    all_voxels_RMSE(ff,:) = RMSE;
                end
                [val, ind] = min(RMSE); disp(val);
            else
                continue
            end
            
            if val <= thres % if the error is within acceptable range, accept it and write it in the list.
                fprintf('Bead Location: %i, Minimum error: %f\n', ind, val);
                fprintf('Peak prom.: %.2f, Peak width: %.2f, RMSE: %.2f.\n', peakProminences(1), peakWidths(1), val);
                
                xlocDetected(i) = ind;
                [vals, inds] = mink(RMSE, spatialLPF);
                [x, y, z] = ind2sub(sizes,inds);
                
                % if new found location is within the best "spatialLPF" guesses, assign the previous location.
                if find(intersect(inds,previousLocs,'sorted'))
                    XYZlocDetected(i,1) = mode(X(x)); XYZlocDetected(i,2) = mode(Y(y)); XYZlocDetected(i,3) = mode(Z(z));
                    %XYZlocDetected(i,1) = int16(mean(X(x))); XYZlocDetected(i,2) = int16(mean(Y(y))); XYZlocDetected(i,3) = int16(mean(Z(z)));
                else
                    % else, take mode of X, Y and Z of best "spatialLPF" guesses.
                    [x, y, z] = ind2sub(sizes,inds);
                    previousLocs = inds;
                    XYZlocDetected(i,1) = mode(X(x)); XYZlocDetected(i,2) = mode(Y(y)); XYZlocDetected(i,3) = mode(Z(z));
                    %XYZlocDetected(i,1) = int16(mean(X(x))); XYZlocDetected(i,2) = int16(mean(Y(y))); XYZlocDetected(i,3) = int16(mean(Z(z)));
                end

            elseif i ~= 1 && singleAllowed ~= 1 % if the error is too high, improve the guess by arranging the z distance comparing it with the highest peak.
                xlocDetected(i) = ind;

                [x, y, z] = ind2sub(sizes,ind);
                XYZlocDetected(i,1) = X(x); XYZlocDetected(i,2) = Y(y); XYZlocDetected(i,3) = Z(z);

                maxSignal = max(signalsExtracted(1,:)); % max signal in the data measured. Found in the first peak.
                medianPeak = median(nonzeros(signalsExtracted(i-1,:))); % median of the previous peak (i-1). Which gives us the hidden peak (i) within.
                zFirstPeak = XYZlocDetected(1,3);

                newZ = int64(sqrt(maxSignal/medianPeak) * zFirstPeak); % A/(z^2) = maxSignal & A/(newZ^2) = medianPeak

                XYZlocDetected(i,3) = newZ;

                fprintf('This error is higher than the threshold (%.2f). Scaling the z guess based on peak ratios.\nNew bead location: %i\n', thres, newZ);
            else
                XYZlocDetected(i,:) = NaN;
                fprintf('Error of the highest peak is higher than the threshold (%.2f). No bead found or the distance is higher than the z limit %i um. Prominence threshold can be adjusted.\n', thres, Z(end));
            end

        end
        %figure(30); hold on; grid; scatter3(XYZlocDetected(:,1),XYZlocDetected(:,2),XYZlocDetected(:,3), 200,'g', 'filled', 'hexagram'); grid; axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;
        if NPeaks ==0
            continue
        %%%%% elseif NPeaks ==1
            %%%%% XYZlocDetected_single_bead_all(ff,:) = XYZlocDetected(1,:);
        else
            XYZlocDetected_single_bead_all(ff,:) = XYZlocDetected(1,:);
            %%%%% XYZlocDetected_single_bead_all(ff,:,2) = XYZlocDetected(2,:);
        end
        disp(XYZlocDetected);

        %figure(31); hold on; grid;
        %scatter3(XYZlocDetected(:,1),XYZlocDetected(:,2),XYZlocDetected(:,3),200,'g','filled','hexagram'); grid;
        %axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;

    end
    
    %scatter3(340,20,40,100,'r','filled','diamond'); grid;
    %scatter3(690,40,50,100,'r','filled','diamond'); grid;
    %scatter3(1005,125,40,100,'r','filled','diamond'); grid;
    %scatter3(1340,90,70,100,'r','filled','diamond'); grid;
    
    %raw2SepRows_quad2020_in_vivo(N_SPAD);
    
    %colors used in the paper
    %7E2F8E
    %D95319
    %77AC30
    
end