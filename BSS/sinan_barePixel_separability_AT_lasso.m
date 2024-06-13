addpath(genpath('/users/syilmaz/MATLAB/20230503_quad2020_in_vivo_bss/20220817_BSS_with_angle_sensitivity'))
addpath(genpath('/users/syilmaz/MATLAB/20230503_quad2020_in_vivo_bss/pca_ica/pca_ica'))

%analyzed fluorescent dye (only green):
    %FluoroSpheres, tau = 6ns, Ex = 505, Em = 515, QY = 0.8, SS = 30nm
    
ncolor = 50;
greencolormapHC = ([54*linspace(0,1,ncolor)' 255*linspace(0,1,ncolor)' zeros(ncolor,1)]/255).^1; 

%% step1: make the D-matrix
deadangle = pi/2;
clusterSize = 16;
xpitch = 24.93; % xpitch = 25.67; for dualShank2018........xpitch = 24.93; for quadShank2020
ypitch = 97.00; % ypitch = 73.22; for dualShank2018........ypitch = 97.00; for quadShank2020
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
[angles, angleSenstivity] = loadAngleSensitivity('angle_sensitivity_data_quad2020/air_sum_data', 1, 3);
angleData = [angles; angleSenstivity]';

%in free space
detectionField = measured_volume(X,Y,Z,angleData, NPixel, DProbe);
maxOfEachPixel = max(detectionField,[],2);

detectionFieldColumnNormalized = detectionField./max(detectionField);
detectionFieldAllNormalized = detectionField./max(max(detectionField));

%% step2: create multiple fluorosphere point sources within the volume and measure the shank response.
%load('angle_sensitivity_data/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x25_y15_z30.mat')
%load('angle_sensitivity_data/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x10_y10_z10.mat')

L=4;
[XYZloc,xloc,x_int] = create_random_beads_3D(X,Y,Z,L,X(1),Y(1),Z(1));

XINTG = reshape(x_int,length(X),length(Y),length(Z));
figure(28); clf;
scatter3(XYZloc(:,1),XYZloc(:,2),XYZloc(:,3),'g', 'o', 'filled');
axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]);

%superposition of these beads' shank response will create the final measurement.
% using detectionField:
linearizedShankResponse = detectionField*x_int;
figure(20); subplot(L+1,1,1);
plot(linearizedShankResponse); xlim([1 NPixel]); title(sprintf('Total Response'));

% if we'd like to see individual beads' effects:
beadResponses = zeros(NPixel, L);
for k=1:L
    beadResponses(:,k) = detectionField(:,xloc(k));
    figure(20); subplot(L+1,1,k+1); plot(beadResponses(:,k)); xlim([1 NPixel]); title(sprintf('Bead Location: X=%i um, Y=%i um, Z=%i um',XYZloc(k,1),XYZloc(k,2),XYZloc(k,3)));
end

%% step 3: add noises and take multiple frames/samples

%load('angle_sensitivity_data/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x25_y15_z30.mat')
%load('angle_sensitivity_data/Plots/simulated_shank/d_pixel_angle_sensitivity_data_resolution_x10_y10_z10.mat')

%add the time-gate to the sources
cntrate = 1e6; %non-dcr counts/sec
%scale by 1e12 for correct usage of imnoise 
cntscaler = cntrate/1e12;
simSIG = 1; % simulate noises for signal data, iterates for the number of samples.
simMAP = 1; % simulate noises for mapping data, iterates for the number of cubic units in the volume.

FPS = 200;
time = 30; %sec
numOfSamples = FPS*time;

if simSIG
    temporal_data = createFramesAddNoises(linearizedShankResponse, numOfSamples);
end

% to create a map of all cubic units in the volume
if simMAP
    mappedDetField = zeros(N, numOfSamples); % [numOfCubicUnits x numOfSamples]
    for n=1:N

        % creating frames with different noise levels for detectionField.
        noised_detField = createFramesAddNoises(detectionField(:,n), numOfSamples);
        
        %fprintf('Mapping continues. n = %i\n', n);
        mappedDetField(n,:) = fastICA(noised_detField, 1, 'kurtosis', 0);
    
    end
end
    
%% step4: reconstruct BSS from SIG_blurred (simulation) or measured data N_SPAD(257:512) (empirical)

%load('angle_sensitivity_data_quad2020/Plots/simulated_shank/detection_field_mapped_200fps_30sec_simSIG_simMAP_resolution_x10_y10_z10.mat')

temporal_data = delta_mean_stim;
[Lhat, PCA_scores, Zica, Zica2, W, Wunnoised] = bss_with_pca_and_ica(temporal_data, 95);

XYZlocDetected = zeros(Lhat,3);
xlocDetected = zeros(Lhat,1);
for b=1:Lhat
    mappedDiff = mappedDetField - Zica2(b);
    RMSE = sqrt(sum(mappedDiff.^2,2));
    [val, ind] = min(RMSE);
    xlocDetected(b) = ind;
    
    [x, y, z] = ind2sub(sizes,ind);
    XYZlocDetected(b,1) = X(x); XYZlocDetected(b,2) = Y(y); XYZlocDetected(b,3) = Z(z);
end

figure(28); hold on; grid;
scatter3(XYZlocDetected(:,1),XYZlocDetected(:,2),XYZlocDetected(:,3),'r', '*');
axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]);

%% step4v2: use low-pass filtering to determine width, diff and relative intensity to find fluorophores.

%load('angle_sensitivity_data_quad2020/Plots/simulated_shank/detection_field_mapped_200fps_30sec_simSIG_simMAP_resolution_x10_y10_z10.mat')

%load('2bead_data/test_data.mat'); linearizedShankResponse = z50u_y400u_xONSPAD(1:128)';
%load('20230423_in_vivo_workdir/1-3_stim_delta_image.mat'); linearizedShankResponse = delta_mean_stim;

%load('20230512_GFP_workdir/savedata/20230603_all_data.mat');

%load('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20230519_10um_bead_video_workdir/savedata/10um_bead_video_day6_3.mat');

% % % % %     XYZlocDetected_single_bead_all = ones(size(N_SPAD_all,2), 3, 2).*(NaN); % two beads are even detected here.
XYZlocDetected_single_bead_all = ones(size(N_SPAD_all,2), 3).*(NaN); % only single bead is allowed to be detected.

    for ff=1+0*400:size(N_SPAD_all,2)
    display(ff);
    %linearizedShankResponse = rollingWindowHotPixRemoval(N_SPAD_all(:,ff),5,3);
    linearizedShankResponse = N_SPAD_all(:,ff);

    %figure(7); plot(linearizedShankResponse); xlim([1 NPixel]); hold on; title('Linearized Shank Response')
    %[autoCorr, w, diffHeight] = autocorrAndDiff(linearizedShankResponse, L);

    smoothData = smooth(linearizedShankResponse,20);
    normalizedSmoothData = smoothData./max(smoothData);

    %figure(1); plot(normalizedSmoothData); title('Normalized Smooth Data'); xlim([1 NPixel]); hold on;
    [peaks, peakLocs, peakWidths, peakProminences] = findpeaks(normalizedSmoothData,'SortStr','descend','MinPeakProminence',0.1);
    NPeaks = length(peaks);
    fprintf('We have %i number of neurons found.\n', NPeaks);
    
    % create a window centered around each peak with width to extract individual peaks
    signalsExtracted = zeros(NPeaks, NPixel);
    signalsExtractedNormalized = zeros(NPeaks, NPixel);
    XYZlocDetected = zeros(NPeaks,3);
    xlocDetected = zeros(NPeaks,1);

    thres = 1.5; % adjust it so that there is at least one peak satisfying this threshold condition.
    for i=1:NPeaks
        window = zeros(NPixel,1);
        if peakLocs(i)-ceil(peakWidths(i)) >= 1 && peakLocs(i)+ceil(peakWidths(i)) <= NPixel
            window(peakLocs(i)-ceil(peakWidths(i)):peakLocs(i)+ceil(peakWidths(i))) = 1;
        elseif peakLocs(i)-ceil(peakWidths(i)) < 1
            window(1:peakLocs(i)+ceil(peakWidths(i))+10) = 1;
        elseif peakLocs(i)+ceil(peakWidths(i)) > NPixel
            window(peakLocs(i)-ceil(peakWidths(i)-10):end) = 1;
        end

        signalsExtracted(i,:) = linearizedShankResponse .* window(1:NPixel);

        %figure(21); subplot(NPeaks,1,i); sgtitle('Extracted Bead Windows'); plot(signalsExtracted(i,:)); xlim([1 NPixel]);hold on;% title(sprintf('Bead Location: X=%i um, Y=%i um, Z=%i um',XYZloc(k,1),XYZloc(k,2),XYZloc(k,3)));

        signalsExtractedNormalized(i,:) = signalsExtracted(i,:) ./ max(signalsExtracted(i,:));
        [B, B_fitinfo] = lasso(detectionFieldColumnNormalized, signalsExtractedNormalized(i,:), 'CV', 3, 'numLambda', 100, 'maxIter', 1e2);
        [vals, inds] = maxk(B(:,B_fitinfo.Index1SE), 5);
        [val, ind] = max(B(:,B_fitinfo.Index1SE));
        
        %diff = detectionFieldColumnNormalized' - signalsExtractedNormalized(i,:);
        %RMSE = sqrt(sum(diff.^2,2));
        %%[val, ind] = min(RMSE);
        %[vals, inds] = mink(RMSE, 20);
        fprintf('Bead Location: %i, Minimum error: %f\n', ind, val);

        if val <= thres % if the error is within acceptable range, accept it and write it in the list.
            xlocDetected(i) = ind;
            % if new found location is within the best 5 guesses, assign the previous location.
            %if find(inds==previousLoc)
                
            %    XYZlocDetected(i,1) = mode(X(x)); XYZlocDetected(i,2) = mode(Y(y)); XYZlocDetected(i,3) = mode(Z(z));
            %else
                [x, y, z] = ind2sub(sizes,inds);
                XYZlocDetected(i,1) = mode(X(x)); XYZlocDetected(i,2) = mode(Y(y)); XYZlocDetected(i,3) = mode(Z(z));
            %end
            
        elseif i ~= 1 % if the error is too high, improve the guess by arranging the z distance comparing it with the highest peak.
            xlocDetected(i) = ind;
            
            [x, y, z] = ind2sub(sizes,ind);
            XYZlocDetected(i,1) = X(x); XYZlocDetected(i,2) = Y(y);% XYZlocDetected(i,3) = Z(z);

            maxSignal = max(signalsExtracted(1,:)); % max signal in the data measured. Found in the first peak.
            medianPeak = median(nonzeros(signalsExtracted(i-1,:))); % median of the previous peak (i-1). Which gives us the hidden peak (i) within.
            zFirstPeak = XYZlocDetected(1,3);

            newZ = int64(sqrt(maxSignal/medianPeak) * zFirstPeak); % A/(z^2) = maxSignal & A/(newZ^2) = medianPeak

            XYZlocDetected(i,3) = newZ;

            fprintf('This error is higher than the threshold (%.2f). Scaling the z guess based on peak ratios.\nNew bead location: %i\n', thres, newZ);
        else
            xlocDetected(i) = ind;

            [x, y, z] = ind2sub(sizes,ind);
            XYZlocDetected(i,1) = NaN; XYZlocDetected(i,2) = NaN; XYZlocDetected(i,3) = NaN;
            fprintf('Error of the first peak is higher than the threshold (%.2f). No bead found or the distance is higher than the z limit %i um. Prominence threshold can be adjusted. Still writing it as a XYZlocDetected.\n', thres, Z(end));
        end

    end
    
    if NPeaks ==0
        continue
    %%%%% elseif NPeaks ==1
        %%%%% XYZlocDetected_single_bead_all(ff,:) = XYZlocDetected(1,:);
    else
        XYZlocDetected_single_bead_all(ff,:) = XYZlocDetected(1,:);
        %%%%% XYZlocDetected_single_bead_all(ff,:,2) = XYZlocDetected(2,:);
    end
    
    previousLoc = inds(1);
    
    %figure(28); hold on; grid;
    %scatter3(XYZlocDetected(:,1),XYZlocDetected(:,2),XYZlocDetected(:,3),'r', '*'); grid;
    %axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;
    end

%% to create a video from a previously analyzed array, through plotting figures

v = VideoWriter('10um_bead_video_day6_3_single_bead_allowed_threshold_1_point_5_modexyz_search40.avi');
v.FrameRate = 400;
open(v);

for ff = 1: size(XYZlocDetected_single_bead_all, 1)
    
    try
        figure(28);
        scatter3(XYZlocDetected_single_bead_all(ff,1,1),XYZlocDetected_single_bead_all(ff,2,1),XYZlocDetected_single_bead_all(ff,3,1), 'g', 'filled', 'hexagram');
        axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]);
        % if you'd like to plot the second bead coming in and out of the view from time to time, uncomment below
        hold on; scatter3(XYZlocDetected_single_bead_all(ff,1,2),XYZlocDetected_single_bead_all(ff,2,2),XYZlocDetected_single_bead_all(ff,3,2), 'g', 'filled', 'pentagram'); grid;
        axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;
        
        grid; frame = getframe(gcf);
        writeVideo(v, frame);
        clf(28,'reset');
    catch
        try
            clf(28,'reset');
            figure(28);
            scatter3(XYZlocDetected_single_bead_all(ff,1,1),XYZlocDetected_single_bead_all(ff,2,1),XYZlocDetected_single_bead_all(ff,3,1), 'g', 'filled', 'hexagram');
            axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;

            grid; frame = getframe(gcf);
            writeVideo(v, frame);
        catch
            clf(28,'reset');
            figure(28); % if not found, assign the previous one here.
            scatter3(NaN,NaN,NaN, 'g', 'filled', 'hexagram');
            axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid;

            grid; frame = getframe(gcf);
            writeVideo(v, frame);
        end
    end
        
end

close(v);