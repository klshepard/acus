%% Measurement: Fluorescent imaging
clear; clc; fclose('all');
cdir = 'C:\Users\opticsroom\Documents\sinan\20230420_quad2020_in_vivo';
addpath(genpath(cdir))
cd(cdir)

% addpath(genpath('C:\Users\opticsroom\Documents\jaebin\instruments'));
addpath(genpath('C:\adriaan\Devices'))
addpath(genpath('C:\adriaan\64 Standard functions') )

%flashfile = 'bitfiles\Nov2020_QP_v3_pll_intclk_sinan_1024pixels.bit'; clkdiv = 40; fvco = 800e6; fspadon = 20e6;
%flashfile = 'bitfiles\Nov2020_QP_v3_pll_extclk_sinan_1024pixels.bit'; clkdiv = 12; fspadon = 80e6; fvco = 960e6;

%NPixel = 1024;
NPixel = 128;

shankNumber = 1;
%flashfile = 'bitfiles\Nov2020_QP_v3_pll_extclk_sinan_128pixels_shank0.bit'; clkdiv = 12; fspadon = 80e6; fvco = 960e6;
flashfile = 'bitfiles\Nov2020_QP_v3_pll_extclk_sinan_128pixels_shank1.bit'; clkdiv = 12; fspadon = 80e6; fvco = 960e6;
%flashfile = 'bitfiles\Nov2020_QP_v3_pll_extclk_sinan_128pixels_shank2.bit'; clkdiv = 12; fspadon = 80e6; fvco = 960e6;
%flashfile = 'bitfiles\Nov2020_QP_v3_pll_extclk_sinan_128pixels_shank3.bit'; clkdiv = 12; fspadon = 80e6; fvco = 960e6;

python_script = 'C:/Users/opticsroom/Documents/sinan/Transferred_from_Confocal/python_2020_QP_sinan/pyspadoledv7_sinan.py';
%exefile = 'exefiles\pyspadoledv7.exe';
%exefile = 'exefiles\pyspadoledv5.exe';

workdir = '20221208_test_workdir';

%% Acquire image
run = 1;
expname = 'imaging_test'; %expname = 'imaging_test'; % if filename == 'imaging_test_1', we delete the rawdata created. Not saved.
rawdir = [workdir '/rawdata/'];
tsdir = [workdir '/timestamp/'];
filename = [expname '_' num2str(run)];

plotOut = 0;
%spadphaselist = 0:4:360;
%for i=1:length(spadphaselist)

spadphase = 180; numPhases = length(spadphase);
spadduty  = 50; %RANGE 1 t o 99 ideally. But for quad2020, 40-99. don't change it. 50 works the best.

%clkdivhex = dec2hex(clkdiv);

numPatterns = 1; % no change to this. This resets everything in Python and then takes a whole another pattern again. Not continuous data.

timepoints = 400*0.025; %500*10; % FPS: 500 Hz How many snaps in total?. 100 MB in total.
usedframes = 40;
usedframesConsidNPixel = usedframes / (1024/NPixel); % we use this because Python code always expects 1024 pixels while determining transfersize.
                                                     % If you see an image with full imager but not individual shanks, this must be the reason.
                                                     % The data sent to the chip via patternFrames is determined by usedframesConsidNPixel. We effectively sum up more frames.
patternFrames = usedframesConsidNPixel * timepoints; % + ignoreFrames = 2; always. Defined in Python code.
framePeriods = 5000; % start from lower. Because this could lead to saturation at the beginning of shank. Which results in very low response at the distal end of the shank where we measure.

sizeKB_per_TimePoint = 2*usedframesConsidNPixel;
sizeMB_total = sizeKB_per_TimePoint*timepoints/1024;

fps = 1/(1/fspadon * framePeriods * usedframes);

restPeriods = 0; soulPeriods = 0; oleddiv = 10;

[N_SPAD, nframes, lambdam, lambdav, ~, fullData, TS] = ...
            program_probe_uLED_SPADv4_sinan(python_script, flashfile, rawdir, filename, tsdir, [], clkdiv, numPatterns, patternFrames, framePeriods, restPeriods, soulPeriods, spadphase, spadphase, spadduty, oleddiv, [], numPhases, fspadon, timepoints, plotOut, NPixel, shankNumber);
       
N_SPAD_hpRem = rollingWindowHotPixRemoval(N_SPAD,5,5);
figure(3); plot(N_SPAD_hpRem); %xlim([1 256]); ylim([0 2000]); title(['Phase ', num2str(spadphaselist(i))]); %hold on; plot(N_SPAD); legend('Hot pix removed', 'Original');
%saveas(gcf,['savephaserun/Phase_', num2str(spadphaselist(i)), '.png']);

%end
%save([workdir '/savedata/' filename '_' num2str(3) '.mat']);

%end

%     framePeriods = 500;
%     restPeriods = 2000;
%     patternFrames = 20;
%     numPatternsInput = 1; %datasize =  2kB*patternFrames*numpatterns
%     numPatterns = 1;
%     K = 16;
%     dutycycle = 20;
%     sumInverseIm = 0;
%     VANLED = 0;
%     IANLED = 0;
%     oleddiv = 500;