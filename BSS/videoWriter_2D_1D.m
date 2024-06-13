
load('20230519_10um_bead_video_workdir/savedata/10um_bead_video_day6_3.mat');

%% 2D
v2 = VideoWriter('20230519_10um_bead_video_workdir/plots/10um_bead_video_day6_3_2D_ymaxed.avi');
v2.FrameRate = 400;
open(v2);

for ff = 1: size(N_SPAD_all,2)
    [physImg] = raw2SepRows_quad2020_in_vivo(N_SPAD_all(:, ff));
    
    frame = getframe(gcf);
    writeVideo(v2, frame);
    
end

close(v2);

%% 1D
NPixel = size(N_SPAD_all, 1);

v1 = VideoWriter('20230519_10um_bead_video_workdir/plots/10um_bead_video_day6_3_1D_ymaxed.avi');
v1.FrameRate = 400;
open(v1);

for ff = 1: size(N_SPAD_all,2)
    figure(4);
    plot(1:NPixel, N_SPAD_all(:, ff),'LineWidth',2); xlim([1,NPixel]); ylim([0, 400]);
    
    frame = getframe(gcf);
    writeVideo(v1, frame);
    
end

close(v1);