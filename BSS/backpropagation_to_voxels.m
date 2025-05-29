fps = 400;
duration = 15; %sec
timepoints = fps*duration;

%% pinv ==> SVD

all_voxels_backpropagated = ones(size(N_SPAD_all_diff,2), size(detectionFieldColumnNormalized,2)).*(NaN);
left_psinv = pinv(detectionFieldColumnNormalized);%, 1e-5);
for i=1:timepoints
    all_voxels_backpropagated(i,:) = left_psinv*N_SPAD_all_diff(:,i);
end

save('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_1_102_400fps_detrended3_filtered_3D_detected_all_voxels_and_backProp.mat','N_SPAD_all_diff', 'XYZlocDetected_single_bead_all', 'all_voxels_RMSE', 'all_voxels_backpropagated', '-v7.3');

%% lsqminnorm ==> QR decomposition + COD [gives almost same results with pinv]
x_not = ones(size(N_SPAD_all_diff,2), size(detectionFieldColumnNormalized,2)).*(NaN);

for i=1:timepoints
    x_not(i,:) = lsqminnorm(detectionFieldColumnNormalized, N_SPAD_all_diff(:,i),1e-2,'warn');
    disp(i);
end
save('/users/syilmaz/MATLAB/20230420_quad2020_in_vivo/20241105_revisions_in_vivo_sparse_GCaMP_workdir/savedata/in_vivo_sst_gcamp6s_animal_1_102_400fps_detrended3_filtered_3D_detected_all_voxels_and_backProp_lsqminnorm_tol1e-2.mat','N_SPAD_all_diff', 'XYZlocDetected_single_bead_all', 'all_voxels_RMSE', 'all_voxels_backpropagated', '-v7.3');

all_voxels_backpropagated = x_not;
%% plot backpropagated voxel vs. time

XYZ_location_center = [70, 30, 150];
spatial_sum = zeros(size(all_voxels_backpropagated, 1), 1);
resolution=10;
figure;
for rb = 0:10:40
    range_border = rb;

    voxel_count=0;
    for i=-range_border:resolution:range_border
        for j=-range_border:resolution:range_border
            for k=-range_border:resolution:range_border
                XYZ_location = [XYZ_location_center(1)+i, XYZ_location_center(2)+j, XYZ_location_center(3)+k];
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
sgtitle(['Backprop. Voxel Coordinates: ' num2str(XYZ_location_center)])
%% 3D heatmap
% % % 
% % % all_voxels_backpropagated_3D = reshape(all_voxels_backpropagated(1193,:),[length(X), length(Y), length(Z)]);
% % % 
% % % figure(28);
% % % for v = 70000:80000%size(all_voxels_backpropagated,2) 
% % %     [x, y, z] = ind2sub(sizes,v);
% % %     
% % %     %scatter3(X(x),Y(y),Z(z), abs(all_voxels_backpropagated(1193,v))*200, 'g', 'filled', 'hexagram'); grid; hold on;
% % %     %axis([X(1) X(end) Y(1) Y(end) 0 Z(end)+1]); view([-10 15]); grid; 
% % % end

%% 2D heatmaps z stack
figure(30);

all_voxels_backpropagated_filtered = detrend(movmean(all_voxels_backpropagated',1),3)';

minimumV = median(all_voxels_backpropagated_filtered, 'all');
maximumV = max(all_voxels_backpropagated_filtered, [], 'all');

for ff = 1193:1193%1:size(all_voxels_backpropagated, 1)
    all_voxels_backpropagated_3D = reshape(all_voxels_backpropagated_filtered(ff,:),[length(X), length(Y), length(Z)]);
    p = 1;
    for zz = 10:length(Z)
        subplot(10,1,p); p=p+1;
        img = all_voxels_backpropagated_3D(:,:,zz)';

        imshow(img, [minimumV, maximumV])
        colormap default;
    end
end

%% 2D heatmaps RMSEs

% % % figure(30);
% % % 
% % % minimumV = min(all_voxels_RMSE,[], 'all', 'omitnan');
% % % maximumV = max(all_voxels_RMSE, [], 'all', 'omitnan');
% % % 
% % % for ff = 50:50%1:size(all_voxels_backpropagated, 1)
% % %     all_voxels_RMSE_3D = reshape(all_voxels_RMSE(ff,:),[length(X), length(Y), length(Z)]);
% % %     p = 1;
% % %     for zz = 10:length(Z)
% % %         subplot(10,1,p); p=p+1;
% % %         img = all_voxels_RMSE_3D(:,:,zz)';
% % % 
% % %         imshow(img, [minimumV, maximumV])
% % %         colormap default; colorbar;
% % %     end
% % % end