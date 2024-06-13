% to plot detection fields
pixel_chosen = 129;
figure;
%extracted_pixel_data = detectionField(pixel_chosen,:);
%extracted_pixel_data = detectionFieldColumnNormalized(pixel_chosen,:);
extracted_pixel_data = detectionFieldAllNormalized(pixel_chosen,:);

threeD_field = reshape(extracted_pixel_data, sizes);

% 19 ==> 200 um
imagesc('XData',X,'YData',Y,'CData',threeD_field(:,:,19)');
colorbar