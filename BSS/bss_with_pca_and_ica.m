function [Lhat, score, Zica, Zica2, W, Wunnoised] = bss_with_pca_and_ica(temporal_data, threshold)
%   Lhat: guess for number of fluorophores.
%   score: columns are principal components, rows are score corresponding to observations.
%
%   measuredData: high dimensional data. (pixel aspect ratio has been changed from [2x128] to [31,8] or [16x16])
%   threshold: necessary variance to stop for looking for more components. (listed in explained)
%

    %% First method to determine number of fluorophores, Lhat. Using only one frame.
    physImg = raw2SepRows_quad2020_in_vivo(temporal_data(:,end)); % Data Format: [2x128] This is the last measurement in temporal_data. Just to show the signal in 2D plane.
    
    % to increase dimensions in our data, we convert the array to a different aspect ratio
    %newAspectRatio = [32, 8];
    newAspectRatio = [16,16];

    higherDimension = zeros(newAspectRatio);
    reshaped = reshape(physImg,size(physImg,1),newAspectRatio(2),[]);
    for i=1:size(reshaped,3)
        higherDimension(2*i-1:2*i,:) = reshaped(:,:,i);
    end

    %[U,S,V] = svd(data_in);
    [coeff,score,latent,tsquared,explained,mu] = pca(higherDimension); % default is SVD based.
    
    % explained has ratios of accumulated eigenvalues to total eigenvalues.
    % representing how much variance is shown with each principal
    % component. Either >95% or >99% is our limit.
    summed = 0;
    for k=1:length(explained)
        summed = summed+explained(k);
        if summed>=threshold %
            fprintf('Converged to %f %%. We have %i number of neurons found.\n', summed, k);
            Lhat = k;
            break
        end
    end
    
    figure(95); clf; scatter(score(:,1),score(:,2)); axis equal; xlabel('1st Principal Component'); ylabel('2nd Principal Component');
    
    %% Second method to determine number of fluorophores, Lhat. Using full temporal_data.
%     [coeff,score,latent,tsquared,explained,mu] = pca(temporal_data);
%     
%     summed = 0;
%     for k=1:length(explained)
%         summed = summed+explained(k);
%         if summed>=threshold %
%             fprintf('Converged to %f %%. We have %i number of neurons found.\n', summed, k);
%             Lhat = k;
%             break
%         end
%     end
%     
%     figure(52); clf; scatter(score(:,1),score(:,2)); axis equal; xlabel('1st Principal Component'); ylabel('2nd Principal Component');
%     
    %% to perform ICA based on the number of sources guessed above.
    
    %model = rica(measuredData, Lhat);
    
    [Zica, W, T, mu] = fastICA(temporal_data,Lhat); % W:demixing matrix, Zica:approximation of underlying sources.
    Wunnoised = W.*(W>abs(0.05));
    Zica2 = Wunnoised * temporal_data;
    
end