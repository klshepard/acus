function [ A ] = measured_volume(X,Y,Z, angleData ,NPixel,DProbe)
%Create a model for angle-specific measurements from point sources in a volume.
%   angleData = [angles, angleSenstivity];

% for every pixel, calculate the power measured from each cubic unit in X,Y,Z volume.

    for nn = 1:NPixel
        for zz=1:length(Z)        

            %center the X-Y to the pixel. Distances to the pixel from all X,Y combinations.
            X0 = (X-DProbe(nn,1));
            Y0 = (Y-DProbe(nn,2));

            % meshgrids
            [Xm,Ym] = meshgrid(X0,Y0);

            %translate to spherical coordinates        
            [~,elM,rM] = cart2sph(Xm(:),Ym(:),Z(zz));
            %azM = reshape(azM,length(Y),length(X));
            elM = reshape(pi/2 - elM,length(Y),length(X)); % radial
            elM_degree = elM .* 180/pi;
            rM  = reshape(rM,length(Y),length(X));

            % we have 2 windowing functions, 1 parallel and 1 perpendicular
            %windowpara = 1+mpara*cos(gammapara*Tvr);
            % perpendicular window 
            %windowperp = 1+0.5*mperp*cos(gammaperp*Pvr+pi);
            % a global lambertian windowing function
            %windowing = cos(elM).^2;
            % And a hard cutoff deadangle window      
            %deadangle_gain = ones(length(Y),length(X));                
            %deadangle_gain(find(abs(elM)>deadangle)) = 0;

            % A scalar for distance (ASPAD = 171um^2)
            scaler = 7.8241e+07 ./ (rM.^2); % determined after a calibration. But we'll still use normalization for both empirical and simulation data to prevent FPS mismatch.
    %        scaler = 1; 
    %        scaler = 1./(4*pi*sqrt(rM)); 
    %        scaler = 1./(4*pi*rM.^2);
    %        scaler = 1./(4*pi*rM);
    %        scaler = 1./(4*pi*rM.^1.5);
    %        scaler = sqrt(171./(4*pi*rM.^2));

            % bare pixel. empirical angle sensitivity effect.
            anglesDetected = zeros(length(elM_degree(:)),1);
            pixAngleResponse = zeros(length(elM_degree(:)),1);
            for e=1:length(elM_degree(:))
                [~, ind] = min(abs(elM_degree(e)-angleData(:,1))); % instead of using 'findnearest'.
                %ind = findnearest(elM_degree(e), angleData(:,1));
                anglesDetected(e) = angleData(ind,1);
                pixAngleResponse(e) = angleData(ind,2);
            end
            
            % angular response of this pixel to each cubic unit in the volume
            pixAngleResponse = reshape(pixAngleResponse,length(Y),length(X));

            % multiply these two effects.
            output(:,:,zz) = (pixAngleResponse .* scaler)';
            %figure(17); imagesc(output(:,:,zz)); colorbar; title(sprintf('Pixel # %i; Z = %i um' , nn, Z(zz)));
            if nn == 129
                %figure(17); clf; imagesc(X,Y,output(:,:,zz)'); set(gca,'YDir','normal'); colorbar; title(sprintf('Pixel # %i; Z = %i um' , nn, Z(zz))); hold on; plot(DProbe(nn,1), DProbe(nn,2), 'x', 'color', 'red');
            end
        end
        A(nn,:) = output(:);
    end
    %A(A<0) = 0;
end