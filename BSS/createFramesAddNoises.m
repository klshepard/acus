function [temporal_data] = createFramesAddNoises(SIG, numOfSamples)
    
    temporal_data = zeros(length(SIG), numOfSamples);
    for s=1:numOfSamples
        %add uniform background, time-gate and apply poisson noise
        BG = 1e2*imnoise(rand(size(SIG))','poisson')';

        %calculate DCR
        integtime = 1; %lets just say 1 second integration time
        DCR = 125; %cps 
        darknoise = 1e7*(ones(size(SIG)) - imnoise(ones(size(SIG))*integtime*DCR,'poisson'));

        %calculate MIMASP signal

        SIG_shot_noise = 1e4*imnoise(SIG./sum(sum(SIG)),'poisson'); % adding green signal photon shot noise.
        SIG_plus_shot_noise = SIG + SIG_shot_noise;
        SIG_blurred = SIG_plus_shot_noise + BG + darknoise; % [256x1] signal with noise components.
        
        temporal_data(:,s) = SIG_blurred;
        
    end

end