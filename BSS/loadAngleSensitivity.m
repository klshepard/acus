function [anglesExtended, normalized_air_int] = loadAngleSensitivity(directory, deltaAngle, interpolation)
    
    % directory: directory to the .mat file containing -90 to 90 angle sensitivity data
    % deltaAngle: steps of angle for interpolation. > 5
    % interpolation: degree of interpolation. Cubic = 3, Linear = 1.
    
    % to extrapolate angle sensitivity data for every angle from -90 to 90.
    angleData = load(directory);
    anglesExtended = -90:deltaAngle:90;

    if interpolation == 3
        air_interpolated = interp1(angleData.angles,angleData.air,anglesExtended,'pchip');
    elseif interpolation == 1
        air_interpolated = interp1(angleData.angles,angleData.air,anglesExtended,'linear');
    end

    % to make it a symmetric waveform to get rid of human-induced errors in the optical setup
    air_interpolated(int64(length(air_interpolated)/2):end) = flip(air_interpolated(1:int64(length(air_interpolated))/2));
    normalized_air_int = air_interpolated ./ max(air_interpolated);
    
    %plot(angleData.angles, angleData.air); hold on;
    %figure(12); plot(anglesExtended, normalized_air_int); xlim([-90 90]); xticks(-90:20:90); xlabel('Angle'); ylabel('Normalized Pixel Response');

end