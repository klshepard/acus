function [autoCorr, lags, w, diffHeight] = autocorrAndDiff(linearizedData, L)
    % unfinished. Decided to go with smoothing the measured signal directly.
    
    [autoCorr,lags] = xcorr(linearizedData, 'coeff');
    figure(2); plot(lags,autoCorr); hold on;
    
    [pks,locs, w, p] = findpeaks(autoCorr,'SortStr','descend','MinPeakProminence',0.1);
    
    pixLoc = zeros(L,1);
    for m=1:L
        [max1, ind1] = max(pks);
        pks(ind1) = -Inf;
        
        pixLoc(m) = ind1;
    end
end