function dat2 = rollingWindowHotPixRemoval(dat1,nwin,nthresh)
%Remove hot pixel with rolling window.
% nwin = rolling window width, should be odd
% nthresh = threshold, ex) 5*median

nsgside = floor(nwin/2);
dat2 = dat1;
for pp = 1:length(dat1)
    win1 = max(1,pp-nsgside);
    win2 = min(length(dat1),pp+nsgside);
    
    tempMed = median(dat1(win1:win2));
    if dat1(pp) > tempMed * nthresh
        dat2(pp) = tempMed; %replace with mean
    end

end