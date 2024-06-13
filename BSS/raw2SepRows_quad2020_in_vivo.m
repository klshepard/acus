function [physImg] = raw2SepRows_quad2020(N_SPAD)
%Transform N_SPAD(512,1) raw spad image to physImg(4,128) physically relevant
%matrix. Raw N_SPAD counts pixels in zigzag format.
%First pixel is bottom-leftmost, representing the physical image.
%N_SPAD   257 258 261 262 ... 510
%         259 260 263 264 ... 512
%         1   2   5   6   ... 254
%         3   4   7   8   ... 256
%physImg  (2,1) (2,2) (2,3) (2,4) ... (2,128)
%         (1,1) (1,2) (1,3) (1,4) ... (1,128)
%To plot physImg in true physical orientation:
%figure; imagesc(physImg); set(gca,'Ydir','normal')


%hot pix removal by rolling window
%N_SPAD_hpRemoved = rollingWindowHotPixRemoval(N_SPAD,5,5);

nshanks = 1; %ceil(length(N_SPAD_hpRemoved)/256); 
physImg = zeros(nshanks*2,64);

for ss = 1:nshanks %shank index
    for pp = 1:128 %each in vivo shank has 128 pixels, 2 rows
        iof4 = mod(pp-1,4); %returns 0,1,2,3
        groupof4 = ceil(pp/4);
        if iof4 < 2 %0 and 1; top row of shank
            toprow = 1;
        else %2 and 3
            toprow = 0;
        end
        newR = (ss-1)*2+1+toprow;
        newC = 1+(iof4+(toprow-1)*2)+(groupof4-1)*2;
        physImg(newR,newC) = N_SPAD((ss-1)*128+pp);
    end
end

figure(802);
subplot(2,1,1); imagesc(physImg); set(gca,'Ydir','normal'); caxis([0 400]);
subplot(2,1,2); plot(physImg','LineWidth',2); xlim([1,length(physImg(1,:))]); legend('Row 1', 'Row 2'); ylim([0, 400]);