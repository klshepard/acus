function [physImg] = raw2SepRows_NASP_only(N_SPAD)
%Transform N_SPAD(256,1) raw spad image to physImg(2,128) physically relevant
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
N_SPAD_hpRemoved = rollingWindowHotPixRemoval(N_SPAD,5,5);

nshanks = ceil(length(N_SPAD_hpRemoved)/256);
physImg = zeros(nshanks*2,128);

for ss = 1:nshanks %shank index
    for pp = 1:256 %each shank has 256 pixels, 2 rows
        iof4 = mod(pp-1,4); %returns 0,1,2,3
        groupof4 = ceil(pp/4);
        if iof4 < 2 %0 and 1; top row of shank
            toprow = 1;
        else %2 and 3
            toprow = 0;
        end
        newR = (ss-1)*2+1+toprow;
        newC = 1+(iof4+(toprow-1)*2)+(groupof4-1)*2;
        physImg(newR,newC) = N_SPAD_hpRemoved((ss-1)*256+pp);
    end
end
physImg = flip(physImg,1); % flipping it so that it complies with our x-y definitions on the image.

figure(72);
subplot(2,1,1); imagesc(physImg(1:2,:)); %set(gca,'Ydir','normal');
subplot(2,1,2); plot(physImg(1:2,:)','LineWidth',1); xlim([1,length(physImg(1,:))]); legend('r1','r2');