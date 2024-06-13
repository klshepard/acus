function [row1,row2] = probePlot2Rows(dat2)

figure(798);
dat2 = rollingWindowHotPixRemoval(dat2,5,3);
subplot(2,1,1); plot(1:256,dat2)
row1 = []; row2 = [];
for pp = 1:256
    if mod(pp,4)==1 || mod(pp,4)==2
        row1 = [row1 dat2(pp)];
    else
        row2 = [row2 dat2(pp)];
    end
end
subplot(2,1,2);
plot(1:128,row1,1:128,row2)