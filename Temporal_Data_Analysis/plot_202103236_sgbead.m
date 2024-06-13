figure; hold on
plot(N_SPAD_10umbead_cont/max(N_SPAD_10umbead_cont))
plot(N_SPAD_10umbead_80mhz/max(N_SPAD_10umbead_80mhz))


%z sweep
chosenpx = 257:512;
figure; hold on
for zz=1:length(dzrange)/2
    lgd = [num2str(dzrange(zz)) char(181) 'm'];
    plot(chosenpx-256,N_SPAD_zsweep(zz,chosenpx),'DisplayName',lgd,'Linewidth',1)
end
legend


%xyz sweep - 1
chosenpx = 257:512;
figure; hold on
xx = find(dxrange==0);
yy = find(dyrange==0);
for zz=1:length(dzrange)
    lgd = [num2str(dzrange(zz)) char(181) 'm'];
    plot(chosenpx-256,squeeze(N_SPAD_xyzsweep(xx,yy,zz,chosenpx)),'DisplayName',lgd,'Linewidth',1)
end
legend

%xyz sweep - 2
chosenpx = 257:512;
figure; hold on
yy = find(dyrange==0);
zz = find(dzrange==20);
for xx=1:4:length(dxrange)
    lgd = ['x' num2str(dxrange(xx)) char(181) 'm '];
    plot(chosenpx-256,squeeze(N_SPAD_xyzsweep(xx,yy,zz,chosenpx)),'DisplayName',lgd,'Linewidth',1)
end
legend
title('X sweep, Y=0, Z=20um')


%xyz sweep - 3
chosenpx = 257:512;
figure; hold on
zz = find(dzrange==20);
for xx=1:length(dxrange)
for yy=1:length(dyrange)
    lgd = ['x' num2str(dxrange(xx)) char(181) 'm ' 'y' num2str(dyrange(yy)) char(181) 'm'];
    plot(chosenpx-256,squeeze(N_SPAD_xyzsweep(xx,yy,zz,chosenpx)),'DisplayName',lgd,'Linewidth',1)
end
end
legend

%xyz sweep - 4. exp03
chosenpx = 257:512;
figure; hold on
for yy=1:2:length(dyrange)
    lgd = [num2str(dyrange(yy)) char(181) 'm'];
    plot(chosenpx-256,squeeze(N_SPAD_xyzsweep(1,yy,1,chosenpx)),'DisplayName',lgd,'Linewidth',1)
end
legend

%xyz sweep - 5. exp04, exp05
chosenpx = 257:512;
figure; hold on
zz = find(dzrange==20);
for xx=1:2:length(dxrange)
    lgd = [num2str(dxrange(xx)) char(181) 'm'];
    plot(chosenpx-256,squeeze(N_SPAD_xyzsweep(xx,1,zz,chosenpx)),'DisplayName',lgd,'Linewidth',1)
end
legend



%0314 y sweep
y000 = find(dyrange==0);
N_SPAD_xyzsweep = N_SPAD_xyzsweep(:,:,:,257:512);
maxaddr = find(N_SPAD_xyzsweep(1,y000,1,:)==max(N_SPAD_xyzsweep(1,y000,1,:)));

% dat = squeeze(N_SPAD_xyzsweep(1,y000,1,:));
% plot(dat)

dat = squeeze(N_SPAD_xyzsweep(1,:,1,maxaddr));
plot(dyrange,dat)

ytheta = atan(dyrange/100)/pi*180;
dist = sqrt(dyrange.^2 + 100^2);
plot(ytheta,dat.*dist.^2)

%0314 y sweep v2
dat = squeeze(N_SPAD_xyzsweep(1,y000,1,:));
plot(1:256,dat)


dat = squeeze(N_SPAD_xyzsweep(1,:,1,:));
dist = sqrt(dyrange.^2 + 100^2);
scatter(dist,dat)
set(gca,'XScale','log','YScale','log')

plot(dyrange,dat(:,119:122))





%0314 x sweep
x000 = find(dxrange==0);
N_SPAD_xyzsweep = N_SPAD_xyzsweep(:,:,:,257:512);
maxaddr = find(N_SPAD_xyzsweep(x000,1,1,:)==max(N_SPAD_xyzsweep(x000,1,1,:)));

dat = squeeze(N_SPAD_xyzsweep(:,1,1,maxaddr))';
plot(dxrange,dat)

dat = squeeze(N_SPAD_xyzsweep(:,1,1,:))';
plot(dxrange,dat)

dat = squeeze(N_SPAD_xyzsweep(:,1,1,:))';
figure; hold on
for pp = 1:256
    if pp<121
        plot(dxrange,dat(pp,:),'r')
    else
        plot(dxrange,dat(pp,:),'b')
    end
end

xtheta = atan(dxrange/100)/pi*180;
dist = sqrt(dxrange.^2 + 100^2);
plot(xtheta,dat.*dist.^2)




%% 0323 y sweep: fiber illuminating from distal end of shank, moved single 10um bead in 25um increments
dat = N_SPAD_xyzsweep(:,:,:,257:512);
dat = dat(:,:,:,180:end);
dat = squeeze(dat);

figure; plot(dat')
xlabel('pixel'); ylabel('photon count')

figure; semilogy(dat')

figure; plot(dyrange,sum(dat,2))
xlabel('x(um)'); ylabel('sum of photon count')


%% 0323 xyz sweep: fiber illuminating from distal end of shank, moved single 10um bead in 25um increments
dat = N_SPAD_xyzsweep(:,:,:,257:512);
%dat = dat(:,:,:,180:end);
dat = squeeze(dat);

% xy distribution at zz=1
zz = 1;
datplot = sum(squeeze(dat(:,:,zz,:)),3);
figure; imagesc(dxrange,dyrange,datplot')
axis equal tight

figure;
clims = sum(dat,4);
clims = [min(clims(:)) max(clims(:))];
for zz = 1:length(dzrange)
    subplot(1,length(dzrange),zz)
    datplot = sum(squeeze(dat(:,:,zz,:)),3);
    imagesc(dxrange,dyrange,datplot')
    axis equal tight
    caxis(clims)
end

% xz distribution at zz=1
yy = 1;
datplot = sum(squeeze(dat(:,yy,:,:)),3);
figure; imagesc(dxrange,dzrange,datplot')
axis equal tight

figure;
clims = sum(dat,4);
clims = [min(clims(:)) max(clims(:))];
nskip = 2;
for yy = 1:nskip:length(dyrange)/nskip
    subplot(1,ceil(length(dyrange)/nskip),yy)
    datplot = sum(squeeze(dat(:,yy,:,:)),3);
    imagesc(dxrange,dzrange,datplot')
    axis equal tight
    caxis(clims)
end

%overlay all raw data
figure; hold on
for xx = 1:length(dxrange)
    for yy = 1:10:length(dyrange)
        plot(squeeze(dat(xx,yy,zz,:)))
    end
end

%% resolution experiment
dat = N_SPAD_xyzsweep(:,:,:,257:512);
% to get rid of one hot pixel and plot the rest of the pixels
dat(:,:,:,35) = [];
% to take first 1 mm of the shank only. The rest of the chip is scratched.
dat= dat(:,:,:,1:150); 
dat = squeeze(dat);

%overlay y sweep in the middle
figure; hold on
legendStrings = strings();
first = 1; 
for yy = 1:6:length(dyrange)
        plot(dat(yy,:))
        
        if first == 1
            legendStrings(1) = strcat('y\_sweep\_', num2str(yy));
            first = 0; 
        else
            legendStrings(end+1) = strcat('y\_sweep\_', num2str(yy));
        end
end
legend(legendStrings)

%% plotting rows separately

dat = N_SPAD_xyzsweep(:,:,:,257:512);
dat = squeeze(dat);
dat2 = dat(17,:);
[r1, r2] = probePlot2Rows(dat2);


for ii = 1:size(dat,1)
    dat3(ii) = dat(ii,pp);
%     [r1, r2] = probePlot2Rows(dat2);
end
figure; plot(dat3)

figure; plot(dat(:,40:2:50))


