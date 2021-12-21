clearvars;
baseFolder = 'D:\SiliconeHalfAndHalfSiliconeSensor2'

% Load static
load([baseFolder,filesep,'Static2.mat']);
baseFolder = 'D:\SiliconeHalfAndHalfSiliconeSensor2'

% Load dynamic
centerVelocity = load([baseFolder,filesep,'CenterVelocity2.mat'],'VELOCITIES');
centerVelocity = centerVelocity.VELOCITIES;
% leftVelocity = load([baseFolder,filesep,'LeftVelocity.mat'],'VELOCITIES');
% leftVelocity = leftVelocity.VELOCITIES;
% rightVelocity = load([baseFolder,filesep,'RightVelocity.mat'],'VELOCITIES');
% rightVelocity = rightVelocity.VELOCITIES;
axisData = load([baseFolder,filesep,'CenterVelocity2.mat'],'xaxis','zaxis');
% axisData = load([baseFolder,filesep,'CenterVelocity2.mat'],'xaxis','zaxis');
xaxis = axisData.xaxis;
% zaxis = axisData.zaxis;
% zaxis = axisData.zaxis_desamp;

endVelocity = centerVelocity;

figure; imagesc(endVelocity(:,:,1));
[~,y] = ginput(1);
y = round(y);
% Set everything below cutoff equal to NaN;
endVelocity(y:end,:,:) = NaN;

frontBar = find(squeeze(endVelocity(1,:,1)) ~= Inf,100,'first')
frontBar = frontBar(1)+1;
endBar = find(squeeze(endVelocity(1,:,1)) ~= Inf,100,'last')
endBar = endBar(end)+1;

frontBar = find(isnan(squeeze(endVelocity(1,1:round(size(endVelocity,2)/2),1))))
FBD = diff(frontBar);
frontBar = find(FBD ~= 1)
frontBar = frontBar(1);
if isempty(frontBar)
    frontBar = find(isnan(squeeze(endVelocity(1,1:round(size(endVelocity,2)/2),1))))
frontBar = frontBar(end);
end

endBar = find(isnan(squeeze(endVelocity(1,round(size(endVelocity,2)/2):end,1))))
EBD = diff(endBar)

endBarLoc = find(EBD ~= 1) + 1;
endBar = endBar(endBarLoc) + round(size(endVelocity,2)/2) - 2;
if isempty(endBar)
    endBar = find(isnan(squeeze(endVelocity(1,round(size(endVelocity,2)/2):end,1))))
    endBar = endBar(1) + round(size(endVelocity,2)/2) -1;
end
figure;
imagesc(endVelocity(:,:,1));
hold on;
vline(frontBar,'g');
vline(endBar,'g');
hold off;

% Load strain
load([baseFolder,filesep,'STRAIN2.mat']);

BScan = mat2gray(abs(IQData(:,:,end)));
BScan = BScan(1:y,:)
figure; imagesc(BScan);
colormap(gray)
caxis([0 0.15])
rho = 1000;

for k = 1:size(dis_cumsum,3)
    displacement_smoothed(:,:,k) = modefilt(squeeze(dis_cumsum(:,:,k)),[33 33]);
    progressbar(k/size(dis_cumsum,3));
end

% % TEST
% TESTVEL = endVelocity;
% TESTVEL(TESTVEL>10) = NaN;
% TESTVEL(TESTVEL ==Inf) = NaN;
% TESTVEL = mean(TESTVEL,3,'omitnan')
% figure; imagesc(TESTVEL); colormap(jet); caxis([0 10])

VELOCITIES = centerVelocity(1:size(centerVelocity,1),:,end);
BSS = abs(b_scan_strain)';
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        speedAtPixel = VELOCITIES(columnIndex,rowIndex);
        strainAtPixel = BSS(columnIndex,rowIndex);
        sigma(columnIndex,rowIndex) = 3*1000*(speedAtPixel^2)*strainAtPixel;
    end
end

% % Try stress at 330
% for depthValue = 10:10:300
% stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
% figure; subplot(2,1,1)
% imagesc(BScan); hold on; hline(depthValue); hold off; 
% subplot(2,1,2);
% plot(stressMM); title(num2str(depthValue))
% end

figure;
for depthValue = 220:size(b_scan_strain)
    subplot(2,1,1)
    title(num2str(depthValue))
    imagesc(sigma);
    hold on; 
    hline(depthValue);
    hold off; 
    subplot(2,1,2)
    plot(movmean(abs(sigma(depthValue,:)),50,'omitnan'));
    pause
end

depthValue = 220
% depthValue = 100
stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
%         stressAtPixel = sigma(columnIndex,rowIndex);
%         stressAtPixel = axialStress(:,rowIndex);
% stressAtPixel = mean(sigma(interface,rowIndex));
stressAtPixel = stressMM(:,rowIndex);
        strainAtPixel = BSS(columnIndex,rowIndex);
        YM(columnIndex,rowIndex) = (abs(stressAtPixel)/abs(strainAtPixel))/1000; % in kPa
    end
end

close all force;
% figure;
% imagesc(sigma); 
% YM = filloutliers(YM,'makima');
% figure; imagesc(xaxis,zaxis(1:y),abs(YM));
figure; imagesc(xaxis,zaxis(1:size(BScan,1)),abs(YM(1:size(BScan,1),:)));

colormap(jet)
caxis([0 100])
title("Young's modulus")
xlabel('Distance (mm)')
ylabel('Depth (mm)');

figure; imagesc(BSS);
figure; imagesc(stressMM);

figure;
subplot(2,1,1)
imagesc(YM); colormap(jet)
caxis([0 100])
subplot(2,1,2)
imagesc(BScan(1:y,:));
set(gca,'colormap',gray)
caxis([0 0.15])
% topPhantom = roipoly;
leftPhantom = roipoly;
rightPhantom = roipoly;

YMTopPhantom = mean(rmoutliers(YM(topPhantom)),'omitnan')
YMLeftPhantom = mean(rmoutliers(YM(leftPhantom)),'omitnan')
YMRightPhantom = mean(rmoutliers(YM(rightPhantom)),'omitnan')
