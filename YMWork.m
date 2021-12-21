%% Solving for stress and YM Calculation
clearvars;
baseFolder = 'E:\HalfHalfAgarDay2NewSensorAuto'

% Load static
load([baseFolder,filesep,'Static2.mat']);
baseFolder = 'E:\HalfHalfAgarDay2NewSensorAuto'

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

% Get last images corresponding to most compressed state
endVelocity(:,:,1) = centerVelocity(:,:,30);
endVelocity(:,:,2) = leftVelocity(:,:,end);
endVelocity(:,:,3) = rightVelocity(:,:,end);

% figure;
% nexttile
% imagesc(endVelocity(1:y,:,1))
% colormap(jet); caxis([0 10])
% nexttile
% imagesc(endVelocity(1:y,:,2))
% colormap(jet); caxis([0 10])
% nexttile
% imagesc(endVelocity(1:y,:,3))
% colormap(jet); caxis([0 10])


% Select cut-off
figure; imagesc(endVelocity(:,:,1));
[~,y] = ginput(1);
y = round(y);
% Set everything below cutoff equal to NaN;
endVelocity(y:end,:,:) = NaN;
figure; imagesc(endVelocity(:,:,1)); figure; imagesc(endVelocity(:,:,2)); figure; imagesc(endVelocity(:,:,3))

% Find end bars
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
% 
% TESTVEL = endVelocity(:,:,1);
% TEST = imresize(TESTVEL(:,frontBar:endBar),[size(clippedBScan,1) size(clippedBScan,2)]);
% figure; imagesc(TEST)

for k = 1:3
    endVelocity(:,:,k) = imresize(endVelocity(:,frontBar:endBar,k), [size(BScan,1) size(BScan,2)]);
end
figure; 
nexttile
imagesc(endVelocity(1:y,:,1))
 colormap(jet); caxis([0 10])
 nexttile
 imagesc(endVelocity(1:y,:,2));
  colormap(jet); caxis([0 10])
  nexttile
  imagesc(endVelocity(1:y,:,3));
   colormap(jet); caxis([0 10])
   
   endVelocity(endVelocity>10) = NaN;
   figure; 
nexttile
imagesc(endVelocity(1:y,:,1))
 colormap(jet); caxis([0 10])
 nexttile
 imagesc(endVelocity(1:y,:,2));
  colormap(jet); caxis([0 10])
  nexttile
  imagesc(endVelocity(1:y,:,3));
   colormap(jet); caxis([0 10])
   
   
% Load strain
load([baseFolder,filesep,'STRAIN2.mat']);

BScan = mat2gray(abs(IQData(:,:,end)));
BScan = BScan(1:y,:)
figure; imagesc(BScan);
colormap(gray)
caxis([0 0.15])
rho = 1000;

% % TEST
% TESTVEL = endVelocity;
% TESTVEL(TESTVEL>10) = NaN;
% TESTVEL(TESTVEL ==Inf) = NaN;
% TESTVEL = mean(TESTVEL,3,'omitnan')
% figure; imagesc(TESTVEL); colormap(jet); caxis([0 10])


clearvars finalImage;
area = size(endVelocity,2)
thirdArea = round(area/3)
finalImage(:,1:thirdArea) = endVelocity(:,1:thirdArea,1);
finalImage(:,thirdArea+1:2*thirdArea) = endVelocity(:,thirdArea+1:2*thirdArea,2);
finalImage(:,2*thirdArea:size(endVelocity,2)) = endVelocity(:,2*thirdArea:size(endVelocity,2),1);
finalImage = finalImage(1:y,:);
finalImage(finalImage == Inf) = NaN;
finalImage(finalImage == -Inf) = NaN;
figure; imagesc(finalImage); colormap(jet); caxis([0 10])

rho = 1000;

figure;
imagesc(BScan);
colormap(gray)
caxis([0 0.15])
topPhantom = roipoly;
leftPhantom = roipoly;
rightPhantom = roipoly;
interface = roipoly;
% 3ρc^2 = σ/ε where rho is mass density, c is wave speed, sigma is stress,
% and epsilon is strain.

% Note: STRAIN is millistrain, so must be multiplied by 1k.

% b_scan_strain = b_scan_strain'*1000;
VELOCITIES = centerVelocity(1:size(b_scan_strain,2),:,:);
BSS = b_scan_strain';
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        speedAtPixel = VELOCITIES(columnIndex,rowIndex);
        strainAtPixel = BSS(columnIndex,rowIndex);
        sigma(columnIndex,rowIndex) = 3*1000*(speedAtPixel^2)*strainAtPixel;
    end
end

% sigma = fillmissing(sigma,'makima')
bottomCut = lowerBound;
for timePoint = 1:size(IQData,3)
    if timePoint >= size(IQData,3)/2;
        bottomCut = 450;
    end
    for counter = 1:size(IQData,2)
        Aline = abs(IQData(1:bottomCut,counter,timePoint));
%         figure(1);
%         plot(Aline);
        peakLoc = find(Aline == max(Aline));
%         hold on; vline(peakLoc);
%         hold off; 
        SDL(:,counter,timePoint) = peakLoc;
    end
    finalSDL(:,timePoint) = movmean(filloutliers(SDL(:,:,timePoint),'makima'),10);
end

    figure('windowState','maximized')

% look at sigma through depth
for k = 1:size(sigma,1)
    subplot(3,1,1)
    imagesc(BScan);
    colormap(gray);
    hold on;
    hline(k); hold off;
    subplot(3,1,2)
    plot(sigma(k,:))
    subplot(3,1,3)
    plot(mean(sigma(1:k,:)));
    pause(0.1)
end

finalSDL = round(finalSDL(:,50));
% stressAtPixel = sigma.*topPhantom;
% figure; imagesc(stressAtPixel)
% stressAtPixel(stressAtPixel==0) = NaN;
% for k = 1:size(stressAtPixel,2)
%     meanStress(:,k) = mean(stressAtPixel(:,k),'omitnan');
%     stdStress(:,k) = std(stressAtPixel(:,k),'omitnan');
% end
% sigma = fillmissing(sigma,'makima');
% sigma = filloutliers(sigma,'makima');
% sigma = sigma;
figure;
subplot(2,1,1)
imagesc(BScan);
subplot(2,1,2);
imagesc(sigma);
colormap(jet)
% caxis([0 0.15])

% Try stress at interface.
for k = 1:size(sigma,2)
stressAtPixel(:,k) = mean(sigma(finalSDL(k)-10:finalSDL(k)+10,k),'omitnan');
end

figure;
plot(stressAtPixel)
% % Try stress in selected ROI 
% interface = roipoly;
% sigmaInterface = sigma.*interface;
% sigmaInterface(sigmaInterface==0) = NaN;
% for k = 1:size(sigmaInterface,2)
% meanStress(:,k) = mean(abs(sigmaInterface(:,k)),'omitnan');
% stdStress(:,k) = std(abs(sigmaInterface(:,k)),'omitnan');
% end
% 
% meanStress = filloutliers(meanStress,'makima')


for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
%         stressAtPixel = sigma(columnIndex,rowIndex);
%         stressAtPixel = axialStress(:,rowIndex);
% stressAtPixel = mean(sigma(interface,rowIndex));
stressAtPixel = meanStress(:,rowIndex);
        strainAtPixel = b_scan_strain(columnIndex,rowIndex);
        YM(columnIndex,rowIndex) = (abs(stressAtPixel)/abs(strainAtPixel))/1000; % in kPa
    end
end

close all force;
% figure;
% imagesc(sigma); 
YM = filloutliers(YM,'makima');
figure; imagesc(xaxis,zaxis(1:y),YM);
colormap(jet)
caxis([0 100])


% Try stress at 330
for depthValue = 10:10:300
stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
figure(1); subplot(2,1,1)
imagesc(BScan); hold on; hline(depthValue); hold off; 
subplot(2,1,2);
plot(stressMM); title(num2str(depthValue))
pause(0.1)
end

depthValue = 130
stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
figure; plot(stressMM)
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
YM = filloutliers(YM,'makima');
figure; imagesc(xaxis,zaxis(1:y),YM);
colormap(jet)
caxis([0 500])


YMTopPhantom = mean(rmoutliers(YM(topPhantom)),'omitnan')
YMLeftPhantom = mean(rmoutliers(YM(leftPhantom)),'omitnan')
YMRightPhantom = mean(rmoutliers(YM(rightPhantom)),'omitnan')

TOP = YM(7,59)
LEFT = YM(291,71)
RIGHT = YM(304,498)


%% Try this

VELOCITIES = finalImage;

speedLeft = mean(rmoutliers(VELOCITIES(leftPhantom)),'omitnan')
speedRight = mean(rmoutliers(VELOCITIES(rightPhantom)),'omitnan')

% Left side
sigma = zeros([size(BScan,1) size(BScan,2)]);
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        if leftPhantom(columnIndex,rowIndex)
%         speedAtPixel = VELOCITIES(columnIndex,rowIndex,end);
        strainAtPixel = b_scan_strain(columnIndex,rowIndex);
        sigma(columnIndex,rowIndex) = 3*1000*(speedLeft^2)*strainAtPixel;
        end
    end
end

% Right side
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        if rightPhantom(columnIndex,rowIndex)
%         speedAtPixel = VELOCITIES(columnIndex,rowIndex,end);
        strainAtPixel = b_scan_strain(columnIndex,rowIndex);
        sigma(columnIndex,rowIndex) = 3*1000*(speedRight^2)*strainAtPixel;
        end
    end
end

YMTracker = [];
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
%         if leftPhantom(columnIndex,rowIndex)
        stressAtPixel = sigma(columnIndex,rowIndex);
        strainAtPixel = b_scan_strain(columnIndex,rowIndex);
%         stressTracker(columnIndex,rowIndex) = stressAtPixel;
%         strainTracker(columnIndex,rowIndex) = strainAtPixel; 
        YM(columnIndex,rowIndex) = (stressAtPixel/strainAtPixel)/1000; % in kPa
%         YMTracker = [ YMTracker (stressAtPixel/strainAtPixel)/1000];
%         figure(1);
%         nexttile(1)
%         imagesc(stressTracker);
%         title(stressAtPixel)
%         nexttile(2)
%         imagesc(strainTracker);
%         title(strainAtPixel)
%         nexttile(3)
%         plot(YMTracker)
%         pause(0.2)
       

%         end
    end
end

close all force;
% figure;
% imagesc(sigma); 
figure; imagesc(xaxis,zaxis(1:y),YM);
caxis([0 50])

YMcut = YM(1:y,:);
 I = uint8(255*mat2gray(BScan(1:y,:), [0 0.15]));
    bg = ind2rgb(I,gray(255));
    bgImg = double(bg);
    alphaFactor = 0.5;
    bgImgAlpha = (1 - alphaFactor) .* bgImg;
    
    II = uint8(255*mat2gray(YMcut, [0 50]));
    II(II==255) = 0;
    im = ind2rgb(II,fireice(255));

    fgImg = double(im);
    fgImgAlpha = alphaFactor .* fgImg;
    
    fusedImg = fgImgAlpha + bgImgAlpha;

    gcf=figure(1);
    imagesc(xaxis,zaxis(1:y),fusedImg)
