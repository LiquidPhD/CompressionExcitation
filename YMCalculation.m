%% Solving for stress and YM Calculation
clearvars;
baseFolder = 'F:\ThreeSiliconePhantomAndSensor'

cutoffFrame = 50;
% Load static
load([baseFolder,filesep,'Static2.mat']);
% Load dynamic
centerVelocity = load([baseFolder,filesep,'CenterVelocity.mat'],'VELOCITIES');
centerVelocity = centerVelocity.VELOCITIES;
leftVelocity = load([baseFolder,filesep,'LeftVelocity.mat'],'VELOCITIES');
leftVelocity = leftVelocity.VELOCITIES;
rightVelocity = load([baseFolder,filesep,'RightVelocity.mat'],'VELOCITIES');
rightVelocity = rightVelocity.VELOCITIES;
axisData = load([baseFolder,filesep,'CenterVelocity.mat'],'xaxis','zaxis');
xaxis = axisData.xaxis;
zaxis = axisData.zaxis;
% Get last images corresponding to most compressed state
endVelocity(:,:,1) = centerVelocity(:,:,end);
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
% 
% TESTVEL = endVelocity(:,:,1);
% TEST = imresize(TESTVEL(:,frontBar:endBar),[size(clippedBScan,1) size(clippedBScan,2)]);
% figure; imagesc(TEST)

for k = 1:3
    endVelocity(:,:,k) = imresize(endVelocity(:,frontBar:endBar,k), [size(clippedBScan,1) size(clippedBScan,2)]);
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
load([baseFolder,filesep,'STRAIN.mat']);

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
% 3ρc^2 = σ/ε where rho is mass density, c is wave speed, sigma is stress,
% and epsilon is strain.

% Note: STRAIN is millistrain, so must be multiplied by 1k.

% b_scan_strain = b_scan_strain'*1000;
VELOCITIES = finalImage;

for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        speedAtPixel = VELOCITIES(columnIndex,rowIndex,end);
        strainAtPixel = b_scan_strain(columnIndex,rowIndex);
        sigma(columnIndex,rowIndex) = 3*1000*(speedAtPixel^2)*strainAtPixel;
    end
end


for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        stressAtPixel = sigma(columnIndex,rowIndex);
        strainAtPixel = b_scan_strain(columnIndex,rowIndex);
        YM(columnIndex,rowIndex) = (stressAtPixel/strainAtPixel)/1000; % in kPa
    end
end

close all force;
% figure;
% imagesc(sigma); 
figure; imagesc(xaxis,zaxis(1:y),YM);
caxis([0 50])

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


 I = uint8(255*mat2gray(BScan(1:y,:), [0 0.15]));
    bg = ind2rgb(I,gray(255));
    bgImg = double(bg);
    alphaFactor = 0.5;
    bgImgAlpha = (1 - alphaFactor) .* bgImg;
    
    II = uint8(255*mat2gray(YM, [0 50]));
    II(II==255) = 0;
    im = ind2rgb(II,fireice(255));

    fgImg = double(im);
    fgImgAlpha = alphaFactor .* fgImg;
    
    fusedImg = fgImgAlpha + bgImgAlpha;

    gcf=figure(1);
    imagesc(xaxis,zaxis(1:y),fusedImg)
