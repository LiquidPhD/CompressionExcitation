%% Work 11/23/21

clearvars;
close all force;
baseFolder = 'E:\GelatinHalfHalf112821';
%
% [lowerBound] = selectLowerBound(folders)
startDepth = 1;
lowerBound = 1395

%% Process static data

% Load static data
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
if ~isempty(rdir([baseFolder,filesep,'combinedIQData.mat']))
    load([baseFolder,filesep,'combinedIQData.mat'])
    folderIndex = 1;
    
else
    [IQData,Parameters] = loadStaticData(folders,lowerBound,length(folders));
folderIndex = 1;
save([baseFolder,filesep,'combinedIQData.mat'],'IQData','Parameters');
end
[vec_phase_diff] = VPD(IQData);
%  [~,~,~,Parameters] = loadDynamicData(folders,folderIndex,lowerBound);

for VPDIndex = 1:size(vec_phase_diff,3)
imageArray = double(abs(vec_phase_diff(:,:,VPDIndex)));
[rows columns numberOfColorBands] = size(imageArray);

k = 1; % Order of the polynomial
windowSize = 41;
verticallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 1);
subplot(2, 2, 2);
imshow(verticallySmoothedImage, [0 255]);
title('Savitzky-Golay filtered in the vertical direction only');

% Apply the Savitzky-Golay filter.
% First apply it in the vertical (row) direction.
k = 1; % Order of the polynomial
windowSize = 21;
horizontallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 2);
subplot(2, 2, 3);
imshow(horizontallySmoothedImage, [0 255]);
title('Savitzky-Golay filtered in the horizontal direction only');

doublySmoothedImage = sgolayfilt(verticallySmoothedImage, k, windowSize, [], 2);
subplot(2, 2, 4);
imshow(doublySmoothedImage, [0 255]);
title('Savitzky-Golay filtered in both directions');

VPD_smoothed(:,:,VPDIndex) = doublySmoothedImage;
end
 
 for k = 1:size(vec_phase_diff,3)
     VPD_unwrapped(:,:,k) = phase_unwrap(VPD_smoothed(:,:,k));
 end
 
 playWaveVideo(VPD_unwrapped)
[particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters);
% particleVelocity = permute(vec_phase_diff,[2 1 3]);
[Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
Nz = 1395;
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);
zaxis = zaxis(startDepth:end);



close all force;
for k = 2:20 %size(TEST,3)
    TEST = VPD_smoothed(:,:,k);
figure; imagesc(TEST)
title(num2str(k))
colormap(jet)
end
autoArrangeFigures


%% Normalize phase diff to surface

figure; himage = imagesc(mat2gray(abs(IQData(:,:,20)))); caxis([0 0.1])
colormap(gray);
[~,y] = ginput(1);

TESTABS = abs(TEST);

figure; plot(mean(TESTABS(round(y),:),1))
TESTFITTING = smooth(mean(TESTABS(round(y),:),1),33);
figure; plot(TESTFITTING)

XFITTING = 1:length(TESTFITTING);

TEST3 = phase_unwrap(TESTABS - TESTFITTING');
% TEST3 = phase_unwrap(TESTABS - TESTLINE);
TEST3 = TEST3(y:end,:);
newMatrix = [];
for lineIterator = 1:size(TESTABS,1)
    if lineIterator < round(y)
        newMatrix(lineIterator,:) = NaN*(1:size(TESTABS,2));
    else
%     newMatrix(lineIterator,:) = phase_unwrap(TESTABS(lineIterator,:)-TESTFITTING');
        newMatrix(lineIterator,:) = TESTABS(lineIterator,:)-TESTFITTING';

    end
end

figure; imagesc(newMatrix); colormap(jet)
newMatrix = newMatrix(round(y):end,:);
figure; imagesc(newMatrix); colormap(jet)
TESTNM = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    newMatrix,100);

% startPoint = 40; 
mouseoverShowDepth(newMatrix)
uiwait;
clearvars surf_dis corr_matrix
for pos=1:size(newMatrix,2)
    surf_dis(pos)=median(TESTABS(startPoint:startPoint+20,pos));
    corr_matrix(:,pos)=TESTABS(:,pos)-surf_dis(pos);
end
    corr_matrix(1:startPoint-1,:)=[];
    
strain=strain_calculation2(corr_matrix,zaxis(2),1,100);
figure; imagesc(abs(strain)/1000); colormap(hot); colorbar; caxis([0,10]);
    
% figure; plot(newMatrix(:,100)); hold on; plot(newMatrix(:,450)); hold off; 
% figure; plot(TEST(:,50)); hold on; plot(TEST(:,450)); hold off; 
% 
% figure; imagesc(TEST2); figure; imagesc(TEST3)

% TEST3(1:selectedSurface,:) = [];
deltaX = diff(xaxis);
deltaX = deltaX(1);
TEST3 = TEST3.*deltaX;
  TEST3STRAIN = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    TEST3,100);  
TESTNM = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    newMatrix,100);
TESTSTRAIN = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    abs(TEST),100); 


Instron = [68.14077 69.76282 63.35117];
Instron = mean(Instron);
stiffInstron = [100.71113 111.63580 120.22536 123.46917 104.89593 108.09592];
stiffInstron = mean(stiffInstron);
softInstron = [84.01424 84.65460 83.87546];
softInstron = mean(softInstron);
figure; himage = imagesc(mat2gray(abs(IQData(171:end,:,20)))); caxis([0 0.1])
colormap(gray);
  
h = impoly(gca);
PosTime = wait(h);

BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));
BW2(BW2>0) = 1;

figure; imagesc(BW2)

TESTLINE = TEST3STRAIN';
TESTLINE(BW2~=1) = NaN;
figure; imagesc(TESTLINE)

TESTLINE = mean(TESTLINE,1,'omitnan');
figure; plot(TESTLINE)

% E = epsilon / sigma
% YM = stress / strain

stressCalculated = Instron .* TESTLINE;
TEST3STRAIN = TEST3STRAIN';
stressCalculated2(1:200) = median(stressCalculated(1:200));
stressCalculated2(201:length(stressCalculated)) = median(stressCalculated(201:length(stressCalculated)));
figure; plot(stressCalculated2)
for k = 1:size( TEST3STRAIN,1)
YMCompression(k,:) = (stressCalculated2 ./ TEST3STRAIN(k,:)');
end
figure; imagesc(abs(YMCompression)./1000)
colormap(jet)
CAX = [quantile(abs(YMCompression(:))./1000,0.25) quantile(abs(YMCompression(:))./1000,0.75)]
caxis(CAX)

disp('Select bad region')
badRegion = roipoly;

leftBlock = roipoly;
rightBlock = roipoly;
YMCompression(badRegion==1)=NaN;
leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlock~=1) = NaN;
rightYM(rightBlock~=1)=NaN;
figure; imagesc(leftYM);
figure; imagesc(rightYM)

YMLeftPhantom = mean(rmoutliers(leftYM),'omitnan')
YMRightPhantom = mean(rmoutliers(rightYM),'omitnan')

YMLeftPhantom = mean(rmoutliers(abs(YMCompression(leftBlock))),'omitnan')
YMRightPhantom = mean(rmoutliers(abs(YMCompression(rightBlock))),'omitnan')

YMLNegRem = YMCompression;
YMLNegRem(YMLNegRem<0) = NaN;
YMRNegRem = YMLNegRem;

YMLeftPhantomNegRem = mean(rmoutliers(abs(YMLNegRem(leftBlock))),'omitnan')
YMRightPhantomNegRem = mean(rmoutliers(abs(YMRNegRem(rightBlock))),'omitnan')


% figure; imagesc(YMCompression)
figure; imagesc(mat2gray(abs(IQData(:,:,20)))); colormap(gray); caxis([0 0.1])
% caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])



leftBlockTop = roipoly;
rightBlockTop = roipoly;

leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlockTop~=1) = NaN;
rightYM(rightBlockTop~=1)=NaN;
figure; imagesc(leftYM);
figure; imagesc(rightYM)

YMLeftPhantom = mean(rmoutliers(leftYM),'omitnan')
YMRightPhantom = mean(rmoutliers(rightYM),'omitnan')

YMLeftPhantom = mean(rmoutliers(abs(YMCompression(leftBlockTop))),'omitnan')
YMRightPhantom = mean(rmoutliers(abs(YMCompression(rightBlockTop))),'omitnan')

YMLNegRem = YMCompression;
YMLNegRem(YMLNegRem<0) = NaN;
YMRNegRem = YMLNegRem;

YMLeftPhantomNegRem = mean(rmoutliers(abs(YMLNegRem(leftBlockTop))),'omitnan')
YMRightPhantomNegRem = mean(rmoutliers(abs(YMRNegRem(rightBlockTop))),'omitnan')

