% Code for static calculation of YM

% Load IQ Data
clearvars -except baseFolder lowerBound Parameters
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
[IQData,VMIQ,vec_phase_diff] = loadStaticData(folders,lowerBound,30);
folderIndex = 1;
save([baseFolder,filesep,'combinedIQData.mat'],'IQData');
 [~,~,~,Parameters] = loadDynamicData(folders,folderIndex,lowerBound);
 axisData = load([baseFolder,filesep,'CenterVelocity2.mat'],'xaxis','zaxis');
xaxis = axisData.xaxis;
zaxis = axisData.zaxis;
% Complex conjugate generation - First frame with all others

for frame = 1:size(IQData,3)-1
    complex_phase_diff_comp(:,:,frame) = IQData(:,:,1).*conj(IQData(:,:,frame+1));
end

% Vector method of complex conj
IQData = IQData(1:lowerBound,:,:);
sdl = ones([1 size(IQData,1)]);
VMIQ = permute(complex_phase_diff_comp,[2 1 3]);
[vector_complex_OCE_data,vec_phase_diff_2] = ...
    vec_meth_snr(VMIQ,sdl,20);
vec_phase_diff_2 = permute(vec_phase_diff_2,[2 1 3]);

for k = 9%:size(vec_phase_diff_2,3)
TEST = vec_phase_diff_2(:,:,k);
TEST2 = phase_unwrap(TEST);
figure;
imagesc(TEST2);
figure;
plot(unwrap(squeeze(TEST2(:,100)),'b')); hold on;
plot(unwrap(squeeze(TEST2(:,450)),'r')); hold off;


% Selection of boundary region and calc of str 

figure; himage = imagesc(mat2gray(abs(IQData(:,:,k)))); caxis([0 0.1])
colormap(gray);
  
h = impoly(gca);
PosTime = wait(h);

BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));
BW2(BW2>0) = 1;

figure; imagesc(BW2)

TESTLINE = TEST2;
TESTLINE(BW2~=1) = NaN;
figure; imagesc(TESTLINE)

TESTLINE = mean(TESTLINE,1,'omitnan');
figure; plot(TESTLINE)

% Normalize to particular depth
TESTABS = abs(TEST);
figure; plot(mean(TESTABS(170:180,:),1))
TESTFITTING = smooth(mean(TESTABS(170:180,:),1),33);
XFITTING = 1:length(TESTFITTING);

TEST3 = phase_unwrap(TESTABS - TESTFITTING');
figure; plot(TEST2(100,:)); hold on; plot(TEST3(100,:)); hold off; 

figure; imagesc(TEST2); figure; imagesc(TEST3)

% Calc strain


  TEST3STRAIN = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    TEST3,100);  
TEST3STRAIN = TEST3STRAIN';
% TEST3STRAINFAST = strain_2D_slash_no_sensor(baseFolder,Parameters.delta_z,...
% abs(TEST3),50,100);

% figure; imagesc(TEST3STRAINFAST)
figure; imagesc(mat2gray(TEST3STRAIN));
caxis([0.5 0.8])

% figure; imagesc(TEST3STRAIN);


% YM Calc using instron data
Instron = [58.81076 54.06557 56.22892 52.13362];
Instron = mean(Instron);
stiffInstron = [281.64155 277.05164 267.3551 273.22116];
stiffInstron = mean(stiffInstron);
softInstron = [140.13218 134.14949 138.51699 135.98369];
softInstron = mean(softInstron);
% figure; himage = imagesc(mat2gray(abs(TEST3STRAIN))); caxis([0 0.5])
% colormap(gray);
%   
% h = impoly(gca);
% PosTime = wait(h);
% 
% BW = createMask(h,himage);
% BW2 = ceil(imgaussfilt(single(BW),3));
% BW2(BW2>0) = 1;

figure; imagesc(BW2)

TESTLINE = TEST3STRAIN;
TESTLINE(BW2~=1) = NaN;
figure; imagesc(TESTLINE)

TESTLINE = mean(TESTLINE,1,'omitnan');
figure; plot(TESTLINE)

% E = epsilon / sigma
% YM = stress / strain

stressCalculated = Instron .* TESTLINE;
TEST3STRAIN = TEST3STRAIN;
stressCalculated2(1:250) = median(stressCalculated(1:250));
stressCalculated2(250:length(stressCalculated)) = median(stressCalculated(250:length(stressCalculated)));
figure; plot(stressCalculated2)
for YMIndex = 1:size( TEST3STRAIN,1)
YMCompression(YMIndex,:) = (stressCalculated2 ./ TEST3STRAIN(YMIndex,:));

end




figure; imagesc(YMCompression)
caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])
errorSoft = 100; errorStiff = 100;
YMCompressionOG = YMCompression;
% while errorSoft > 10 || errorStiff > 10 
close all force; 
I = uint8(255*mat2gray(abs(IQData(:,:,k)),[0 1e8]));
figure; imagesc(I)
% caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])
disp('Select bad region')
badRegion = roipoly;


YMCompression = YMCompressionOG;
YMCompression(badRegion==1) = NaN;

[leftBlock,xLeft,yLeft] = roipoly;
[rightBlock,xRight,yRight] = roipoly;

leftCoords = [xLeft yLeft];
rightCoords = [xRight yRight];


[leftLabeled, ~] = bwlabel(leftBlock ~= 0);
measurements = regionprops(leftLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsLeft = allCentroids(1);
yCentroidsLeft = allCentroids(2);

[rightLabeled, ~] = bwlabel(rightBlock ~= 0);
measurements = regionprops(rightLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsRight = allCentroids(1:2:end);
yCentroidsRight = allCentroids(2:2:end);

leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlock~=1) = NaN;
rightYM(rightBlock~=1)=NaN;
% figure; imagesc(leftYM);
% figure; imagesc(rightYM)

YMLeftPhantom = mean(rmoutliers(leftYM),'omitnan');
YMRightPhantom = mean(rmoutliers(rightYM),'omitnan');

YMLeftPhantom = mean(rmoutliers(abs(YMCompression(leftBlock))),'omitnan');
YMRightPhantom = mean(rmoutliers(abs(YMCompression(rightBlock))),'omitnan');

YMLNegRem = YMCompression;
YMLNegRem(YMLNegRem<0) = NaN;
YMRNegRem = YMLNegRem;

YMLeftPhantomNegRem = mean(rmoutliers(abs(YMLNegRem(leftBlock))),'omitnan');
YMRightPhantomNegRem = mean(rmoutliers(abs(YMRNegRem(rightBlock))),'omitnan');
errorStiff = 100*((stiffInstron - YMLeftPhantomNegRem) / stiffInstron);
errorSoft = 100*((softInstron - YMRightPhantomNegRem) / softInstron);

disp(['Soft/Right ',' Instron: ',num2str(softInstron),' Calc: ',num2str(YMRightPhantomNegRem),' Error: ', num2str(errorSoft),'%' ])
disp(['Stiff/Left ',' Instron: ',num2str(stiffInstron),' Calc: ',num2str(YMLeftPhantomNegRem),' Error: ', num2str(errorStiff),'%'])

close all force; 


 I = uint8(255*mat2gray(abs(IQData(:,:,k)),[0 1e8]));
    bg = ind2rgb(I,gray(255));
%     bgImg = double(bg);
%     alphaFactor = 0.5;
%     bgImgAlpha = (1 - alphaFactor) .* bgImg;
%     
%     overlayImage = YMCompression;
%     overlayImage(leftBlock~=1 & rightBlock~=1) = NaN;
%     figure; imagesc(overlayImage)
%     II = uint8(255*mat2gray(overlayImage, [0 255]));
%     II(II==255) = 0;
%     im = ind2rgb(II,jet(255));
% 
%     fgImg = double(im);
%     fgImgAlpha = alphaFactor .* fgImg;
%     
%     fusedImg = fgImgAlpha + bgImgAlpha;
% 
%     gcf=figure(1);
% %     imagesc(xaxis,zaxis(1:y),fusedImg)
%     imagesc(fusedImg)
%     
%      I = uint8(255*mat2gray(abs(IQData(:,:,k+1)),[0 1e8]));
%     bg = ind2rgb(I,gray(255));
%     bgImg = double(bg);
%     alphaFactor = 0.5;
%     bgImgAlpha = (1 - alphaFactor) .* bgImg;
%     
%     overlayImage = YMCompression;
%     overlayImage(leftBlock~=1 & rightBlock~=1) = NaN;
%    overlayImage(~isnan(overlayImage)) = 1;
%    figure; imagesc(overlayImage)
%     overlayImage(isnan(overlayImage))=0;
%    figure; imagesc(labeloverlay(bg,overlayImage))
%    
%    
%    
%    
%    
%     I = uint8(255*mat2gray(abs(IQData(:,:,k)),[0 1e8]));
%     bg = ind2rgb(I,gray(255));
%     bgImg = double(bg);
%    figure; imagesc(I)
% % caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])
% 
% YMCompressionOG = YMCompression;
% YMCompression(badRegion==1) = NaN;
% 
% leftBlock = roipoly;
% 
% leftYM = YMCompression;
% rightYM = YMCompression;
% leftYM(leftBlock~=1) = NaN;
% rightYM(rightBlock~=1)=NaN;
% figure; imagesc(leftYM);
% figure; imagesc(rightYM)
% 
% YMLeftPhantom = mean(rmoutliers(leftYM),'omitnan');
% YMRightPhantom = mean(rmoutliers(rightYM),'omitnan');
% 
% YMLeftPhantom = mean(rmoutliers(abs(YMCompression(leftBlock))),'omitnan');
% YMRightPhantom = mean(rmoutliers(abs(YMCompression(rightBlock))),'omitnan');
% 
% YMLNegRem = YMCompression;
% YMLNegRem(YMLNegRem<0) = NaN;
% YMRNegRem = YMLNegRem;
% 
% YMLeftPhantomNegRem = mean(rmoutliers(abs(YMLNegRem(leftBlock))),'omitnan');
% YMRightPhantomNegRem = mean(rmoutliers(abs(YMRNegRem(rightBlock))),'omitnan');
% errorStiff = 100*((stiffInstron - YMLeftPhantomNegRem) / stiffInstron);
% errorSoft = 100*((softInstron - YMRightPhantomNegRem) / softInstron);
% 
% disp(['Soft/Right ',' Instron: ',num2str(softInstron),' Calc: ',num2str(YMRightPhantomNegRem),' Error: ', num2str(errorSoft),'%' ])
% disp(['Stiff/Left ',' Instron: ',num2str(stiffInstron),' Calc: ',num2str(YMLeftPhantomNegRem),' Error: ', num2str(errorStiff),'%'])

overlayImage = YMCompression;
overlayImage(leftBlock~=1 & rightBlock~=1) = NaN;
overlayImage(~isnan(overlayImage)) = 1;
% figure; imagesc(overlayImage)
overlayImage(isnan(overlayImage))=0;
figure; 
    t = tiledlayout(1,1);
        ax1 = axes(t);

imagesc(xaxis,zaxis(1:size(I,1)),labeloverlay(bg,overlayImage))

 hold on;
    h1 = drawpolygon('Position',leftCoords,'InteractionsAllowed','none','FaceAlpha',0);
    h1.Color = 'yellow'
    h2 = drawpolygon('Position',rightCoords,'InteractionsAllowed','none','FaceAlpha',0);
    h2.Color = 'yellow'
    
      
    hold off;
    ax2 = axes(t)
    xlim(ax2,[1 length(xaxis)]);
    ylim(ax2,[1 size(I,1)]);
    ax2.YDir = 'reverse'
    ax2.YAxis.Visible = 'off'
    ax2.XAxis.Visible = 'off';
    ax2.Color = 'none';
    %         hline(depthValue,'w')
    text(ax2,xCentroidsLeft,yCentroidsLeft,num2str(YMLeftPhantom),'FontSize',20, ...
        'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
    text(ax2,xCentroidsRight,yCentroidsRight,num2str(YMRightPhantom),'FontSize',20, ...
        'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);

    
disp(['Soft/Right ',' Instron: ',num2str(softInstron),' Calc: ',num2str(YMRightPhantomNegRem),' Error: ', num2str(abs(errorSoft)),'%' ])
disp(['Stiff/Left ',' Instron: ',num2str(stiffInstron),' Calc: ',num2str(YMLeftPhantomNegRem),' Error: ', num2str(abs(errorStiff)),'%'])
% end
    
export_fig([baseFolder,filesep,'StaticYM',filesep,'YMFrame',num2str(k),'.png'],'-native')

leftStiffness(:,k) = YMLeftPhantomNegRem;
rightStiffness(:,k) = YMRightPhantomNegRem;
errorLeft(:,k) = errorStiff;
errorRight(:,k) = errorSoft;
% save([baseFolder,filesep,'StaticYM',filesep,'Frame',num2str(k),'workspace.mat'],'-v7.3')
end