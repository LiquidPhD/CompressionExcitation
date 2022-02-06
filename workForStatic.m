
figure; imagesc(YMCompression)
caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])

YMCompressionOG = YMCompression;
YMCompression(badRegion==1) = NaN;

leftBlock = roipoly;

leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlock~=1) = NaN;
rightYM(rightBlock~=1)=NaN;
figure; imagesc(leftYM);
figure; imagesc(rightYM)

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



 I = uint8(255*mat2gray(abs(IQData(:,:,k)),[0 1e8]));
    bg = ind2rgb(I,gray(255));
    bgImg = double(bg);
    alphaFactor = 0.5;
    bgImgAlpha = (1 - alphaFactor) .* bgImg;
    
    overlayImage = YMCompression;
    overlayImage(leftBlock~=1 & rightBlock~=1) = NaN;
    figure; imagesc(overlayImage)
    II = uint8(255*mat2gray(overlayImage, [0 255]));
    II(II==255) = 0;
    im = ind2rgb(II,jet(255));

    fgImg = double(im);
    fgImgAlpha = alphaFactor .* fgImg;
    
    fusedImg = fgImgAlpha + bgImgAlpha;

    gcf=figure(1);
%     imagesc(xaxis,zaxis(1:y),fusedImg)
    imagesc(fusedImg)
    
     I = uint8(255*mat2gray(abs(IQData(:,:,k+1)),[0 1e8]));
    bg = ind2rgb(I,gray(255));
    bgImg = double(bg);
    alphaFactor = 0.5;
    bgImgAlpha = (1 - alphaFactor) .* bgImg;
    
    overlayImage = YMCompression;
    overlayImage(leftBlock~=1 & rightBlock~=1) = NaN;
   overlayImage(~isnan(overlayImage)) = 1;
   figure; imagesc(overlayImage)
    overlayImage(isnan(overlayImage))=0;
   figure; imagesc(labeloverlay(bg,overlayImage))
   
   
   
   
   
    I = uint8(255*mat2gray(abs(IQData(:,:,k)),[0 1e8]));
    bg = ind2rgb(I,gray(255));
    bgImg = double(bg);
   figure; imagesc(I)
% caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])

YMCompressionOG = YMCompression;
YMCompression(badRegion==1) = NaN;

leftBlock = roipoly;

leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlock~=1) = NaN;
rightYM(rightBlock~=1)=NaN;
figure; imagesc(leftYM);
figure; imagesc(rightYM)

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

overlayImage = YMCompression;
overlayImage(leftBlock~=1 & rightBlock~=1) = NaN;
overlayImage(~isnan(overlayImage)) = 1;
figure; imagesc(overlayImage)
overlayImage(isnan(overlayImage))=0;
figure; imagesc(labeloverlay(bg,overlayImage))

