close all force;
minKernelSize = 10;
maxKernelSize = 50;
kernelCounter = 1;
delete(gcp('nocreate'));
parpool(12)
numberOfKernels = maxKernelSize-minKernelSize+1;
for kernelIndex = minKernelSize:maxKernelSize
    N_radius = kernelIndex;
    [TOF_speed(:,:,kernelCounter)] = runTOFparforTEST(particleVelocity,Parameters,kernelIndex,position,N_radius);
    % May need to adjust the weight equation below due to kernelIndex no
    % longer starting at 1.
    w(:,kernelCounter) = (2*(numberOfKernels - kernelCounter + 1)) / (numberOfKernels * (numberOfKernels + 1));
    kernelCounter = kernelCounter+1;
end
TOF_speed(isnan(TOF_speed)) =0;
for k = 1:size(TOF_speed,3)
    TOF_speed_MKWA(:,:,k) = w(:,k).*TOF_speed(:,:,k);
end
TOF_speed_MKWA = sum(TOF_speed_MKWA,3)
figure; imagesc(TOF_speed_MKWA)
caxis([0 8])
colormap(jet)
for k = 1:2:size(TOF_speed,3)
    figure; imagesc(TOF_speed(:,:,k))
    caxis([0 8])
    colormap(jet)
end
playWaveVideo(TOF_speed,0.2,jet,[0 8])