function [particleVelocity] = directionalFiltering(particleVelocity)

% Perform directional filtering
parfor imageDepth = 1:size(particleVelocity,1)
    TEST = squeeze(particleVelocity(imageDepth,:,:));
    Y = fft2(TEST);
    %         FFTFILTERBOT = imresize(FILT,[size(Y,1) size(Y,2)])
    %         FFTFILTERBOT(FFTFILTERBOT<0.1) = 0;
    %         FFTFILTERBOT(:,round(size(FFTFILTERBOT,2)/2):end) = 0;
    FFTFILTERRIGHT = ones([size(Y,1) size(Y,2)]);
    FFTFILTERRIGHT(1:round(size(FFTFILTERRIGHT,1)/2),1:round(size(FFTFILTERRIGHT,2)/2)) = 0;
    FFTFILTERRIGHT(round(size(FFTFILTERRIGHT,1)/2)+1:end,round(size(FFTFILTERRIGHT,2)/2)+1:end) = 0;
    FFTFILTERLEFT = fliplr(FFTFILTERRIGHT);
    NEWTOP = FFTFILTERLEFT .* Y;
    NEWBOT = FFTFILTERRIGHT .* Y;
    newleftifft = ifft2(ifftshift(NEWTOP));
    newrightifft = ifft2(ifftshift(NEWBOT));
    newPVLEFT(imageDepth,:,:) = newleftifft;
    newPVRIGHT(imageDepth,:,:) = newrightifft;
end

particleVelocityCombined = [newPVLEFT(:,1:round(size(newPVLEFT,2)/2),:) newPVRIGHT(:,round(size(newPVRIGHT,2)/2)+1:end,:)];
particleVelocityCombined2 = particleVelocityCombined.*2;
%    playWaveVideo(abs(particleVelocityCombined))
particleVelocity = abs(particleVelocityCombined);