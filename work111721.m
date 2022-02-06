%% NEED THESE STEPS.
%% Normalize to surface by subtraction


vec_phase_diff = permute(vec_phase_diff,[2 1 3]);
for VPDIndex = 1:size(vec_phase_diff,3)
imageArray = double(vec_phase_diff(:,:,VPDIndex));
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

playWaveVideo(VPD_smoothed)

playWaveVideo(VPD_smoothed,0.1,fireice,0.05)

    figure('WindowState','maximized')
for k = 1:size(vec_phase_diff,3)
    subplot(1,2,1)
    imagesc(vec_phase_diff(:,:,k));
    colormap(fireice);
    caxis([-0.05 0.05])
    subplot(1,2,2)
    imagesc(VPD_smoothed(:,:,k));
    colormap(fireice);
    caxis([-0.05 0.05])
    pause(0.2)
end

for k = 1:size(VPD_smoothed,3)
TEST = mat2gray(VPD_smoothed(:,:,k));
figure; histogram(TEST)
minVal = 0.12
maxVal = 0.62
numDivisions = 10
sizeOfDivision = (maxVal-minVal)/numDivisions

TEST(TEST<minVal+sizeOfDivision) = 0.2;
TEST(TEST<minVal+2*sizeOfDivision & TEST > minVal+sizeOfDivision) = 0.4;
TEST(TEST<minVal+3*sizeOfDivision & TEST > minVal+2*sizeOfDivision) = 0.6;
TEST(TEST<minVal+4*sizeOfDivision & TEST > minVal+3*sizeOfDivision) = 0.8;
TEST(TEST<minVal+5*sizeOfDivision & TEST > minVal+4*sizeOfDivision) = 1;
TEST(TEST<minVal+5*sizeOfDivision & TEST > minVal+4*sizeOfDivision) = 1;
TEST(TEST<minVal+5*sizeOfDivision & TEST > minVal+4*sizeOfDivision) = 1;
TEST(TEST<minVal+5*sizeOfDivision & TEST > minVal+4*sizeOfDivision) = 1;
TEST(TEST<minVal+5*sizeOfDivision & TEST > minVal+4*sizeOfDivision) = 1;
TEST(TEST<minVal+5*sizeOfDivision & TEST > minVal+4*sizeOfDivision) = 1;

reducedColors(:,:,k) = TEST;
end


%% Histogram equalization before this.

p = twomodegauss(0.5,0.2,0.0001,0.0001,1,0.01,0); % Actually a 1-mode Gaussian...
for dispIndex = 1:size(VPD_smoothed,3)
    tempImage = mat2gray(VPD_smoothed(:,:,dispIndex));
    g = histeq(tempImage,p);
    VPDGauss(:,:,dispIndex) = g;
end


playWaveVideo(VPDGauss,0.2,fireice)
for k = 1:size(VPDGauss,3)-1
TEST(:,:,k) = movvar(VPDGauss(:,:,k+1)-VPDGauss(:,:,k),25);
end
playWaveVideo(TEST)

figure; 
for k = 1:size(TEST,3)
    imagesc(TEST(:,:,k))
%     pause(0.1)
    intermed = TEST(:,:,k);
    intermed(intermed>0.01) = NaN;
    TEST2(:,:,k) = intermed;
end

% playWaveVideo(TEST2)

TEST3 = VPDGauss;
for k = 1:size(TEST3,3)-3
    intermed = TEST3(:,:,k);
    intermed(isnan(TEST2(:,:,k))) = 0;
    TEST3(:,:,k) = intermed;
end
playWaveVideo(TEST3)




p = twomodegauss(0.5,0.2,0.0001,0.0001,1,0.01,0); % Actually a 1-mode Gaussian...
for dispIndex = 1:size(vec_phase_diff,3)
    tempImage = mat2gray(vec_phase_diff(:,:,dispIndex));
    g = histeq(tempImage,p);
    VPDGauss(:,:,dispIndex) = g;
end

for VPDIndex = 1:size(vec_phase_diff,3)
imageArray = double(VPDGauss(:,:,VPDIndex));
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

playWaveVideo(VPD_smoothed)



 figure;
            imagesc(vec_phase_diff(:,:,10));
            [x,y] = ginput(1)
            
            waveform_left = squeeze(vec_phase_diff(round(y),round(x),:))

  Ts = diff(taxis);
  Ts = Ts(1);
            Line = waveform_left;
            FFT = fft(Line,2^12);
            freq = linspace(0,1,2^12)/Ts;
            
            Time = [0:Ts:(length(waveform_left)-1)*Ts]*1e3;
            
            fig = figure;
            plot(freq,abs(FFT))
            grid on
            ylabel('Magnitude (Arb.)');
            xlabel('Frequency (Hz)');
            title('FFT of the signal');
            axis([0 4000 0 max(abs(FFT))*1.1])
            h = imrect(gca,[1e3-200,-200,400,(max(abs(FFT))*1.1)+400]);
            position = wait(h);
            b = fir1(50,[position(1)*2*Ts (position(1)+position(3))*2*Ts]);
            
            for k = 1:size(vec_phase_diff,3)
                for l = 1:size(vec_phase_diff,1)
               newProfiles(l,:,k) = filter(b,1,vec_phase_diff(l,:,k));     
                end
            end
playWaveVideo(newProfiles,0.2,fireice,0.05)


for VPDIndex = 1:size(vec_phase_diff,3)
imageArray = double(newProfiles(:,:,VPDIndex));
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

playWaveVideo(VPD_smoothed,0.2,fireice,0.05)

p = twomodegauss(0.5,0.2,0.0001,0.0001,1,0.01,0); % Actually a 1-mode Gaussian...
for dispIndex = 1:size(VPD_smoothed,3)
    tempImage = mat2gray(VPD_smoothed(:,:,dispIndex));
    g = histeq(tempImage,p);
    VPD_smoothed(:,:,dispIndex) = g;
end

figure;
    for k = 1:40 %size(particleVelocity,3)
        imagesc(xaxis,zaxis,abs(VPD_smoothed(:,:,k)))
        xlabel('Distance (mm)');
        ylabel('Distance (mm');
%         caxis([-1e-7 1e-7])
% caxis([-0.1 0.1])
% caxis([-0.02 0.02])
        colormap(fireice);
        colorbar;
        title({'Particle Velocity';[num2str((folderIndex-1)*0.1),' mm compressed']})
        if k == 1
            gif([baseFolder,filesep,'particleVelocity.gif'],'DelayTime', 0.1)
        else
            gif
        end
        
    end
