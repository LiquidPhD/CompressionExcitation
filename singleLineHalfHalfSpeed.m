% Code to get avg SW speed from left and right side of half/half phantom
clearvars;
baseFolder = 'D:\020322AgarHalfAndHalf';
%
foldersToProc = rdir([baseFolder,filesep,'Dynamic\*\**\*speedProcessing*']);
load([baseFolder,filesep,'labels.mat'])
for k = 1:size(labels,3)
    label = labels(:,:,k);
    sensor = zeros([size(label,1) size(label,2)]);
    sensor(label== 'Label1') = 1;
    left = zeros([size(label,1) size(label,2)]);
    left(label== 'Label2') = 2;
    right = zeros([size(label,1) size(label,2)]);
    right(label== 'Label3') = 3;
    sensorMaskStack(:,:,k) = logical(sensor);
    leftMaskStack(:,:,k) = logical(left);
    rightMaskStack(:,:,k) = logical(right);
end

sensorMaskStack = sensorMaskStack(1:750,:,:);
leftMaskStack = leftMaskStack(1:750,:,:);
rightMaskStack = rightMaskStack(1:750,:,:);
workingIndices = [10:75]
clearvars leftSpeedAvg rightSpeedAvg YMLeft YMRight
for imageIterator = 1:length(foldersToProc)
load(foldersToProc(imageIterator).name,'vec_phase_diff','Parameters')

 dispLeft = zeros([size(vec_phase_diff,1) size(vec_phase_diff,2) workingIndices(end)-workingIndices(1)+1]);
    dispRight = dispLeft;
    
     for imageDepth = 1:size(vec_phase_diff,1)
            TEST = squeeze(vec_phase_diff(imageDepth,:,workingIndices));
            Y = fft2(TEST);
            %         FFTFILTERBOT = imresize(FILT,[size(Y,1) size(Y,2)])
            %         FFTFILTERBOT(FFTFILTERBOT<0.1) = 0;
            %         FFTFILTERBOT(:,round(size(FFTFILTERBOT,2)/2):end) = 0;
            FFTFILTERBOT = ones([size(Y,1) size(Y,2)]);
            FFTFILTERBOT(1:round(size(FFTFILTERBOT,1)/2),1:round(size(FFTFILTERBOT,2)/2)) = 0;
            FFTFILTERBOT(round(size(FFTFILTERBOT,1)/2)+1:end,round(size(FFTFILTERBOT,2)/2)+1:end) = 0;
            FFTFILTERTOP = fliplr(FFTFILTERBOT);
            NEWTOP = FFTFILTERTOP .* Y;
            NEWBOT = FFTFILTERBOT .* Y;
            newtopifft = ifft2(ifftshift(NEWTOP));
            newbotifft = ifft2(ifftshift(NEWBOT));
            dispLeft(imageDepth,:,:) = abs(newtopifft);
            dispRight(imageDepth,:,:) = abs(newbotifft);
            disp([num2str(imageDepth),'/',num2str(size(vec_phase_diff,1))])
     end
        [VelocityLeft] = getSpeedSingleLineFixedDepth(dispLeft,Parameters,1,450)
                [VelocityRight] = getSpeedSingleLineFixedDepth(dispRight,Parameters,2,450)

% [VelocityLeft,VelocityRight] = getSpeedSingleLineFixedDepth(vec_phase_diff,Parameters,1,[10:75],450);
leftSpeed(:,imageIterator) = VelocityLeft;
rightSpeed(:,imageIterator) = VelocityRight;

YMLeft(:,imageIterator) = 3*1000*VelocityLeft^2/10000
YMRight(:,imageIterator) = 3*1000*VelocityRight^2/10000

end

figure; plot(leftSpeedAvg); hold on; plot(rightSpeedAvg); hold off;
figure; plot(YMLeft); hold on; plot(YMRight); hold off;
save([baseFolder,filesep,'singleLineSpeeds.mat'],'leftSpeed','rightSpeed','YMLeft','YMRight')