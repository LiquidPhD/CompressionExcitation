% Code to get avg SW speed from left and right side of half/half phantom
clearvars;
baseFolder = 'D:\021722GelatinBlueYellow';
cutoffDepth = 750;
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

sensorMaskStack = sensorMaskStack(1:cutoffDepth,:,:);
leftMaskStack = leftMaskStack(1:cutoffDepth,:,:);
rightMaskStack = rightMaskStack(1:cutoffDepth,:,:);


for imageIterator = 1:length(foldersToProc)
load(foldersToProc(imageIterator).name,'TOF_speed','corrCoeff')
if size(corrCoeff,2) < size(TOF_speed,2)
    sizeAdjustment = size(TOF_speed,2) - size(corrCoeff,2);
    corrCoeff = [corrCoeff zeros([size(TOF_speed,1) sizeAdjustment])];
end
if imageIterator == 1
    figure; imagesc(TOF_speed);
    caxis([0 10])
    [pulseBoundarySelectors,~] = ginput(2);
    close all force;
    pulseBoundarySelectors = round(pulseBoundarySelectors);
end
TOF_speed = TOF_speed(1:cutoffDepth,:);
    TOF_speed(:,pulseBoundarySelectors(1):pulseBoundarySelectors(2)) = NaN;
leftSpeed = TOF_speed.*leftMaskStack(:,:,imageIterator);
rightSpeed = TOF_speed.*rightMaskStack(:,:,imageIterator);

leftSpeed(leftSpeed==0) = NaN;
rightSpeed(rightSpeed==0) = NaN;
% for k = 1:size(leftSpeed,1)
%     for l = 1:size(leftSpeed,2)
%         if corrCoeff(k,l) > 0.7
%             leftSpeed(k,l) = leftSpeed(k,l);
%             rightSpeed(k,l) = rightSpeed(k,l);
%         else
%             leftSpeed(k,l) = NaN;
%             rightSpeed(k,l) = NaN;
%         end
%     end
% end
leftSpeedAvg(:,imageIterator) = mean(leftSpeed,'all','omitnan');
rightSpeedAvg(:,imageIterator) = mean(rightSpeed,'all','omitnan');

YMFromSpeedLeft(:,imageIterator) = 3*1000*leftSpeedAvg(:,imageIterator)^2;
YMFromSpeedRight(:,imageIterator) = 3*1000*rightSpeedAvg(:,imageIterator)^2;
end

save([baseFolder,filesep,'speedAndYMFromSamples.mat'],'leftSpeedAvg','rightSpeedAvg','YMFromSpeedLeft','YMFromSpeedRight')