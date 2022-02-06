 %% YM Calc
%% Justin Rippy
%% 10/04/2021
% Code to compare dynamic and static methods for YM calculation.

clearvars;
BF = 'E:\HalfHalfAgarDay2NewSensorAuto';
baseFolder = BF;

%% Dynamic
% Calculate YM using E 3ρc^2

% Load data
centerVelocity = load([baseFolder,filesep,'CenterVelocity2.mat'],'VELOCITIES');
centerVelocity = centerVelocity.VELOCITIES;

axisData = load([baseFolder,filesep,'CenterVelocity2.mat'],'xaxis','zaxis');
xaxis = axisData.xaxis;
zaxis = axisData.zaxis;

figure(1);
imagesc(xaxis,zaxis(1:size(centerVelocity,1)),centerVelocity(:,:,1));
colormap(jet);
[leftPhantom,xLeft,yLeft] = roipoly;
[rightPhantom,xRight,yRight] = roipoly;


leftCoords = [xLeft yLeft];
rightCoords = [xRight yRight];

[leftLabeled, ~] = bwlabel(leftPhantom ~= 0);
measurements = regionprops(leftLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsLeft = allCentroids(1);
yCentroidsLeft = allCentroids(2);

[rightLabeled, ~] = bwlabel(rightPhantom ~= 0);
measurements = regionprops(rightLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsRight = allCentroids(1:2:end);
yCentroidsRight = allCentroids(2:2:end);


% for k = 1:size(centerVelocity,3)
%     velocityFrame = centerVelocity(:,:,k);
%     leftAvgVelocity(:,:,k) = mean(velocityFrame(leftPhantom),'omitnan');
%     rightAvgVelocity(:,:,k) = mean(velocityFrame(rightPhantom),'omitnan');
%     
%     YMMap(:,:,k) = 3*1000*(velocityFrame).^2; % In Pascals
% end
% figure;
% mkdir([baseFolder,filesep,'DynamicYM']);
% for k = 1:size(YMMap,3)
%     YM = YMMap(:,:,k)/1000; % Gives kPa
%     t = tiledlayout(1,1);
%     ax1 = axes(t);
%     imagesc(xaxis,zaxis(1:size(YM,1)),abs(YM));
%     
%     colormap(jet)
%     caxis([0 300])
%     colorbar;
%     title({"Young's modulus";['Frame: ',num2str(k)]})
%     xlabel('Distance (mm)')
%     ylabel('Depth (mm)');
%     
%     YMLeftPhantom = mean(rmoutliers(YM(leftPhantom)),'omitnan');
%     YMRightPhantom = mean(rmoutliers(YM(rightPhantom)),'omitnan');
%     hold on;
%     h1 = drawpolygon('Position',leftCoords,'InteractionsAllowed','none','FaceAlpha',0);
%     h1.Color = 'yellow'
%     h2 = drawpolygon('Position',rightCoords,'InteractionsAllowed','none','FaceAlpha',0);
%     h2.Color = 'yellow'
%     
%     hold off;
%     ax2 = axes(t)
%     xlim(ax2,[1 length(xaxis)]);
%     ylim(ax2,[1 size(centerVelocity,1)]);
%     ax2.YDir = 'reverse'
%     ax2.YAxis.Visible = 'off'
%     ax2.XAxis.Visible = 'off';
%     ax2.Color = 'none';
%     %         hline(depthValue,'w')
%     text(ax2,xCentroidsLeft,yCentroidsLeft,num2str(YMLeftPhantom),'FontSize',20, ...
%         'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
%     text(ax2,xCentroidsRight,yCentroidsRight,num2str(YMRightPhantom),'FontSize',20, ...
%         'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
%     export_fig([baseFolder,filesep,'DynamicYM',filesep,'YMFrame',num2str(k),'.png'],'-native')
%     
%     
% end


% Circshift and moving selected region up over time
circshiftMultiplier = 1
for k = 1:size(centerVelocity,3)
    velocityFrame = centerVelocity(:,:,k);
    LP = leftPhantom;
    RP = rightPhantom;
    if mod(k,2) == 0
    LP = circshift(leftPhantom,-1*circshiftMultiplier,1);
    RP = circshift(rightPhantom,-1*circshiftMultiplier,1);
    end
    leftAvgVelocity(:,:,k) = mean(velocityFrame(LP),'omitnan');
    rightAvgVelocity(:,:,k) = mean(velocityFrame(RP),'omitnan');
    
    YMMap(:,:,k) = 3*1000*(velocityFrame).^2; % In Pascals
end

figure;
mkdir([baseFolder,filesep,'DynamicYM']);
deltaZ = diff(zaxis);
deltaZ = deltaZ(1);
for k = 1:size(YMMap,3)
    YM = YMMap(:,:,k)/1000; % Gives kPa
    t = tiledlayout(1,1);
    ax1 = axes(t);
    imagesc(xaxis,zaxis(1:size(YM,1)),abs(YM));
    if mod(k,2) == 0
    leftCoords =[leftCoords(:,1) leftCoords(:,2)-(1*circshiftMultiplier)*deltaZ]
    rightCoords =[rightCoords(:,1) rightCoords(:,2)-(1*circshiftMultiplier)*deltaZ]
    end
    colormap(jet)
    caxis([0 300])
    colorbar;
    title({"Young's modulus";['Compression: ',num2str((k-1)*0.1),' mm']})
    xlabel('Distance (mm)')
    ylabel('Depth (mm)');
    
    YMLeftPhantom = mean(rmoutliers(YM(leftPhantom)),'omitnan');
    YMRightPhantom = mean(rmoutliers(YM(rightPhantom)),'omitnan');
    hold on;
    h1 = drawpolygon('Position',leftCoords,'InteractionsAllowed','none','FaceAlpha',0);
    h1.Color = 'yellow'
    h2 = drawpolygon('Position',rightCoords,'InteractionsAllowed','none','FaceAlpha',0);
    h2.Color = 'yellow'
    
    hold off;
    ax2 = axes(t)
    xlim(ax2,[1 length(xaxis)]);
    ylim(ax2,[1 size(centerVelocity,1)]);
    ax2.YDir = 'reverse'
    ax2.YAxis.Visible = 'off'
    ax2.XAxis.Visible = 'off';
    ax2.Color = 'none';
    %         hline(depthValue,'w')
    text(ax2,xCentroidsLeft,yCentroidsLeft,num2str(YMLeftPhantom),'FontSize',20, ...
        'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
    text(ax2,xCentroidsRight,yCentroidsRight,num2str(YMRightPhantom),'FontSize',20, ...
        'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
    export_fig([baseFolder,filesep,'DynamicYM',filesep,'YMFrame',num2str(k),'.png'],'-native')
    
    if k == 1
        gif([baseFolder,filesep,'DynamicYM',filesep,'YMGIF.gif'],0.2)
    else
        gif
    end
    
end

%% Static portion

clearvars -except BF lowerBound;
baseFolder = BF;

load([baseFolder,filesep,'STRAIN2.mat']);
load([baseFolder,filesep,'combinedIQData.mat']);

% Instron values in kPa
instronValueStiff = mean([242.35154,224.50661,223.05130]); 
instronValueSoft = mean([175.65058,171.56215,172.42562]);
instronValueSensor = mean([36.60783,35.17744,38.64111]);

% Use sensor YM to calculate strain for bottom samples.
%  E = σ/ε
% figure; imagesc(mat2gray(abs(IQData(:,:,2))));
% colormap(gray);
% caxis([0 0.2])
% [sensorRegion,sensorX,sensorY] = roipoly;
figure;
for k = 1:size(BScanStrainAll,3)
    
himage= imagesc(mat2gray(abs(IQData(:,:,k))));
colormap(gray);
caxis([0 0.3])
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if k == 1
% caxis([i_thresh_low i_thresh_high])


h = impoly(gca);
PosTime = wait(h);
else
    h = impoly(gca,PosTime)
    PosTime = wait(h);
    
end
BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));
BW2(BW2>0) = 1;
close all force; 
strainSensorRegion = BScanStrainAllFast(:,:,k);
strainSensorRegion = strainSensorRegion.*BW(1:size(strainSensorRegion,1),:,:);
strainSensorRegion(strainSensorRegion==0) = NaN;
figure; imagesc(strainSensorRegion)

% Get mean sensor strain across sensor layer
strainSensorValues(:,k) = mean(strainSensorRegion,1,'omitnan');
figure; plot(strainSensorValues(:,k))
pause
end
% Get sigma
for k = 1:size(

% Calculate E for bottom layers using this strain

