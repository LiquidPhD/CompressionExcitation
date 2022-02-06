% Either load label data or make it
try
    load([baseFolder,filesep,'labels.mat'])
catch
volumeSegmenter
uiwait;
save([baseFolder,filesep,'labels.mat'],'labels')
end


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

VELOCITIES = VELOCITIES(:,:,1:40);

for imageIterator = 1:size(VELOCITIES,3)
    YMTotalSmoothed = (VELOCITIES(:,:,imageIterator).^2)*3*1000; % in Pa
    
    
% MODIFY FOR LEFT AND RIGHT PHANTOMS.    
speed_sensor = VELOCITIES(:,:,imageIterator).*sensorMaskStack(1:size(VELOCITIES,1),:,imageIterator);
% speed_sensor(speed_sensor==0) = NaN;
speed_sensor(isnan(speed_sensor))=0;
CC = CC(1:size(speed_sensor,1),:,:);
CCCurrent = zeros([size(speed_sensor,1) size(speed_sensor,2)]);
CCCurrent(:,21:end) = squeeze(CC(:,:,imageIterator));
for k = 1:size(speed_sensor,1)
    for l = 1:size(speed_sensor,2)
        if CCCurrent(k,l) > 0.7
            speed_sensorCC(k,l) = speed_sensor(k,l);
        else
            speed_sensorCC(k,l) = NaN;
        end
    end
end


    YMLeftAvg3 = mean(YMTotalSmoothed.*leftMaskStack(1:size(YMTotalSmoothed,1),:,imageIterator),'all','omitnan');
YMRightAvg3 = mean(YMTotalSmoothed.*rightMaskStack(1:size(YMTotalSmoothed,1),:,imageIterator),'all','omitnan');

  
close all force; 
figure;
ax1 = axes;

TotalMask = 2*leftMaskStack(:,:,imageIterator) + rightMaskStack(:,:,imageIterator);

[leftLabeled, ~] = bwlabel(leftMaskStack(:,:,imageIterator) ~= 0);
measurements = regionprops(leftLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsLeft = allCentroids(1);
yCentroidsLeft = allCentroids(2);

[rightLabeled, ~] = bwlabel(rightMaskStack(:,:,imageIterator) ~= 0);
measurements = regionprops(rightLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsRight = allCentroids(1:2:end);
yCentroidsRight = allCentroids(2:2:end);


TEST = mat2gray(abs(IQData(1:size(VELOCITIES,1),:,imageIterator)));
% TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
TEST = 255*TEST;
imagesc(ax1,xaxis,zaxis(1:size(VELOCITIES,1)),TEST(1:size(VELOCITIES,1),:));
ylabel('Depth (mm)')
xlabel('Distance (mm)')
colormap(ax1,'gray');
caxis([0 25])
ax2 = axes;
imagesc(ax2,TotalMask(1:size(VELOCITIES,1),:),'alphadata',TotalMask(1:size(VELOCITIES,1),:))
% imagesc(ax2,YMTotalSmoothed(1:670,:),'alphadata',YMTotalSmoothed(1:670,:)>0);
% colormap(ax2,'jet');
% TEST2 = YMTotalSmoothed(1:670,:)
% TEST2(TEST2==0) = NaN;
% CLow = quantile(TEST2,0.25,'all');
% CHigh = quantile(TEST2,0.75,'all');
% caxis(ax2,[0 1000]);
% caxis(ax2,[CLow CHigh])
  text(ax2,xCentroidsLeft,yCentroidsLeft,num2str(YMLeftAvg3/1000),'FontSize',20, ...
        'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
    text(ax2,xCentroidsRight,yCentroidsRight,num2str(YMRightAvg3/1000),'FontSize',20, ...
        'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
% h = colorbar;
% ylabel(h,'YM (kPa)')
title(ax1,"Mean Young's modulus using shear wave")
