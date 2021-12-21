close all force; 
for filterSizeIterator = 10:10:100
opticFlow=opticalFlowFarneback; %% the method I used, I only tried LK and Farneback so far, and LK didn't show shit
opticFlow.NumPyramidLevels=4;
opticFlow.NeighborhoodSize=7;
opticFlow.FilterSize=filterSizeIterator;


images = mat2gray(abs(IQData));
flow_Vx=zeros(size(images));
flow_Vy=flow_Vx;
flow_O=flow_Vx;
flow_M=flow_Vx;
% 
% h = figure;
% movegui(h);
% hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
% hPlot = axes(hViewPanel);

% figure;
% imagesc(images(:,:,10),[0.55,0.8]); colormap(gray);
% figure;
% imagesc(imgaussfilt(images(:,:,10),3),[0.55,0.8]); colormap(gray);


for frame=1:size(images,3) %% I believe you have to use gray scale images.
   

    image=images(:,:,frame);
    flow = estimateFlow(opticFlow,image); %% This function gives you flow velocity in x and y, phase angle, and magnitude. I'm storing them all.
    
    flow_Vx(:,:,frame)=flow.Vx;
    flow_Vy(:,:,frame)=flow.Vy;
    flow_O(:,:,frame)=flow.Orientation;
    flow_M(:,:,frame)=flow.Magnitude;
%     
%     imshow(image,[0.5,.8])
%     title([num2str(frame),' out of ',num2str(size(images,3))]);
%     hold on
%     plot(flow,'DecimationFactor',[5 5],'ScaleFactor',2,'Parent',hPlot);
%     hold off
%     pause(0.1)
end

T = cumsum(flow_Vy,3);
figure; imagesc(T(:,:,end));
title(['Filter Size: ',num2str(filterSizeIterator)])
end
autoArrangeFigures


sensor_dl=[];
intensities = flow_Vy(:,:,end)
figure
himage= imagesc(intensities(:,:,end));
% caxis([i_thresh_low i_thresh_high])
colormap(gray);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

h = impoly(gca);
PosTime = wait(h);
BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));

for i = 1:size(BW,2)
    idx = find(BW2(:,i) == 1);
    sensor_dl(:,i) = idx(1);
end

figure
imagesc(intensities(:,:,1))
caxis([i_thresh_low i_thresh_high])
colormap(gray)
hold on
plot(smooth(sensor_dl+5),'LineWidth',3);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
hold off; 
% sensor_dl=smooth(sensor_dl+5);

T = cumsum(flow_Vy,3);
figure; imagesc(T(:,:,end));

hold on
plot(smooth(sensor_dl+5),'LineWidth',3);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
hold off; 
% 
% figure; 
% for k = 1:size(IQData,3)
%     imagesc(mat2gray(abs(IQData(:,:,k))));
%     colormap(gray)
%     caxis([0 0.2])
%     pause
% end
TEST = mean(T(sensor_dl:sensor_dl+20,:,end),1)



for k = 1:size(flow_Vy,3)
    
    sensor_dl=[];
intensities = flow_Vy(:,:,k)
figure
himage= imagesc(images(:,:,k));
% caxis([i_thresh_low i_thresh_high])
colormap(gray);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

if k == 1
    h = impoly(gca);
    PosTime = wait(h);
else
   h = impoly(gca,PosTime) 
   PosTime = wait(h);
end
BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));

for i = 1:size(BW,2)
    idx = find(BW2(:,i) == 1);
    sensor_dl(:,i) = idx(1);
end
sensors(:,k) = sensor_dl;
% figure(1);
% plot(sensor_dl);
% pause
close all force;
end


for k = 1:size(flow_Vy,3)
    figure(1);
    imagesc(images(:,:,k))
    hold on;
    plot(squeeze(sensors(:,k)),'g');
    hold off;
    pause
end

for k = 1:size(sensors,2)
    figure(1)
    subplot(2,1,1)
    imagesc(mat2gray(images(:,:,k)))
    caxis([0 0.3])
    colormap(gray)
    hold on;
    plot(sensors(:,k))
    hold off;
    subplot(2,1,2)
    plot(mean(flow_Vy(sensors(:,k):sensors(:,k)+10,:,k),1));
    pause
end


agarYM = [58.7117,59.94864,61.07395];
agarYM = mean(agarYM);
instronYM = agarYM; % kPa

deltaZ = Parameters.delta_z % in m

figure
himage= imagesc(images(:,:,k));
% caxis([i_thresh_low i_thresh_high])
colormap(gray);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

initialPhantomSize = sensors(:,1)
for k = 1:size(sensors,2)-1
    diffSensors(:,k) = sensors(:,k+1)-sensors(:,k);
end

% Using E = σ/ε = F/A / ΔL/L0

% instronYM = sigma / (diffSensors(:,k) / initialPhantomSize)
for k = 1:size(diffSensors,2)
sigma(:,k) = (instronYM*1000) * ((abs(squeeze(diffSensors(:,k)))*deltaZ) ./ (initialPhantomSize*deltaZ));
end
sigmaCS = cumsum(sigma,2);
figure; plot(sigmaCS(:,end))


for k = 1:size(b_scan_strain,1)
    YM(:,k) = squeeze(abs(sigmaCS(:,end))) ./ squeeze(abs(b_scan_strain(k,:))/1000)';
end
figure; 
imagesc(YM)
caxis([0 100])

