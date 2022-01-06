clearvars;
baseFolder = 'E:\GelatinHalfHalf20211221';
startDepth = 1;
lowerBound = 1395

% Load velocity data
load([baseFolder,filesep,'Dynamic',filesep,'CenterVelocity.mat'])
% Load static data
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
folders = natsortfiles(folders);
if ~isempty(rdir([baseFolder,filesep,'combinedIQData.mat']))
    load([baseFolder,filesep,'combinedIQData.mat'])
    
temp = load([folders(folderIndex).folder,filesep,'Parameter.mat']);
Parameter = temp.Parameter;
Parameters = Parameter;
clearvars temp Parameter;
xAxis = Parameters.Trans.ElementPos(:,1);

% For old data backwards-compatibility
try
    Parameters.delta_z = Parameters.PData(3).PDelta(3)*Parameters.Trans.lambda;
    Parameters.delta_x = Parameters.PData(3).PDelta(1)*Parameters.Trans.lambda;
catch
    Parameters.delta_z = Parameters.PData(1).PDelta(3)*Parameters.Trans.lambda;
    Parameters.delta_x = Parameters.PData(1).PDelta(1)*Parameters.Trans.lambda;
end
Parameters.delta_t = 100*1e-6;

timeValue = 1e-4; % in seconds; 100 microseconds
tAxis = 0:timeValue:size(IQData,3);

AxialLengthToTest = 5;
try
    AxialResolution = Parameters.PData(3).PDelta(3);
catch
    AxialResolution = Parameters.PData(1).PDelta(3);
end
Parameters.M = ceil(AxialLengthToTest/AxialResolution);
Parameters.N = 3;
Parameters.fc = Parameters.Trans.frequency*1e6;
Parameters.c=1500;
    folderIndex = 1;
else
[IQData,Parameters] = loadStaticData(folders,lowerBound,length(folders));
folderIndex = 1;
save([baseFolder,filesep,'combinedIQData.mat'],'IQData','Parameters');
% load('StrainImages.mat')
% load('combinedIQData.mat')
end

% Take out known bad areas
% IQData(:,:,28) = []; IQData(:,:,17) = []; 

% Calc vec phase diff
[vec_phase_diff] = VPD(IQData);

% 2D Loupas
sdl = ones([1 size(IQData,2)]);
[~,Loupas_phase_shift] = ...
    Loupas_estimator_USE(IQData, sdl);
Loupas_phase_shift = permute(Loupas_phase_shift,[2 1 3]);

% Use frames 1:26
Loupas_phase_shift = Loupas_phase_shift(:,:,1:26);
% Sizes for further processing

[Nz,Nx,Nt]= size(Loupas_phase_shift);         % The dimensions of data in the z axis, the x axis and time.
% Nz = 1395;
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);
% zaxis = zaxis(startDepth:end);
images = mat2gray(abs(IQData));

% Either load label data or make it
try
    load([baseFolder,filesep,'labels.mat'])
catch
volumeSegmenter
save([baseFolder,filesep,'labels.mat'],'labels')
end

for imageIterator = 2:size(Loupas_phase_shift,3)
% imageIterator = 2

imageToProcess = Loupas_phase_shift(:,:,imageIterator);

% Separate out labels

for k = 1:size(labels,3)
    label = labels(:,:,k);
    sensor = zeros([size(label,1) size(label,2)]);
    sensor(label== 'Label1') = 1;
    left = zeros([size(label,1) size(label,2)]);
    left(label== 'Label2') = 2;
    right = zeros([size(label,1) size(label,2)]);
    right(label== 'Label3') = 3;
    sensorMaskStack(:,:,k) = sensor;
    leftMaskStack(:,:,k) = left;
    rightMaskStack(:,:,k) = right;
end


vec_sensor = imageToProcess.*sensorMaskStack(:,:,imageIterator);
vec_left = imageToProcess.*leftMaskStack(:,:,imageIterator);
vec_right = imageToProcess.*rightMaskStack(:,:,imageIterator);
vec_whole = vec_sensor+vec_left+vec_right;
% figure; imagesc(vec_whole); colormap(jet)

% Get top and bottom tracking
LBL = sensorMaskStack(:,:,imageIterator);
for pos=1:size(LBL,2)
    sensor_bdl(pos)=find(LBL(:,pos)==1,1,'last');
end


for pos=1:size(LBL,2)
    sensor_sdl(pos)=find(LBL(:,pos)==1,1)+10;
end

mask2=zeros(size(LBL(:,:)));
for pos=1:size(LBL,2)
    mask2(sensor_bdl(pos):end,pos)=1;
end

mask=zeros(size(LBL(:,:)));
for pos=1:size(LBL,2)
    mask(sensor_sdl(pos):sensor_bdl(pos),pos)=1;
end


unwrapped_whole = phase_unwrap(vec_whole);
for pos = 1:size(vec_whole,2)
    surf_dis(pos) = nanmean(unwrapped_whole(sensor_sdl(pos)+5:sensor_sdl(pos)+20,pos));
end
[xData, yData] = prepareCurveData( [], surf_dis );
    
    % Set up fittype and options.
    ft = fittype( 'fourier3' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    surf_dis2 = feval(fitresult,1:length(surf_dis));
    figure; plot(surf_dis2); hold on; plot(surf_dis,'g'); hold off; 
corrected_whole = unwrapped_whole-surf_dis2';
figure; imagesc(corrected_whole);

corrected_sensor = corrected_whole.*sensorMaskStack(:,:,imageIterator);
corrected_left = corrected_whole.*leftMaskStack(:,:,imageIterator);
corrected_right = corrected_whole.*rightMaskStack(:,:,imageIterator);
corrected_whole = corrected_sensor+corrected_left+corrected_right;



x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
corrected_whole(isnan(corrected_whole)) = 0;
vec_image = corrected_whole';
BScanWidth = size(vec_image,1);
sensor_sdl2 = ones([1 size(vec_image,2)])
sensor_bdl2 = ones([1 size(vec_image,2)]).*size(vec_image,2);
how_many_depths = size(vec_image,2);
fitted_displacement = zeros([size(vec_image,2) size(vec_image,1)]);
for x = 1:BScanWidth%size(displacement,1)
    if x < 1
        x_start = 1;
    else
        x_start = x;
    end
    
    if x > BScanWidth - x_kern
        x_end = BScanWidth;
    else
        x_end = x+x_kern;
    end
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl2(x):sensor_bdl2(1,x)),1);
    weights = logical(x_disp);

    % find nans
%     NaNChecker = isnan(x_disp)
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'cubicspline' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    opts.Weights = weights;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x) = feval(fitresult,1:length(x_disp));
%     figure(1);
%     subplot(3,1,1);
%     imagesc(vec_image'); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
%     %     title(pathname,'Interpreter','none');
%     
%     subplot(3,1,2);
%     %imagesc(displacement'); colormap(jet);
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x));
%     title(x);
%     axis tight;
%     drawnow;
%     subplot(3,1,3);
%     
%     imagesc(fitted_displacement); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
% % %         title(pathname,'Interpreter','none');
    
end

corrected_whole = fitted_displacement;
figure; imagesc(fitted_displacement);


corrected_sensor = abs(corrected_whole).*sensorMaskStack(:,:,imageIterator);

corrected_left = abs(corrected_whole).*leftMaskStack(:,:,imageIterator);
% corrected_left(corrected_left==0) = NaN;
corrected_right= abs(corrected_whole).*rightMaskStack(:,:,imageIterator);
figure; imagesc(corrected_sensor);
figure; imagesc(corrected_left);
figure; imagesc(corrected_right)

speed_sensor = VELOCITIES(:,:,imageIterator).*sensorMaskStack(1:size(VELOCITIES,1),:,imageIterator);
% speed_sensor(speed_sensor==0) = NaN;
speed_sensor(isnan(speed_sensor))=0;

x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
vec_image = speed_sensor';
BScanWidth = size(vec_image,1);
sensor_sdl3 = ones([1 size(vec_image,2)])
sensor_bdl3 = ones([1 size(vec_image,2)]).*size(vec_image,2);
how_many_depths = size(vec_image,2);

clearvars fitted_displacement
for x = 1:BScanWidth%size(displacement,1)
    if x < 1
        x_start = 1;
    else
        x_start = x;
    end
    
    if x > BScanWidth - x_kern
        x_end = BScanWidth;
    else
        x_end = x+x_kern;
    end
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl3(x):sensor_bdl3(1,x)),1);
    weights = logical(x_disp);

    % find nans
    NaNChecker = isnan(x_disp)
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'cubicinterp' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    opts.Weights = weights;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl3(1,x):sensor_bdl3(1,x),x) = feval(fitresult,1:length(x_disp));
%     figure(1);
%     subplot(3,1,1);
%     imagesc(vec_image'); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
%     %     title(pathname,'Interpreter','none');
%     
%     subplot(3,1,2);
%     %imagesc(displacement'); colormap(jet);
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl3(1,x):sensor_bdl3(1,x),x));
%     title(x);
%     axis tight;
%     drawnow;
%     subplot(3,1,3);
%     
%     imagesc(fitted_displacement); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
% % %         title(pathname,'Interpreter','none');
    
end


speed_sensor = fitted_displacement;
speed_sensor = speed_sensor.*sensorMaskStack(1:size(speed_sensor,1),:,imageIterator);

% figure; imagesc(speed_sensor)
% [~,y] = ginput(1)
% Shear wave YM: 3*rho*c^2
YM_sensor = 3*1000*speed_sensor.^2;
YMMean = mean(YM_sensor(sensor_bdl-50:sensor_bdl,:),1,'omitnan');
figure; plot(YMMean)



 [xData, yData] = prepareCurveData( [], YMMean );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    YMFit = feval(fitresult,1:length(YMMean));
    figure; plot(YMFit); hold on; plot(YMMean); hold off; 
    YMMean = YMFit;
    
    
strain_sensor=strain_calculation2(corrected_sensor,zaxis(2),1,10);
strain_left=strain_calculation2(corrected_left,zaxis(2),1,10);
strain_right=strain_calculation2(corrected_right,zaxis(2),1,10);
strain_oneshot = strain_calculation2(corrected_whole,zaxis(2),1,10);
strain_whole = strain_sensor+strain_left+strain_right;
close all force; 
figure;
imagesc(xaxis,zaxis,abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,75]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain')
export_fig([baseFolder,filesep,'StrainMap',num2str(imageIterator),'.png'])
save([baseFolder,filesep,'StrainMap',num2str(imageIterator),'.fig'])


stressSensor = abs(YMMean'.*strain_sensor); stressSensor = mean(stressSensor,1,'omitnan');figure; plot(stressSensor)
stressSensor = mean(stressSensor,1,'omitnan');
figure; plot(stressSensor)
% stressSensor = fliplr(stressSensor)
SS = ones([size(strain_whole,1) size(strain_whole,2)]);
for k = 1:size(SS,2)
    stressSensor2 = stressSensor.*SS(:,k);
end


YMLeft = abs(stressSensor2./(abs(strain_left)/1000));
YMLeft(isinf(YMLeft)) = 0;

YMRight = abs(stressSensor2./(abs(strain_right)/1000));
YMRight(isinf(YMRight)) = 0;
YMTotal = YMLeft + YMRight;

figure; imagesc(YMTotal/1000); colormap(jet); caxis([0 500])
% figure; imagesc(YMLeft); caxis([0 600]);
% 
% colormap(jet)
% 
% disp('Select good region')
% goodLeft = roipoly;
% 
% YMLeftAvg = mean(rmoutliers(YMLeft(goodLeft)),'omitnan')
% YMLeft(YMLeft==0) = NaN;
% YMLeftAvg2 = mean(rmoutliers(YMLeft),'all','omitnan')
% 
% YMRight = abs(stressSensor2./(abs(strain_right)/1000));
% YMRight(isinf(YMRight)) = 0;
% figure; imagesc(YMRight); caxis([0 100]);
% 
% colormap(jet)
% disp('Select good region')
% goodRight = roipoly;
% 
% YMRightAvg = mean(rmoutliers(YMRight(goodRight)),'omitnan')
% YMRight(YMRight==0) = NaN;
% YMRightAvg2 = mean(rmoutliers(YMRight),'all','omitnan')
% 
% YMTotal = YMLeft + YMRight;
% YMTotal = rmoutliers(YMTotal);
YMLeftAvg2 = mean(YMLeft/1000,'all','omitnan')
YMRightAvg2 = mean(YMRight/1000,'all','omitnan')

x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
% corrected_whole(isnan(corrected_whole)) = 0;
vec_image = YMTotal';
BScanWidth = size(vec_image,1);
sensor_sdl2 = ones([1 size(vec_image,2)]);
sensor_bdl2 = ones([1 size(vec_image,2)]).*size(vec_image,2);
how_many_depths = size(vec_image,2);
fitted_displacement = zeros([size(vec_image,2) size(vec_image,1)]);
for x = 1:BScanWidth%size(displacement,1)
    if x < 1
        x_start = 1;
    else
        x_start = x;
    end
    
    if x > BScanWidth - x_kern
        x_end = BScanWidth;
    else
        x_end = x+x_kern;
    end
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl2(x):sensor_bdl2(1,x)),1);
    weights = logical(x_disp);

    % find nans
%     NaNChecker = isnan(x_disp)
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'cubicspline' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    opts.Weights = weights;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x) = feval(fitresult,1:length(x_disp));
%     figure(1);
%     subplot(3,1,1);
%     imagesc(vec_image'); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
%     %     title(pathname,'Interpreter','none');
%     
%     subplot(3,1,2);
%     %imagesc(displacement'); colormap(jet);
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x));
%     title(x);
%     axis tight;
%     drawnow;
%     subplot(3,1,3);
%     
%     imagesc(fitted_displacement); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
% % %         title(pathname,'Interpreter','none');
    
end

YMTotalSmoothed = fitted_displacement/1000;

YMLeftAvg3 = mean(YMLeft.*leftMaskStack(:,:,imageIterator),'all','omitnan')
YMRightAvg3 = mean(YMRight.*rightMaskStack(:,:,imageIterator),'all','omitnan')

figure; imagesc(YMTotalSmoothed); caxis([0 750])
figure; imagesc(YMTotal); caxis([0 750])

close all force; 
figure;
ax1 = axes;

TEST = mat2gray(abs(IQData(:,:,imageIterator)));
% TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
TEST = 255*TEST;
imagesc(ax1,xaxis,zaxis,TEST);
ylabel('Depth (mm)')
xlabel('Distance (mm)')
colormap(ax1,'gray');
caxis([0 25])
ax2 = axes;
imagesc(ax2,YMTotalSmoothed,'alphadata',YMTotalSmoothed>0);
colormap(ax2,'jet');
caxis(ax2,[0 1000]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
colorbar;
title(ax1,"Young's modulus using shear wave")

export_fig([baseFolder,filesep,'YMOverlaySW',num2str(imageIterator),'.png'])
save([baseFolder,filesep,'YMOverlaySW',num2str(imageIterator),'.fig'])

sensorInstron = [163.65177
200.94859
214.18276
213.10991
214.29698
213.9483
215.44237
];
sensorInstron = mean(sensorInstron);
stressSensor = sensorInstron.*strain_sensor;
stressSensor = mean(stressSensor,1,'omitnan');
figure; plot(stressSensor)
% stressSensor = fliplr(stressSensor)
SS = ones([size(strain_whole,1) size(strain_whole,2)]);
for k = 1:size(SS,2)
    stressSensor2 = stressSensor.*SS(:,k);
end



YMLeft = abs(stressSensor2./(abs(strain_left)/1000));
YMLeft(isinf(YMLeft)) = 0;

YMRight = abs(stressSensor2./(abs(strain_right)/1000));
YMRight(isinf(YMRight)) = 0;
YMTotal = YMLeft + YMRight;

figure; imagesc(YMTotal); colormap(jet); caxis([0 500])
% figure; imagesc(YMLeft); caxis([0 600]);
% 
% colormap(jet)
% 
% disp('Select good region')
% goodLeft = roipoly;
% 
% YMLeftAvg = mean(rmoutliers(YMLeft(goodLeft)),'omitnan')
% YMLeft(YMLeft==0) = NaN;
% YMLeftAvg2 = mean(rmoutliers(YMLeft),'all','omitnan')
% 
% YMRight = abs(stressSensor2./(abs(strain_right)/1000));
% YMRight(isinf(YMRight)) = 0;
% figure; imagesc(YMRight); caxis([0 100]);
% 
% colormap(jet)
% disp('Select good region')
% goodRight = roipoly;
% 
% YMRightAvg = mean(rmoutliers(YMRight(goodRight)),'omitnan')
% YMRight(YMRight==0) = NaN;
% YMRightAvg2 = mean(rmoutliers(YMRight),'all','omitnan')
% 
% YMTotal = YMLeft + YMRight;
% YMTotal = rmoutliers(YMTotal);
YMLeftAvg2 = mean(YMLeft,'all','omitnan')
YMRightAvg2 = mean(YMRight,'all','omitnan')

x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
% corrected_whole(isnan(corrected_whole)) = 0;
vec_image = YMTotal';
BScanWidth = size(vec_image,1);
sensor_sdl2 = ones([1 size(vec_image,2)]);
sensor_bdl2 = ones([1 size(vec_image,2)]).*size(vec_image,2);
how_many_depths = size(vec_image,2);
fitted_displacement = zeros([size(vec_image,2) size(vec_image,1)]);
for x = 1:BScanWidth%size(displacement,1)
    if x < 1
        x_start = 1;
    else
        x_start = x;
    end
    
    if x > BScanWidth - x_kern
        x_end = BScanWidth;
    else
        x_end = x+x_kern;
    end
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl2(x):sensor_bdl2(1,x)),1);
    weights = logical(x_disp);

    % find nans
%     NaNChecker = isnan(x_disp)
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'cubicspline' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    opts.Weights = weights;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x) = feval(fitresult,1:length(x_disp));
%     figure(1);
%     subplot(3,1,1);
%     imagesc(vec_image'); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
%     %     title(pathname,'Interpreter','none');
%     
%     subplot(3,1,2);
%     %imagesc(displacement'); colormap(jet);
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x));
%     title(x);
%     axis tight;
%     drawnow;
%     subplot(3,1,3);
%     
%     imagesc(fitted_displacement); colormap(jet);
%     hold on;
%     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
%     hold off;
%     colorbar;
% % %         title(pathname,'Interpreter','none');
    
end

YMTotalSmoothed = fitted_displacement;

YMLeftAvg3 = mean(YMLeft.*leftMaskStack(:,:,imageIterator),'all','omitnan')
YMRightAvg3 = mean(YMRight.*rightMaskStack(:,:,imageIterator),'all','omitnan')

figure; imagesc(YMTotalSmoothed); caxis([0 750])
figure; imagesc(YMTotal); caxis([0 750])


figure;
ax1 = axes;
TEST = mat2gray(abs(IQData(:,:,imageIterator)));
% TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
TEST = 255*TEST;
imagesc(TEST);
colormap(ax1,'gray');
caxis([0 25])
ax2 = axes;
imagesc(ax2,YMTotalSmoothed,'alphadata',YMTotalSmoothed>0);
colormap(ax2,'jet');
caxis(ax2,[0 500]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
colorbar;
title(ax1,"Young's modulus using sensor instron value")

export_fig([baseFolder,filesep,'YMOverlayInstron',num2str(imageIterator),'.png'])
save([baseFolder,filesep,'YMOverlayInstron',num2str(imageIterator),'.fig'])

end