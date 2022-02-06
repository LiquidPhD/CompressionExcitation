clearvars;
close all force;
sensorInstron = [47.91195
44.8882
39.04582
40.04896
39.25353
38.7425
39.2035
    ]
% sensorInstron = [32.10241
% 30.76959
% 28.32373
% 27.39786
% ]
sensorInstron = mean(sensorInstron);
blueInstron = [89.57214
83.90152
84.22232
82.16962
80.83222
79.64514

    ]
blueInstron = mean(blueInstron)
redInstron = [376.10455
348.21764
342.70138
339.93964
338.0063
327.22941
319.80587

    ]
redInstron = mean(redInstron)


baseFolder = 'E:\ThinAgarHalfHalf121021';
%
% [lowerBound] = selectLowerBound(folders)
startDepth = 1;
lowerBound = 1395

%% Process static data

% Load static data
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
if ~isempty(rdir([baseFolder,filesep,'combinedIQData.mat']))
    load([baseFolder,filesep,'combinedIQData.mat'])
    folderIndex = 1;
else
[IQData,Parameters] = loadStaticData(folders,lowerBound,length(folders));
folderIndex = 1;
save([baseFolder,filesep,'combinedIQData.mat'],'IQData','Parameters');
% load('StrainImages.mat')
% load('combinedIQData.mat')
end
[vec_phase_diff] = VPD(IQData);
[particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters);
% particleVelocity = permute(vec_phase_diff,[2 1 3]);


sdl = ones([1 size(IQData,2)])
[~,Loupas_phase_shift] = ...
    Loupas_estimator_USE(IQData, sdl);
Loupas_phase_shift = permute(Loupas_phase_shift,[2 1 3]);

[Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
% Nz = 1395;
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);
% zaxis = zaxis(startDepth:end);
images = mat2gray(abs(IQData));

try
    load('labels.mat')
catch
volumeSegmenter
end

imageIterator = 2
% load phase difference
vec_image=(vec_phase_diff(:,:,imageIterator));
% vec_image = particleVelocity(:,:,imageIterator)
% Segment image into parts


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


vec_sensor = vec_image.*sensorMaskStack(:,:,imageIterator);
vec_left = vec_image.*leftMaskStack(:,:,imageIterator);
vec_right = vec_image.*rightMaskStack(:,:,imageIterator);

LBL = modefilt(sensorMaskStack(:,:,imageIterator),[51 51]);
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


vec_image_sensor=vec_image.*mask;
vec_image_sample=vec_image.*mask2;

unwrapped_sensor=phase_unwrap(vec_image_sensor);

for pos=1:size(vec_image_sensor,2)
    
    surf_dis(pos)=nanmean(unwrapped_sensor(sensor_sdl(pos)+5:sensor_sdl(pos)+20,pos));
    corrected_sensor(:,pos)=(unwrapped_sensor(:,pos)-surf_dis(pos)).*mask(:,pos);
end
    

unwrapped_sensor(sensor_sdl(pos)+5:sensor_sdl(pos)+20,pos);

figure;
plot(unwrapped_sensor(:,124));
figure;
plot(corrected_sensor(:,124));


figure;
imagesc(corrected_sensor);
figure;
imagesc(unwrapped_sensor);

% corrected_sensor(corrected_sensor==0) = NaN;

unwrapped_sample=phase_unwrap(vec_image_sample);

for pos=1:size(vec_image_sample,2)
    surf_dis(pos)=nanmedian(unwrapped_sample(sensor_bdl(pos):sensor_bdl(pos)+10,pos));
    corrected_sample(:,pos)=(unwrapped_sample(:,pos)-surf_dis(pos)).*mask2(:,pos);
end

imagesc(corrected_sample);
% corrected_sample(corrected_sample==0) = NaN;

corrected_left = corrected_sample.*leftMaskStack(:,:,imageIterator);
% corrected_left(corrected_left==0) = NaN;
corrected_right= corrected_sample.*rightMaskStack(:,:,imageIterator);
% corrected_right(corrected_right==0) = NaN;
figure; imagesc(corrected_left); figure; imagesc(corrected_right)

close all force;

x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
vec_image = abs(corrected_sensor)';
BScanWidth = size(vec_image,1);
sensor_sdl = ones([1 size(vec_image,2)])
sensor_bdl = ones([1 size(vec_image,2)]).*size(vec_image,2);
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
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl(x):sensor_bdl(1,x)),1);
    % find nans
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl(1,x):sensor_bdl(1,x),x) = feval(fitresult,1:length(x_disp));
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
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl(1,x):sensor_bdl(1,x),x));
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
    %     title(pathname,'Interpreter','none');
%     disp(num2str(x/BScanWidth))
    
end

corrected_sensor = fitted_displacement;

x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
vec_image = corrected_left';
BScanWidth = size(vec_image,1);
sensor_sdl = ones([1 size(vec_image,2)])
sensor_bdl = ones([1 size(vec_image,2)]).*size(vec_image,2);
how_many_depths = size(vec_image,2);

fitted_displacement = zeros([size(vec_image,2) size(vec_image,1)])
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
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl(x):sensor_bdl(1,x)),1);
    % find nans
    NaNChecker = isnan(x_disp)
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl(1,x):sensor_bdl(1,x),x) = feval(fitresult,1:length(x_disp));
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
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl(1,x):sensor_bdl(1,x),x));
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
    %     title(pathname,'Interpreter','none');
    
end

corrected_left = fitted_displacement;

x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
vec_image = corrected_right';
BScanWidth = size(vec_image,1);
sensor_sdl = ones([1 size(vec_image,2)])
sensor_bdl = ones([1 size(vec_image,2)]).*size(vec_image,2);
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
    
    x_disp = median(vec_image( x_start:x_end,sensor_sdl(x):sensor_bdl(1,x)),1);
    % find nans
    NaNChecker = isnan(x_disp)
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl(1,x):sensor_bdl(1,x),x) = feval(fitresult,1:length(x_disp));
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
%     plot(xData, yData, xData, fitted_displacement(sensor_sdl(1,x):sensor_bdl(1,x),x));
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
    %     title(pathname,'Interpreter','none');
    
end

corrected_right = fitted_displacement;

figure; imagesc(corrected_sensor);
figure; imagesc(corrected_left);
figure; imagesc(corrected_right)

corrected_sensor = abs(corrected_sensor).*sensorMaskStack(:,:,imageIterator);

corrected_left = abs(corrected_left).*leftMaskStack(:,:,imageIterator);
% corrected_left(corrected_left==0) = NaN;
corrected_right= abs(corrected_right).*rightMaskStack(:,:,imageIterator);
% corrected_right(corrected_right==0) = NaN;
figure; imagesc(corrected_sensor);
figure; imagesc(corrected_left);
figure; imagesc(corrected_right)
% 
% corrected_sensor(corrected_sensor==0) = NaN;
% corrected_left(corrected_left==0) = NaN;
% corrected_right(corrected_right==0) = NaN;

corrected_sensor2 = medfilt2(corrected_sensor,[100 100]);
figure; imagesc(corrected_sensor2)
corrected_left2 = medfilt2(corrected_left,[100 100]);
figure; imagesc(corrected_left2)
corrected_right2 = medfilt2(corrected_right,[100 100]);
figure; imagesc(corrected_right2);

strain_sensor=strain_calculation2(corrected_sensor,zaxis(2),1,100);
strain_left=strain_calculation2(corrected_left,zaxis(2),1,100);
strain_right=strain_calculation2(corrected_right,zaxis(2),1,100);

strain_sensor2=strain_calculation2(corrected_sensor2,zaxis(2),1,100);
strain_left2=strain_calculation2(corrected_left2,zaxis(2),1,100);
strain_right2=strain_calculation2(corrected_right2,zaxis(2),1,100);

strain_sensor(isnan(strain_sensor)) = 0;
strain_left(isnan(strain_left)) = 0;
strain_right(isnan(strain_right)) = 0;

strain_sensor2(isnan(strain_sensor2)) = 0;
strain_left2(isnan(strain_left2)) = 0;
strain_right2(isnan(strain_right2)) = 0;
strain_whole2 = strain_sensor2+strain_left2+strain_right2;
figure;
imagesc(xaxis,zaxis,abs(strain_whole2)/1000); colormap(hot); colorbar; caxis([0,1]);
title('med filtered')

% strain_sensor = strain_sensor';
% strain_left = strain_left';
% strain_right = strain_right';
% strain_whole=strain_sensor.*sensorMaskStack(:,:,imageIterator)+strain_left.*double(logical(leftMaskStack(:,:,imageIterator)))+strain_right.*double(logical(rightMaskStack(:,:,imageIterator)));
strain_whole = strain_sensor+strain_left+strain_right;
figure;
% strain_whole(isnan(strain_whole)) = 0;
imagesc(xaxis,zaxis,abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,1]);


figure;
% strain_whole(isnan(strain_whole)) = 0;
TEST = abs(strain_whole)/1000;
TEST(TEST<0.1) = NaN;
imagesc(xaxis,zaxis,TEST); colormap(hot); colorbar; caxis([0,1]);


figure; imagesc(strain_sensor); figure; imagesc(strain_left); figure; imagesc(strain_right)


strain_sensor2 = strain_sensor;
figure; imagesc(abs(strain_sensor2))
colormap(jet)
% CAX = [quantile(abs(strain_sensor2)./1000,0.25) quantile(abs(strain_sensor2)./1000,0.75)]
% caxis(CAX)

disp('Select boundary region')
boundary = roipoly;
strain_sensor2(boundary==1) = NaN;
strain_sensor2(strain_sensor==0) = NaN;

stressSensor = sensorInstron.*strain_sensor2;
stressSensor = mean(stressSensor,1,'omitnan');
figure; plot(stressSensor)
% stressSensor = fliplr(stressSensor)
SS = ones([size(strain_whole,1) size(strain_whole,2)]);
for k = 1:size(SS,2)
    stressSensor2 = stressSensor.*SS(:,k);
end


YMLeft = stressSensor2./(abs(strain_left)/1000);
YMLeft(isinf(YMLeft)) = 0;
figure; imagesc(YMLeft/1000); caxis([0 100]);

colormap(jet)
% CAX = [quantile(abs(strain_sensor2)./1000,0.25) quantile(abs(strain_sensor2)./1000,0.75)]
% caxis(CAX)

disp('Select good region')
goodLeft = roipoly;
% 
% leftMask = logical(leftMaskStack(:,:,imageIterator));
% YMLeftAvg = mean(rmoutliers(YMLeft(leftMask)),'omitnan');
YMLeftAvg = mean(rmoutliers(YMLeft(goodLeft)),'omitnan')


YMRight = stressSensor2./(abs(strain_right)/1000);
YMRight(isinf(YMRight)) = 0;
figure; imagesc(YMRight/1000); caxis([0 100]);

colormap(jet)
% CAX = [quantile(abs(strain_sensor2)./1000,0.25) quantile(abs(strain_sensor2)./1000,0.75)]
% caxis(CAX)

disp('Select good region')
goodRight = roipoly;
% 
% rightMask = logical(rightMaskStack(:,:,imageIterator));
YMRightAvg = mean(rmoutliers(YMRight(goodRight)),'omitnan');
YMLeftAvg = mean(rmoutliers(YMLeft(goodLeft)),'omitnan')

close all force;
disp(num2str(YMLeftAvg/1000))
disp(num2str(YMRightAvg/1000))


YMWhole = YMLeft + YMRight; 
figure; imagesc(YMWhole/1000); caxis([0 25]); colormap(jet)

YMTotal = stressSensor2./(abs(strain_whole)/1000);
YMTotal(isinf(YMTotal)) = 0;
YMTotal(isnan(YMTotal)) = 0;
figure; imagesc(YMTotal); colormap(jet); caxis([0 10000])


