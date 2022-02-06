for imageIterator = 2:size(Loupas_phase_shift,3)
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
    x_disp(isnan(x_disp)) = 0;
    x_disp(isinf(x_disp)) = 0;

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
YM_sensor = 3*1000*speed_sensor.^2; % in Pa
YMMean = mean(YM_sensor(sensor_bdl-50:sensor_bdl,:),1,'omitnan');
figure; plot(YMMean)
YMMean(1:length(YMMean)/2) = mean(YMMean(50:150));
YMMean(length(YMMean)/2+1:length(YMMean)) = mean(YMMean(400:500));
figure; plot(YMMean)


%  [xData, yData] = prepareCurveData( [], YMMean );
%     
%     % Set up fittype and options.
%     ft = fittype( 'fourier8' );
%     opts = fitoptions( 'Method', 'LinearLeastSquares' );
%     opts.Robust = 'Bisquare';
%     
%     % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     YMFit = feval(fitresult,1:length(YMMean));
%     figure; plot(YMFit); hold on; plot(YMMean); hold off; 
%     YMMean = YMFit;
    
strain_kernel = 10;    
strain_sensor=strain_calculation2(corrected_sensor,zaxis(2),1,strain_kernel);
strain_left=strain_calculation2(corrected_left,zaxis(2),1,strain_kernel);
strain_right=strain_calculation2(corrected_right,zaxis(2),1,strain_kernel);
strain_oneshot = strain_calculation2(corrected_whole,zaxis(2),1,strain_kernel);
strain_whole = strain_sensor+strain_left+strain_right;
close all force; 
TEST = abs(strain_whole);
TEST(TEST==0) = NaN;
CLow = quantile(TEST,0.25,'all');
if CLow < 0 
    CLow = 0;
end
CHigh = quantile(TEST,0.75,'all');
if isnan(CHigh)
    CHigh = 10;
end
figure;
imagesc(xaxis,zaxis(1:670),abs(strain_whole(1:670,:))); colormap(hot); colorbar; 
caxis([CLow CHigh])
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain')
export_fig([baseFolder,filesep,'StrainMapNew',num2str(imageIterator),'.png'])
save([baseFolder,filesep,'StrainMapNew',num2str(imageIterator),'.fig'])
strainStack(:,:,imageIterator) = abs(strain_whole);

stressSensor = abs(YMMean.*(strain_sensor));
stressSensor(stressSensor==0) = NaN;
stressSensor = mean(stressSensor,1,'omitnan');
figure; plot(stressSensor)

% Try: Uniform stress in L and R sides.

stressSensor2(1:round(length(stressSensor)/2)) = mean(stressSensor(145:210));
stressSensor2(round(length(stressSensor)/2)+1:length(stressSensor)) = mean(stressSensor(400:500));
figure; plot(stressSensor2)
% stressSensor = fliplr(stressSensor)
SS = ones([size(strain_whole,1) size(strain_whole,2)]);
for k = 1:size(SS,2)
    stressSensor3 = stressSensor.*SS(:,k);
end


% YM in Pa, not kPa
YMLeft = abs(stressSensor3)./(abs(strain_left));
YMLeft(isinf(YMLeft)) = 0;

YMRight = abs(stressSensor3)./(abs(strain_right));
YMRight(isinf(YMRight)) = 0;
YMTotal = YMLeft + YMRight;

figure; imagesc(YMTotal/1000); colormap(jet); caxis([0 100]) % in kPa
figure; imagesc(YMTotal)

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

YMLeftAvg3 = mean(YMTotalSmoothed.*leftMaskStack(:,:,imageIterator),'all','omitnan')
YMRightAvg3 = mean(YMTotalSmoothed.*rightMaskStack(:,:,imageIterator),'all','omitnan')

figure; imagesc(YMTotalSmoothed/1000); caxis([0 500])
figure; imagesc(YMTotal/1000); caxis([0 750])

% Convert to kPa
YMTotalSmoothed = YMTotalSmoothed/1000;
close all force; 
figure;
ax1 = axes;

TEST = mat2gray(abs(IQData(1:670,:,imageIterator)));
% TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
TEST = 255*TEST;
imagesc(ax1,xaxis,zaxis(1:670),TEST(1:670,:));
ylabel('Depth (mm)')
xlabel('Distance (mm)')
colormap(ax1,'gray');
caxis([0 25])
ax2 = axes;
imagesc(ax2,YMTotalSmoothed(1:670,:),'alphadata',YMTotalSmoothed(1:670,:)>0);
colormap(ax2,'jet');
TEST2 = YMTotalSmoothed(1:670,:)
TEST2(TEST2==0) = NaN;
CLow = quantile(TEST2,0.25,'all');
CHigh = quantile(TEST2,0.75,'all');
% caxis(ax2,[0 1000]);
caxis(ax2,[CLow CHigh])
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
h = colorbar;
ylabel(h,'YM (kPa)')
title(ax1,"Young's modulus using shear wave")

export_fig([baseFolder,filesep,'YMOverlaySWNew',num2str(imageIterator),'.png'])
save([baseFolder,filesep,'YMOverlaySWNew',num2str(imageIterator),'.fig'])
YMSW(:,:,imageIterator) = YMTotalSmoothed;


close all force; 
figure;
ax1 = axes;

TotalMask = leftMaskStack(:,:,imageIterator) + rightMaskStack(:,:,imageIterator);

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


TEST = mat2gray(abs(IQData(1:670,:,imageIterator)));
% TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
TEST = 255*TEST;
imagesc(ax1,xaxis,zaxis(1:670),TEST(1:670,:));
ylabel('Depth (mm)')
xlabel('Distance (mm)')
colormap(ax1,'gray');
caxis([0 25])
ax2 = axes;
imagesc(ax2,TotalMask(1:670,:),'alphadata',TotalMask(1:670,:))
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


export_fig([baseFolder,filesep,'YMAvgOverlay',num2str(imageIterator),'.png'])
save([baseFolder,filesep,'YMAvgOverlay',num2str(imageIterator),'.fig'])
% sensorInstron = [47.91195
% 44.8882
% 39.04582
% 40.04896
% 39.25353
% 38.7425
% 39.2035
% ];
% sensorInstron = mean(sensorInstron);
% stressSensor = sensorInstron.*strain_sensor;
% stressSensor = mean(stressSensor,1,'omitnan');
% figure; plot(stressSensor)
% % stressSensor = fliplr(stressSensor)
% SS = ones([size(strain_whole,1) size(strain_whole,2)]);
% for k = 1:size(SS,2)
%     stressSensor2 = stressSensor.*SS(:,k);
% end
% 
% 
% 
% YMLeft = abs(stressSensor2./(abs(strain_left)/1000));
% YMLeft(isinf(YMLeft)) = 0;
% 
% YMRight = abs(stressSensor2./(abs(strain_right)/1000));
% YMRight(isinf(YMRight)) = 0;
% YMTotal = YMLeft + YMRight;
% YMTotal = YMTotal / 1000 ; %in kPa
% figure; imagesc(YMTotal); colormap(jet); caxis([0 500])
% figure; imagesc(YMLeft); caxis([0 600]);
% % 
% % colormap(jet)
% % 
% % disp('Select good region')
% % goodLeft = roipoly;
% % 
% % YMLeftAvg = mean(rmoutliers(YMLeft(goodLeft)),'omitnan')
% % YMLeft(YMLeft==0) = NaN;
% % YMLeftAvg2 = mean(rmoutliers(YMLeft),'all','omitnan')
% % 
% % YMRight = abs(stressSensor2./(abs(strain_right)/1000));
% % YMRight(isinf(YMRight)) = 0;
% % figure; imagesc(YMRight); caxis([0 100]);
% % 
% % colormap(jet)
% % disp('Select good region')
% % goodRight = roipoly;
% % 
% % YMRightAvg = mean(rmoutliers(YMRight(goodRight)),'omitnan')
% % YMRight(YMRight==0) = NaN;
% % YMRightAvg2 = mean(rmoutliers(YMRight),'all','omitnan')
% % 
% % YMTotal = YMLeft + YMRight;
% % YMTotal = rmoutliers(YMTotal);
% YMLeftAvg2 = mean(YMLeft,'all','omitnan')
% YMRightAvg2 = mean(YMRight,'all','omitnan')
% 
% x_kern = 20;
% % dis_cumsum = dis_cumsum(:,:,end);
% % dis_cumsum = dis_cumsum';
% % corrected_whole(isnan(corrected_whole)) = 0;
% vec_image = YMTotal';
% BScanWidth = size(vec_image,1);
% sensor_sdl2 = ones([1 size(vec_image,2)]);
% sensor_bdl2 = ones([1 size(vec_image,2)]).*size(vec_image,2);
% how_many_depths = size(vec_image,2);
% fitted_displacement = zeros([size(vec_image,2) size(vec_image,1)]);
% for x = 1:BScanWidth%size(displacement,1)
%     if x < 1
%         x_start = 1;
%     else
%         x_start = x;
%     end
%     
%     if x > BScanWidth - x_kern
%         x_end = BScanWidth;
%     else
%         x_end = x+x_kern;
%     end
%     
%     x_disp = median(vec_image(x_start:x_end,sensor_sdl2(x):sensor_bdl2(1,x)),1);
%     weights = logical(x_disp);
% 
%     % find nans
% %     NaNChecker = isnan(x_disp)
%     [xData, yData] = prepareCurveData( [], x_disp );
%     
%     % Set up fittype and options.
%     ft = fittype( 'cubicspline' );
%     opts = fitoptions( 'Method', 'LinearLeastSquares' );
%     opts.Robust = 'Bisquare';
%     opts.Weights = weights;
%     
%     % Fit model to data.
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x) = feval(fitresult,1:length(x_disp));
% %     figure(1);
% %     subplot(3,1,1);
% %     imagesc(vec_image'); colormap(jet);
% %     hold on;
% %     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
% %     hold off;
% %     colorbar;
% %     %     title(pathname,'Interpreter','none');
% %     
% %     subplot(3,1,2);
% %     %imagesc(displacement'); colormap(jet);
% %     plot(xData, yData, xData, fitted_displacement(sensor_sdl2(1,x):sensor_bdl2(1,x),x));
% %     title(x);
% %     axis tight;
% %     drawnow;
% %     subplot(3,1,3);
% %     
% %     imagesc(fitted_displacement); colormap(jet);
% %     hold on;
% %     plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
% %     hold off;
% %     colorbar;
% % % %         title(pathname,'Interpreter','none');
%     
% end
% 
% YMTotalSmoothed = fitted_displacement;
% 
% YMLeftAvg3 = mean(YMLeft.*leftMaskStack(:,:,imageIterator),'all','omitnan')
% YMRightAvg3 = mean(YMRight.*rightMaskStack(:,:,imageIterator),'all','omitnan')
% 
% figure; imagesc(YMTotalSmoothed); caxis([0 750])
% figure; imagesc(YMTotal); caxis([0 750])
% 
% % 
% % figure;
% % ax1 = axes;
% % TEST = mat2gray(abs(IQData(1:670,:,imageIterator)));
% % % TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
% % TEST = 255*TEST;
% % imagesc(TEST);
% % colormap(ax1,'gray');
% % caxis([0 25])
% % ax2 = axes;
% % imagesc(ax2,YMTotalSmoothed(1:670,:),'alphadata',YMTotalSmoothed(1:670,:)>0);
% % colormap(ax2,'jet');
% % caxis(ax2,[0 500]);
% % ax2.Visible = 'off';
% % linkprop([ax1 ax2],'Position');
% % colorbar;
% 
% close all force; 
% figure;
% ax1 = axes;
% 
% TEST = mat2gray(abs(IQData(1:670,:,imageIterator)));
% % TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
% TEST = 255*TEST;
% imagesc(ax1,xaxis,zaxis(1:670),TEST(1:670,:));
% ylabel('Depth (mm)')
% xlabel('Distance (mm)')
% colormap(ax1,'gray');
% caxis([0 25])
% ax2 = axes;
% imagesc(ax2,YMTotalSmoothed(1:670,:),'alphadata',YMTotalSmoothed(1:670,:)>0);
% colormap(ax2,'jet');
% TEST2 = YMTotalSmoothed(1:670,:)
% TEST2(TEST2==0) = NaN;
% CLow = quantile(TEST2,0.25,'all');
% CHigh = quantile(TEST2,0.75,'all');
% % caxis(ax2,[0 1000]);
% caxis(ax2,[CLow CHigh])
% ax2.Visible = 'off';
% linkprop([ax1 ax2],'Position');
% h = colorbar;
% ylabel(h,'YM (kPa)')
% title(ax1,"Young's modulus using sensor instron value")
% 
% export_fig([baseFolder,filesep,'YMOverlayInstron',num2str(imageIterator),'.png'])
% save([baseFolder,filesep,'YMOverlayInstron',num2str(imageIterator),'.fig'])
% YMI(:,:,imageIterator) = YMTotalSmoothed;
end