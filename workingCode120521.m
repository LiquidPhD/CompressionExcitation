% baseFolder = 'C:\Users\justi\Desktop\CompressionAndExcitation\HalfHalfGelatin112821'
clearvars;
close all force;
baseFolder = 'E:\ThinAgar121221';
%
% [lowerBound] = selectLowerBound(folders)
startDepth = 1;
lowerBound = 1395

%% Process static data

% Load static data
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
% if ~isempty(rdir([baseFolder,filesep,'combinedIQData.mat']))
%     load([baseFolder,filesep,'combinedIQData.mat'])
%     folderIndex = 1;    
% else
    [IQData,Parameters] = loadStaticData(folders,lowerBound,length(folders));
folderIndex = 1;
save([baseFolder,filesep,'combinedIQData.mat'],'IQData','Parameters');
% end
[vec_phase_diff] = VPD(IQData);
[particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters);
% particleVelocity = permute(vec_phase_diff,[2 1 3]);
[Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
Nz = 1395;
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);
zaxis = zaxis(startDepth:end);
images = mat2gray(abs(IQData));

volumeSegmenter



%% Run all the rest above before running below separately
for imageIterator = 2:20

vec_image=abs(vec_phase_diff(:,:,imageIterator));




x_kern = 20;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
vec_image = vec_image';
BScanWidth = size(vec_image,1);
sensor_sdl = ones([1 size(vec_image,2)])
bdl = ones([1 size(vec_image,2)]).*size(vec_image,2);
how_many_depths = size(vec_image,2);

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
    
    x_disp = median(vec_image(x_start:x_end,sensor_sdl(x):bdl(1,x)),1);
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_sdl(1,x):bdl(1,x),x) = feval(fitresult,1:length(x_disp));
    figure(1);
    subplot(3,1,1);
    imagesc(vec_image'); colormap(jet);
    hold on;
    plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
    hold off;
    colorbar;
%     title(pathname,'Interpreter','none');
    
    subplot(3,1,2);
    %imagesc(displacement'); colormap(jet);
    plot(xData, yData, xData, fitted_displacement(sensor_sdl(1,x):bdl(1,x),x));
    title(x);
    axis tight;
    drawnow;   
    subplot(3,1,3);

    imagesc(fitted_displacement); colormap(jet);
    hold on;
    plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
    hold off;
    colorbar;
%     title(pathname,'Interpreter','none');
    
end

vec_image = fitted_displacement;


imagesc(vec_image); colormap(jet);

images=20*log10(abs(IQData));

LBL = modefilt(labels(:,:,imageIterator),[51 51]);
for pos=1:size(LBL,2)
    sensor_bdl(pos)=find(LBL(:,pos)==1,1,'last');
end


for pos=1:size(LBL,2)
    sensor_sdl(pos)=find(LBL(:,pos)==1,1)+20;
end

mask2=zeros(size(LBL(:,:)));
for pos=1:size(LBL,2);
    mask2(sensor_bdl(pos):end,pos)=1;
end

mask=zeros(size(LBL(:,:)));
for pos=1:size(LBL,2)
    mask(sensor_sdl(pos):sensor_bdl(pos),pos)=1;
end

% 
% figure; imagesc(unwrapped_sample); 
% figure; imagesc(abs(unwrapped_sample));
% figure; imagesc(images(:,:,2)); colormap(gray);

vec_image_sensor=vec_image.*mask;
vec_image_sample=vec_image.*mask2;

unwrapped_sensor=phase_unwrap(vec_image_sensor);

for pos=1:size(vec_image_sensor,2)
    surf_dis(pos)=nanmedian(unwrapped_sensor(sensor_sdl(pos)+5:sensor_sdl(pos)+15,pos));
    corr_sensor(:,pos)=(unwrapped_sensor(:,pos)-surf_dis(pos)).*mask(:,pos);
end


imagesc(corr_sensor);

unwrapped_sample=phase_unwrap(vec_image_sample);

for pos=1:size(vec_image_sample,2)
    surf_dis(pos)=nanmedian(unwrapped_sample(sensor_bdl(pos):sensor_bdl(pos)+10,pos));
    corr_sample(:,pos)=(unwrapped_sample(:,pos)-surf_dis(pos)).*mask2(:,pos);
end

strain_sensor=strain_calculation2(corr_sensor,zaxis(2),1,100);
figure;
imagesc(abs(strain_sensor)/1000); colormap(hot); colorbar; caxis([0,0.5]);

strain_sample=strain_calculation2(corr_sample,zaxis(2),1,100);
figure;
imagesc(abs(strain_sample)/1000); colormap(hot); colorbar; caxis([0,0.5]);

strain_whole=strain_sensor.*mask+strain_sample.*mask2;
close all force;
figure;
imagesc(xaxis,zaxis,abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,0.5]);
export_fig([baseFolder,filesep,'AbsSmoothedStrain',num2str(imageIterator),'.png'])
strainImages(:,:,imageIterator) = strain_whole;
end
save([baseFolder,filesep,'StrainImages.mat'],'strainImages','labels','xaxis','zaxis')
notify('Done.')