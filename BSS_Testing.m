
for BSS_Window_size = 10:10:100
b_scan_strain = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    dis_cumsum(:,:,end),BSS_Window_size);

VELOCITIES = centerVelocity(1:size(b_scan_strain,2),:,:);
BSS = b_scan_strain';
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
        speedAtPixel = VELOCITIES(columnIndex,rowIndex);
        strainAtPixel = BSS(columnIndex,rowIndex);
        sigma(columnIndex,rowIndex) = 3*1000*(speedAtPixel^2)*strainAtPixel;
    end
end

% Try stress at 330
for depthValue = 10:10:300
stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
figure; subplot(2,1,1)
imagesc(BScan); hold on; hline(depthValue); hold off; 
subplot(2,1,2);
plot(stressMM); title(num2str(depthValue))
end

depthValue = 200
stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
for columnIndex = 1:size(VELOCITIES,1)
    for rowIndex = 1:size(VELOCITIES,2)
%         stressAtPixel = sigma(columnIndex,rowIndex);
%         stressAtPixel = axialStress(:,rowIndex);
% stressAtPixel = mean(sigma(interface,rowIndex));
stressAtPixel = stressMM(:,rowIndex);
        strainAtPixel = BSS(columnIndex,rowIndex);
        YM(columnIndex,rowIndex) = (abs(stressAtPixel)/abs(strainAtPixel))/1000; % in kPa
    end
end

close all force;
% figure;
% imagesc(sigma); 
YM = filloutliers(YM,'makima');
fig = figure; imagesc(xaxis,zaxis(1:y),YM);
colormap(jet)
caxis([0 100])
saveas(fig,['I:\SiliconeHalfAndHalfSiliconeSensor2\BSS_Testing',filesep,'BSS_Window',num2str(BSS_Window_size),'.fig'])
end

figures = rdir('I:\SiliconeHalfAndHalfSiliconeSensor2\BSS_Testing\*.fig*')

for k = 1:length(figures)
    open(figures(k).name)
    a = get(gca,'Children')
    dataFromPlots(:,:,k) = a.CData;
end

for k = 1:length(figures)
    w(:,k) = (2*(length(figures) - k + 1)) / (length(figures) * (length(figures) + 1));
end

for k = 1:size(dataFromPlots,3)
    weightedImage(:,:,k) = dataFromPlots(:,:,k).*w(:,k);
end

weightedImage = sum(weightedImage,3);
figure;
imagesc(abs(weightedImage))
caxis([0 50])

figure; imagesc(xaxis,zaxis(1:y),abs(filloutliers(weightedImage,'linear')));
caxis([0 100])
colormap(jet)
