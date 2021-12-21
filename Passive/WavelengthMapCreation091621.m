clearvars;
% folders1 = rdir(['G:\GelatinPhantoms092420\USE\Used\*\**\*Param*'])
% folders2 = rdir(['D:\Ultrasound112720\*\**\*Param*'])
% folders3 = rdir(['D:\UltrasoundData\SKIN_ELASTOGRAPHY\*\**\*Param*'])
% folders4 = rdir(['D:\RippySkin\*\**\*Param*'])
% folders = [folders1; folders2; folders3; folders4];
% folders = rdir(['D:\RippySkin\*\**\*Param*'])
% folders = rdir(['D:\GelatinPhantoms20210122\*\**\*Param*'])
folders = rdir(['D:\SiliconeHalfAndHalfSiliconeSensor2\Static\*\**\*Param*'])
progressbar(0)
for folderIndex = 1:length(folders)
    try
clearvars -except folders folderIndex
processUSEForWorkspaceWavelengthUse(folderIndex,folders)
load([folders(folderIndex).folder,filesep,'wavelengthWorkspace.mat'],'displacement','IQData','delta_z','xaxis','zaxis','taxis')
tic
progressbar(0);
for k = 1:size(displacement,3)
    displacement_smoothed(:,:,k) = modefilt(squeeze(displacement(:,:,k)),[33 33]);
    progressbar(k/size(displacement,3));
end
toc
win_size = 3;
% playWaveVideo(displacement_smoothed)
displacement = displacement_smoothed(:,:,:);
figure; imagesc(squeeze(abs(IQData(:,:,1))))
sdl = ones(1,size(displacement,2));
sdl = sdl.*50;
bdl = ones(1,size(displacement,2));
bdl = bdl.*size(displacement,1);
[tr_disp] = time_reversal(displacement,sdl,bdl,win_size);
% figure; imagesc(tr_disp);

[strain_images] = strain_calculation(displacement,delta_z,1,100);
figure; imagesc(strain_images); colormap(fireice)
% playWaveVideo(strain_images)

% for k = 1:size(displacement,3)
%     figure(1)
%     imagesc(squeeze(displacement(:,:,k)));
%     colormap(jet)
%     figure(2)
%     imagesc(squeeze(strain_images(:,:,k)));
%     colormap(jet)
%     autoArrangeFigures
%     pause(0.1)
% end
[tr_strain] = time_reversal(strain_images,sdl,bdl,win_size);
% figure; imagesc(tr_strain);

wavelength = 2*pi*sqrt(tr_disp./tr_strain);
pause
close all force; 
[frequency,FT_plot] = plot_spectrum(squeeze(displacement(200,200,:)),10000)
[selectedFreq,~] = ginput(1);
close all force; 
figure; imagesc(xaxis,zaxis,abs(wavelength)*selectedFreq/1000)
colormap(jet)
% caxis([1e-5 1e-4])
caxis([0 20])

export_fig([folders(folderIndex).folder,filesep,'WavelengthMap.png'],'-png','-native')
% saveas(gcf,[folders(folderIndex).folder,filesep,'WavelengthMap.fig']);
% saveas(gcf,[folders(folderIndex).folder,filesep,'WavelengthMap.png']);
progressbar(folderIndex/length(folders))
    save([folders(folderIndex).folder,filesep,'WavelengthMapWorkspace.mat'])

    catch
        progressbar(folderIndex/length(folders))

    end 
end