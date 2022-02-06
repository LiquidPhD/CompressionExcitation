folders = rdir('C:\Users\justi\Desktop\Passive\*\**\Parameter.mat') 
for k= 1:length(folders)
processUSEForWorkspace(k,folders)
end

for folderIndex = 1:length(folders)
load([folders(folderIndex).folder,filesep,'workspace.mat'],'displacement','IQData','delta_z')
tic
for k = 1:size(displacement,3)
    displacement_smoothed(:,:,k) = modefilt(squeeze(displacement(:,:,k)),[33 33]);
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

[strain_images] = strain_calculation(displacement,delta_z,1,60);
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
figure; imagesc(abs(wavelength))
colormap(jet)
caxis([1e-5 1e-4])
pause
end