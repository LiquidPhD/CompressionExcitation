figure; imagesc(squeeze(displacement(:,:,1)));
% 
% for k = 1:150
% showWavePlot(displacement,200,200-k)
% pause(0.1)
% end
[freq,FT_plot] = plot_spectrum(squeeze(displacement(200,200,:)),10000)

figure; 
imagesc(xaxis,zaxis,abs(wavelength*683)); colormap(jet)
% displacement = displacement*1000;
displacementCS = cumsum(displacement,3);
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

wavelength = 2*pi*sqrt(tr_disp./(tr_strain./1000));
[frequency,FT_plot] = plot_spectrum(squeeze(displacement(200,200,:)),10000)
close all force; 
figure; imagesc(xaxis,zaxis,abs(wavelength)*256)
colormap(jet)
caxis([1e-5 1e-4])