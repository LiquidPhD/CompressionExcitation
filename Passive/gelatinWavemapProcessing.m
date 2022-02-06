clearvars;
% folders1 = rdir(['G:\GelatinPhantoms092420\USE\Used\*\**\*Param*'])
% folders2 = rdir(['D:\Ultrasound112720\*\**\*Param*'])
% folders3 = rdir(['D:\UltrasoundData\SKIN_ELASTOGRAPHY\*\**\*Param*'])
% folders4 = rdir(['D:\RippySkin\*\**\*Param*'])
% folders = [folders1; folders2; folders3; folders4];
% folders = rdir(['D:\RippySkin\*\**\*Param*'])
folders = rdir(['D:\GelatinPhantoms20210122\*\**\*Param*'])
progressbar(0)
for folderIndex = 10:length(folders)
    try
clearvars -except folders folderIndex
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
continueFlag = 0;
locationFlag = 0; 
locationNumber = 1; 
while locationFlag == 0
while continueFlag == 0
close all force; 
playWaveVideo(displacement,0.05,jet)
figure; imagesc(displacement(:,:,5));
[x,y] = ginput(1);
close all force;
x = round(x);
y = round(y);
[frequency,FT_plot] = plot_spectrum(squeeze(displacement(y,x,:)),10000)

answer = questdlg('Good?', ...
	'CHOOSE MORTAL', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        continueFlag = 1
    case 'No'

end
end


answer = inputdlg('NumberOfPeaks')
close all force;
for peakSelection = 1:str2num(cell2mat(answer))
[frequency,FT_plot,P1] = plot_spectrum(squeeze(displacement(y,x,:)),10000)

[selectedFreq,~] = ginput(1);
close all force;
figure;
subplot(2,1,1)
imagesc(displacement(:,:,5))
hold on;
scatter(x,y,'r*'); hold off;
subplot(2,1,2);
plot(frequency,P1)
hold on;
vline(selectedFreq); hold off; 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
export_fig([folders(folderIndex).folder,filesep,'Location',num2str(locationNumber),'AndFFT',num2str(peakSelection),'.png'],'-png','-native')
savefig([folders(folderIndex).folder,filesep,'Location',num2str(locationNumber),'AndFFT',num2str(peakSelection),'.fig'])
close all force; 
figure; imagesc(xaxis,zaxis,abs(wavelength)*selectedFreq)
colormap(jet)
% caxis([1e-5 1e-4])
if max(max(abs(wavelength)*selectedFreq)) > 20
    caxis([0 20])
else
    caxis([0 max(max(abs(wavelength)*selectedFreq))])
end
% caxis([0 20])
colorbar;
export_fig([folders(folderIndex).folder,filesep,'WavelengthMapLocation',num2str(locationNumber),'Peak',num2str(peakSelection),'.png'],'-png','-native')
savefig([folders(folderIndex).folder,filesep,'WavelengthMapLocation',num2str(locationNumber),'Peak',num2str(peakSelection),'.fig'])
% saveas(gcf,[folders(folderIndex).folder,filesep,'WavelengthMap.fig']);
% saveas(gcf,[folders(folderIndex).folder,filesep,'WavelengthMap.png']);
end
answer2 = questdlg('Another Location?', ...
	'CHOOSE MORTAL', ...
	'Yes','No','No');
% Handle response
switch answer2
    case 'Yes'
    locationNumber = locationNumber+1;
    continueFlag = 0; 
    case 'No'
    locationFlag = 1;
end

end
progressbar(folderIndex/length(folders))
    save([folders(folderIndex).folder,filesep,'WavelengthMapWorkspace.mat'])

    catch
        progressbar(folderIndex/length(folders))

    end 
end