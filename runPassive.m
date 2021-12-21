function [wavelength,selectedFreq,frequencyMap] = runPassive(particleVelocity,BScan,Parameters,folderIndex)

displacement = particleVelocity(:,:,20:end);
win_size = 10;
sdl = ones(1,size(displacement,2));
sdl = sdl.*50;
bdl = ones(1,size(displacement,2));
bdl = bdl.*size(displacement,1);
[tr_disp] = time_reversal_FFT(displacement,sdl,bdl,win_size);
% cumsumDisplacement = cumsum(displacement,3);
[strain_images] = strain_calculation(displacement,Parameters.delta_z,1,5);
% [strain_images] = strain_3D_robust_no_sensor2(delta_x,clippedBScan(:,:,1),...
%     displacement,10)
% strain_images = strain_images';
[tr_strain] = time_reversal_FFT(strain_images,sdl,bdl,win_size);

wavelength = 2*pi*sqrt(tr_disp./tr_strain);
close all force;
if folderIndex == 1
    try
        notify('Select frequency!')
    catch
    end
    %     figure; imagesc(displacement(:,:,3));
    %     [x,y] = ginput(1);
%     try
%     finalPosition = evalin('base','finalPosition')
%     catch
    disp_z_smooth_desamp_interp = displacement;
    delta_t_interp = Parameters.delta_t;
    mouseoverFFT(BScan,delta_t_interp)
    uiwait
%     end
end
finalPosition = evalin('base','finalPosition')
waveform_left = squeeze(particleVelocity(round(finalPosition(1,2)),round(finalPosition(1,1)),:))
maskFlag = 0;
% Fernando filtering
if folderIndex == 1
    Ts = Parameters.delta_t;
    Line = waveform_left;
    FFT = fft(Line,2^12);
    freq = linspace(0,1,2^12)/Ts;
    
    Time = [0:Ts:(length(waveform_left)-1)*Ts]*1e3;
    
    fig = figure;
    plot(freq,abs(FFT))
    grid on
    ylabel('Magnitude (Arb.)');
    xlabel('Frequency (Hz)');
    title('FFT of the signal');
    axis([0 4000 0 max(abs(FFT))*1.1])
    [selectedFreq,~] = ginput(1)
    selectedFreq = round(selectedFreq);
end
mask = ones([size(particleVelocity,1) size(particleVelocity,2)]);
[frequencyMap] = freqMap(displacement,delta_t_interp,mask,maskFlag);
figure;
[Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
imagesc(xaxis,zaxis,abs(wavelength).*selectedFreq)

colormap(jet);
colorbar;
% caxis([1e-5 1e-4])
caxis([0 20])