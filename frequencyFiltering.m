function [frequencyFilteredPV] = frequencyFiltering(particleVelocity,BScan,Parameters,xaxis)


displacement = particleVelocity(:,:,20:110);
delta_t_interp = Parameters.delta_t;
assignin('base','displacement',displacement)
mouseoverFFT(BScan(:,:,1),delta_t_interp)
uiwait
finalPosition = evalin('base','finalPosition')
waveform_left = squeeze(particleVelocity(round(finalPosition(1,2)),round(finalPosition(1,1)),:))

Ts = delta_t_interp;
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
%         [selectedFreq,~] = ginput(1)
%         selectedFreq = round(selectedFreq);
h = imrect(gca,[1e3-200,-200,400,(max(abs(FFT))*1.1)+400]);
position = wait(h);
b = fir1(100,[position(1)*2*Ts (position(1)+position(3))*2*Ts]);
delay = mean(grpdelay(b));


frequencyFilteredPV=zeros(size(displacement,1),size(displacement,2),size(displacement,3));
for pos=1:size(displacement,2)
    clc; disp([num2str(round(pos/length(xaxis)*100)),'%']);
    
    for depth=1:size(displacement,1)
        
        signal=squeeze(displacement(depth,pos,:));
        if isnan(signal(1))
            signal = zeros(length(signal),1);
        end
        
        signal_tmp = filter(b,1,signal);
        env = abs(hilbert(signal_tmp));
        signal_tmp = signal_tmp./env;
        
        frequencyFilteredPV(depth,pos,:) = signal_tmp;
    end
    
    
end

% playWaveVideo(filtered_motion)