function mouseoverFFT(BScan,delta_t)
% mouseoverFFT Shows FFT when you mouse over an image and allows you to
% select a particular point.
%
% Output is called "finalPosition" and is directly assigned to base
% workspace. 
%
% You must have variables "displacement" and "delta_t" in your workspace
% prior to calling this function. 
%
% Justin Rippy
% 06/28/21


% Generate unique figure showing BScan
fig = figure(555);
fig.Units = 'normalized';
fig.OuterPosition = [0 0 0.5 1];
% figure('Units','normalized','OuterPosition',[0 0 0.5 1])
imagesc(mat2gray(BScan(:,:,1)));
caxis([0 0.1]); colormap(gray);
ustring = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','clicked object')

% This sets the properties for the mouseover and click.
set(fig,'WindowButtonmotionFcn',{@position,ustring});
set(fig,'WindowButtonDownFcn',{@closeWindowAndReturn})

% Generate unique figure showing FFT
fig2 = figure(777);
fig2.Units = 'normalized';
fig2.OuterPosition = [0.5 0 0.5 1];
waveform = squeeze(BScan(100,100,:));
Ts = delta_t;
Line = waveform;
FFT = fft(Line,2^12);
freq = linspace(0,1,2^12)/Ts;
Time = [0:Ts:(length(waveform)-1)*Ts]*1e3;

plot(freq,abs(FFT))
grid on
ylabel('Magnitude (Arb.)');
xlabel('Frequency (Hz)');
title('FFT of the signal');
axis([0 4000 0 max(abs(FFT))*1.1])
% autoArrangeFigures;
evalin('base','firstClick = 0');

fig3 = figure(666);
fig3.Units = 'normalized';
fig3.OuterPosition = [1 0 0.5 1];
end

function position(hobject,event,ustring)

% Evaluate whether or not the figure has been initially clicked or not
firstClick = evalin('base','firstClick');

% Click inside the plot to select it
if firstClick == 1
    figure(555)
    axesHandle = gca;
    pt = get(axesHandle, 'CurrentPoint');
    % Show coordinates
    set(ustring,'String',num2str(pt));
    
    displacement = evalin('base','particleVelocity;');
    delta_t = evalin('base','Parameters.delta_t;');
    
    waveform = squeeze(displacement(round(pt(1,2)),round(pt(1,1)),:));
    
    
    Ts = delta_t;
    Line = waveform;
    FFT = fft(Line,2^12);
    freq = linspace(0,1,2^12)/Ts;
    
    Time = [0:Ts:(length(waveform)-1)*Ts]*1e3;
    figure(777)
    plot(freq,abs(FFT))
    grid on
    ylabel('Magnitude (Arb.)');
    xlabel('Frequency (Hz)');
    title('FFT of the signal');
    axis([0 4000 0 max(abs(FFT))*1.1])
    fig3 = figure(666);
fig3.Units = 'normalized';
fig3.OuterPosition = [1 0 0.5 1];
plot(waveform)
end

end

function closeWindowAndReturn(hobject,event)
firstClick = evalin('base','firstClick');

if firstClick == 0
    evalin('base','firstClick=1');
else
    figure(555)
    axesHandle = gca;
    pt = get(axesHandle, 'CurrentPoint');
    close all force;
    pt = round(pt(1,1:2));
    assignin('base','finalPosition',pt)
    
    
end
end