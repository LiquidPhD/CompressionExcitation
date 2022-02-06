function [YMNegNotRemoved,YMNegRem,Block,points] = selectYMRegion(I,YM)
% Generate unique figure showing BScan
YMOG = YM;
fig = figure; 
fig.WindowState = 'maximized'; 

imagesc(I)


% YMCompression = YMCompressionOG;

    h = impoly(gca)
    PosTime = wait(h);
    points = h.getPosition
    Block = h.createMask;

finishFlag = 0; 
while finishFlag == 0
YM(Block~=1) = NaN;
% figure; imagesc(leftYM);
% figure; imagesc(rightYM)

YMNegNotRemoved = mean(rmoutliers(abs(YM(Block))),'omitnan');

YMLNegRem = YM;
YMLNegRem(YMLNegRem<0) = NaN;

YMNegRem = mean(rmoutliers(abs(YMLNegRem(Block))),'omitnan');

% 
% % Evaluate whether or not the figure has been initially clicked or not
% firstClick = evalin('base','firstClick');
% 
% % Click inside the plot to select it
%     figure(555)
%     axesHandle = gca;
%     pt = get(axesHandle, 'CurrentPoint');
%     % Show coordinates
%     set(ustring,'String',num2str(pt));
%     
%     displacement = evalin('base','particleVelocity;');
%     delta_t = evalin('base','Parameters.delta_t;');
%     
%     waveform = squeeze(displacement(round(pt(1,2)),round(pt(1,1)),:));
%     
%     
%     Ts = delta_t;
%     Line = waveform;
%     FFT = fft(Line,2^12);
%     freq = linspace(0,1,2^12)/Ts;
%     
%     Time = [0:Ts:(length(waveform)-1)*Ts]*1e3;
%     figure(777)
%     plot(freq,abs(FFT))
%     grid on
%     ylabel('Magnitude (Arb.)');
%     xlabel('Frequency (Hz)');
%     title('FFT of the signal');
%     axis([0 4000 0 max(abs(FFT))*1.1])
%     fig3 = figure(666);
% fig3.Units = 'normalized';
% fig3.OuterPosition = [1 0 0.5 1];
% plot(waveform)
clf(fig) 
fig = figure(1);
fig.Units = 'normalized';
fig.OuterPosition = [0 0 0.5 1];
fig.CloseRequestFcn = 'closereq'

imagesc(I)

title({num2str(YMNegNotRemoved);[num2str(YMNegRem)]})
drawnow;
% YMCompression = YMCompressionOG;
answer = questdlg('Good?', ...
	'CHOOSE', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        finishFlag = 1;
    case 'No'
        YM = YMOG;
        h = impoly(gca,PosTime)
    PosTime = wait(h);
    points = h.getPosition
    Block = h.createMask;
    assignin('base','Block',Block)
    assignin('base','points',points)
end

 
    
    
end

