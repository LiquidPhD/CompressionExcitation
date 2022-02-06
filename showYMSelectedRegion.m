function showYMSelectedRegion(I,YMCompression)
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
fig = figure(999);
fig.Units = 'normalized';
fig.OuterPosition = [0 0 0.5 1];


imagesc(I)


% YMCompression = YMCompressionOG;

    h = impoly(gca,[3.77764976958525,418.068513119534;267.300691244240,456.202623906706;275.167050691244,1134.98979591837;5.08870967741935,1131.17638483965;3.77764976958525,418.068513119534])
%     PosTime = wait(h);
    points = h.getPosition
    Block = h.createMask;
    assignin('base','Block',Block)
    assignin('base','points',points)
    set(fig,'WindowButtonDownFcn',@calculateYM);
set(fig,'CloseRequestFcn',@my_closereq)
assignin('base','PosTime',PosTime)
end 

function calculateYM(hobject,event,YMCompression,Block)
YM = evalin('base','YMCompression');
Block = evalin('base','Block')
YM(Block~=1) = NaN;
% figure; imagesc(leftYM);
% figure; imagesc(rightYM)

YMNegNotRemoved = mean(rmoutliers(abs(YM(Block))),'omitnan');

YMLNegRem = YM;
YMLNegRem(YMLNegRem<0) = NaN;

YMNegRem = mean(rmoutliers(abs(YMLNegRem(Block))),'omitnan');
title({num2str(YMNegNotRemoved);[num2str(YMNegRem)]})
drawnow;
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

 h = impoly(gca)
%     PosTime = wait(h);
    points = h.getPosition
    Block = h.createMask;
    assignin('base','Block',Block)
    assignin('base','points',points)
end


function my_closereq(src,event)
% Close request function 
% to display a question dialog box 
selection = questdlg('Close This Figure?', ...
    'Close Request Function', ...
    'Yes','No','Yes'); 
switch selection 
    case 'Yes'
        delete(gcf)
    case 'No'
        return 
end
end
% 
% function closeWindowAndReturn(hobject,event)
% firstClick = evalin('base','firstClick');
% 
% if firstClick == 0
%     evalin('base','firstClick=1');
% else
%     figure(555)
%     axesHandle = gca;
%     pt = get(axesHandle, 'CurrentPoint');
%     close all force;
%     pt = round(pt(1,1:2));
%     assignin('base','finalPosition',pt)
%     
%     
% end
% end