load([baseFolder,filesep,'stressStrainandYM.mat'],'YMSensorStack','strainStack','sensorStressStack',...
    'YMStack','YMLeftStack','YMRightStack')

close all force; 
for k = 1:size(YMSensorStack,3)
    figure(1)
    h1 = scatter(k,YMSensorStack(:,100,k),'r')
    hold on;
    h2 = scatter(k,YMSensorStack(:,400,k),'b')
    
end
title('YM in sensor vs compression #')
xlabel('Compression #')
ylabel('YM (kPa')
b = [h1 h2]
legend(b,'YM Left side','YM Right side')

for k = 1:size(sensorStressStack,3)
    figure(1);
    plot(sensorStressStack(:,:,k))
    pause(0.1)
    hold on; 
end
    

xaxis = 1:40;
xaxis = xaxis.*10 % axis in um

close all force; 
for k = 2:40%size(YMLeftStack,3)
    figure(1)
    h1 = scatter(k,YMLeftStack(:,:,k),'r')
    hold on;
    h2 = scatter(k,YMRightStack(:,:,k),'b')
 
end
title('YM vs compression #')
xlabel('Compression #')
ylabel('YM (kPa)')
b = [h1 h2]
legend(b,'YM Left side','YM Right side')


close all force; 
for k = 2:40%size(YMLeftStack,3)
    figure(1)
    h1 = scatter(k,YMLeftStack(:,:,k),'r')
    hold on;
%     h2 = scatter(k,YMRightStack(:,:,k),'b')
 
end
title('YM vs compression #, Stiff side')
xlabel('Compression #')
ylabel('YM (kPa)')

close all force; 
for k = 2:40%size(YMLeftStack,3)
    figure(1)
%     h1 = scatter(k,YMLeftStack(:,:,k),'r')
    
    h2 = scatter(k,YMRightStack(:,:,k),'b')
 hold on;
end
title('YM vs compression #')
xlabel('Compression #')
ylabel('YM (kPa)')

% Same graphs but with um in x axis
% 
% close all force; 
% for k = 2:40%size(YMLeftStack,3)
%     figure(1)
%     h1 = scatter(xaxis(k),YMLeftStack(:,:,k),'r','filled')
%     hold on;
%     h2 = scatter(xaxis(k),YMRightStack(:,:,k),'b','filled')
%  
% end
% title('YM vs compression')
% xlabel('Distance compressed (um)')
% xlim([0 405])
% ylim([0 50])
% ylabel('YM (kPa)')
% b = [h1 h2]
% legend(b,'YM Left side','YM Right side')


close all force; 

h1 = scatter(xaxis(2:40),squeeze(YMLeftStack(:,:,2:40)),'r','filled');
hold on;
h2 = scatter(xaxis(2:40),squeeze(YMRightStack(:,:,2:40)),'b','filled');
h = lsline;
title('YM vs compression')
xlabel('Distance compressed (um)')
xlim([0 405])
ylim([0 50])
ylabel('YM (kPa)')
b = [h1 h2]
legend(b,'YM Left side','YM Right side')

close all force; 
% for k = 2:40%size(YMLeftStack,3)
%     figure(1)
%     h1 = scatter(xaxis(k),YMLeftStack(:,:,k),'r','filled')
%     hold on;
% %     h2 = scatter(k,YMRightStack(:,:,k),'b')
%  
% end
scatter(xaxis(2:40),squeeze(YMLeftStack(:,:,2:40)),'r','filled')
h = lsline;
h.Color = 'r';
title('YM vs compression, stiff side')
xlabel('Distance compressed (um)')
xlim([0 405])
ylim([0 50])
ylabel('YM (kPa)')
% 
close all force; 
% for k = 2:40%size(YMLeftStack,3)
%     figure(1)
% %     h1 = scatter(k,YMLeftStack(:,:,k),'r')
%     
%     h2 = scatter(xaxis(k),YMRightStack(:,:,k),'b','filled')
%  hold on;
% end
scatter(xaxis(2:40),squeeze(YMRightStack(:,:,2:40)),'b','filled');
h = lsline;
h.Color = 'b';
title('YM vs compression, soft side')
xlabel('Distance compressed (um)')
xlim([0 405])
ylim([0 10])
ylabel('YM (kPa)')