figure; 
h1 = plot(STRAINm*100,yprime/1000); hold on; h2 = scatter(STRAINm(1:length(YMLeftStack))*100,YMLeftStack,'filled')
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM Instron vs computed values, stiff side')
xlim([0 4])
ylim([0 50])
b = [h1 h2];
legend(b,'Instron','Computed','Location','SE')




figure; 
h1 = plot(STRAINm*100,yprime/1000); hold on; h2 = plot(STRAINm(1:length(YMRightStack))*100,YMRightStack)
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM Instron vs computed values, soft side')
xlim([0 4])
b = [h1 h2];
legend(b,'Instron','Computed','Location','SE')

DTNA = squeeze(cumsum(diff(thicknessNotAvg(:,80,:)))/8.02);
figure; scatter(1:length(DTNA),DTNA)

figure; 
h1 = plot(STRAINm*100,yprime/1000); hold on; h2 = scatter(abs(DTNA(3:length(YMLeftStack)))*100,YMLeftStack(3:end),'filled')
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM Instron vs computed values, stiff side')
xlim([0 4])
ylim([0 50])
b = [h1 h2];
legend(b,'Instron','Computed','Location','SE')



DTNA = squeeze(cumsum(diff(thicknessNotAvg(:,430,:)))/8.58);
figure; scatter(1:length(DTNA),DTNA)


figure; 
h1 = plot(STRAINm*100,yprime/1000); hold on; h2 = scatter(abs(DTNA(1:length(YMRightStack)))*100,YMRightStack,'filled')
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM Instron vs computed values, soft side')
xlim([0 4])
b = [h1 h2];
legend(b,'Instron','Computed','Location','SE')
ylim([0 25])
