figure;
scatter(1:40,squeeze(thicknessNotAvg(:,81,:)))
title('Thickness change with compression #')
ylabel('Thickness (mm)')
xlabel('Compression #')
title('Thickness change with compression #, stiff side')
xlim([0 26])

labelpoints(1,8.0152,'8.02','position','E','buffer',0.3)
labelpoints(25,7.81257,'7.81','position','N','buffer',0.1)

figure;
scatter(1:40,squeeze(thicknessNotAvg(:,430,:)))
title('Thickness change with compression #')
ylabel('Thickness (mm)')
xlabel('Compression #')
title('Thickness change with compression #, soft side')
xlim([0 26])

labelpoints(1,8.57807,'8.58','position','E','buffer',0.3)
labelpoints(25,8.28538,'8.29','position','N','buffer',0.1)

DTNA = squeeze(cumsum(diff(thicknessNotAvg(:,400,:)))/8.58);
figure; scatter(1:length(DTNA),DTNA)


figure; 
h1 = plot(STRAINm*100,yprime/1000); hold on; h2 = plot(abs(DTNA(1:length(YMRightStack)))*100,YMRightStack)
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM Instron vs computed values, soft side')
xlim([0 4])
b = [h1 h2];
legend(b,'Instron','Computed','Location','SE')


figure
for k = 284:size(thicknessNotAvg,2)
    DTNA = squeeze(cumsum(diff(thicknessNotAvg(:,k,:)))/8.58);
h1 = plot(STRAINm(4:end)*100,yprime(4:end)/1000); hold on; h2 = scatter(abs(DTNA(1:length(YMRightStack)))*100,YMRightStack); hold off; 
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title(num2str(k))
xlim([0 4])
b = [h1 h2];
legend(b,'Instron','Computed','Location','SE')
pause
end