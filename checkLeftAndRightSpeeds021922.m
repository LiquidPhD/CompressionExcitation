figure; 
h1 = scatter(1:length(YMFromSpeedLeft),YMFromSpeedLeft/1000); hold on; 
h2 = scatter(1:length(YMFromSpeedLeft),YMFromSpeedRight/1000); 
h3 = scatter(1:length(YMLeftStack),YMLeftStack); 
h4 = scatter(1:length(YMLeftStack),YMRightStack); 
hold off;
b = [h1 h2 h3 h4];
legend(b,'YM from speed left','YM from speed right','YM calc left','YM calc right')

figure; 
h1 = scatter(1:length(YMFromSpeedLeft),YMFromSpeedLeft/1000); hold on; 
h3 = scatter(1:length(YMLeftStack),YMLeftStack); 
hold off;
b = [h1  h3 ];
legend(b,'YM from speed left','YM calc left')


figure; 
h2 = scatter(1:length(YMFromSpeedLeft),YMFromSpeedRight/1000); hold on; 
h4 = scatter(1:length(YMLeftStack),YMRightStack); 
hold off;
b = [ h2 h4 ];
legend(b,'YM from speed right','YM calc right')

SpeedX = 1:length(YMFromSpeedLeft);
CalcX = 1:length(YMLeftStack);

figure;
plot((leftSpeedAvg).^2.*3);
hold on;
plot(rightSpeedAvg); hold off; 