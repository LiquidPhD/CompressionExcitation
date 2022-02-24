clearvars thicknessNotAvgL thicknessNotAvgR thicknessNotAvg

for l = 1:size(leftMaskStack,3)
for k = 1:size(firstFrame,2)
    try
    topIndices(:,k,l) = find(leftMaskStack(:,k,l)==1,1,'first');
    bottomIndices(:,k,l) = find(leftMaskStack(:,k,l)==1,1,'last');
    catch
        topIndices(:,k,l) = NaN;
        bottomIndices(:,k,l) = NaN;
    end
end
end

for k = 1:size(leftMaskStack,3)
   topMean = round(topIndices(:,:,k))
      botMean = round(bottomIndices(:,:,k))
    thicknessNotAvg(:,:,k) = (botMean-topMean)*deltaZ
end

thicknessNotAvgL = thicknessNotAvg;
figure;
scatter(1:size(thicknessNotAvgL,3),squeeze(thicknessNotAvgL(:,81,:)))


for l = 1:size(leftMaskStack,3)
for k = 1:size(firstFrame,2)
    try
    topIndices(:,k,l) = find(rightMaskStack(:,k,l)==1,1,'first');
    bottomIndices(:,k,l) = find(rightMaskStack(:,k,l)==1,1,'last');
    catch
        topIndices(:,k,l) = NaN;
        bottomIndices(:,k,l) = NaN;
    end
end
end

for k = 1:size(leftMaskStack,3)
   topMean = round(topIndices(:,:,k))
      botMean = round(bottomIndices(:,:,k))
    thicknessNotAvg(:,:,k) = (botMean-topMean)*deltaZ
end

thicknessNotAvgR = thicknessNotAvg;
figure;
scatter(1:size(thicknessNotAvgR,3),squeeze(thicknessNotAvgR(:,430,:)))


thicknessNotAvgL(isnan(thicknessNotAvgL)) = 0;
thicknessNotAvgR(isnan(thicknessNotAvgR)) = 0;
thicknessNotAvg = squeeze(thicknessNotAvgL+thicknessNotAvgR);


leftThickness = squeeze(thicknessNotAvgL(:,81,:))*1000;
rightThickness = squeeze(thicknessNotAvgR(:,430,:))*1000;

leftStrain = abs(cumsum(diff(leftThickness))/leftThickness(1));
rightStrain = abs(cumsum(diff(rightThickness))/rightThickness(1));
% Looking at soft vs stiff
baseFolder = 'D:\020322AgarHalfAndHalf';
load([baseFolder,filesep,'stressStrainandYM.mat'],'YMSensorStack','strainStack','sensorStressStack',...
    'YMStack','YMLeftStack','YMRightStack')

% clearvars;
CSVFiles = natsortfiles(rdir('D:\020322AgarHalfAndHalf\Instron\agar_phantom_tealyellow020322.is_comp_RawData\*.csv'))
sensor = 1:6;
soft = 13:17
stiff = 8:11
counter = 1;
clearvars Extension Load
for k = stiff
opts = delimitedTextImportOptions("NumVariables", 3);
    
    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Var1", "Extension", "Load"];
    opts.SelectedVariableNames = ["Extension", "Load"];
    opts.VariableTypes = ["string", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
    
    % Import the data
    SpecimenRawData1 = readtable(CSVFiles(k).name, opts);
    
    Extension(counter,:) = table2array(SpecimenRawData1(1:120,1));
    Load(counter,:) = table2array(SpecimenRawData1(1:120,2));
    counter = counter+1;
end

figure;
for k = 1:size(Extension,1)
    plot(Extension(k,:),Load(k,:));
    hold on;
end
hold off; 


figure;
for k = 1:size(Extension,1)
    plot(Extension(k,:),Load(k,:));
    hold on;
end
hold off; 


Extension = mean(Extension); % mm
Load = mean(Load); % N 
figure; plot(Extension,Load)
d = 48.3   
h = 9.36
r = d/2;
SA = pi*r^2 % mm^2
SAm = SA*1e-6
hm = h*1e-3
STRESS = Load/SA; % N / mm^2 is MPa, 1x10^6 Pa
STRAIN = Extension/h; % mm / mm is  unitless
figure; plot(STRAIN*100,STRESS)
ylabel('Stress (MPa)')
xlabel('Strain (%)')
title('Stiff stress vs strain')
% Units: N/m^2

STRESSm = Load/SAm
STRAINm = Extension/h; % should be the same whether changing to m since it is a ratio...

figure; plot(STRAINm,STRESSm); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain')
title('Stiff stress vs strain')

figure; plot(STRAINm*100,STRESSm); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain (%)')
title('Stiff phantom stress vs strain, Instron')
% Get slopes
% From cftool: 2.412e+06*x.^2+1.095e+05*x-380.3 if in Pa or
% 2.412*x.^2+0.1095*x-0.0003803 if using MPa
% Derivative would be... 2*2.412e+06*x+1.095e+05 or 2*2.412*x+0.1095
% y = 2.412e+06*STRAINm.^2+1.095e+05*STRAINm-380.3;
% yprime = 2*2.412e+06*STRAINm+1.095e+05;

% f(x) = p1*x^3 + p2*x^2 + p3*x + p4
%        y = 5.099e+06*STRAINm.^3 + 1.054e+05*STRAINm.^2 + 7.271e+04*STRAINm + 40.51;
%        yprime = 3*5.099e+06*STRAINm.^2+2*1.054e+05*STRAINm+7.271e+04;

       
       y = 4.101e+05.*STRAINm.^2+1.404+04.*STRAINm+207.9;
       yprime1 = 2*4.101e+05*STRAINm+1.404e+04;
figure; plot(y)
figure; plot(yprime1)

sensor = 1:6;
soft = 13:17
stiff = 8:11
counter = 1;
clearvars Extension Load
for k = [8,15]
opts = delimitedTextImportOptions("NumVariables", 3);
    
    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Var1", "Extension", "Load"];
    opts.SelectedVariableNames = ["Extension", "Load"];
    opts.VariableTypes = ["string", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
    
    % Import the data
    SpecimenRawData1 = readtable(CSVFiles(k).name, opts);
    
    Extension(counter,:) = table2array(SpecimenRawData1(1:120,1));
    Load(counter,:) = table2array(SpecimenRawData1(1:120,2));
    counter = counter+1;
end

figure;
for k = 1:size(Extension,1)
    plot(Extension(k,:),Load(k,:));
    hold on;
end
hold off; 


figure;
for k = 1:size(Extension,1)
    plot(Extension(k,:),Load(k,:));
    hold on;
end
hold off; 

% Stiff
k = 1
figure; plot(Extension(k,:),Load(k,:))
d = 48.21   
h = 11.19
r = d/2;
SA = pi*r^2 % mm^2
SAm = SA*1e-6
hm = h*1e-3
STRESS = Load(k,:)/SA; % N / mm^2 is MPa, 1x10^6 Pa
STRAIN = Extension(k,:)/h; % mm / mm is  unitless
figure; plot(STRAIN*100,STRESS)
ylabel('Stress (MPa)')
xlabel('Strain (%)')
title('Stiff stress vs strain')
% Units: N/m^2

STRESSm1 = Load(k,:)/SAm
STRAINm1 = Extension(k,:)/h; % should be the same whether changing to m since it is a ratio...

figure; plot(STRAINm1,STRESSm1); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain')
title('Stiff stress vs strain')

figure; plot(STRAINm1*100,STRESSm1); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain (%)')
title('Agar phantom stress vs strain, Instron')

figure; plot(STRAINm1(1:end-1)*100,diff(STRESSm1))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm1(7:end)*100,STRESSm1(7:end))
title('Agar phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm1(7:end-1)*100,diff(STRESSm1(7:end)))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

% Soft
k = 2
figure; plot(Extension(k,:),Load(k,:))
d = 48.3   
h = 9.36
r = d/2;
SA = pi*r^2 % mm^2
SAm = SA*1e-6
hm = h*1e-3
STRESS = Load(k,:)/SA; % N / mm^2 is MPa, 1x10^6 Pa
STRAIN = Extension(k,:)/h; % mm / mm is  unitless
figure; plot(STRAIN*100,STRESS)
ylabel('Stress (MPa)')
xlabel('Strain (%)')
title('Stiff stress vs strain')
% Units: N/m^2

STRESSm2 = Load(k,:)/SAm
STRAINm2 = Extension(k,:)/h; % should be the same whether changing to m since it is a ratio...




%% Getting fit and deriv with poly 3...
[xData, yData] = prepareCurveData( STRAINm2(1:38), STRESSm2(1:38) );

% Set up fittype and options.
ft = fittype( 'poly3' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf -Inf -Inf 0];
opts.Upper = [Inf Inf Inf 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'STRESSm2 vs. STRAINm2', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'STRAINm2', 'Interpreter', 'none' );
ylabel( 'STRESSm2', 'Interpreter', 'none' );
grid on


coeff1 = fitresult.p1
coeff2 = fitresult.p2
coeff3 = fitresult.p3
coeff4 = fitresult.p4

y = coeff1*STRAINm2.^3 + coeff2*STRAINm2.^2 + coeff3*STRAINm2+coeff4
yprime2= 3*coeff1*STRAINm2.^2 + 2*coeff2*STRAINm2+coeff3;

figure; 
h1 = scatter(leftStrain*100, YMLeft(1:end-1),'filled','b');
hold on;
h2 = scatter(leftStrain*100,YMLeftStack,'filled','r');
h3 = scatter(STRAINm*100,yprime1/1000,'filled','y');
hold off;
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM comparison, stiff side')
xlim([0 4])
ylim([0 50])
b = [h1 h2 h3];
legend(b,'Shear wave','Computed','Instron Fit','Location','SE')

figure; 
h1 = scatter(rightStrain*100, YMRight(1:end-1),'filled','b');
hold on;
h2 = scatter(rightStrain*100,YMRightStack,'filled','r');
h3 = scatter(STRAINm2*100,yprime2/1000,'filled','y');
hold off;
ylabel("Young's modulus (kPa)")
xlabel('Strain (%)')
title('YM comparison, soft side')
xlim([0 4])
ylim([0 20])    
b = [h1 h2 h3];
legend(b,'Shear wave','Computed','Instron Fit','Location','SE')

% figure; 
% h1 = scatter(STRAINm*100,yprime1/1000,'filled','r'); hold on; 
% h2 = scatter(STRAINm2*100,yprime2/1000,'filled','b');
