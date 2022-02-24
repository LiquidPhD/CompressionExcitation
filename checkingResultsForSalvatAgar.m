% Gelatin testing
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
for k = 8
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

STRESSm = Load(k,:)/SAm
STRAINm = Extension(k,:)/h; % should be the same whether changing to m since it is a ratio...

figure; plot(STRAINm,STRESSm); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain')
title('Stiff stress vs strain')

figure; plot(STRAINm*100,STRESSm); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain (%)')
title('Agar phantom stress vs strain, Instron')

figure; plot(STRAINm(1:end-1)*100,diff(STRESSm))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm(7:end)*100,STRESSm(7:end))
title('Agar phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')


% Getting fit and derivative with poly1...
[xData, yData] = prepareCurveData( STRAINm(7:end), STRESSm(7:end) );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'STRESSm vs. STRAINm', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'STRAINm', 'Interpreter', 'none' );
ylabel( 'STRESSm', 'Interpreter', 'none' );
grid on

coeff1 = fitresult.p1
coeff2 = fitresult.p2


y = coeff1*STRAINm + coeff2;
yprime = coeff1*ones([length(y)-1 1]);

y1 = coeff1*STRAINm + coeff2;
yprime1 = coeff1*ones([length(y)-1 1]);

figure; plot(y); hold on; plot(yprime); hold off; 
figure; plot(STRAINm,y); hold on; plot(STRAINm,STRESSm); hold off;

figure; plot(STRAINm(1:end-1),squeeze(yprime)/1000)


figure; plot(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
title('Derivative of stress vs strain %, gelatin')
ylabel('Stress (Pa)')
xlabel('Strain (%)')
hold on; 
plot(STRAINm(1:end-1)*100,yprime/1000);
hold off;

figure; 
h1 = plot(STRAINm(7:end)*100,STRESSm(7:end)); hold on; 
h2 = plot(STRAINm(7:end)*100,y(7:end)); hold off; 
ylim([0 10000])

xlabel('Strain (%)')
ylabel('Stress (Pa)')
title('Agar Instron vs poly1 fit')
b = [h1 h2];
legend(b,'Instron','Poly1 fit','Location','NW')
text(8,3000,['SSE: ',num2str(gof.sse)])
text(8,2500,['R^2: ',num2str(gof.rsquare)])
text(8,2000,['dfe: ',num2str(gof.dfe)])
text(8,1500,['Adjusted R^2: ',num2str(gof.adjrsquare)])
text(8,1000,['RMSE: ',num2str(gof.rmse)])

figure; scatter(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
h1 = lsline;
title('Derivative of stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')
hold on; 
h2 = plot(STRAINm(1:end-1)*100,yprime/1000);
hold off;
b = [h1 h2]
legend(b,'Least squares line','Deriv of poly1 fit','Location','NW')


% Getting fit and derivative with poly2...
[xData, yData] = prepareCurveData( STRAINm(7:end), STRESSm(7:end) );

% Set up fittype and options.
ft = fittype( 'poly2' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf -Inf 0];
opts.Upper = [Inf Inf 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'STRESSm vs. STRAINm', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'STRAINm', 'Interpreter', 'none' );
ylabel( 'STRESSm', 'Interpreter', 'none' );
grid on

coeff1 = fitresult.p1
coeff2 = fitresult.p2
coeff3 = fitresult.p3

y = coeff1*STRAINm.^2 + coeff2*STRAINm + coeff3
yprime = 2*coeff1*STRAINm + coeff2

y2 = coeff1*STRAINm.^2 + coeff2*STRAINm + coeff3
yprime2 = 2*coeff1*STRAINm + coeff2
yprime2 = yprime2(1:end-1);

figure; plot(y); hold on; plot(yprime); hold off; 
figure; plot(STRAINm,y); hold on; plot(STRAINm,STRESSm); hold off;

figure; plot(STRAINm,yprime/1000)


figure; 
h1 = plot(STRAINm(7:end)*100,STRESSm(7:end)); hold on; 
h2 = plot(STRAINm(7:end)*100,y(7:end)); hold off; 
xlabel('Strain (%)')
ylabel('Stress (Pa)')
title('Agar Instron vs poly2 fit')
b = [h1 h2];
legend(b,'Instron','Poly2 fit','Location','NW')
text(8,3000,['SSE: ',num2str(gof.sse)])
text(8,2500,['R^2: ',num2str(gof.rsquare)])
text(8,2000,['dfe: ',num2str(gof.dfe)])
text(8,1500,['Adjusted R^2: ',num2str(gof.adjrsquare)])
text(8,1000,['RMSE: ',num2str(gof.rmse)])

figure; scatter(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
h1 = lsline;
title('Derivative of stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')
hold on; 
h2 = plot(STRAINm(1:end-1)*100,yprime(1:end-1)/1000);
hold off;
b = [h1 h2]
legend(b,'Least squares line','Deriv of poly2 fit','Location','NW')

%% Getting fit and deriv with poly 3...
[xData, yData] = prepareCurveData( STRAINm(7:end), STRESSm(7:end) );

% Set up fittype and options.
ft = fittype( 'poly3' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'STRESSm vs. STRAINm', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'STRAINm', 'Interpreter', 'none' );
ylabel( 'STRESSm', 'Interpreter', 'none' );
grid on



coeff1 = fitresult.p1
coeff2 = fitresult.p2
coeff3 = fitresult.p3
coeff4 = fitresult.p4

y = coeff1*STRAINm.^3 + coeff2*STRAINm.^2 + coeff3*STRAINm+coeff4
yprime = 3*coeff1*STRAINm.^2 + 2*coeff2*STRAINm+coeff3;

y3 = coeff1*STRAINm.^3 + coeff2*STRAINm.^2 + coeff3*STRAINm+coeff4
yprime3 = 3*coeff1*STRAINm.^2 + 2*coeff2*STRAINm+coeff3;

yprime3 = yprime3(1:end-1);
figure; plot(y); hold on; plot(yprime); hold off; 
figure; plot(STRAINm,y); hold on; plot(STRAINm,STRESSm); hold off;

figure; plot(STRAINm,yprime/1000)


figure; plot(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
title('Derivative of stress vs strain %, gelatin')
ylabel('Stress (Pa)')
xlabel('Strain (%)')
hold on; 
plot(STRAINm(1:end-1)*100,yprime(1:end-1)/1000);
hold off;

figure; scatter(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
lsline;
title('Derivative of stress vs strain %, gelatin')
ylabel('Stress (Pa)')
xlabel('Strain (%)')
hold on; 
plot(STRAINm(1:end-1)*100,yprime(1:end-1)/1000);
hold off;


figure; 
h1 = plot(STRAINm(7:end)*100,STRESSm(7:end)); hold on; 
h2 = plot(STRAINm(7:end)*100,y(7:end)); hold off; 
xlabel('Strain (%)')
ylabel('Stress (Pa)')
title('Agar Instron vs poly3 fit')
b = [h1 h2];
legend(b,'Instron','Poly3 fit','Location','NW')
text(8,3000,['SSE: ',num2str(gof.sse)])
text(8,2500,['R^2: ',num2str(gof.rsquare)])
text(8,2000,['dfe: ',num2str(gof.dfe)])
text(8,1500,['Adjusted R^2: ',num2str(gof.adjrsquare)])
text(8,1000,['RMSE: ',num2str(gof.rmse)])

figure; scatter(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
h1 = lsline;
title('Derivative of stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')
hold on; 
h2 = plot(STRAINm(1:end-1)*100,yprime(1:end-1)/1000);
hold off;
b = [h1 h2]
legend(b,'Least squares line','Deriv of poly3 fit','Location','NW')


% Overlay with all fittings

figure; scatter(STRAINm(7:end-1)*100,diff(STRESSm(7:end)))
h1 = lsline;
title('Derivative of stress vs strain %, gelatin')
ylabel('Stress (Pa)')
xlabel('Strain (%)')
hold on; 
h2 = plot(STRAINm(1:end-1)*100,yprime1(1:end)/1000);
h3 = plot(STRAINm(1:end-1)*100,yprime2(1:end)/1000);
h4 = plot(STRAINm(1:end-1)*100,yprime3(1:end)/1000);
hold off;
b = [h1 h2 h3 h4];
legend(b,'Least squares line','Poly1 fit','Poly2 fit','Poly3 fit')


% Looking at soft vs stiff

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

figure; plot(STRAINm2,STRESSm2); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain')
title('Soft stress vs strain')

figure; plot(STRAINm2*100,STRESSm2); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain (%)')
title('Agar phantom stress vs strain, Instron')

figure; plot(STRAINm2(1:end-1)*100,diff(STRESSm2))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(1:end-19)*100,STRESSm2(20:end))
title('Agar phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(1:end-20)*100,diff(STRESSm2(20:end)))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

% Plot both
figure; h1 = plot(STRAINm1*100,STRESSm1); hold on;
h2 = plot(STRAINm2*100,STRESSm2); % Uses Pa
hold off;  
ylabel('Stress (Pa)')
xlabel('Strain (%)')
b = [h1 h2];
legend(b,'Stiff','Soft','Location','NW')
title('Agar phantom stress vs strain, Instron')

figure; h1 = plot(STRAINm1(1:end-1)*100,diff(STRESSm1)); hold on; 
h2 = plot(STRAINm2(1:end-1)*100,diff(STRESSm2))
hold off;
legend(b,'Stiff','Soft','Location','NW')
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(7:end)*100,STRESSm2(7:end))
title('Agar phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(7:end-1)*100,diff(STRESSm2(7:end)))
title('Derivative of Instron stress vs strain %, agar')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')


figure;
h1 = plot(diff(STRESSm1(7:end))); hold on; h2 = plot(diff(STRESSm2(4:end))); hold off; 
title('Derivative of Instron stress vs data point')
ylabel('\DeltaStress (Pa)')
xlabel('Compression sampling point')


figure;
h1 = plot(STRAINm2(1:end-7)*100,diff(STRESSm1(7:end))./diff(STRESSm2(4:end-3)));
title('Ratio of stiff YM to soft using Instron data, agar')
ylabel('stiff / soft')
xlabel('Strain (%)')


YMSubtraction = abs(diff(STRAINm1(1:length(STRAINm2)))-diff(STRAINm2));
figure; plot(YMSubtraction)

YMDivision = diff(STRAINm1(1:length(STRAINm2)))./diff(STRAINm2);
figure; plot(YMDivision)
