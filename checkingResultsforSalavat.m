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
for k = 1:size(Extension,1)

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
title('Stiff phantom stress vs strain, Instron')

figure; plot(STRAINm(1:end-1)*100,diff(STRESSm))
title('Derivative of stress vs strain %, stiff Agar')
ylabel('Stress (Pa)')
xlabel('Strain (%)')


% Getting fit and derivative with poly2...
[xData, yData] = prepareCurveData( STRAINm, STRESSm );

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

figure; plot(y); hold on; plot(yprime); hold off; 
figure; plot(STRAINm,y); hold on; plot(STRAINm,STRESSm); hold off;

figure; plot(STRAINm,yprime/1000)


figure; plot(STRAINm(1:end-1)*100,diff(STRESSm))
title('Derivative of stress vs strain %, stiff Agar')
ylabel('Stress (Pa)')
xlabel('Strain (%)')
hold on; 
plot(STRAINm(1:end-1)*100,yprime(1:end-1)/1000);
hold off;

%% Getting fit and deriv with poly 3...
[xData, yData] = prepareCurveData( STRAINm, STRESSm );

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

figure; plot(y); hold on; plot(yprime); hold off; 
figure; plot(STRAINm,y); hold on; plot(STRAINm,STRESSm); hold off;

figure; plot(STRAINm,yprime/1000)


figure; plot(STRAINm(1:end-1)*100,diff(STRESSm))
title('Derivative of stress vs strain %, stiff Agar')
ylabel('Stress (Pa)')
xlabel('Strain (%)')
hold on; 
plot(STRAINm(1:end-1)*100,yprime(1:end-1)/1000);
hold off;


