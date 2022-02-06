% clearvars;
CSVFiles = natsortfiles(rdir('D:\AgarHalfHalf01142022\agar_phantom_RIP_9.is_comp_RawData\*.csv'))
sensor = 1:6;
soft = 7:12
stiff = 1:6
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

Extension = mean(Extension); % mm
Load = mean(Load); % N 
figure; plot(Extension,Load)
d = 32.52 
h = 9.2
r = d/2;
SA = pi*r^2 % mm^2
SAm = SA*1e-6
hm = h*1e-3
STRESS = Load/SA;
STRAIN = Extension/h;
figure; plot(STRAIN*100,STRESS)
ylabel('Stress')
xlabel('Strain (%)')
title('Stiff stress vs strain')


xaxisPCT = (xaxis./4000)*100

close all force; 
scatter(xaxisPCT(2:40),squeeze(YMLeftStack(:,:,2:40)),'r','filled')
h = lsline;
h.Color = 'r';
title('YM vs compression, stiff side')
xlabel('Strain (%)')

ylabel('YM (kPa)')

% Fit from CFTool

y = 0.0002115*xaxis.^2+0.001061*xaxis-0.0004074
% deriv = 
yprime = 2*0.0002115*xaxis+0.001061

figure; plot(xaxis,y); hold on; plot(STRAIN*100,STRESS); hold off;
figure; plot(xaxis,yprime)
