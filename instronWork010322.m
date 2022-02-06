% clearvars;
CSVFiles = natsortfiles(rdir('G:\AgarHalfHalf01142022\agar_phantom_RIP_9.is_comp_RawData\*.csv'))
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
d = 33.66   
h = 10.17
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
title('Stiff stress vs strain')
% Get slopes
% From cftool: 2.412e+06*x.^2+1.095e+05*x-380.3 if in Pa or
% 2.412*x.^2+0.1095*x-0.0003803 if using MPa
% Derivative would be... 2*2.412e+06*x+1.095e+05 or 2*2.412*x+0.1095
y = 2.412e+06*STRAINm.^2+1.095e+05*STRAINm-380.3;
yprime = 2*2.412e+06*STRAINm+1.095e+05;
figure; plot(y)
figure; plot(yprime)

figure; 
plot(STRAINm,STRESSm); hold on;
plot(STRAINm,y); hold off; 

figure; 
plot(STRAINm*100,yprime)
ylabel('Change in stress (Pa)')
xlabel('Strain (%)')
title('Stiff stress vs strain')
% STRAIN = STRAIN(62:end);
% STRESS = STRESS(62:end);
% figure; plot(STRAIN,STRESS)
% 
% firstIndex = 70
% secondIndex = 40
% (STRESS(firstIndex)-STRESS(secondIndex)) / (STRAIN(firstIndex)-STRAIN(secondIndex))
% 
% % Units: N/m^2
% 
% STRESSm = Load/SAm
% STRAINm = Extension/hm;
% 
% (STRESSm(firstIndex)-STRESSm(secondIndex)) / (STRAINm(firstIndex)-STRAINm(secondIndex))
% 
% % Strain pct
% figure; plot(STRAIN,STRESS)
% 
% % Sliding window slopes
% windowLength = 10;
% for k = windowLength:length(STRAIN)
%     slopeInstron(k-(windowLength-1)) = (STRESSm(k)-STRESSm(k-(windowLength-1))) / (STRAINm(k)-STRAINm(k-(windowLength-1)))
% end
% 
% figure; plot(slopeInstron)
% 
% STRESSm(/STRAINm

