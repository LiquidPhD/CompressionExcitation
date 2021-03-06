baseFolder = 'D:\021722GelatinBlueYellow';
load([baseFolder,filesep,'stressStrainandYM.mat'],'YMSensorStack','strainStack','sensorStressStack',...
    'YMStack','YMLeftStack','YMRightStack')

% clearvars;
CSVFiles = natsortfiles(rdir('D:\021722GelatinBlueYellow\Instron\GelatinBlueYellow021722.is_comp_RawData\*.csv'))


sensor = 6:10;
soft = 1:5
stiff = 11:15
counter = 1;
clearvars Extension Load
for k = [13,2]
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
d = 48.58   
h = 11.45
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

STRESSm1 = STRESSm1(26:end);
STRAINm1 = STRAINm1(26:end);
figure; plot(STRAINm1,STRESSm1); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain')
title('Stiff stress vs strain')

figure; plot(STRAINm1*100,STRESSm1); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain (%)')
title('Gelatin phantom stress vs strain, Instron')

figure; plot(STRAINm1(1:end-1)*100,diff(STRESSm1))
title('Derivative of Instron stress vs strain %, gelatin')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm1*100,STRESSm1)
title('Gelatin phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm1(1:end-1)*100,diff(STRESSm1(:)))
title('Derivative of Instron stress vs strain %, gelatin')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

% Soft
k = 2
figure; plot(Extension(k,:),Load(k,:))
d = 48.03   
h = 10.07
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

STRESSm2 = STRESSm2(5:end);
STRAINm2 = STRAINm2(5:end);
figure; plot(STRAINm2,STRESSm2); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain')
title('Soft stress vs strain')

figure; plot(STRAINm2*100,STRESSm2); % Uses Pa
ylabel('Stress (Pa)')
xlabel('Strain (%)')
title('gelatin phantom stress vs strain, Instron')

figure; plot(STRAINm2(1:end-1)*100,diff(STRESSm2))
title('Derivative of Instron stress vs strain %, gelatin')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(1:end-19)*100,STRESSm2(20:end))
title('gelatin phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(1:end-20)*100,diff(STRESSm2(20:end)))
title('Derivative of Instron stress vs strain %, gelatin')
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
title('gelatin phantom stress vs strain, Instron')

figure; h1 = plot(STRAINm1(1:end-1)*100,diff(STRESSm1)); hold on; 
h2 = plot(STRAINm2(1:end-1)*100,diff(STRESSm2))
hold off;
b = [h1 h2]
legend(b,'Stiff','Soft','Location','NW')
title('Derivative of Instron stress vs strain %, gelatin')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(7:end)*100,STRESSm2(7:end))
title('gelatin phantom stress vs strain %, Instron')
ylabel('Stress (Pa)')
xlabel('Strain (%)')

figure; plot(STRAINm2(7:end-1)*100,diff(STRESSm2(7:end)))
title('Derivative of Instron stress vs strain %, gelatin')
ylabel('\DeltaStress (Pa)')
xlabel('Strain (%)')


figure;
h1 = plot(diff(STRESSm1(7:end))); hold on; h2 = plot(diff(STRESSm2(4:end))); hold off; 
title('Derivative of Instron stress vs data point')
ylabel('\DeltaStress (Pa)')
xlabel('Compression sampling point')


figure;
h1 = plot(STRAINm2(1:end-7)*100,diff(STRESSm1(7:end))./diff(STRESSm2(4:end-3)));
title('Ratio of stiff YM to soft using Instron data, gelatin')
ylabel('stiff / soft')
xlabel('Strain (%)')



% Check change in thickness of phantom



deltaZ = Parameters.delta_z;

comboStack = leftMaskStack+rightMaskStack;
for l = 1:size(comboStack,3)
for k = 1:size(comboStack,2)
    try
    topIndices(:,k,l) = find(comboStack(:,k,l)==1,1,'first');
    bottomIndices(:,k,l) = find(comboStack(:,k,l)==1,1,'last');
    catch
        topIndices(:,k,l) = NaN;
        bottomIndices(:,k,l) = NaN;
    end
end
end

for k = 1:size(comboStack,3)
   topMean = round(topIndices(:,:,k))
      botMean = round(bottomIndices(:,:,k))
    thicknessNotAvg(:,:,k) = (botMean-topMean)*deltaZ*1000
end

scatter(1:40,squeeze(thicknessNotAvg(:,81,:)))

% for k = 1:40
%     figure(1)
%     plot(thicknessNotAvg(:,:,k));
%     hold on; 
%     pause
% end

figure;
scatter(1:size(thicknessNotAvg,3),squeeze(thicknessNotAvg(:,124,:))); hold on; 
scatter(1:size(thicknessNotAvg,3),squeeze(thicknessNotAvg(:,509,:))); hold off; 

stiffDiff = abs(thicknessNotAvg(:,124,end)-thicknessNotAvg(:,124,1))
softDiff = abs(thicknessNotAvg(:,509,end)-thicknessNotAvg(:,509,1))

figure(1)
for k = 1:size(comboStack,3)
 
    ST_Diff(:,k) = abs(thicknessNotAvg(:,124,1)-thicknessNotAvg(:,124,k))/thicknessNotAvg(:,124,1);
    SO_Diff(:,k) = abs(thicknessNotAvg(:,509,1)-thicknessNotAvg(:,509,k))/thicknessNotAvg(:,509,1);
end
figure;
h1 = scatter(1:k,ST_Diff);
yyaxis right;
h2 = scatter(1:k,SO_Diff);
b = [h1 h2]
legend(b,'Stiff','Soft')


figure;
scatter(1:40,squeeze(thicknessNotAvg(:,124,:))); hold on; 
scatter(1:40,squeeze(thicknessNotAvg(:,509,:))); hold off; 
for k = 1:40
    figure(1)
    imagesc(images(:,:,k))
    colormap(gray)
    caxis([0 0.1])
    pause
end



figure; h1 = plot(STRAINm1(26:end-1)*100,diff(STRESSm1(26:end))); hold on;
h2 = scatter(ST_Diff(10:end-14)*100,YMLeftStack(18:end-4),'filled'); 
h3 = hline(mean(YMLeftStack(18:end-4)));hold off; 
title("Young's modulus vs strain %, Stiff side of phantom")
ylabel("YM (kPa)")
xlabel('Strain (%)')
b = [h1 h2 h3];
legend(b,'Diff of Instron','Calculated value', 'Mean of calculated value')


figure; h1 = plot(STRAINm2(26:end-13)*100,diff(STRESSm2(26:end-12))); hold on;
h2 = scatter(SO_Diff(10:end-18)*100,YMRightStack(18:end-8),'filled');
h3 = hline(mean(YMRightStack(18:end-4)));hold off; 

title("Young's modulus vs strain %, Soft side of phantom")
ylabel("YM (kPa)")
xlabel('Strain (%)')
b = [h1 h2 h3];
legend(b,'Diff of Instron','Calculated value', 'Mean of calculated value')
