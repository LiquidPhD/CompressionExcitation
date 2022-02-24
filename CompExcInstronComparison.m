clearvars;
baseFolder = 'D:\021722GelatinBlueYellow';
load([baseFolder,filesep,'stressStrainandYM.mat'],'YMSensorStack','strainStack','sensorStressStack',...
    'YMStack','YMLeftStack','YMRightStack')

% clearvars;
CSVFiles = natsortfiles(rdir('D:\021722GelatinBlueYellow\Instron\GelatinBlueYellow021722.is_comp_RawData\*.csv'))
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
       yprime = 2*4.101e+05*STRAINm+1.404e+04;
figure; plot(y)
figure; plot(yprime)

figure; 
plot(STRAINm,STRESSm); hold on;
plot(STRAINm,y); hold off; 

figure; 
plot(STRAINm*100,yprime)
ylabel('Change in stress (Pa)')
xlabel('Strain (%)')
title('Stiff stress vs strain derivative, Instron')

sensorThickness = 4.05 % mm

compressionAxis = (1:39)*0.01 % compression in mm
strainCompressionAxis = compressionAxis./sensorThickness % Unitless

figure; h1 = plot(strainCompressionAxis(1:length(YMLeftStack))*100,squeeze(YMLeftStack))
hold on; h2 = plot(STRAINm*100,yprime/1000)
b = [h1 h2]
legend(b,'YM from method','YM from Instron')

figure; subplot(2,1,1); plot(strainCompressionAxis*100,squeeze(YMLeftStack(:,:,2:40))); subplot(2,1,2); plot(STRAINm*100,yprime/1000); 
figure; h1 = plot(strainCompressionAxis*100,squeeze(YMLeftStack(:,:,2:40))); ylabel("Young's modulus (kPa)"); ylim([100 700])
yyaxis right; h2 = plot(STRAINm*100,yprime/1000); ylabel("Young's modulus (kPa)"); ylim([100 700])
xlabel('Compression (%)')
title('Calculated data vs Instron')
b = [h1 h2];
legend(b,'Calculated YM','Instron YM')

figure; plot(strainCompressionAxis*100,squeeze(YMLeftStack(:,:,2:40))); hold on; plot(STRAINm*100,yprime/1000); 

figure; plot(yprime(34:74)/1000)
hold on; plot(squeeze(YMLeftStack(:,:,2:40)))
STRAINm(34:74)

figure; h1 = plot(STRAINm(34:74),yprime(34:74)/1000)
hold on; h2 = plot(STRAINm(34:74),squeeze(YMLeftStack(:,:,2:42)))
title('Calculate data vs Instron')
xlabel('Strain')
ylabel("Young's modulus (kPa)");
b = [h1 h2];
legend(b,'Instron YM','Calculated YM')


%% Check compression % using stiff phantom measurements...

firstFrame = leftMaskStack(:,:,1);
figure; imagesc(firstFrame)
for k = 1:size(firstFrame,2)
    try
    topIndex(:,k) = find(firstFrame(:,k)==1,1,'first');
    bottomIndex(:,k) = find(firstFrame(:,k)==1,1,'last');
    catch
        topIndex(:,k) = NaN;
        bottomIndex(:,k) = NaN;
    end
end

firstFrameTopMean = round(mean(topIndex,'omitnan'));
firstFrameBotMean = round(mean(bottomIndex,'omitnan'));
deltaZ = Parameters.delta_z;
stiffThickness = (firstFrameBotMean-firstFrameTopMean)*deltaZ*1000 % in mm

for l = 1:40
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

for k = 1:40
   topMean = round(mean(topIndices(:,:,k),'omitnan'));
      botMean = round(mean(bottomIndices(:,:,k),'omitnan'));
    thickness(:,k) = (botMean-topMean)*deltaZ
end

for l = 1:40
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

for k = 1:40
   topMean = round(topIndices(:,:,k))
      botMean = round(bottomIndices(:,:,k))
    thicknessNotAvg(:,:,k) = (botMean-topMean)*deltaZ
end

scatter(1:40,squeeze(thicknessNotAvg(:,81,:)))

for k = 1:40
    figure(1)
    plot(thicknessNotAvg(:,:,k));
    hold on; 
    pause
end

figure;
scatter(1:40,squeeze(thicknessNotAvg(:,124,:))); hold on; 
scatter(1:40,squeeze(thicknessNotAvg(:,509,:))); hold off; 

stiffDiff = abs(thicknessNotAvg(:,124,end)-thicknessNotAvg(:,124,1))
softDiff = abs(thicknessNotAvg(:,509,end)-thicknessNotAvg(:,509,1))

figure(1)
for k = 1:40
 
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


% Try with relabeled data: labelsForThickness

for k = 1:size(labelsForThickness,3)
    label = labelsForThickness(:,:,k);
    stiff = zeros([size(label,1) size(label,2)]);
    stiff(label== 'Stiff') = 1;
    soft = zeros([size(label,1) size(label,2)]);
    soft(label== 'Soft') = 2;
    stiffStack(:,:,k) = logical(stiff);
    softStack(:,:,k) = logical(soft);

end


firstFrame = stiffStack(:,:,1);
figure; imagesc(firstFrame)
for k = 1:size(firstFrame,2)
    try
    topIndex(:,k) = find(firstFrame(:,k)==1,1,'first');
    bottomIndex(:,k) = find(firstFrame(:,k)==1,1,'last');
    catch
        topIndex(:,k) = NaN;
        bottomIndex(:,k) = NaN;
    end
end

firstFrameTopMean = round(mean(topIndex,'omitnan'));
firstFrameBotMean = round(mean(bottomIndex,'omitnan'));
deltaZ = zaxis(2)
stiffThickness = (firstFrameBotMean-firstFrameTopMean)*deltaZ % in mm

for l = 1:40
for k = 1:size(leftMaskStack,2)
    try
    topIndices(:,k,l) = find(stiffStack(:,k,l)==1,1,'first');
    bottomIndices(:,k,l) = find(stiffStack(:,k,l)==1,1,'last');
    catch
        topIndices(:,k,l) = NaN;
        bottomIndices(:,k,l) = NaN;
    end
%     figure(1)
%     imagesc(images(:,:,l));
%     colormap(gray)
%     caxis([0 0.1]);
%     hold on;
%     hline(topIndices(:,k,l))
%     hline(bottomIndices(:,k,l))
%     hold off;
%     pause
end
end
% 
% for k = 1:40
%    topMean = round(mean(topIndices(:,:,k),'omitnan'));
%       botMean = round(mean(bottomIndices(:,:,k),'omitnan'));
%     thickness(:,k) = (botMean-topMean)*deltaZ
% end

for l = 1:40
for k = 1:size(topIndices,2)
   topMean = round(mean(topIndices(:,k,l),'omitnan'));
      botMean = round(mean(bottomIndices(:,k,l),'omitnan'));
      
    thickness2(:,k,l) = bottomIndices(:,k,l)-topIndices(:,k,l);
end
end
% 
% for k = 1:size(thickness2,3)
%     figure(1);
%     plot(squeeze(thickness2(:,:,k)))
%     pause
% end

for k = 1:size(thickness2,3)
    avgThicknessStiff(:,k) = mean(squeeze(thickness2(:,50:250,k)),'omitnan');
end
    
    
% for k = 1:40
%     figure(1)
%     imagesc(images(:,:,k))
%     colormap(gray)
%     caxis([0 0.1])
%     pause
% end

% Try with relabeled data: labelsForThickness on soft side
for k = 1:size(labelsForThickness,3)
    label = labelsForThickness(:,:,k);
    stiff = zeros([size(label,1) size(label,2)]);
    stiff(label== 'Stiff') = 1;
    soft = zeros([size(label,1) size(label,2)]);
    soft(label== 'Soft') = 2;
    stiffStack(:,:,k) = logical(stiff);
    softStack(:,:,k) = logical(soft);

end


firstFrame = softStack(:,:,1);
figure; imagesc(firstFrame)
for k = 1:size(firstFrame,2)
    try
    topIndex(:,k) = find(firstFrame(:,k)==1,1,'first');
    bottomIndex(:,k) = find(firstFrame(:,k)==1,1,'last');
    catch
        topIndex(:,k) = NaN;
        bottomIndex(:,k) = NaN;
    end
end

firstFrameTopMean = round(mean(topIndex,'omitnan'));
firstFrameBotMean = round(mean(bottomIndex,'omitnan'));
deltaZ = zaxis(2)
softThickness = (firstFrameBotMean-firstFrameTopMean)*deltaZ % in mm

for l = 1:40
for k = 1:size(softStack,2)
    try
    topIndices(:,k,l) = find(softStack(:,k,l)==1,1,'first');
    bottomIndices(:,k,l) = find(softStack(:,k,l)==1,1,'last');
%      figure(1)
%     imagesc(images(:,:,l));
%     colormap(gray)
%     caxis([0 0.1]);
%     hold on;
%     hline(topIndices(:,k,l))
%     hline(bottomIndices(:,k,l))
%     hold off;
%     pause
    catch
        topIndices(:,k,l) = NaN;
        bottomIndices(:,k,l) = NaN;
    end
   
end
end
% 
% for k = 1:40
%    topMean = round(mean(topIndices(:,:,k),'omitnan'));
%       botMean = round(mean(bottomIndices(:,:,k),'omitnan'));
%     thickness(:,k) = (botMean-topMean)*deltaZ
% end

for l = 1:40
for k = 1:size(topIndices,2)
   topMean = round(mean(topIndices(:,k,l),'omitnan'));
      botMean = round(mean(bottomIndices(:,k,l),'omitnan'));
      
    thickness2(:,k,l) = bottomIndices(:,k,l)-topIndices(:,k,l);
end
end
% 
% for k = 1:size(thickness2,3)
%     figure(1);
%     plot(squeeze(thickness2(:,:,k)))
%     pause
% end

for k = 1:size(thickness2,3)
    avgThicknessSoft(:,k) = mean(squeeze(thickness2(:,350:550,k)),'omitnan');
end

figure; h1 = plot(avgThicknessStiff); hold on; h2 = plot(avgThicknessSoft); hold off;
b = [h1 h2]
legend(b,'Stiff','Soft')


% Try with sensor

firstFrame = sensorMaskStack(:,:,1);
figure; imagesc(firstFrame)
for k = 1:size(firstFrame,2)
    try
    topIndex(:,k) = find(firstFrame(:,k)==1,1,'first');
    bottomIndex(:,k) = find(firstFrame(:,k)==1,1,'last');
    catch
        topIndex(:,k) = NaN;
        bottomIndex(:,k) = NaN;
    end
end

firstFrameTopMean = round(mean(topIndex,'omitnan'));
firstFrameBotMean = round(mean(bottomIndex,'omitnan'));
deltaZ = zaxis(2)
sensorThickness = (firstFrameBotMean-firstFrameTopMean)*deltaZ % in mm

for l = 1:40
for k = 1:size(sensorMaskStack,2)
    try
    topIndices(:,k,l) = find(sensorMaskStack(:,k,l)==1,1,'first');
    bottomIndices(:,k,l) = find(sensorMaskStack(:,k,l)==1,1,'last');
%      figure(1)
%     imagesc(images(:,:,l));
%     colormap(gray)
%     caxis([0 0.1]);
%     hold on;
%     hline(topIndices(:,k,l))
%     hline(bottomIndices(:,k,l))
%     hold off;
%     pause
    catch
        topIndices(:,k,l) = NaN;
        bottomIndices(:,k,l) = NaN;
    end
   
end
end
% 
% for k = 1:40
%    topMean = round(mean(topIndices(:,:,k),'omitnan'));
%       botMean = round(mean(bottomIndices(:,:,k),'omitnan'));
%     thickness(:,k) = (botMean-topMean)*deltaZ
% end

for l = 1:40
for k = 1:size(topIndices,2)
   topMean = round(mean(topIndices(:,k,l),'omitnan'));
      botMean = round(mean(bottomIndices(:,k,l),'omitnan'));
      
    thicknessSensor(:,k,l) = bottomIndices(:,k,l)-topIndices(:,k,l);
end
end

for k = 1:size(thicknessSensor,3)
    figure(1)
    plot(thicknessSensor(:,:,k))
    pause
end
% 
% for k = 1:size(thickness2,3)
%     figure(1);
%     plot(squeeze(thickness2(:,:,k)))
%     pause
% end

for k = 1:size(thickness2,3)
    avgThicknessSensorStiffSide(:,k) = mean(squeeze(thicknessSensor(:,50:250,k)),'omitnan');
    avgThicknessSensorSoftSide(:,k) = mean(squeeze(thicknessSensor(:,350:550,k)),'omitnan');

end

figure; h1 = plot(avgThicknessStiff); hold on; h2 = plot(avgThicknessSoft); h3 = plot(avgThicknessSensorStiffSide); h4 = plot(avgThicknessSensorSoftSide); hold off;
b = [h1 h2 h3 h4]
legend(b,'Stiff','Soft','Sensor (stiff)', 'Sensor (soft)')

figure; h3 = plot(avgThicknessSensorStiffSide); hold on; h4 = plot(avgThicknessSensorSoftSide); hold off;
b = [h3 h4]
legend(b,'Sensor (stiff)', 'Sensor (soft)')

avgThicknessStiffMM = avgThicknessSensorStiffSide*deltaZ;
avgThicknessSoftMM = avgThicknessSensorSoftSide*deltaZ;

figure; h3 = plot(avgThicknessStiffMM); hold on; h4 = plot(avgThicknessSoftMM); hold off;
b = [h3 h4]
legend(b,'Sensor (stiff)', 'Sensor (soft)')

 % See how much sensor changes with each compression

for k = 1:length(avgThicknessStiffMM)
    stiffThicknessChange(:,k) = avgThicknessStiffMM(:,1) - avgThicknessStiffMM(:,k);
    softThicknessChange(:,k) = avgThicknessSoftMM(:,1) - avgThicknessSoftMM(:,k);
end

figure; h1 = plot(stiffThicknessChange); hold on; h2 = plot(softThicknessChange);
b = [h1 h2]
legend(b,'Stiff','Soft')

% Knowing that each compression was actually 0.01 mm
stiffAmtChange = (1:40)*0.01-stiffThicknessChange
softAmtChange = (1:40)*0.01-softThicknessChange

xaxisForCompression = (1:40)*0.01;
figure; 
h1 = plot(xaxisForCompression,stiffAmtChange); hold on; h2 = plot(xaxisForCompression,softAmtChange); hold off;
b = [h1 h2]; legend(b,'Stiff','Soft')
xlabel('Compression (mm)')
ylabel('Change in thickness')


figure; 
h1 = plot(xaxisForCompression,stiffAmtChange./xaxisForCompression); hold on; h2 = plot(xaxisForCompression,softAmtChange./xaxisForCompression); hold off;
b = [h1 h2]; legend(b,'Stiff','Soft')
xlabel('Compression (mm)')
ylabel('Change in thickness (%)')

