clearvars;
baseFolder = 'F:\AgarPhantomHalfHalfAgarSensor'
load([baseFolder,filesep,'BScanStrainAllSteps.mat']);
load([baseFolder,filesep,'Static2.mat']);
baseFolder = 'F:\AgarPhantomHalfHalfAgarSensor'
centerVelocity = load([baseFolder,filesep,'CenterVelocity2.mat'],'VELOCITIES');
centerVelocity = centerVelocity.VELOCITIES;
axisData = load([baseFolder,filesep,'CenterVelocity2.mat'],'xaxis','zaxis');
xaxis = axisData.xaxis;
zaxis = axisData.zaxis;

figure(1);
imagesc(xaxis,zaxis(1:size(IQData,1)),mat2gray(abs(IQData(:,:,1)))); caxis([0 0.1]);
colormap(gray);
[leftPhantom,xLeft,yLeft] = roipoly;
[rightPhantom,xRight,yRight] = roipoly;

leftCoords = [xLeft yLeft];
rightCoords = [xRight yRight];

[leftLabeled, ~] = bwlabel(leftPhantom ~= 0);
measurements = regionprops(leftLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsLeft = allCentroids(1);
yCentroidsLeft = allCentroids(2);

[rightLabeled, ~] = bwlabel(rightPhantom ~= 0);
measurements = regionprops(rightLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsRight = allCentroids(1:2:end);
yCentroidsRight = allCentroids(2:2:end);

for k = 1:20 % :size(BScanStrainAll,3)
    VELOCITIES = centerVelocity(:,:,k+1);
    BSS = abs(squeeze(BScanStrainAll(:,:,k)))';
    VELOCITIES = VELOCITIES(1:size(BSS,1),:);
    for columnIndex = 1:size(VELOCITIES,1)
        for rowIndex = 1:size(VELOCITIES,2)
            speedAtPixel = VELOCITIES(columnIndex,rowIndex);
            strainAtPixel = BSS(columnIndex,rowIndex);
            sigma(columnIndex,rowIndex) = 3*1000*(speedAtPixel^2)*strainAtPixel;
        end
    end
     DY = diff(zaxis);
        DY = DY(1);
        leftCoords = [xLeft yLeft-(k-1)*5*DY];
rightCoords = [xRight yRight-(k-1)*5*DY];

[leftLabeled, ~] = bwlabel(leftPhantom ~= 0);
measurements = regionprops(leftLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsLeft = allCentroids(1);
yCentroidsLeft = allCentroids(2);

[rightLabeled, ~] = bwlabel(rightPhantom ~= 0);
measurements = regionprops(rightLabeled, 'Centroid');
allCentroids = [measurements.Centroid];
xCentroidsRight = allCentroids(1:2:end);
yCentroidsRight = allCentroids(2:2:end);
    for depthValue = 10:10:450
        stressMM = movmean(abs(sigma(depthValue,:)),50,'omitnan');
        for columnIndex = 1:size(VELOCITIES,1)
            for rowIndex = 1:size(VELOCITIES,2)
                %         stressAtPixel = sigma(columnIndex,rowIndex);
                %         stressAtPixel = axialStress(:,rowIndex);
                % stressAtPixel = mean(sigma(interface,rowIndex));
                stressAtPixel = stressMM(:,rowIndex);
                strainAtPixel = BSS(columnIndex,rowIndex);
                YM(columnIndex,rowIndex) = (abs(stressAtPixel)/abs(strainAtPixel))/1000; % in kPa
            end
        end
       


        figure;
        t = tiledlayout(1,1);
        ax1 = axes(t);
        imagesc(xaxis,zaxis(1:size(YM,1)),abs(YM));
        
        colormap(jet)
        caxis([0 100])
        colorbar;
        title({"Young's modulus";['Frame: ',num2str(k),' StressDepth: ',num2str(depthValue)]})
        xlabel('Distance (mm)')
        ylabel('Depth (mm)');
        
        YMLeftPhantom = mean(rmoutliers(YM(leftPhantom)),'omitnan')
        YMRightPhantom = mean(rmoutliers(YM(rightPhantom)),'omitnan')
        hold on;
        h1 = drawpolygon('Position',leftCoords,'InteractionsAllowed','none','FaceAlpha',0);
        h1.Color = 'yellow'
        h2 = drawpolygon('Position',rightCoords,'InteractionsAllowed','none','FaceAlpha',0);
        h2.Color = 'yellow'
        
        hold off;
        ax2 = axes(t)
        xlim(ax2,[1 length(xaxis)]);
        ylim(ax2,[1 size(IQData,1)]);
        ax2.YDir = 'reverse'
        ax2.YAxis.Visible = 'off'
        ax2.XAxis.Visible = 'off';
        ax2.Color = 'none';
        hline(depthValue,'w')
        text(ax2,xCentroidsLeft,yCentroidsLeft,num2str(YMLeftPhantom),'FontSize',20, ...
            'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
        text(ax2,xCentroidsRight,yCentroidsRight,num2str(YMRightPhantom),'FontSize',20, ...
            'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
        export_fig([baseFolder,filesep,'ParameterTesting',filesep,'YMFrame',num2str(k),'Depth',num2str(depthValue),'.png'],'-native')
        
        close all force;
        
       figure;
        t = tiledlayout(1,1);
        ax1 = axes(t);
        imagesc(xaxis,zaxis(1:size(YM,1)),mat2gray(abs(IQData(1:size(YM,1),:,k))));
        colormap(gray)
        caxis([0 0.1])
        colorbar;
        title({"Young's modulus";['Frame: ',num2str(k),' StressDepth: ',num2str(depthValue)]})
        xlabel('Distance (mm)')
        ylabel('Depth (mm)');
        
        YMLeftPhantom = mean(rmoutliers(YM(leftPhantom)),'omitnan')
        YMRightPhantom = mean(rmoutliers(YM(rightPhantom)),'omitnan')
        hold on;
        h1 = drawpolygon('Position',leftCoords,'InteractionsAllowed','none','FaceAlpha',0);
        h1.Color = 'yellow'
        h2 = drawpolygon('Position',rightCoords,'InteractionsAllowed','none','FaceAlpha',0);
        h2.Color = 'yellow'
        hold off;
        ax2 = axes(t)
        xlim(ax2,[1 length(xaxis)]);
        ylim(ax2,[1 size(IQData,1)]);
        ax2.YDir = 'reverse'
        ax2.YAxis.Visible = 'off'
        ax2.XAxis.Visible = 'off';
        ax2.Color = 'none';
        text(ax2,xCentroidsLeft,yCentroidsLeft,num2str(YMLeftPhantom),'FontSize',20, ...
            'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
        text(ax2,xCentroidsRight,yCentroidsRight,num2str(YMRightPhantom),'FontSize',20, ...
            'HorizontalAlignment','center','BackgroundColor','black','Color','white','Margin',1);
        hold on;
                hline(depthValue,'g')
                hold off; 

                export_fig([baseFolder,filesep,'ParameterTesting',filesep,'BScanFrame',num2str(k),'Depth',num2str(depthValue),'.png'],'-native')

                close all force
                
                figure;
                plot(xaxis,stressMM)
                axis padded;
                xlabel('Distance (mm)')
                ylabel('Stress')
                title(['Stress plot at depth: ',num2str(depthValue)])
                export_fig([baseFolder,filesep,'ParameterTesting',filesep,'StressPlotFrame',num2str(k),'Depth',num2str(depthValue),'.png'],'-native')
close all force;
leftYM(:,depthValue,k) = YMLeftPhantom;
rightYM(:,depthValue,k) = YMRightPhantom;
    end
end

save([baseFolder,filesep,'leftAndRightYM.mat'],'leftYM','rightYM')