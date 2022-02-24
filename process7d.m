clearvars;
% filesToProcess =
% folders = rdir('G:\ChickenBreast\*\**\*Param*')
folders = rdir('C:\Users\Verasonics\Documents\Vantage-4.3.0-2009141600\7d_20210422\*\**\*displacementProcessing*')
% folders = rdir('D:\Justin\StaticDynamic040821\GelatinSilicone\Dynamic\*\**\*Param*')
% depth = 25;
% folders = rdir('C:\Users\Verasonics\Documents\Vantage-4.3.0-2009141600\dynamic7d_8PercentGelatinPhantom_Aquaflex_25-Mar-2021_73302\*Param*')
for folderIndex = 1:length(folders)
   load(folders(folderIndex).name)
    
    [Nz,Nx,Nt]= size(displacement);         % The dimensions of data in the z axis, the x axis and time.
    zaxis = linspace(0,(Nz-1)*delta_z,Nz)*1e3;                      %(mm) Aixial axis.
    xaxis = linspace(-(Nx-1)/2*delta_x,(Nx-1)/2*delta_x,Nx)*1e3;    %(mm) Lateral axis.
    taxis = linspace(0,(Nt-1)*delta_t,Nt);
    
%     filename = [folders(folderIndex).folder,filesep,'threeViews.gif'];
%     pauseAmount = 0.2;
%     figure;
%     for gifIndex = 1:size(displacement,3)
%         subplot(3,1,1)
%         imagesc(xaxis,zaxis,im2uint8((displacement(:,:,gifIndex))))
%         subplot(3,1,2)
%         imagesc(xaxis,zaxis,im2uint8(abs(newDispLEFT(:,:,gifIndex))))
%         subplot(3,1,3)
%         imagesc(xaxis,zaxis,im2uint8(abs(newDispRIGHT(:,:,gifIndex))))
%         colormap(fireice)
%         pause(0.2)
%         if gifIndex == 1
%             gif(filename,'DelayTime',pauseAmount)
%         else
%             gif
%         end
%     end
    
    % figure; imagesc(abs(squeeze(newDispLEFT(275,:,:))))
    % colormap(fireice)
    % figure; imagesc(abs(squeeze(newDispRIGHT(275,:,:))))
    % colormap(fireice)
    
%     figure; imagesc(displacement(:,:,50))
    
    for sideIndex = 1% :2
        loopExit=0
%          roi(2) = 290
            roi(4) = 25;
        while loopExit == 0
            if sideIndex == 1
%                 figure; imagesc(abs(newDispLEFT(:,:,10)))
playWaveVideo(abs(newDispLEFT))
                % figure; h = imagesc(abs(clippedBScan(:,:,1))); alpha = abs(newDispLEFT(:,:,10)); set(h,'alphaData',alpha*0.9); colormap(jet)
                [~,depth] = ginput(1);
                roi(2) = round(depth);
                close all force;
                disp2 = abs(newDispLEFT);
                roi(1) = 1
                roi(3) = round(size(disp2,2)/2);
                
            else
                figure; imagesc(abs(newDispRIGHT(:,:,10)))
                % figure; h = imagesc(abs(clippedBScan(:,:,1))); alpha = abs(newDispRIGHT(:,:,10)); set(h,'alphaData',alpha*0.9); colormap(jet)
                [~,depth] = ginput(1);
                roi(2) = round(depth);
                disp2 = abs(newDispRIGHT);
                roi(1) = round(size(disp2,2)/2)+1;
                roi(3) = round(size(disp2,2)/2)-1
                
            end
            xRegion = round(roi(1)):round(roi(1)+roi(3))-1;
            yRegion = round(roi(2)):round(roi(2)+roi(4));
            % yRegion = 380-15:380+15;
            
            windowXAxis = xaxis(xRegion)
            % ROI Scaling
            xScale = xaxis(2)-xaxis(1) % mm
            zScale = zaxis(2)-zaxis(1) % mm
            roiScaledZ = roi(2)*zScale
            %     roiScaledX = roi(1)*xScale
            roiScaledX = xaxis(round(roi(1))) % Because 0 is in middle of scale
            roiXDistance = roi(3)*xScale;
            roiZDistance = roi(4)*zScale;
            roiScaled = [roiScaledX roiScaledZ roiXDistance roiZDistance]
            
            
            clearvars spacetime
            % Step 2: Get waves in that region and convert to displacement/time map
            for k = 1:size(disp2,3)
                spacetime(:,k) = mean(squeeze(disp2(yRegion,xRegion,k)));
            end
            if roiScaled(1) <0
                spacetime = flipud(spacetime)
            end
            
            figure; imagesc(spacetime)
            % if sideIndex == 1
            spacetime = flip(spacetime,1);
            
            Spacetime = medfilt2(spacetime,[2 2],'Symmetric');
            Spacetime = Spacetime;
            figure;
            himage= imagesc(Spacetime)
            h = impoly(gca);
            PosTime = wait(h);
            pos = getPosition(h);
            title('Select time region to be calculated');
            BW = createMask(h,himage);
            
            Line_t = sum(BW,2);
            idx_t1 = find(Line_t==0);
            Jump = find(diff(idx_t1)>2);
            t_ini = idx_t1(Jump(1));
            t_end = idx_t1(Jump(1)+1);
            
            
            Spacetime_Mask = Spacetime.*BW;
            count = 1
            clearvars Tdata Xdata
            for i = 1:size(Spacetime_Mask,2)
                [val,idx] = max(Spacetime_Mask(:,i));
                if val~=0 & (isnan(val)==0) & (isempty(idx)==0)
                    Tdata(count) = i;
                    Xdata(count) = idx;
                    count = count+1;
                end
            end
            %
            %
            % count = 1;
            % clearvars xCount tCount xData tData
            % for i = 1: size(Spacetime_Mask,2)
            %     [val,idx] = max(Spacetime_Mask(:,i));
            %     if val~=0 & (isnan(val)==0) & (isempty(idx)==0)
            %         Xdata(count) = Xaxis(i);
            %         xCount(count) = i;
            %         Tdata(count) = Time1(idx);
            %         tCount(count) = idx;
            %         count = count +1;
            %     end
            % end
            
            waveLocation = Xdata(5:end-4);
            timeLocation = Tdata(5:end-4);
            figure; imagesc(Spacetime_Mask); hold on; plot(Tdata(5:end-4),Xdata(5:end-4),'r*'); hold off;
            
            waveLocationExistingValues = waveLocation(~isnan(waveLocation))
            velocityFit = fit(taxis(timeLocation)',-1*xaxis(waveLocationExistingValues)','poly1');
            if sideIndex == 1
                VelocityLeft(:,folderIndex) = abs(velocityFit.p1/1000)
            else
                VelocityRight(:,folderIndex) = abs(velocityFit.p1/1000)
                
            end
            
            
            [xData,yData] = prepareCurveData([],waveLocationExistingValues)
            lengthXData = length(xData)
            waveLocationExistingValues(waveLocationExistingValues>length(windowXAxis)) = []
            yData = windowXAxis(waveLocationExistingValues)
            for k = lengthXData+1:length(waveLocation)
                xData(k) = NaN;
                yData(k) = NaN;
            end
            plot(taxis(1:length(yData)),yData)
            close all force;
            
%             %% Without image showing location taken
%             figure;
%             tiledlayout(1,1);
%             %             tiledlayout(2,1);
%             
%             % Tile 1
%             ax1 = nexttile
%             %             %         imagesc(xaxis,zaxis,disp2(:,:,k));
%             %             imagesc(xaxis,zaxis,mat2gray(abs(clippedBScan(:,:,1))))
%             %             hold on;
%             %             rectangle('Position',roiScaled,'EdgeColor','r')
%             %             hold off;
%             %             colormap(ax1,gray)
%             %             caxis([0 0.2])
%             %             title('Selected Region')
%             %             xlabel('Distance (mm)')
%             %             ylabel('Distance (mm)')
%             %             % Tile 2
%             %             ax2 =  nexttile
%             %         imagesc(roiScaledX:roiScaledX+roiXDistance,roiScaledZ:roiScaledZ+roiZDistance,plottableStuff)
%             %         imagesc(taxis,xaxis(xRegion),plottableStuff)
%             %         hold on;
%             % %         plot(xaxis(waveLocationExistingValues),'b*');
%             %         plot(taxis,yData,'r*')
%             % %         plot(velocityFit)
%             %         plot(velocityFit,yData,taxis)
%             % % plot(velocityFit,taxisinterpolated,xaxis(xRegion))
%             %         hold off;
            velocityFit = fit(taxis(timeLocation)',...
                windowXAxis(waveLocationExistingValues)','poly1');
%             imagesc(taxis,xaxis(xRegion),Spacetime);
%             xl = xlim
%             yl = ylim
%             hold on;
%             plot(taxis(timeLocation),yData,'r*')
%             ylabel('Distance (mm)')
%             yyaxis right;
%             plot(velocityFit)
%             % if roiScaled(1)<0
%             ylim(yl)
%             set(gca,'Ydir','reverse')
%             % end
%             legend('hide')
%             ax = gca;
%             ax.YAxis(2).Visible = 'off'
%             xlabel('Time (s)')
%             
%             
%             
%             %             colormap(ax2,parula)
%             colormap(ax1,parula)
            if sideIndex == 1
                title(['Velocity: ',num2str(VelocityLeft(:,folderIndex)),' m/s'])
            else
                title(['Velocity: ',num2str(VelocityRight(:,folderIndex)),' m/s'])
            end
% %             titleToUse = regexp(folder,filesep)
% %             digit = regexp(folder(titleToUse(end)+1:end),'\d')
% %             if digit(2) == 2
% %                 numbers = folder(titleToUse(end)+1:titleToUse(end)+2)
% %             else
% %                 numbers = folder(titleToUse(end)+1)
% %             end
% %             phantomNum = regexp(folder,'Num')
% %             titleToUse = [numbers,'% Gelatin Phantom #',folder(phantomNum+3)];
% With location info
  figure;
%             tiledlayout(1,1);
%                         tiledlayout(2,1);
            
            % Tile 1
%             ax1 = nexttile
                        %         imagesc(xaxis,zaxis,disp2(:,:,k));
%                         imagesc(xaxis,zaxis,mat2gray(abs(clippedBScan(:,:,1))))
%                         hold on;
%                         rectangle('Position',roiScaled,'EdgeColor','r')
%                         hold off;
%                         colormap(ax1,gray)
%                         caxis([0 0.2])
%                         title('Selected Region')
%                         xlabel('Distance (mm)')
%                         ylabel('Distance (mm)')
%                         % Tile 2
%                         ax2 =  nexttile
%                     imagesc(roiScaledX:roiScaledX+roiXDistance,roiScaledZ:roiScaledZ+roiZDistance,spacetime)
                    imagesc(taxis,xaxis(xRegion),spacetime)
                    hold on;
            %         plot(xaxis(waveLocationExistingValues),'b*');
                    plot(taxis(timeLocation),yData,'r*')
            %         plot(velocityFit)
                    plot(velocityFit,yData,taxis(timeLocation))
            % plot(velocityFit,taxisinterpolated,xaxis(xRegion))
                    hold off;
%             velocityFit = fit(taxis(timeLocation)',...
%                 windowXAxis(waveLocationExistingValues)','poly1');
%             imagesc(taxis,xaxis(xRegion),Spacetime);
%             xl = xlim
%             yl = ylim
%             hold on;
%             plot(taxis(timeLocation),yData,'r*')
%             ylabel('Distance (mm)')
%             yyaxis right;
%             plot(velocityFit)
%             % if roiScaled(1)<0
%             ylim(yl)
%             set(gca,'Ydir','reverse')
%             % end
%             legend('hide')
%             ax = gca;
%             ax.YAxis(2).Visible = 'off'
%             xlabel('Time (s)')
            
            
            
%                         colormap(ax2,parula)
%             colormap(ax1,gray)
            if sideIndex == 1
                title(['Velocity: ',num2str(VelocityLeft(:,folderIndex)),' m/s'])
            else
                title(['Velocity: ',num2str(VelocityRight(:,folderIndex)),' m/s'])
            end
            folderName = folders(folderIndex).folder;
            percentFind = regexp(folderName,'pct')
            underscoreFind = regexp(folderName,'7d_')
            sgtitle([folderName(underscoreFind+3:percentFind-1),'% Gelatin phantom'])
            
            %         pause
            
            if sideIndex == 1
                export_fig([folders(folderIndex).folder,filesep,'DistanceVsTimeLeft.png'])
            else
                export_fig([folders(folderIndex).folder,filesep,'DistanceVsTimeRight.png'])
            end
            figure;
            uitable('Data',VelocityLeft')
            folders(folderIndex).folder
            answer = questdlg('Good?')
            switch answer
                case 'Yes'
                    loopExit = 1;
                case 'No'
                case 'cancel'
                    break
            end
        end
    end
%     save([folders(folderIndex).folder,filesep,'workspace.mat'])
end
save('C:\Users\Verasonics\Documents\Vantage-4.3.0-2009141600\7d_20210422\Velocities.mat','VelocityLeft')
% save('D:\Justin\StaticDynamic040521\Dynamic\SpeedsAt75pixelsindepth.mat','VelocityLeft')