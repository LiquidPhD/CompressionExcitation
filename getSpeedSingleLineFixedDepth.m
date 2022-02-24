function [Velocity] = getSpeedSingleLineFixedDepth(displacement,Parameters,sideIndexSelector,depth)
        %% Data loading.
        % load backward; % [ADJUSTABLE]
        disp_z = displacement;  % The first frame is normally very noisy.
        M = Parameters.M; N = Parameters.N;
        % Basic parameters.
        density = 1000;                   %(kg/m^3) Mass density of the medium.       [ADJUSTABLE]
        [Nz,Nx,Nt]= size(disp_z);         % The dimensions of data in the z axis, the x axis and time.
        zaxis = linspace((M/2)*Parameters.delta_z,(Nz-M/2)*Parameters.delta_z,Nz)*1e3;        % (mm) Aixial axis.
        xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    % (mm) Lateral axis.
        taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);                          % (s) Time axis.
        
     
        disp_z_smooth_desamp = disp_z;
        zaxis_desamp = zaxis;
      
        % 1.3: Data interpolation in time domain.
        % Motivation: The 2018-07-04 signal processing project shows that temporal
        %             interpolation can eliminate stairs of the TOF profiles.
        N_interp = 5;       % Ratio of interpolation. [ADJUSTABLE]
        taxis_interp = linspace(min(taxis),max(taxis),Nt*N_interp);
        delta_t_interp = min(diff(taxis_interp));
        tic();
        disp_z_smooth_desamp_interp = zeros([size(disp_z_smooth_desamp,1) size(disp_z_smooth_desamp,2) N_interp*size(disp_z_smooth_desamp,3)]);
        % for ii = 1:Nz_desamp
%         clearvars disp_z_smooth_desamp_interp
        for ii = 1:size(disp_z_smooth_desamp,1)
            signal = squeeze(disp_z_smooth_desamp(ii,:,:));
            signal_interp = interp1(taxis,signal',taxis_interp,'spline');
            disp_z_smooth_desamp_interp(ii,:,:) = signal_interp';
        end
        toc();
    
    
    for sideIndex = sideIndexSelector  % :2
        loopExit=0
%          roi(2) = 290
            roi(4) = 25;
        while loopExit == 0
            if sideIndex == 1
%                 figure; imagesc(abs(newDispLEFT(:,:,10)))
% playWaveVideo(abs(disp_z_smooth_desamp_interp(:,:,newWorkingIndices)),0.01)
                % figure; h = imagesc(abs(clippedBScan(:,:,1))); alpha = abs(newDispLEFT(:,:,10)); set(h,'alphaData',alpha*0.9); colormap(jet)
%                 [~,depth] = ginput(1);
                roi(2) = round(depth);
                close all force;
                disp2 = abs(disp_z_smooth_desamp_interp);
                roi(1) = 1
                roi(3) = round(size(disp_z_smooth_desamp_interp,2)/2);
                
            else
%                 figure; imagesc(abs(newDispRIGHT(:,:,10)))
                % figure; h = imagesc(abs(clippedBScan(:,:,1))); alpha = abs(newDispRIGHT(:,:,10)); set(h,'alphaData',alpha*0.9); colormap(jet)
%                 [~,depth] = ginput(1);
                roi(2) = round(depth);
%                 disp2 = abs(newDispRIGHT);
                disp2 = abs(disp_z_smooth_desamp_interp);

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
                spacetime(:,k) = mean(squeeze(disp2(yRegion-50:yRegion+50,xRegion,k)));
            end
            if roiScaled(1) <0
                spacetime = flipud(spacetime)
            end
            
            figure; imagesc(spacetime)
            % if sideIndex == 1
            spacetime = flip(spacetime,1);
            
            Spacetime = medfilt2(spacetime,[2 2],'Symmetric');
            Spacetime = Spacetime;
            try
                if sideIndex == 1
            pos1 = evalin('base','pos1')
             figure('WindowState','maximized');
            himage= imagesc(Spacetime)
            h = impoly(gca,pos1);
            PosTime = wait(h);
            pos1 = getPosition(h);
            title('Select time region to be calculated');
            assignin('base','pos1',pos1)
                else
                    pos2 = evalin('base','pos2')
             figure('WindowState','maximized');
            himage= imagesc(Spacetime)
            h = impoly(gca,pos2);
            PosTime = wait(h);
            pos2 = getPosition(h);
            title('Select time region to be calculated');
            assignin('base','pos2',pos2)
                end
            catch
             figure('WindowState','maximized');
            himage= imagesc(Spacetime)
            h = impoly(gca);
            PosTime = wait(h);
            
            if sideIndex == 1
                pos1 = getPosition(h);
            title('Select time region to be calculated');
            assignin('base','pos1',pos1)
            else
                pos2 = getPosition(h);
            title('Select time region to be calculated');
                assignin('base','pos2',pos2)
            end
            end
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
            
            
            waveLocation = Xdata(5:end-4);
            timeLocation = Tdata(5:end-4);
            figure; imagesc(Spacetime_Mask); hold on; plot(Tdata(5:end-4),Xdata(5:end-4),'r*'); hold off;
            
            waveLocationExistingValues = waveLocation(~isnan(waveLocation))
            velocityFit = fit(taxis_interp(timeLocation)',-1*xaxis(waveLocationExistingValues)','poly1');
                Velocity = abs(velocityFit.p1/1000)
          if sideIndex == 1
              VelocityLeft = Velocity;
          else
              VelocityRight = Velocity;
          end
            
            
            [xData,yData] = prepareCurveData([],waveLocationExistingValues)
            lengthXData = length(xData)
            waveLocationExistingValues(waveLocationExistingValues>length(windowXAxis)) = []
            yData = windowXAxis(waveLocationExistingValues)
            for k = lengthXData+1:length(waveLocation)
                xData(k) = NaN;
                yData(k) = NaN;
            end
            plot(taxis_interp(1:length(yData)),yData)
            close all force;

            velocityFit = fit(taxis_interp(timeLocation)',...
                windowXAxis(waveLocationExistingValues)','poly1');

                title(['Velocity: ',num2str(Velocity),' m/s'])
           

  figure;

                    imagesc(taxis_interp,xaxis(xRegion),spacetime)
                    hold on;
                    plot(taxis_interp(timeLocation),yData,'r*')
                    plot(velocityFit,yData,taxis_interp(timeLocation))
                    hold off;
              ylabel('Distance (mm)')
              xlabel('Time (s)')
                title(['Velocity: ',num2str(Velocity),' m/s'])
%                 export_fig([folders(folderIndex).folder,filesep,'DistanceVsTimeFitting.png'])

            figure;
          
            uitable('Data',Velocity')
            
%             folders(folderIndex).folder
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
end