clearvars;
% folders1 = rdir(['G:\GelatinPhantoms092420\USE\Used\*\**\*Param*'])
% folders2 = rdir(['D:\Ultrasound112720\*\**\*Param*'])
% folders3 = rdir(['D:\UltrasoundData\SKIN_ELASTOGRAPHY\*\**\*Param*'])
% folders4 = rdir(['D:\RippySkin\*\**\*Param*'])
% folders = [folders1; folders2; folders3; folders4];
% folders = rdir(['D:\RippySkin\*\**\*Param*'])
% folders = rdir(['D:\GelatinPhantoms20210122\HalfAndHalf\*\**\*wavelengthMapWorkspace*'])
folders = rdir(['D:\comboElastography\*\**\*Param*'])
progressbar(0)
for folderIndex = 1:length(folders)
    
    clearvars -except folders folderIndex
    load([folders(folderIndex).folder,filesep,'wavelengthMapWorkspace.mat'],'displacement','IQData','delta_z','xaxis','zaxis','taxis','wavelength')
    load([folders(folderIndex).folder,filesep,'Parameter.mat'])
    
    %% Parameters for speed mapping using region
    
    xAxis = Parameter.Trans.ElementPos(:,1);
    if length(Parameter.PData) == 3
        delta_z = Parameter.PData(3).PDelta(3)*Parameter.Trans.lambda;
        delta_x = Parameter.PData(3).PDelta(1)*Parameter.Trans.lambda;
        AxialResolution = Parameter.PData(3).PDelta(3);
    else
        delta_z = Parameter.PData.PDelta(3)*Parameter.Trans.lambda;
        delta_x = Parameter.PData.PDelta(1)*Parameter.Trans.lambda;
        AxialResolution = Parameter.PData.PDelta(3);
    end
    delta_t = 100*1e-6;
    
    topClip = 1;
    bottomClip = 600;
    AxialLengthToTest = 5;
    M = ceil(AxialLengthToTest/AxialResolution);
    N = 3;
    fc = Parameter.Trans.frequency*1e6
    c=1500;
    IQDataClip = IQData(topClip:bottomClip,:,:);
    
    ROI_dim_z = size(IQDataClip,1)-M;
    ROI_dim_x = size(IQDataClip,2);
    ROI_dim_t = size(IQDataClip,3)-N;
    xaxisleft = xaxis(1:length(xaxis)/2);
    xaxisright = xaxis(length(xaxis)/2+1:end);
    
    parfor imageDepth = 1:size(displacement,1)
        TEST = squeeze(displacement(imageDepth,:,:));
        Y = fft2(TEST);
        %         FFTFILTERBOT = imresize(FILT,[size(Y,1) size(Y,2)])
        %         FFTFILTERBOT(FFTFILTERBOT<0.1) = 0;
        %         FFTFILTERBOT(:,round(size(FFTFILTERBOT,2)/2):end) = 0;
        FFTFILTERRIGHT = ones([size(Y,1) size(Y,2)]);
        FFTFILTERRIGHT(1:round(size(FFTFILTERRIGHT,1)/2),1:round(size(FFTFILTERRIGHT,2)/2)) = 0;
        FFTFILTERRIGHT(round(size(FFTFILTERRIGHT,1)/2)+1:end,round(size(FFTFILTERRIGHT,2)/2)+1:end) = 0;
        FFTFILTERLEFT = fliplr(FFTFILTERRIGHT);
        NEWTOP = FFTFILTERLEFT .* Y;
        NEWBOT = FFTFILTERRIGHT .* Y;
        newleftifft = ifft2(ifftshift(NEWTOP));
        newrightifft = ifft2(ifftshift(NEWBOT));
        newDispLEFT(imageDepth,:,:) = newleftifft;
        newDispRIGHT(imageDepth,:,:) = newrightifft;
    end
    %%
    
    
    continueFlag = 0;
    locationFlag = 0;
    locationNumber = 1;
    
    dispLeft = abs(displacement(:,1:round(size(displacement,2)/2),:));
    BScanLeft = fliplr(abs(IQData(:,1:round(size(displacement,2)/2),:)));
    dispLeft = fliplr(dispLeft);
    dispRight = abs(displacement(:,round(size(displacement,2)/2)+1:end,:));
    BScanRight =abs(IQData(:,round(size(displacement,2)/2)+1:end,:));
    
    while locationFlag == 0
        while continueFlag == 0
            
            
            close all force;
            playWaveVideo(displacement,0.05,jet)
            figure; imagesc(displacement(:,:,5));
            [x,y] = ginput(1);
            close all force;
            x = round(x);
            y = round(y);
            [frequency,FT_plot] = plot_spectrum(squeeze(displacement(y,x,:)),10000)
            
            answer = questdlg('Good?', ...
                'CHOOSE MORTAL', ...
                'Yes','No','No');
            % Handle response
            switch answer
                case 'Yes'
                    continueFlag = 1
                case 'No'
                    
            end
            
            
            
            answer = inputdlg('NumberOfPeaks')
            close all force;
            for peakSelection = 1:str2num(cell2mat(answer))
                [frequency,FT_plot,P1] = plot_spectrum(squeeze(displacement(y,x,:)),10000)
                
                [selectedFreq,~] = ginput(1);
                close all force;
                figure;
                subplot(2,1,1)
                imagesc(displacement(:,:,5))
                hold on;
                scatter(x,y,'r*'); hold off;
                subplot(2,1,2);
                plot(frequency,P1)
                hold on;
                vline(selectedFreq); hold off;
                title('Single-Sided Amplitude Spectrum of X(t)')
                xlabel('f (Hz)')
                ylabel('|P1(f)|')
                export_fig([folders(folderIndex).folder,filesep,'Location',num2str(locationNumber),'AndFFT',num2str(peakSelection),'.png'],'-png','-native')
                savefig([folders(folderIndex).folder,filesep,'Location',num2str(locationNumber),'AndFFT',num2str(peakSelection),'.fig'])
                close all force;
                figure; imagesc(xaxis,zaxis,abs(wavelength)*selectedFreq)
                colormap(jet)
                % caxis([1e-5 1e-4])
                if max(max(abs(wavelength)*selectedFreq)) > 20
                    caxis([0 20])
                else
                    caxis([0 max(max(abs(wavelength)*selectedFreq))])
                end
                % caxis([0 20])
                colorbar;
                export_fig([folders(folderIndex).folder,filesep,'WavelengthMapLocation',num2str(locationNumber),'Peak',num2str(peakSelection),'.png'],'-png','-native')
                savefig([folders(folderIndex).folder,filesep,'WavelengthMapLocation',num2str(locationNumber),'Peak',num2str(peakSelection),'.fig'])
                % saveas(gcf,[folders(folderIndex).folder,filesep,'WavelengthMap.fig']);
                % saveas(gcf,[folders(folderIndex).folder,filesep,'WavelengthMap.png']);
            end
            
            %% Speed original way
            %     figure;
            % imagesc(mat2gray(BScanLeft(:,:,1)));
            % caxis([0 0.7])
            % colormap(gray)
            % hold on; hline(y,'r:','Depth selected'); hold off
            %% LEFT SIDE
            try
                clearvars plottableStuff
            catch
            end
            for k = 1:size(dispLeft,3)
                plottableStuff(:,k,1) = abs(mean(squeeze(dispLeft(y:y+20,:,k))));
                %     plottableStuff(:,k,2) = abs(mean(squeeze(newDispRIGHT(250:280,:,k))));
            end
            Fs = diff(taxis);
            Fs = Fs(1);
            Fs = 1/Fs;
            
            Spacetime = plottableStuff;
            
            figure;
            himage= imagesc(Spacetime)
            h = impoly(gca);
            PosTime = wait(h);
            pos = getPosition(h);
            title('Select time region to be calculated');
            BW = createMask(h,himage);
            savedPosLeft = pos;
            Line_t = sum(BW,2);
            idx_t1 = find(Line_t==0);
            Jump = find(diff(idx_t1)>2);
            t_ini = idx_t1(Jump(1));
            t_end = idx_t1(Jump(1)+1);
            
            Spacetime_Mask = Spacetime.*BW;
            clearvars t_delay
            xPosFlag = 0
            xPosCounter = 1;
            while xPosFlag == 0
                ref_profile = detrend(squeeze(Spacetime_Mask(xPosCounter,:)));
                if mean(ref_profile) == 0
                    xPosCounter = xPosCounter+1;
                else
                    xPosFlag = 1;
                end
            end
            
            for xPos = xPosCounter:size(Spacetime_Mask,1)
                ref_profile = detrend(squeeze(Spacetime_Mask(xPosCounter,:)));
                xcorr_profile = detrend(squeeze(Spacetime_Mask(xPos,:)));
                x_corr = xcorr(ref_profile,xcorr_profile,'coeff');
                [~,max_lag] = max(x_corr);
                t_delay(xPos) = max_lag-2;
                %         disp(['Depth: ',num2str(depth),' Sample: ',num2str(xPos),'/',num2str(size(dispLeft,2))])
                
                
                figure(1);
                subplot(2,2,1);
                plot(ref_profile);
                subplot(2,2,2);
                plot(xcorr_profile);
                subplot(2,2,3);
                plot(x_corr);
                subplot(2,2,4);
                plot(t_delay(:));
                pause(0.01);
            end
            
            
            figure;
            himage= plot(t_delay);
            h = imrect(gca,[25,0,...
                50,max(t_delay)+20]);
            delayTime = wait(h);
            delayPos = getPosition(h);
            Cut_Lat_ini = round(delayTime(1));
            Cut_Lat_end = round(delayTime(1))+round(delayTime(3));
            
            
            deltaX = diff(xaxis);
            deltaX = deltaX(1);
            xaxis = (1:size(t_delay,2))*deltaX
            xaxis = xaxis(Cut_Lat_ini:Cut_Lat_end);
            t_delay_test = squeeze(t_delay(1,Cut_Lat_ini:Cut_Lat_end))
            deltaT = diff(taxis);
            deltaT = deltaT(1);
            t_delay_test = t_delay_test*deltaT;
            t_delay_test(t_delay_test<=0) = 0;
            t_delay_testNZ = nonzeros(t_delay_test);
            
            %% Fit: 'untitled fit 1'.
            % [xData, yData] = prepareCurveData( xaxis, t_delay_test );
            [xData, yData] = prepareCurveData(xaxis(1:length(t_delay_testNZ)),t_delay_testNZ);
            
            
            
            % Set up fittype and options.
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            opts.Robust = 'Bisquare';
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            
            % Plot fit with data.
            figure( 'Name', 'untitled fit 1' );
            h = plot( fitresult, xData, yData );
            legend( h, 't_delay_test vs. xaxis', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
            % Label axes
            xlabel( 'xaxis', 'Interpreter', 'none' );
            ylabel( 't_delay_test', 'Interpreter', 'none' );
            grid on
            
            velocitiesLeft(locationNumber) = abs(1/fitresult.p1)/1000
            
            figure;
            imagesc(xaxisleft,zaxis,mat2gray(BScanLeft(50:end,:,1)));
            caxis([0 0.7])
            colormap(gray)
            hold on; hline(zaxis(y),'g','Selected Depth'); hold off
            title(['Speed at depth: ',num2str(abs(1/fitresult.p1)/1000,3),' m/s'])
            export_fig([folders(folderIndex).folder,filesep,'SpeedAtLeftLocation',num2str(locationNumber),'Peak.png'],'-png','-native')
            savefig([folders(folderIndex).folder,filesep,'SpeedAtLeftLocation',num2str(locationNumber),'Peak.fig'])
            %% RIGHT SIDE
            clearvars plottableStuff
            for k = 1:size(dispRight,3)
                plottableStuff(:,k,1) = abs(mean(squeeze(dispRight(y:y+20,:,k))));
                %     plottableStuff(:,k,2) = abs(mean(squeeze(newDispRIGHT(250:280,:,k))));
            end
            Fs = diff(taxis);
            Fs = Fs(1);
            Fs = 1/Fs;
            
            Spacetime = plottableStuff;
            
            figure;
            himage= imagesc(Spacetime)
            h = impoly(gca);
            PosTime = wait(h);
            pos = getPosition(h);
            title('Select time region to be calculated');
            BW = createMask(h,himage);
            savedPosLeft = pos;
            Line_t = sum(BW,2);
            idx_t1 = find(Line_t==0);
            Jump = find(diff(idx_t1)>2);
            t_ini = idx_t1(Jump(1));
            t_end = idx_t1(Jump(1)+1);
            
            Spacetime_Mask = Spacetime.*BW;
            clearvars t_delay
            xPosFlag = 0
            xPosCounter = 1;
            while xPosFlag == 0
                ref_profile = detrend(squeeze(Spacetime_Mask(xPosCounter,:)));
                if mean(ref_profile) == 0
                    xPosCounter = xPosCounter+1;
                else
                    xPosFlag = 1;
                end
            end
            
            for xPos = xPosCounter:size(Spacetime_Mask,1)
                ref_profile = detrend(squeeze(Spacetime_Mask(xPosCounter,:)));
                xcorr_profile = detrend(squeeze(Spacetime_Mask(xPos,:)));
                x_corr = xcorr(ref_profile,xcorr_profile,'coeff');
                [~,max_lag] = max(x_corr);
                t_delay(xPos) = max_lag-2;
                %         disp(['Depth: ',num2str(depth),' Sample: ',num2str(xPos),'/',num2str(size(dispLeft,2))])
                
                
                figure(1);
                subplot(2,2,1);
                plot(ref_profile);
                subplot(2,2,2);
                plot(xcorr_profile);
                subplot(2,2,3);
                plot(x_corr);
                subplot(2,2,4);
                plot(t_delay(:));
                pause(0.01);
            end
            
            
            figure;
            himage= plot(t_delay);
            h = imrect(gca,[25,0,...
                50,max(t_delay)+20]);
            delayTime = wait(h);
            delayPos = getPosition(h);
            Cut_Lat_ini = round(delayTime(1));
            Cut_Lat_end = round(delayTime(1))+round(delayTime(3));
            
            
            deltaX = diff(xaxis);
            deltaX = deltaX(1);
            xaxis = (1:size(t_delay,2))*deltaX
            xaxis = xaxis(Cut_Lat_ini:Cut_Lat_end);
            t_delay_test = squeeze(t_delay(1,Cut_Lat_ini:Cut_Lat_end))
            deltaT = diff(taxis);
            deltaT = deltaT(1);
            t_delay_test = t_delay_test*deltaT;
            t_delay_test(t_delay_test<=0) = 0;
            t_delay_testNZ = nonzeros(t_delay_test);
            
            %% Fit: 'untitled fit 1'.
            % [xData, yData] = prepareCurveData( xaxis, t_delay_test );
            [xData, yData] = prepareCurveData(xaxis(1:length(t_delay_testNZ)),t_delay_testNZ);
            
            
            
            % Set up fittype and options.
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            opts.Robust = 'Bisquare';
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            
            % Plot fit with data.
            figure( 'Name', 'untitled fit 1' );
            h = plot( fitresult, xData, yData );
            legend( h, 't_delay_test vs. xaxis', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
            % Label axes
            xlabel( 'xaxis', 'Interpreter', 'none' );
            ylabel( 't_delay_test', 'Interpreter', 'none' );
            grid on
            
            velocitiesRight(locationNumber) = abs(1/fitresult.p1)/1000
            
            figure;
            imagesc(xaxisright,zaxis,mat2gray(BScanRight(50:end,:,1)));
            caxis([0 0.7])
            colormap(gray)
            hold on; hline(zaxis(y),'g','Selected Depth'); hold off
            title(['Speed at depth: ',num2str(abs(1/fitresult.p1)/1000,3),' m/s'])
            export_fig([folders(folderIndex).folder,filesep,'SpeedAtRightLocation',num2str(locationNumber),'Peak.png'],'-png','-native')
            savefig([folders(folderIndex).folder,filesep,'SpeedAtRightLocation',num2str(locationNumber),'Peak.fig'])
            answer2 = questdlg('Another Location?', ...
                'CHOOSE MORTAL', ...
                'Yes','No','No');
            % Handle response
            switch answer2
                case 'Yes'
                    locationNumber = locationNumber+1;
                    continueFlag = 0;
                case 'No'
                    locationFlag = 1;
            end
            
        end
        progressbar(folderIndex/length(folders))
        
    end
end