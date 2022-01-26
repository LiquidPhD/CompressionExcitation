debugMode = 0;

disp_z = Loupas_phase_shift(:,:,2:90);  % The first frame is normally very noisy.
        
        % Basic parameters.
        density = 1000;                   %(kg/m^3) Mass density of the medium.       [ADJUSTABLE]
        [Nz,Nx,Nt]= size(disp_z);         % The dimensions of data in the z axis, the x axis and time.
        zaxis = linspace((Parameters.M/2)*Parameters.delta_z,(Nz-Parameters.M/2)*Parameters.delta_z,Nz)*1e3;        % (mm) Aixial axis.
        xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    % (mm) Lateral axis.
        taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);                          % (s) Time axis.
        
        % % Checking.
        % disp_min = min(min(min(disp_z))); disp_max = max(max(max(disp_z)));
        % figure
        % for ii = 1:Nt
        %     imagesc(disp_z(:,:,ii),[disp_min disp_max]/150)
        %     title(ii)
        %     pause(0.15)
        % 	colormap('jet')
        % end
        
        %% I: PRE-PROCESSING
        
        % 1.1: Data smoothing using 'sgolay' filter.
        disp_z_smooth = zeros(Nz,Nx,Nt);
        tic();
        for ii=1:Nt
            data = disp_z(:,:,ii);
            data_s = sgolayfilt(data,1,21,[],2);     % [ADJUSTABLE]
            disp_z_smooth(:,:,ii) = data_s;
        end
        toc();
        
        % % Checking.
        % disp_min = min(min(min(disp_z_smooth))); disp_max = max(max(max(disp_z_smooth)));
        % figure
        % for ii = 1:Nt
        %     imagesc(disp_z_smooth(:,:,ii),[disp_min disp_max]/50)
        %     title(ii)
        %     pause(0.25)
        % 	colormap('jet')
        % end
        
        % % 1.2: Axial desampling.
        % % Motivation: Normally, the axial resolution is much higher than the
        % % lateral resolution. For elasticity map, it is not nessary to have a high
        % % axial/lateral resolution ratio.
        % ratio_resolution = 1; % Define a ratio of axial/lateral resolution. [ADJUSTABLE]
        % N_desamp = ceil((delta_x/delta_z)/ratio_resolution);
        % Nz_desamp = ceil(Nz/N_desamp);
        % disp_z_smooth_desamp = zeros(Nz_desamp,Nx,Nt);
        % zaxis_desamp = zeros(1,Nz_desamp);
        % tic();
        % for ii = 1:Nz_desamp
        %     index = (ii-1)*N_desamp+1;
        %     disp_z_smooth_desamp(ii,:,:) = disp_z_smooth(index,:,:);
        %     zaxis_desamp(ii) = zaxis(index);
        % end
        % toc();
        disp_z_smooth_desamp = disp_z_smooth;
        zaxis_desamp = zaxis;
        % % Checking.
        % disp_min = min(min(min(disp_z_smooth_desamp)));
        % disp_max = max(max(max(disp_z_smooth_desamp)));
        % figure
        % for ii = 1:Nt
        %     imagesc(disp_z_smooth_desamp(:,:,ii),[disp_min disp_max]/10)
        %     title(ii)
        %     pause(0.1)
        % 	colormap('jet')
        % end
        
        % 1.3: Data interpolation in time domain.
        % Motivation: The 2018-07-04 signal processing project shows that temporal
        %             interpolation can eliminate stairs of the TOF profiles.
        N_interp = 5;       % Ratio of interpolation. [ADJUSTABLE]
        taxis_interp = linspace(min(taxis),max(taxis),Nt*N_interp);
        delta_t_interp = min(diff(taxis_interp));
        tic();
        disp_z_smooth_desamp_interp = zeros([size(disp_z_smooth_desamp,1) size(disp_z_smooth_desamp,2) size(disp_z_smooth_desamp,3).*N_interp]);
        % for ii = 1:Nz_desamp
        for ii = 1:size(disp_z_smooth_desamp,1)
            signal = squeeze(disp_z_smooth_desamp(ii,:,:));
            signal_interp = interp1(taxis,signal',taxis_interp,'spline');
            disp_z_smooth_desamp_interp(ii,:,:) = signal_interp';
        end
        toc();
       
        
        
        
        % TOF = zeros(Nz_desamp,Nx);
        TOF = zeros(size(disp_z_smooth_desamp_interp,1),Nx)
        N_radius = 25; % Half distance for choosing reference waveform. [ADJUSTABLE]
        tic();
        
%         if folderIndex == 1 %|| ~exist('b','var')
%             % for ii = 1:Nz_desamp
% %             figure;
% %             imagesc(disp_z_smooth_desamp_interp(:,:,10));
% %             [x,y] = ginput(1)
% %             
% %             waveform_left = squeeze(disp_z_smooth_desamp_interp(round(y),round(x),:))
% %             
% %             
% %             
% %             % Fernando filtering
% %             
%             Ts = delta_t_interp;
% %             Line = waveform_left;
% %             FFT = fft(Line,2^12);
% %             freq = linspace(0,1,2^12)/Ts;
% %             
% %             Time = [0:Ts:(length(waveform_left)-1)*Ts]*1e3;
% %             
% %             fig = figure;
% %             plot(freq,abs(FFT))
% %             grid on
% %             ylabel('Magnitude (Arb.)');
% %             xlabel('Frequency (Hz)');
% %             title('FFT of the signal');
% %             axis([0 4000 0 max(abs(FFT))*1.1])
% %             h = imrect(gca,[1e3-200,-200,400,(max(abs(FFT))*1.1)+400]);
% %             position = wait(h);
% position = [25.8064516129029,-200.010743994973,1478.34101382489,401.228396758631]
%             b = fir1(50,[position(1)*2*Ts (position(1)+position(3))*2*Ts]);
% 
%             delay = mean(grpdelay(b));
%             assignin('base','b',b);
%         else
%             b = evalin('base','b');
%             delay = mean(grpdelay(b));
%         end
        
        
        progressbar('ii','jj')
        for ii = 1:size(disp_z_smooth_desamp_interp,1)
            counter = 1;
            for jj = N_radius+1:Nx-N_radius
                waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,jj-N_radius,:));
%                 waveform_left = [filter(b,1,waveform_left(delay+1:end)); zeros(1,delay)'];
                waveform_left = detrend(waveform_left);
                waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,jj+N_radius,:));
%                 waveform_right = [filter(b,1,waveform_right(delay+1:end)); zeros(1,delay)'];
                waveform_right = detrend(waveform_right);
                
                [acor,lag] = xcorr(waveform_left,waveform_right);
                [pks, locs] = findpeaks(acor,'SortStr','descend','NPeaks',3);
                peaks = lag(locs);
                
                
                if folderIndex <= 50 % Center excitation
                    if jj < round(size(particleVelocity,2)/2)
                        newPeaks = abs(peaks(find(peaks > 10)));
                        if ~isempty(newPeaks)
                            newPeaks = newPeaks(1);
                            newLag = newPeaks;
                        else
                            newLag = 0;
                        end
                    else
                        newPeaks = abs(peaks(find(peaks < -10)));
                        if ~isempty(newPeaks)
                            newPeaks = newPeaks(1);
                            newLag = newPeaks;
                        else
                            newLag = 0;
                        end
                    end
                elseif 50 < folderIndex && folderIndex <= 100 % Left excitation
                    if jj < 133
                        newPeaks = abs(peaks(find(peaks > 10)));
                        if ~isempty(newPeaks)
                            newPeaks = newPeaks(1);
                            newLag = newPeaks;
                        else
                            newLag = 0;
                        end
                    else
                        newPeaks = abs(peaks(find(peaks < -10)));
                        if ~isempty(newPeaks)
                            newPeaks = newPeaks(1);
                            newLag = newPeaks;
                        else
                            newLag = 0;
                        end
                    end
                else % Right excitation
                    if jj < 446
                        newPeaks = abs(peaks(find(peaks > 10)));
                        if ~isempty(newPeaks)
                            newPeaks = newPeaks(1);
                            newLag = newPeaks;
                        else
                            newLag = 0;
                        end
                    else
                        newPeaks = abs(peaks(find(peaks < -10)));
                        if ~isempty(newPeaks)
                            newPeaks = newPeaks(1);
                            newLag = newPeaks;
                        else
                            newLag = 0;
                        end
                    end
                end
                %                     [~,I] = min(acor);
                %                 [~,I] = max(acor);
                %                 end
                %                 TOF(ii,jj) = abs(lag(I)*delta_t_interp);
                TOF(ii,jj) = abs(newLag*delta_t_interp);
                
                progressbar([],counter/length(N_radius+1:Nx-N_radius))
                counter = counter+ 1;
            end
            progressbar(ii/size(disp_z_smooth_desamp_interp,1))
        end
        toc();
        TOF_speed = 2*N_radius*Parameters.delta_x./TOF;
        elasticity = density*TOF_speed.^2/1000;
        
        TOF_speed(TOF_speed == Inf) = NaN;
        TOF_speed(TOF_speed == -Inf) = NaN;
        
        
        
        close all force;
        
        
        if debugMode == 1
            notify('Debug!')
            figure; imagesc(TOF_speed); colormap(jet)
            caxis([0 10])
            left = roipoly
            right = roipoly
            mean(mean(TOF_speed(left),'omitnan'))
            mean(mean(TOF_speed(right),'omitnan'))
            close all force;
            pause
        end

        figure; imagesc(xaxis,zaxis_desamp,TOF_speed,[0,10]);axis equal tight
        title('Shear speed map');ylabel('Axial'); xlabel('Lateral');colormap('jet');
        h = colorbar; xlabel(h,'(m/s)','FontSize',14)