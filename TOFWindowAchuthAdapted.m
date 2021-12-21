TOF = NaN.*zeros(size(disp_z_smooth_desamp_interp,1),Nx)
N_radius = 90; % Half distance for choosing reference waveform. [ADJUSTABLE]
tic();
progressbar('ii','jj')
for ii = 1:size(disp_z_smooth_desamp_interp,1)
    counter = 1;
    for jj = 1:size(disp_z_smooth_desamp_interp,2);
        if jj < N_radius
            x_start = 1;
        else
            x_start = jj;
        end
        
        if jj > Nx-N_radius
            x_end = size(disp_z_smooth_desamp_interp,2);
        else
            x_end = jj + N_radius;
        end
        
        % Split the data 
        leftSide = disp_z_smooth_desamp_interp(:,1:round(size(disp_z_smooth_desamp_interp,2)/2),:);
        leftSide = fliplr(leftSide);
        rightSide = disp_z_smooth_desamp_interp(:,round(size(disp_z_smooth_desamp_interp,2)/2)+1:end,:);
        if jj < round(size(particleVelocity,2)/2)
     
        waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,x_start,:));
        waveform_left = [filter(b,1,waveform_left(delay+1:end)); zeros(1,delay)'];
        waveform_left = detrend(waveform_left);
        waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,x_end,:));
        waveform_right = [filter(b,1,waveform_right(delay+1:end)); zeros(1,delay)'];
        waveform_right = detrend(waveform_right);
        
        [acor,lag] = xcorr(waveform_left,waveform_right);
        [pks, locs] = findpeaks(acor,'SortStr','descend','NPeaks',3);
        peaks = lag(locs);
        
        
    
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

   TOF_speed = 2*N_radius*Parameters.delta_x./TOF;
   figure; imagesc(TOF_speed)
   title(['Sliding Window including Ends, N_radius: ',num2str(N_radius)])
   caxis([0 10])
   
   TOF2 = NaN.*zeros(size(disp_z_smooth_desamp_interp,1),Nx)
    jjMin = N_radius+1;
        jjMax = Nx-N_radius
%         progressbar('ii','jj')
        parfor ii = 1:size(disp_z_smooth_desamp_interp,1)
            counter = 1;
            for jj = jjMin:jjMax
                
                waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,jj-N_radius,:));
                waveform_left = [filter(b,1,waveform_left(delay+1:end)); zeros(1,delay)'];
                waveform_left = detrend(waveform_left);
                waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,jj+N_radius,:));
                waveform_right = [filter(b,1,waveform_right(delay+1:end)); zeros(1,delay)'];
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
                TOF2(ii,jj) = abs(newLag*delta_t_interp);
                
%                 progressbar([],counter/length(N_radius+1:Nx-N_radius))
                counter = counter+ 1;
            end
%             progressbar(ii/size(disp_z_smooth_desamp_interp,1))
        end
        toc();
       TOF_speed2 = 2*N_radius*Parameters.delta_x./TOF2;
   figure; imagesc(TOF_speed2)
   title(['Sliding Window Original, N_radius: ',num2str(N_radius)])
   caxis([0 10])
   
   
      TOF3 = NaN.*zeros(size(disp_z_smooth_desamp_interp,1),Nx)
    jjMin = N_radius+1;
        jjMax = Nx-N_radius
%         progressbar('ii','jj')
        parfor ii = 1:size(disp_z_smooth_desamp_interp,1)
            counter = 1;
            for jj = jjMin:jjMax
                
                waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,jj-N_radius,:));
                waveform_left = [filter(b,1,waveform_left(delay+1:end)); zeros(1,delay)'];
                waveform_left = detrend(waveform_left);
                waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,jj+N_radius,:));
                waveform_right = [filter(b,1,waveform_right(delay+1:end)); zeros(1,delay)'];
                waveform_right = detrend(waveform_right);
                
                [acor,lag] = xcorr(waveform_left,waveform_right);
                [pks, locs] = findpeaks(acor,'SortStr','descend','NPeaks',3);
                peaks = lag(locs);
                
                
                    newPeaks = abs(peaks(find(peaks > 10)));
                    if ~isempty(newPeaks)
                        newPeaks = newPeaks(1);
                        newLag = newPeaks;
                    else
                        newLag = 0;
                    end
                   
                %                     [~,I] = min(acor);
                %                 [~,I] = max(acor);
                %                 end
                %                 TOF(ii,jj) = abs(lag(I)*delta_t_interp);
                TOF3(ii,jj) = abs(newLag*delta_t_interp);
                
%                 progressbar([],counter/length(N_radius+1:Nx-N_radius))
                counter = counter+ 1;
            end
%             progressbar(ii/size(disp_z_smooth_desamp_interp,1))
        end
        toc();
       TOF_speed3 = 2*N_radius*Parameters.delta_x./TOF3;
   figure; imagesc(TOF_speed3)
   title(['Sliding Window Single Pass, N_radius: ',num2str(N_radius)])
   caxis([0 10])