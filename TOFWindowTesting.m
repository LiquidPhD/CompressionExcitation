TOF = NaN.*zeros(size(disp_z_smooth_desamp_interp,1),Nx)
N_radius = 25; % Half distance for choosing reference waveform. [ADJUSTABLE]
tic();
for ii = 1:size(disp_z_smooth_desamp_interp,1)
    counter = 1;
    for jj = 1:size(disp_z_smooth_desamp_interp,2);
        if jj < N_radius
            
            if jj == 1
                TOF_speed(ii,jj) = NaN;
            else
                waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,1,:));
                waveform_left = [filter(b,1,waveform_left(delay+1:end)); zeros(1,delay)'];
                waveform_left = detrend(waveform_left);
                waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,jj+1,:));
                waveform_right = [filter(b,1,waveform_right(delay+1:end)); zeros(1,delay)'];
                waveform_right = detrend(waveform_right);
                
                [acor,lag] = xcorr(waveform_left,waveform_right);
                [pks, locs] = findpeaks(acor,'SortStr','descend','NPeaks',3);
                peaks = lag(locs);
            end
            
            disp(num2str(jj))
        elseif jj > Nx-N_radius
            
            if jj == size(disp_z_smooth_desamp_interp,2)
                TOF_speed(ii,jj) = NaN;
            end
        else
            pause
            waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,jj+1,:));
            waveform_left = [filter(b,1,waveform_left(delay+1:end)); zeros(1,delay)'];
            waveform_left = detrend(waveform_left);
            waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,end,:));
            waveform_right = [filter(b,1,waveform_right(delay+1:end)); zeros(1,delay)'];
            waveform_right = detrend(waveform_right);
            
            [acor,lag] = xcorr(waveform_left,waveform_right);
            [pks, locs] = findpeaks(acor,'SortStr','descend','NPeaks',3);
            peaks = lag(locs);
            
        end
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