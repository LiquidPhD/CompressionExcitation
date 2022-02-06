function [sdl,intensities,pathname] = ...
    get_sdl_USE(i_thresh_low,i_thresh_high,surf_thresh,how_many_depths,...
    surf_window)

    % select directory with the phase files
    pathname = uigetdir(start_path,'Select directory with the phase files');
    filename = dir(fullfile(pathname,'*.txt'));
    loops = round(length(filename));
    
    disp(['Processing ',pathname]);
    disp([num2str(loops),' files']);
    
    total_num_xs = str2num(filename(end).name(6:10));
    
    temp = load(fullfile(pathname,filename(1).name));
    if length(temp) < how_many_depths
        how_many_depths = length(temp);
    end
    clear temp;
    
    % initialize variables
    sdl = zeros(1,total_num_xs);
    intensities = zeros(total_num_xs, how_many_depths);
    
    % create a progress bar
    close all force;
%     wait_bar = waitbar(0,'Processing');
    
    % data processing loop
    for counter = 1:loops

        % update the wait bar
%         wait_bar = waitbar(counter/loops,wait_bar,...
%             ['Processing File ', num2str(counter), ' of ', num2str(loops)]);

        % cet the file path
        file_path = filename(counter).name;

        if isempty(regexp(file_path,'IntensityAscan','once'))
            continue % move to next file if it is not a phase profile file
        end

        % get the frame number from the file name
        frame_number = str2double(file_path(6:10));
        
        % load the A scan intensity profile
        aline = load(fullfile(pathname,file_path))';
        aline = aline(1:how_many_depths);

        % generate the binary mask to be used on the data to deteremine which
        % depth levels to use/have a strong enough signal
        smooth_aline = medfilt1(aline,2);
        
        figure(1);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        subplot(2,1,1);
        plot(1:how_many_depths,aline,...
            1:how_many_depths,ones(1,how_many_depths)*(mean(aline)+surf_thresh));
        title(frame_number);
        set(gca,'FontSize',24);
        legend('Aline','Intensity Threshold','Location','NorthEast');
        ylim([i_thresh_low,i_thresh_high]);

        i_pass_1 = find(smooth_aline >= round(mean(aline)+surf_thresh),1);
        
        if ~isempty(i_pass_1)
            hold on;
            plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
            hold off;
        end

        % locate surface for the first frame if threshold too high
        if (counter == 1) && (isempty(i_pass_1))
            i_pass_1 = find(smooth_aline >= round(mean(aline)+(surf_thresh-3)),1);
            
            figure(1);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            subplot(2,1,1);
            plot(1:how_many_depths,aline,...
                1:how_many_depths,ones(1,how_many_depths)*(mean(aline)+(surf_thresh-3)));
            title(frame_number);
            set(gca,'FontSize',24);
            legend('Aline','Intensity Threshold','Location','NorthEast');
            ylim([i_thresh_low,i_thresh_high]);
            
            if ~isempty(i_pass_1)
                hold on;
                plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                hold off;
            end
            
            if isempty(i_pass_1)
                i_pass_1 = find(smooth_aline >= round(mean(aline)+(surf_thresh-6)),1);
                
                figure(1);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                subplot(2,1,1);
                plot(1:how_many_depths,aline,...
                    1:how_many_depths,ones(1,how_many_depths)*(mean(aline)+(surf_thresh-6)));
                title(frame_number);
                set(gca,'FontSize',24);
                legend('Aline','Intensity Threshold','Location','NorthEast');
                ylim([i_thresh_low,i_thresh_high]);
                
                if ~isempty(i_pass_1)
                    hold on;
                    plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                    hold off;
                end
            
                    if isempty(i_pass_1)
                        i_pass_1 = find(smooth_aline >= round(mean(aline)+(surf_thresh-9)),1);
                        
                        figure(1);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        subplot(2,1,1);
                        plot(1:how_many_depths,aline,...
                            1:how_many_depths,ones(1,how_many_depths)*(mean(aline)+(surf_thresh-6)));
                        title(frame_number);
                        set(gca,'FontSize',24);
                        legend('Aline','Intensity Threshold','Location','NorthEast');
                        ylim([i_thresh_low,i_thresh_high]);
                        
                        if ~isempty(i_pass_1)
                            hold on;
                            plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                            hold off;
                        end
                        
                        if isempty(i_pass_1)
                            [~, i_pass_1] = max(smooth_aline);
                        
                            figure(1);
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            subplot(2,1,1);
                            plot(1:how_many_depths,aline,...
                                1:how_many_depths,ones(1,how_many_depths)*(mean(aline)+(surf_thresh-6)));
                            title(frame_number);
                            set(gca,'FontSize',24);
                            legend('Aline','Intensity Threshold','Location','NorthEast');
                            ylim([i_thresh_low,i_thresh_high]);
                        
                            if ~isempty(i_pass_1)
                                hold on;
                                plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                                hold off;
                            end
                        end
                    end
            end
        end

        % locate surface for all other frames if threshold is too high
        if (counter > 1) && (isempty(i_pass_1))
            i_pass_1 = find(smooth_aline >= round(mean(aline)+surf_thresh),1);
            if isempty(i_pass_1) ||...
                    (abs(i_pass_1 - sdl(frame_number - 1)) > surf_window/2)
                temp_b = aline(sdl(frame_number-1)-(surf_window):...
                    sdl(frame_number-1)+surf_window);
                i_pass_1 = find(medfilt1(temp_b) >= round(mean(aline) + (surf_thresh-3)),1) + ...
                    sdl(frame_number-1)-(surf_window+1);
                
                figure(1);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                subplot(2,1,1);
                plot(1:how_many_depths,aline,...
                    1:how_many_depths,ones(1,how_many_depths)*(mean(aline)+(surf_thresh-3)));
                title(frame_number);
                set(gca,'FontSize',24);
                legend('Aline','Intensity Threshold','Location','NorthEast');
                ylim([i_thresh_low,i_thresh_high]);
                
                if ~isempty(i_pass_1)
                    hold on;
                    plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                    hold off;
                end
                
                if isempty(i_pass_1) ||...
                    (abs(i_pass_1 - sdl(frame_number - 1)) > surf_window)
                    temp_b = aline(sdl(frame_number-1)-surf_window:...
                        sdl(frame_number-1)+surf_window);
                    [~, i_pass_1] = max(medfilt1(temp_b));
                    i_pass_1 = i_pass_1 + sdl(frame_number-1)-(surf_window+1);
                    
                    figure(1);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    subplot(2,1,1);
                    plot(1:how_many_depths,aline);
                    title(frame_number);
                    set(gca,'FontSize',24);
                    ylim([i_thresh_low,i_thresh_high]);
                
                    if ~isempty(i_pass_1)
                        hold on;
                        plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                        hold off;
                    end
                end
            end
        end

        % if cannot find the surface of the first location
        if (counter > 1) && (sdl(frame_number-1)) ~= 0
            if (abs(i_pass_1 - sdl(frame_number-1)) > surf_window/2)
                temp_b = aline(sdl(frame_number-1)-surf_window:...
                    sdl(frame_number-1)+(surf_window/2));
                i_pass_1 = find(medfilt1(temp_b) >= round(mean(aline) + (surf_thresh-3)),1) + ...
                    sdl(frame_number-1)-(surf_window+1);
                if isempty(i_pass_1)
                    i_pass_1 = sdl(frame_number-1);
                end
                
                figure(1);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                subplot(2,1,1);
                plot(1:how_many_depths,aline);
                title(frame_number);
                set(gca,'FontSize',24);
                ylim([i_thresh_low,i_thresh_high]);
                
                if ~isempty(i_pass_1)
                    hold on;
                    plot(i_pass_1,aline(i_pass_1),'xr','MarkerSize',26);
                    hold off;
                end
                
            end
            if isempty(i_pass_1)
                i_pass_1 = sdl(frame_number-1);
            end
        end

        if (find(smooth_aline >= i_thresh_high))
            [~, i_pass_2] = max(aline); % if saturated
            i_pass_1 = i_pass_2(1);
        end

        i_pass = i_pass_1+1;% shift surface because of blur
        sdl(frame_number) = i_pass;

        if length(aline) < how_many_depths
            how_many_depths = length(aline);
        end
        
        hold off;
        intensities(frame_number,1:how_many_depths) = aline(1:how_many_depths);
        
        subplot(2,1,2);
        imagesc(intensities',[i_thresh_low, i_thresh_high]);
        xlim([0 total_num_xs]);
        title({'Structural Image and Surface Depth Levels';pathname},...
            'interpreter','none','FontSize',28);
        colormap(gray);
        hold on;
        plot(sdl,'.r');
        hold off;
        pause(0.01);
        
    end
    
    close all force;

    intensities = intensities';
    
end