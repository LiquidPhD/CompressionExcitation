function b_scan_strain = strain_2D_robust_no_sensor2(pixel_res,...
    displacement,z_kern)%,bdl,sensor_dl) %,vec_kern,z_kern,bdl,sensor_dl)

delete(gcp('nocreate'));
poolcfg = gcp('nocreate');
if isempty(poolcfg)
    parpool;
end
displacement = displacement';
% intensities = intensities';

[bscanwidth,how_many_depths] = size(displacement);
bdl = ones([1 bscanwidth])*how_many_depths;
sensor_dl = ones([1 bscanwidth]);
z_kern = floor(z_kern/2)*2;
new_disp = zeros(size(displacement));
for x = 1:bscanwidth
%     new_disp(x,:) = displacement(x,:) - median(displacement(x,1:17));
%     new_disp(x,bdl(x)+1:end) = NaN;
end
new_disp = displacement;
x_kern = 50;
newdisp2 = zeros(size(new_disp));
for x = 1:bscanwidth
    if x-(x_kern/2) < 1
        x_start = 1;
    else
        x_start = x-(x_kern/2);
    end
    if x+(x_kern/2) > bscanwidth
        x_end = bscanwidth;
    else
        x_end = x +(x_kern/2);
    end
    newdisp2(x,:) = mean(new_disp(x_start:x_end,:),1);
end
%% robust linear fitting of the strain
close all force;
b_scan_strain = NaN(size(new_disp));

wait_bar = waitbar(0,'Robust Linear Fit');
for x = 1:bscanwidth
    
    parfor z1 = sensor_dl(x):bdl(x)
        
        if (z1 - z_kern/2) < 2
            z_start = 2;
        else
            z_start = z1 - z_kern/2;
        end
        
        if (z1 + z_kern/2) > bdl(x)
            z_end = bdl(x);
        else
            z_end = z1 + z_kern/2;
        end
        
        
        %     z_fit_size = z_end-z_start+1;
        pix_dists = (z_start:z_end)*pixel_res;
        sub_displacement = newdisp2(x,z_start:z_end);
%         sub_intensities = intensities(x,z_start:z_end);
        
        % 32s for loop, 10s parfor loop
        warning('off','all');
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'Bisquare';
%         sub_x_ints =  sub_intensities+abs(min(sub_intensities))+1;
%         [x_out,y_out, opts.Weights] = prepareCurveData(pix_dists(:),sub_displacement(:),sub_x_ints(:));
                [x_out,y_out] = prepareCurveData(pix_dists(:),sub_displacement(:));

        [fitresult, gof] = fit(x_out(:), y_out(:), ft, opts );
        warning('on','all');
        x_strain(z1) = fitresult.p1;
        
    end
    b_scan_strain(x,:) = x_strain;
%     
%     figure(1);
%     set(gcf,'units','normalized','outerposition',[0,0.5,0.5,0.5]);
%     imagesc(mat2gray(intensities')); colormap(gray); caxis([0,0.1]);
%     hold on;
%     plot(ones(1,bdl(x)).*x,(1:bdl(x)),'-y','LineWidth',2);
%     hold off;
%     title(pathname,'interpreter','none');
    
    figure(2);
    set(gcf,'units','normalized','outerposition',[0.5,0.5,0.5,0.5]);
    imagesc(new_disp'); colormap(jet); title('Displacement');
    colorbar; %caxis([0,1]);
    
    figure(3);
    set(gcf,'units','normalized','outerposition',[0,0,0.5,0.5]);
    imagesc(b_scan_strain'); colormap(hot); title('Strain');
    
    wait_bar = waitbar(x/bscanwidth,wait_bar,['Processing ',...
        num2str(x/bscanwidth*100),'%']);
end
close all force;
b_scan_strain = b_scan_strain*1000;

x_scan = 10;
x_ticks = 0:bscanwidth/4:bscanwidth;
x_tick_labels = x_ticks.*(x_scan/bscanwidth);
y_ticks = 0:how_many_depths/4:how_many_depths;
y_tick_labels = round((y_ticks*pixel_res/1000)*10)/10;

figure(4);
set(gcf,'units','normalized','outerposition',[0,0,1,0.5]);
imagesc(b_scan_strain');
cbar1 = colorbar; colormap(hot); caxis([0,1]);
ylabel(cbar1,'Strain (m\epsilon)','Rotation',270,'Position',[5,0.5]);
set(gca,'XTick',x_ticks);
set(gca,'XTickLabels',x_tick_labels);
xlabel('X (mm)');
set(gca,'YTick',y_ticks);
set(gca,'YTickLabels',y_tick_labels);
ylabel('Depth (mm)');
set(gca,'FontSize',20);


% title({'Strain';pathname},'Interpreter','none','FontSize',12);
% print(gcf,fullfile(pathname,'strain',...
%     ['strain_robust_',num2str(vec_kern),'_',num2str(z_kern),'.jpg']),'-djpeg');
% pause(0.01); close all force;
% 
% save(fullfile(pathname,'strain',...
%     ['strain_robust_',num2str(vec_kern),'_',num2str(z_kern),'.mat']),...
%     'b_scan_strain','-v7.3');