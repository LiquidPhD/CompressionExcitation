for k = 1:size(NF,3)
   NF2(:,:,k) = phase_unwrap(angle(NF(:,:,k))); 
end

x_kern = 25;

% Smooth NF2 and run this through slow strain code
for sectionIndex = 1:size(NF2,3)
%     dis_cumsum = vec_phase_diff(:,:,sectionIndex);
    NF3 = NF2(:,:,sectionIndex)';
    BScanWidth = size(NF3,1);
    sensor_dl = ones([1 size(NF3,2)]);
bdl = ones([1 size(NF3,2)]).*size(NF3,2);
how_many_depths = size(NF3,2);
for x = 1:BScanWidth%size(displacement,1)
    if x < 1
        x_start = 1;
    else
        x_start = x;
    end
    
    if x > BScanWidth - x_kern
        x_end = BScanWidth;
    else
        x_end = x+x_kern;
    end
    
    x_disp = median(NF3(x_start:x_end,sensor_dl(x):bdl(1,x)),1);
    [xData, yData] = prepareCurveData( [], x_disp );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitted_displacement(sensor_dl(1,x):bdl(1,x),x) = feval(fitresult,1:length(x_disp));
    figure(1);
    subplot(3,1,1);
    imagesc(NF3'); colormap(jet);
    hold on;
    plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
    hold off;
    colorbar;
%     title(pathname,'Interpreter','none');
    
    subplot(3,1,2);
    %imagesc(displacement'); colormap(jet);
    plot(xData, yData, xData, fitted_displacement(sensor_dl(1,x):bdl(1,x),x));
    title(x);
    axis tight;
    drawnow;   
    subplot(3,1,3);

    imagesc(fitted_displacement); colormap(jet);
    hold on;
    plot(ones(1,how_many_depths).*x,1:how_many_depths,'-k','LineWidth',2)
    hold off;
    colorbar;
%     title(pathname,'Interpreter','none');
    
end
FDAll(:,:,sectionIndex) = fitted_displacement;
end

for k = 1 :size(FDAll,3)
  BScanStrainAll(:,:,k) = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    FDAll(:,:,k),10);  
  BScanStrainAllAbs(:,:,k) = strain_2D_robust_no_sensor2(Parameters.delta_z,...
      abs(FDAll(:,:,k)),10);
end

save([baseFolder,filesep,'VectorPlusSmoothingPlusStrainxkern25.mat'],'-v7.3')