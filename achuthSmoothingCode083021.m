
x_kern = 50;
% dis_cumsum = dis_cumsum(:,:,end);
% dis_cumsum = dis_cumsum';
BScanWidth = size(dis_cumsum,1);
sensor_dl = ones([1 size(dis_cumsum,2)])
bdl = ones([1 size(dis_cumsum,2)]).*size(dis_cumsum,2);
how_many_depths = size(dis_cumsum,2);

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
    
    x_disp = median(dis_cumsum(x_start:x_end,sensor_dl(x):bdl(1,x)),1);
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
    imagesc(dis_cumsum'); colormap(jet);
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
