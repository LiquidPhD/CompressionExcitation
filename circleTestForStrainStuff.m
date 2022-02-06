
corrected_whole = unwrapped_whole-surf_dis2';
figure; imagesc(corrected_whole)

paper = zeros(size(corrected_whole,1),size(corrected_whole,2)); %create empty array
[y x] = size(paper); %define y,x as size of array
r = 80; %define radius of a circle
for i=1:y
    for j=1:x
        if ((i-y/2)^2)+((j-x/2)^2)<(r^2);  %define origin is at the center
            paper(i,j) = 10;  %define array inside the circle eq. = 1
        end
    end
end
imshow(paper);  %show image
paperMask = paper;
paperMask(paperMask>0) = 1;
% corrected_whole(paper>0) = 0;
corrected_whole = corrected_whole+paper;
figure; imagesc(corrected_whole)

corrected_sensor = corrected_whole.*sensorMaskStack(:,:,imageIterator);
corrected_inclusion = corrected_whole.*paperMask;
corrected_bottom = corrected_whole-corrected_sensor-corrected_inclusion;
corrected_bottom(corrected_bottom<0) = 0; 
figure; imagesc(corrected_sensor); figure; imagesc(corrected_inclusion); figure; imagesc(corrected_bottom)
% 
corrected_sensor(corrected_sensor==0) = NaN;
corrected_inclusion(corrected_inclusion==0) = NaN;
corrected_bottom(corrected_bottom==0) = NaN;
figure; imagesc(corrected_sensor); figure; imagesc(corrected_inclusion); figure; imagesc(corrected_bottom)

strain_whole = strain_calculation2(corrected_whole,zaxis(2),1,10);
strain_sensor=strain_calculation2(corrected_sensor,zaxis(2),1,10);
strain_inclusion=strain_calculation2(corrected_inclusion,zaxis(2),1,10);
strain_bottom=strain_calculation2(corrected_bottom,zaxis(2),1,10);
strain_sensor(isnan(strain_sensor)) = 0;
strain_inclusion(isnan(strain_inclusion)) = 0;
strain_bottom(isnan(strain_bottom)) = 0;
strain_added = strain_sensor+strain_inclusion+strain_bottom;



close all force; 
figure;
imagesc(xaxis,zaxis,abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,5]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain done in one shot')

figure;
imagesc(xaxis,zaxis,abs(strain_added)/1000); colormap(hot); colorbar; caxis([0,5]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain parts added')


  STRAINTESTWHOLE = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    strain_whole,100);  
STRAINTEST = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    TEST3,100);  


%% Try without circle...
bottomMask = sensorMaskStack(:,:,imageIterator)
figure; imagesc(bottomMask)
bottomMask(bottomMask==1) = -1; bottomMask(bottomMask==0) = 1; bottomMask(bottomMask==-1) = 0; figure; imagesc(bottomMask)
corrected_sensor = corrected_whole.*sensorMaskStack(:,:,imageIterator);
corrected_bottom = corrected_whole.*bottomMask;
figure; imagesc(corrected_sensor); figure; imagesc(corrected_bottom)

strain_whole = strain_calculation2(corrected_whole,zaxis(2),1,10);
strain_sensor=strain_calculation2(corrected_sensor,zaxis(2),1,10);
strain_bottom=strain_calculation2(corrected_bottom,zaxis(2),1,10);
strain_added = strain_sensor+strain_bottom;

close all force; 
figure;
imagesc(xaxis,zaxis,abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,5]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain done in one shot')

figure;
imagesc(xaxis,zaxis,abs(strain_added)/1000); colormap(hot); colorbar; caxis([0,5]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain parts added')

700/5
569/5

% Split into rectangular windows...
% 568 px: 8 windows.
% 700 px: 7 windows. 
% Windows are 100 x 71 px
windowY = 100;
windowX = 71;
corrected_whole = corrected_whole(:,1:568);
xCounter = 1
yCounter = 1;
for windowY = 1:100:700
    for windowX = 1:71:568
    temp = strain_calculation2(corrected_whole(windowY:windowY+100-1,windowX:windowX+71-1),zaxis(2),1,10);
    TESTWindows(windowY:windowY+100-1,windowX:windowX+71-1) = temp;
    end
end

figure;
imagesc(xaxis,zaxis,abs(TESTWindows)/1000); colormap(hot); colorbar; caxis([0,5]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain parts added')

figure;
imagesc(xaxis,zaxis,abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,5]);
xlabel('Distance (mm)')
ylabel('Depth (mm)')
title('Strain done in one shot')