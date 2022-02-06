function createDynamicGIFsSmoothed(folders,CAXIS,CMAP)


% if isempty(regexp(folders(1).folder,'Center'))
    numberOfPositions = length(folders)
    for folderIndex = 1:numberOfPositions
        load([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'xaxis','zaxis','TOF_speed')
        VELOCITIES(:,:,folderIndex) = TOF_speed;
    end
    
    for VPDIndex = 1:size(VELOCITIES,3)
imageArray = double(VELOCITIES(:,:,VPDIndex));
[rows columns numberOfColorBands] = size(imageArray);

k = 1; % Order of the polynomial
windowSize = 21;
verticallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 1);
subplot(2, 2, 2);
imshow(verticallySmoothedImage, [0 255]);
title('Savitzky-Golay filtered in the vertical direction only');

% Apply the Savitzky-Golay filter.
% First apply it in the vertical (row) direction.
k = 1; % Order of the polynomial
windowSize = 11;
horizontallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 2);
subplot(2, 2, 3);
imshow(horizontallySmoothedImage, [0 255]);
title('Savitzky-Golay filtered in the horizontal direction only');

doublySmoothedImage = sgolayfilt(verticallySmoothedImage, k, windowSize, [], 2);
subplot(2, 2, 4);
imshow(doublySmoothedImage, [0 255]);
title('Savitzky-Golay filtered in both directions');

VELOCITIES(:,:,VPDIndex) = doublySmoothedImage;
end
    % playWaveVideo(VELOCITIES)
    folderToSave = folders(folderIndex).folder
    rootFolder = regexp(folderToSave,filesep)
    rootFolder = folderToSave(1:rootFolder(end))
    figure;
    for frame = 1:numberOfPositions
        imagesc(xaxis,zaxis,VELOCITIES(:,:,frame))
        h = colorbar
        colormap(CMAP)
        caxis(CAXIS)
        xlabel('Distance (mm)')
        ylabel('Depth (mm)')
        ylabel(h, 'Speed (m/s)')
        if frame == 1
            gif([rootFolder,filesep,'dynamicSmoothed.gif'],'DelayTime', 0.2)
        else
            gif
        end
    end
    save([rootFolder,filesep,'CenterVelocitySmoothed.mat'])
% else
%     numberOfPositions = length(folders)/3
%     clearvars VELOCITIES
%     for folderIndex = 1:numberOfPositions
%         load([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'xaxis','zaxis','TOF_speed')
%         VELOCITIES(:,:,folderIndex) = TOF_speed;
%     end
%     
%     % playWaveVideo(VELOCITIES)
%     folderToSave = folders(folderIndex).folder
%     rootFolder = regexp(folderToSave,filesep)
%     rootFolder = folderToSave(1:rootFolder(end))
%     figure;
%     for frame = 1:numberOfPositions
%         imagesc(xaxis,zaxis,VELOCITIES(:,:,frame))
%         colorbar
%         colormap(jet)
%         caxis([1 5])
%         if frame == 1
%             gif([rootFolder,filesep,'dynamicCenter.gif'],'DelayTime', 0.1)
%         else
%             gif
%         end
%     end
%     save([rootFolder,filesep,'CenterVelocity.mat'])
%     clearvars VELOCITIES
%     counter = 1;
%     for folderIndex = numberOfPositions+1:2*numberOfPositions
%         load([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'xaxis','zaxis','TOF_speed')
%         VELOCITIES(:,:,counter) = TOF_speed;
%         counter = counter+1;
%     end
%     
%     % playWaveVideo(VELOCITIES)
%     folderToSave = folders(folderIndex).folder
%     rootFolder = regexp(folderToSave,filesep)
%     rootFolder = folderToSave(1:rootFolder(end))
%     figure;
%     for frame = 1:numberOfPositions
%         imagesc(xaxis,zaxis,VELOCITIES(:,:,frame))
%         colorbar
%         colormap(jet)
%         caxis([1 5])
%         if frame == 1
%             gif([rootFolder,filesep,'dynamicLeft.gif'],'DelayTime', 0.1)
%         else
%             gif
%         end
%     end
%     save([rootFolder,filesep,'LeftVelocity.mat'])
%     
%     clearvars VELOCITIES
%     counter = 1;
%     for folderIndex = 2*numberOfPositions+1:3*numberOfPositions
%         load([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'xaxis','zaxis','TOF_speed')
%         VELOCITIES(:,:,counter) = TOF_speed;
%         counter = counter+1;
%     end
%     
%     % playWaveVideo(VELOCITIES)
%     folderToSave = folders(folderIndex).folder
%     rootFolder = regexp(folderToSave,filesep)
%     rootFolder = folderToSave(1:rootFolder(end))
%     figure;
%     for frame = 1:numberOfPositions
%         imagesc(xaxis,zaxis,VELOCITIES(:,:,frame))
%         colorbar
%         colormap(jet)
%         caxis([1 5])
%         if frame == 1
%             gif([rootFolder,filesep,'dynamicRight.gif'],'DelayTime', 0.1)
%         else
%             gif
%         end
%     end
%     save([rootFolder,filesep,'RightVelocity.mat'])
%     
% end