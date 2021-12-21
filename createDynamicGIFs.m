function createDynamicGIFs(folders,CAXIS,CMAP)


% if isempty(regexp(folders(1).folder,'Center'))
    numberOfPositions = length(folders)
    folders = natsortfiles(folders);
    for folderIndex = 1:numberOfPositions
        load([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'xaxis','zaxis','TOF_speed')
        TOF_speed = imresize(TOF_speed,[length(zaxis) length(xaxis)]);
        VELOCITIES(:,:,folderIndex) = TOF_speed;
    end
    
    figure; imagesc(VELOCITIES(:,:,1))
    caxis(CAXIS);
    colormap(CMAP)
    [~,y] = ginput(1);
    VELOCITIES = VELOCITIES(1:round(y),:,:);
    % playWaveVideo(VELOCITIES)
    folderToSave = folders(folderIndex).folder
    rootFolder = regexp(folderToSave,filesep)
    rootFolder = folderToSave(1:rootFolder(end))
    figure;
    for frame = 1:numberOfPositions
        imagesc(xaxis,zaxis(1:round(y)),VELOCITIES(:,:,frame))
        h = colorbar
        colormap(CMAP)
        caxis(CAXIS)
        xlabel('Distance (mm)')
        ylabel('Depth (mm)')
        ylabel(h, 'Speed (m/s)')
        if frame == 1
            gif([rootFolder,filesep,'dynamic.gif'],'DelayTime', 0.2)
        else
            gif
        end
    end
    save([rootFolder,filesep,'CenterVelocity.mat'])
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