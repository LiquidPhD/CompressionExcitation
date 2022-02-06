clearvars;
baseFolder = 'E:\ThinAgar121221';
%
folders = rdir([baseFolder,filesep,'Dynamic\*\**\*Param*']);
% [lowerBound] = selectLowerBound(folders)
lowerBound = 1395
%% Process dynamic data
folders = natsortfiles(folders);
for folderIndex = 1:length(folders)
    if ~isempty(rdir([folders(folderIndex).folder,filesep,'speedProcessing.mat']))
        clearvars -except folders lowerBound baseFolder folderIndex
        TEST = matfile([folders(folderIndex).folder,filesep,'speedProcessing.mat']);
         if size(TEST.TOF_speed,1) == 1395
        continue
        end
    load([folders(folderIndex).folder,filesep,'speedProcessing.mat'])
   
    else
    % Load dynamic data
    [IQData,VMIQ,vec_phase_diff,Parameters] = loadDynamicData(folders,folderIndex,lowerBound);
    
    % Get particle velocity
%     [particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters);
    BScan = mat2gray(abs(IQData));
    % % Directional filtering
    % [particleVelocity] = directionalFiltering(particleVelocity);
%     particleVelocity = vec_phase_diff;
    
    sdl = ones([1 size(IQData,2)])
[~,Loupas_phase_shift] = ...
    Loupas_estimator_USE(IQData, sdl);
Loupas_phase_shift = permute(Loupas_phase_shift,[2 1 3]);
particleVelocity = Loupas_phase_shift;
    % Get sizes and set axes
    end
    [Nz,Nx,Nt]= size(Loupas_phase_shift);         % The dimensions of data in the z axis, the x axis and time.
    zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
    xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
    taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);
    
    %     % % Frequency Filtering
    %     [frequencyFilteredPV] = frequencyFiltering(particleVelocity,BSc an,Parameters,xaxis);
    %     playWaveVideo(frequencyFilteredPV)
    %
    %
    %     % Passive algorithm
    %     [wavelength,selectedFreq,frequencyMap] = runPassive(particleVelocity,BScan,Parameters,folderIndex)
    % TOF algorithm
    position = [0 0 0 0]
    N_radius = 30
    
    [TOF_speed] = TOF_noFilter(Loupas_phase_shift,Parameters,folderIndex,position,N_radius)
%     [TOF_speed] = runTOF(double(vec_phase_diff),Parameters,folderIndex,position,N_radius)
    pause(10)
    save([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'-v7.3')
    
    figure; imagesc(xaxis,zaxis,TOF_speed,[0,15]);axis equal tight
    title('Shear speed map');ylabel('Axial'); xlabel('Lateral');colormap('jet');
    h = colorbar; xlabel(h,'(m/s)','FontSize',14)
    export_fig([folders(folderIndex).folder,filesep,'ShearSpeedMap.png'])
    close all force;
    
    figure; imagesc(xaxis,zaxis,mat2gray(BScan(:,:,1))); axis equal tight;
    caxis([0 0.3])
    title('Bscan');
    ylabel('Axial'); xlabel('Lateral');colormap('gray');
    export_fig([folders(folderIndex).folder,filesep,'BScan.png'])
    close all force;
    
    figure;
    for k = 1:40% size(particleVelocity,3)
        imagesc(xaxis,zaxis,particleVelocity(:,:,k))
        xlabel('Distance (mm)');
        ylabel('Distance (mm');
%         caxis([-1e-7 1e-7])
caxis([-0.1 0.1])
        colormap(fireice);
        colorbar;
        title({'Particle Velocity';[num2str((folderIndex-1)*0.1),' mm compressed']})
        if k == 1
            gif([baseFolder,filesep,'particleVelocity.gif'],'DelayTime', 0.1)
        else
            gif
        end
        
    end
    
end


createDynamicGIFs(folders,[0 15],'jet')