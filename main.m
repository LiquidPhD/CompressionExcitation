% clearvars;
baseFolder = 'I:\SiliconeHalfAndHalfSiliconeSensor2';
%
folders = rdir([baseFolder,filesep,'Dynamic\*\**\*Param*']);
% [lowerBound] = selectLowerBound(folders)
lowerBound = 1018
%% Process dynamic data

for folderIndex = 1:length(folders)
    
    % Load dynamic data
    [IQData,VMIQ,vec_phase_diff,Parameters] = loadDynamicData(folders,folderIndex,lowerBound);
    
    % Get particle velocity
    [particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters);
    
    % % Directional filtering
    % [particleVelocity] = directionalFiltering(particleVelocity);
    
    % Get sizes and set axes
    [Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
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
    N_radius = 50
    [TOF_speed] = runTOF(particleVelocity,Parameters,folderIndex,position,N_radius)
    
    save([folders(folderIndex).folder,filesep,'speedProcessing.mat'],'-v7.3')
    
    figure; imagesc(xaxis,zaxis,TOF_speed,[0,10]);axis equal tight
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
    for k = 1:size(particleVelocity,3)
        imagesc(xaxis,zaxis,particleVelocity(:,:,k))
        xlabel('Distance (mm)');
        ylabel('Distance (mm');
        caxis([-1e-7 1e-7])
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


createDynamicGIFs(folders)
%% Process static data

% Load static data
clearvars -except baseFolder lowerBound Parameters
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
[IQData] = loadStaticData(folders,lowerBound);
folderIndex = 1;
 [~,Parameters] = loadDynamicData(folders,folderIndex,lowerBound);

[particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters);

[Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);

for frame = 1:size(IQData,3)-1
    complex_phase_diff_comp(:,:,frame) = IQData(:,:,frame).*conj(IQData(:,:,frame+1));
end

cutoffFrame = 50;
figure;
for frame = 1:cutoffFrame-1%size(complex_phase_diff_comp,3)
    imagesc(angle(complex_phase_diff_comp(:,:,frame)))
    colorbar
    colormap(jet)
    if frame == 1
        gif([baseFolder,filesep,'complexconj2.gif'],'DelayTime', 0.1)
    else
        gif
    end
    
end

figure;
for frame = 1:cutoffFrame-3 %size(particleVelocity,3)
    imagesc(xaxis,zaxis,particleVelocity(:,:,frame))
    xlabel('Distance (mm)')
    ylabel('Distance (mm)')
    colorbar
    caxis([-3e-5 3e-5])
    colormap(jet)
    if frame == 1
        gif([baseFolder,filesep,'PVel2.gif'],'DelayTime', 0.1)
    else
        gif
    end
end




BScanWidth = size(particleVelocity,2);
for timeIndex = 1:cutoffFrame-3%size(particleVelocity,3)
    for pos=1:BScanWidth
        
        %temp=round(sensor_dl(:,pos));
        temp2=particleVelocity(:,pos,timeIndex);
        %     if temp>0
        temp3=temp2(1:21);
        temp4=median(temp3);
        temp5=temp2-temp4;
        
        displacement_corr(:,pos,timeIndex)=temp5;
    end
end
%
%     figure;
%     for k = 1:size(displacement,3)
%     subplot(2,1,1)
%     imagesc(displacement(:,:,k));
%     subplot(2,1,2)
%     imagesc(displacement_corr(:,:,k));
%     pause(0.2)
% end
% caxisMax = max(max(max(displacement_corr)));
% caxisMin = min(min(min(displacement_corr)));
% [CAXIS] = selectCaxis3D(displacement_corr);

dis_cumsum=cumsum(abs(displacement_corr),3);
% figure;
% for k = 1:size(displacement_corr,3)
%     imagesc(dis_cumsum(:,:,k))
%     caxis([0,1e-3])
%     colorbar
%     colormap(jet)
%     pause(0.2)
% end


figure;
for frame = 1:size(dis_cumsum,3)
    imagesc(xaxis,zaxis,dis_cumsum(:,:,frame))
    xlabel('Distance (mm)')
    ylabel('Distance (mm)')
    colorbar
    %     caxis([0 0.4e-3])
    colormap(jet)
    if frame == 1
        gif([baseFolder,filesep,'displacement2.gif'],'DelayTime', 0.1)
    else
        gif
    end
end


save([folders(1).folder,filesep,'Static2.mat'])
%     ppm = ParforProgressbar(size(displacement,3),'parpool','local')
%
% parfor k = 1:size(displacement,3)
%     b_scan_strain = strain_2D_robust_no_sensor2([baseFolder,filesep,'Static'],delta_z,clippedBScan(:,:,1),...
%     displacement(:,:,k),100);
% STRAIN(:,:,k) = b_scan_strain;
% ppm.increment;
% end
% delete(ppm);

b_scan_strain2 = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    dis_cumsum(:,:,end),75);


figure;
imagesc(xaxis,zaxis,b_scan_strain');
colormap(hot)
xlabel('Distance (mm)')
ylabel('Distance (mm)')
export_fig([baseFolder,filesep,'Strain2.png'])


save([baseFolder,filesep,'STRAIN2.mat'],'b_scan_strain')

