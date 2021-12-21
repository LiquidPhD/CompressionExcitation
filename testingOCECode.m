clearvars;
baseFolder = 'F:\AgarPhantomHalfHalfAgarSensor';
%
folders = rdir([baseFolder,filesep,'Dynamic\*\**\*Param*']);
[lowerBound] = selectLowerBound(folders)

folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
[IQData,VMIQ,vec_phase_diff] = loadStaticData(folders,lowerBound,20);
folderIndex = 1;
 [~,~,~,Parameters] = loadDynamicData(folders,folderIndex,lowerBound);

[particleVelocityL,BScan] = calculateParticleVelocity(IQData,Parameters);
particleVelocity = permute(vec_phase_diff,[2 1 3]);
[Nz,Nx,Nt]= size(particleVelocity);         % The dimensions of data in the z axis, the x axis and time.
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);

deltaX = diff(xaxis);
deltaX = deltaX(1);
displacement = cumsum(particleVelocity,3);
b_scan_strain = strain_2D_slash_no_sensor(baseFolder,deltaX,...
   displacement(:,:,10),50,50);

figure;
TEST = imgaussfilt(b_scan_strain,'FilterSize',3); 
% imagesc(imgaussfilt(b_scan_strain),'FilterSize',3); 
imagesc(TEST)
colormap(jet)

for k = 1:size(particleVelocity,3)
    BSS(:,:,k) = strain_2D_slash_no_sensor(baseFolder,deltaX,...
   displacement(:,:,k),50,50);
    PVS(:,:,k) = strain_2D_slash_no_sensor(baseFolder,deltaX,...
   particleVelocity(:,:,k),50,50);
end


