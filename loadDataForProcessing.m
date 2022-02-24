function [Parameters,IQData,vec_phase_diff,Loupas_phase_shift] = loadDataForProcessing(baseFolder,lowerBound)
% Load velocity data
load([baseFolder,filesep,'Dynamic',filesep,'CenterVelocity.mat'])
% Load static data
folders = rdir([baseFolder,filesep,'Static\*\**\*Param*']);
folders = natsortfiles(folders);
if ~isempty(rdir([baseFolder,filesep,'combinedIQData.mat']))
    load([baseFolder,filesep,'combinedIQData.mat'])
    
temp = load([folders(folderIndex).folder,filesep,'Parameter.mat']);
Parameter = temp.Parameter;
Parameters = Parameter;
clearvars temp Parameter;
xAxis = Parameters.Trans.ElementPos(:,1);

% For old data backwards-compatibility
try
    Parameters.delta_z = Parameters.PData(3).PDelta(3)*Parameters.Trans.lambda;
    Parameters.delta_x = Parameters.PData(3).PDelta(1)*Parameters.Trans.lambda;
catch
    Parameters.delta_z = Parameters.PData(1).PDelta(3)*Parameters.Trans.lambda;
    Parameters.delta_x = Parameters.PData(1).PDelta(1)*Parameters.Trans.lambda;
end
Parameters.delta_t = 100*1e-6;

timeValue = 1e-4; % in seconds; 100 microseconds
tAxis = 0:timeValue:size(IQData,3);

AxialLengthToTest = 5;
try
    AxialResolution = Parameters.PData(3).PDelta(3);
catch
    AxialResolution = Parameters.PData(1).PDelta(3);
end
Parameters.M = ceil(AxialLengthToTest/AxialResolution);
Parameters.N = 3;
Parameters.fc = Parameters.Trans.frequency*1e6;
Parameters.c=1500;
    folderIndex = 1;
else
[IQData,Parameters] = loadStaticDataGaussFilt(folders,lowerBound,length(folders));
folderIndex = 1;
save([baseFolder,filesep,'combinedIQData.mat'],'IQData','Parameters');
% load('StrainImages.mat')
% load('combinedIQData.mat')
end

% Take out known bad areas
% IQData(:,:,28) = []; IQData(:,:,17) = []; 
% IQData = IQData(:,:,1:26);
% IQData = IQData(:,:,10:end);
% VELOCITIES = VELOCITIES(:,:,10:end);

% Calc vec phase diff
[vec_phase_diff,VMIQ] = VPD(IQData);

% 2D Loupas
sdl = ones([1 size(IQData,2)]);
[~,Loupas_phase_shift] = ...
    Loupas_estimator_USE(IQData, sdl);
Loupas_phase_shift = permute(Loupas_phase_shift,[2 1 3]);

displacementData = cumsum(Loupas_phase_shift,3);

% Use frames 1:26
% Loupas_phase_shift = Loupas_phase_shift(:,:,1:26);
% Sizes for further processing

[Nz,Nx,Nt]= size(vec_phase_diff);         % The dimensions of data in the z axis, the x axis and time.
% Nz = 1395;
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt);
% zaxis = zaxis(startDepth:end);
images = mat2gray(abs(IQData));
figFlag = 0;
pulseSelectorFlag = 0 ;