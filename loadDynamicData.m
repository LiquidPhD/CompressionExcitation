function [IQData,VMIQ,vec_phase_diff,Parameters] = loadDynamicData(folders,folderIndex,lowerBound)

% Load raw data
load([folders(folderIndex).folder,filesep,'IData.mat']);
load([folders(folderIndex).folder,filesep,'QData.mat']);
temp = load([folders(folderIndex).folder,filesep,'Parameter.mat']);
Parameter = temp.Parameter;
Parameters = Parameter;
try
    IData = squeeze(IData);
    QData = squeeze(QData);
catch
    IData = squeeze(IBuffer);
    QData = squeeze(QBuffer);
end

% Convert to complex data
IQData = complex(IData,QData);
xAxis = Parameter.Trans.ElementPos(:,1);

% For old data backwards-compatibility
try
    Parameters.delta_z = Parameter.PData(3).PDelta(3)*Parameter.Trans.lambda;
    Parameters.delta_x = Parameter.PData(3).PDelta(1)*Parameter.Trans.lambda;
catch
    Parameters.delta_z = Parameter.PData(1).PDelta(3)*Parameter.Trans.lambda;
    Parameters.delta_x = Parameter.PData(1).PDelta(1)*Parameter.Trans.lambda;
end
Parameters.delta_t = 100*1e-6;

timeValue = 1e-4; % in seconds; 100 microseconds
tAxis = 0:timeValue:size(IData,3);

AxialLengthToTest = 5;
try
    AxialResolution = Parameter.PData(3).PDelta(3);
catch
    AxialResolution = Parameter.PData(1).PDelta(3);
end
Parameters.M = ceil(AxialLengthToTest/AxialResolution);
Parameters.N = 3;
Parameters.fc = Parameter.Trans.frequency*1e6;
Parameters.c=1500;


% Vector method for profile reconstruction
IQData = IQData(1:lowerBound,:,:);
sdl = ones([1 size(IQData,1)]);
VMIQ = permute(IQData,[2 1 3]);
[vector_complex_OCE_data,vec_phase_diff] = ...
    vec_meth_snr(VMIQ,sdl,20);
VMIQ = permute(vector_complex_OCE_data,[2 1 3]);
vec_phase_diff = permute(vec_phase_diff,[2 1 3]);
end