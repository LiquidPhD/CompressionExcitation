%% Load workspace or process ultrasound elastography data
function processUSEForWorkspaceWavelengthUseMultipleAngles(folderIndex,folders,nAng)
%  if isempty((rdir([folders(folderIndex).folder,filesep,'wavelengthWorkspace.mat'])))
folder = folders(folderIndex).folder;
load([folders(folderIndex).folder,filesep,'IData.mat'])
load([folders(folderIndex).folder,filesep,'QData.mat'])
load([folders(folderIndex).folder,filesep,'Parameter.mat'])
try
    IData = squeeze(IData);
    QData = squeeze(QData);
catch
    IData = squeeze(IBuffer);
    QData = squeeze(QBuffer);
end
IQData = complex(IData,QData);
numIQcmpnd = floor(size(IQData,3)/nAng);  % The time-domain image number = Total number of frames.

IQDataCmpnd = zeros(size(IQData,1),size(IQData,2),numIQcmpnd);
for ii = 1:numIQcmpnd
    IQDataCmpnd(:,:,ii) = mean(IQData(:,:,(ii-1)*nAng+1:ii*nAng),3);
end

xAxis = Parameter.Trans.ElementPos(:,1);
if length(Parameter.PData) == 3
    delta_z = Parameter.PData(3).PDelta(3)*Parameter.Trans.lambda;
    delta_x = Parameter.PData(3).PDelta(1)*Parameter.Trans.lambda;
    AxialResolution = Parameter.PData(3).PDelta(3);
else
    delta_z = Parameter.PData.PDelta(3)*Parameter.Trans.lambda;
    delta_x = Parameter.PData.PDelta(1)*Parameter.Trans.lambda;
    AxialResolution = Parameter.PData.PDelta(3);
end
delta_t = 100*1e-6;

timeValue = 1e-4; % in seconds; 100 microseconds
tAxis = 0:timeValue:size(IData,3);
% yAxis = 0:269.6016*Parameter.Trans.lambda*1000
% figure;
% for k = 1:size(IQData,3)
%     imagesc(abs(IQData(:,:,k)))
%     pause(0.2)
% end

% figure; imagesc(abs(IQData(:,:,1)));
% [~,y] = ginput(1)
% close all force;
% topClip = 1;
% bottomClip = round(y(1));
topClip = 1;
bottomClip = size(IQData,1);
AxialLengthToTest = 5;
M = ceil(AxialLengthToTest/AxialResolution);
N = 3;
fc = Parameter.Trans.frequency*1e6
c=1500;
IQDataClip = IQDataCmpnd(topClip:bottomClip,:,:);

ROI_dim_z = size(IQDataClip,1)-M;
ROI_dim_x = size(IQDataClip,2);
ROI_dim_t = size(IQDataClip,3)-N;

IQ = permute(IQDataClip, [1 3 2] );    % For more convenient access!!!
% Pay special attention to its affect on the sumation calculation.



%% PARFOR LOOP
tic
realIQ = real(IQ);
imagIQ = imag(IQ);
parfor_progress(ROI_dim_z);
BScan = abs(IQData(:,:,1));

%         SGO = 3;
%         SGFL = 11;
%         for sgolayfiltindex = 1:size(realIQ,3)
%             realIQ(:,:,sgolayfiltindex) = sgolayfilt(double(realIQ(:,:,sgolayfiltindex)),SGO,SGFL);
%             imagIQ(:,:,sgolayfiltindex) = sgolayfilt(double(imagIQ(:,:,sgolayfiltindex)),SGO,SGFL);
%         end
clearvars displacement clippedBScan
parfor_progress(ROI_dim_z);
parfor ii_z = 1:ROI_dim_z
    %     fprintf('Axial......%d/%d\n',ii_z,ROI_dim_z);
    for ii_t = 1:ROI_dim_t
        I = realIQ( ii_z:ii_z+M-1, ii_t:ii_t+N-1, 1:ROI_dim_x );
        Q = imagIQ( ii_z:ii_z+M-1, ii_t:ii_t+N-1, 1:ROI_dim_x );
        uu = Q(1:M,1:N-1,:).*I(1:M,2:N,:) - I(1:M,1:N-1,:).*Q(1:M,2:N,:);
        ud = I(1:M,1:N-1,:).*I(1:M,2:N,:) + Q(1:M,1:N-1,:).*Q(1:M,2:N,:);
        du = Q(1:M-1,1:N,:).*I(2:M,1:N,:) - I(1:M-1,1:N,:).*Q(2:M,1:N,:);
        dd = I(1:M-1,1:N,:).*I(2:M,1:N,:) + Q(1:M-1,1:N,:).*Q(2:M,1:N,:);
        uu = sum(sum(uu));
        ud = sum(sum(ud));
        du = sum(sum(du));
        dd = sum(sum(dd));
        displacement(ii_z,:,ii_t) = c/(4*pi*fc)*atan(uu./ud)./(1+atan(du./dd)/(2*pi));
        clippedBScan(ii_z,:,ii_t) = BScan( ii_z,:,1 );
    end
    parfor_progress;
end
parfor_progress(0);
toc

[Nz,Nx,Nt]= size(displacement);         % The dimensions of data in the z axis, the x axis and time.
zaxis = linspace(0,(Nz-1)*delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*delta_x,(Nx-1)/2*delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*delta_t,Nt);

save([folders(folderIndex).folder,filesep,'wavelengthWorkspace.mat'],'displacement','IQData','delta_z','xaxis','zaxis','taxis','clippedBScan');
%  else
%         % Do nothing, since we will load the file in the next step.
%  end
end