function [IQData,Parameters,startDepth] = loadStaticDataGaussFilt(folders,lowerBound,cutoffFolder)

% Preallocate matrix for speed
IQData = complex(zeros([lowerBound 569 cutoffFolder]),0);
IData = zeros([1395 569 120]);
QData = IData;
% Load data, average over time since static, and load into IQData
tic
for folderIndex = 1:cutoffFolder
       fprintf('Folder......%d/%d\n',folderIndex,length(folders));
    X = load([folders(folderIndex).folder,filesep,'IData.mat']);
    Y = load([folders(folderIndex).folder,filesep,'QData.mat']);
%     load([folders(folderIndex).folder,filesep,'Parameter.mat'])
    try
        IData = squeeze(X.IData);
        QData = squeeze(Y.QData);
    catch
        IData = squeeze(X.IBuffer);
        QData = squeeze(Y.QBuffer);
    end
    IData = imgaussfilt(IData);
    QData = imgaussfilt(QData);
    IQIntermed = mean(complex(IData,QData),3);
    IQData(:,:,folderIndex) = IQIntermed(1:lowerBound,:,:);
% 
% for loopSize = 1:size(IData,3)
% k = 1; % Order of the polynomial
% windowSize = 11;
% verticallySmoothedImage = sgolayfilt(IData(:,:,loopSize), k, windowSize, [], 1);
% 
% % Apply the Savitzky-Golay filter.
% % First apply it in the vertical (row) direction.
% k = 1; % Order of the polynomial
% windowSize = 11;
% % horizontallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 2);
% 
% doublySmoothedImage = sgolayfilt(verticallySmoothedImage, k, windowSize, [], 2);
% % subplot(2, 2, 4);
% % imshow(doublySmoothedImage, [0 255]);
% % title('Savitzky-Golay filtered in both directions');
% 
% % VPD_smoothed(:,:,VPDIndex) = doublySmoothedImage;
% IData(:,:,loopSize) = doublySmoothedImage;
% 
% k = 1; % Order of the polynomial
% windowSize = 21;
% verticallySmoothedImage = sgolayfilt(QData(:,:,loopSize), k, windowSize, [], 1);
% 
% % Apply the Savitzky-Golay filter.
% % First apply it in the vertical (row) direction.
% k = 1; % Order of the polynomial
% windowSize = 21;
% % horizontallySmoothedImage = sgolayfilt(imageArray, k, windowSize, [], 2);
% 
% doublySmoothedImage = sgolayfilt(verticallySmoothedImage, k, windowSize, [], 2);
% % subplot(2, 2, 4);
% % imshow(doublySmoothedImage, [0 255]);
% % title('Savitzky-Golay filtered in both directions');
% 
% % VPD_smoothed(:,:,VPDIndex) = doublySmoothedImage;
% QData(:,:,loopSize) = doublySmoothedImage;
% end
% 
%     IQIntermed = mean(complex(IData,QData),3);
%     IQDataFiltered(:,:,folderIndex) = IQIntermed(1:lowerBound,:,:);

%     if folderIndex == 1
%         figure; imagesc(mat2gray(abs(IQData(:,:,1)))); colormap(gray); caxis([0 0.1])
% [~,y] = ginput(1)
% close all force;
%     end
end
y = 1;

% load([folders(folderIndex).folder,filesep,'Parameter.mat']);
% figure; imagesc(IQData(:,:,1));
% [~,y] = ginput(1)
% IQData = IQData(round(y):lowerBound,:,:);
% % IQData = IQData(1:lowerBound,:,:);
% sdl = ones([1 size(IQData,1)]);
% VMIQ = permute(IQData,[2 1 3]);
% [vector_complex_OCE_data,vec_phase_diff] = ...
%     vec_meth_snr(VMIQ,sdl,20);
% VMIQ = permute(vector_complex_OCE_data,[2 1 3]);
% vec_phase_diff = permute(vec_phase_diff,[2 1 3]);
% toc

temp = load([folders(folderIndex).folder,filesep,'Parameter.mat']);
Parameter = temp.Parameter;
Parameters = Parameter;

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
startDepth = y;
end