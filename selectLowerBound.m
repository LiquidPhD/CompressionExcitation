function [lowerBound] = selectLowerBound(folders)

% Can choose folder index; select 1st since that's how I use it
folderIndex = 1;
% Load data
folder = folders(folderIndex).folder;
load([folders(folderIndex).folder,filesep,'IData.mat'])
load([folders(folderIndex).folder,filesep,'QData.mat'])
load([folders(folderIndex).folder,filesep,'Parameter.mat'])

% Try/Catch/End for backward compatibility
try
    IData = squeeze(IData);
    QData = squeeze(QData);
catch
    IData = squeeze(IBuffer);
    QData = squeeze(QBuffer);
end
IQData = complex(IData,QData);

% Select bottom depth
% notify('Select depth...')
figure;
imagesc(mat2gray(abs(IQData(:,:,1))));
caxis([0 0.15]);
title('Select lower bound')
colormap(gray);
[~,lowerBound] = ginput(1);
lowerBound = round(lowerBound);
close all force;

end