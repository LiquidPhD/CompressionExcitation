% Processing for compression/excitation stuff. Rewrite of a few files.
% Justin Rippy
% 02/19/22


clearvars;
baseFolder = 'D:\021722GelatinBlueYellow';
lowerBound = 1395
figFlag = 0;
pulseSelectorFlag = 0 ;
printExtraImages = 0;
% Load data for processing
[Parameters,IQData,vec_phase_diff,Loupas_phase_shift] = ...
    loadDataForProcessing(baseFolder,lowerBound)
displacementData = cumsum(Loupas_phase_shift,3);

% Calculate sizes
[Nz,Nx,Nt]= size(vec_phase_diff);         % The dimensions of data in the z axis, the x axis and time.
% Nz = 1395;
zaxis = linspace(0,(Nz-1)*Parameters.delta_z,Nz)*1e3;                      %(mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*Parameters.delta_x,(Nx-1)/2*Parameters.delta_x,Nx)*1e3;    %(mm) Lateral axis.
taxis = linspace(0,(Nt-1)*Parameters.delta_t,Nt); % time axis

% Get images for labels
images = mat2gray(abs(IQData));


% Either load label data or make it
try
    load([baseFolder,filesep,'labels.mat'])
catch
    volumeSegmenter
    mydlg = warndlg('Close when done...', 'Close when done');
    waitfor(mydlg);
    save([baseFolder,filesep,'labels.mat'],'labels')
end

for imageIterator = 2:size(displacementData,3)
    % Get frame to process
    imageToProcess = displacementData(:,:,imageIterator);
    
    % Print B-scan if desired
    if printExtraImages == 1
        close all force;
        figure;
        ax1 = axes;
        
        TEST = mat2gray(abs(IQData(1:size(VELOCITIES,1),:,imageIterator)));
        % TEST(TEST>0.1) = 0.1; TEST(TEST<0) = 0; figure; imagesc(TEST)
        TEST = 255*TEST;
        imagesc(ax1,xaxis,zaxis(1:size(VELOCITIES,1)),TEST(1:size(VELOCITIES,1),:));
        ylabel('Depth (mm)')
        xlabel('Distance (mm)')
        colormap(ax1,'gray');
        caxis([0 25])
        % h = colorbar;
        % ylabel(h,'YM (kPa)')
        % title(ax1,"Young's modulus using shear wave")
        filename = [baseFolder,filesep,'BScan',num2str(imageIterator)];
        figSave(filename,'.png',figFlag)
        
        
    end
    
    % Print wrapped displacement for diagnostic purposes.
    if printExtraImages == 1
        figure;
        imagesc(xaxis,zaxis(1:size(VELOCITIES,1)),imageToProcess(1:size(VELOCITIES,1),:))
        ylabel('Depth (mm)')
        xlabel('Distance (mm)')
        title("Wrapped Displacement")
        colorbar;
        filename = [baseFolder,filesep,'Disp',num2str(imageIterator)];
        figSave(filename,'.png',figFlag)
    end
    
    
% Separate out labels
for k = 1:size(labels,3)
    label = labels(:,:,k);
    sensor = zeros([size(label,1) size(label,2)]);
    sensor(label== 'Label1') = 1;
    left = zeros([size(label,1) size(label,2)]);
    left(label== 'Label2') = 2;
    right = zeros([size(label,1) size(label,2)]);
    right(label== 'Label3') = 3;
    sensorMaskStack(:,:,k) = logical(sensor);
    leftMaskStack(:,:,k) = logical(left);
    rightMaskStack(:,:,k) = logical(right);
end


vec_sensor = imageToProcess.*sensorMaskStack(:,:,imageIterator);
vec_left = imageToProcess.*leftMaskStack(:,:,imageIterator);
vec_right = imageToProcess.*rightMaskStack(:,:,imageIterator);
vec_whole = vec_sensor+vec_left+vec_right;
    
end
