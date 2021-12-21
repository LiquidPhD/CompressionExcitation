%% Load workspace or process ultrasound elastography data
function processUSEForWorkspace(folderIndex,folders)
 if isempty((rdir([folders(folderIndex).folder,filesep,'workspace.mat'])))
        folder = folders(folderIndex).folder;
        load([folders(folderIndex).folder,filesep,'IData.mat'])
        load([folders(folderIndex).folder,filesep,'QData.mat'])
        load([folders(folderIndex).folder,filesep,'Parameter.mat'])
        IData = squeeze(IData);
        QData = squeeze(QData);
        IQData = complex(IData,QData);
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
        IQDataClip = IQData(topClip:bottomClip,:,:);
        
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
        
        SGO = 3;
        SGFL = 11;
        for sgolayfiltindex = 1:size(realIQ,3)
            realIQ(:,:,sgolayfiltindex) = sgolayfilt(double(realIQ(:,:,sgolayfiltindex)),SGO,SGFL);
            imagIQ(:,:,sgolayfiltindex) = sgolayfilt(double(imagIQ(:,:,sgolayfiltindex)),SGO,SGFL);
        end
        clearvars displacement clippedBScan
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
        
        % figure; imagesc(clippedBScan(:,:,1));
        % colormap(gray)
        %
        % figure;
        % for k = 1:size(disp,3)
        %     imagesc(disp(:,:,k))
        %     colormap(fireice)
        %     pause(0.2)
        % end
        %
        
        
        p = twomodegauss(0.5,0.2,0.0001,0.0001,1,0.01,0) % Actually a 1-mode Gaussian...
        for dispIndex = 1:size(displacement,3)
            TEST = mat2gray(displacement(:,:,dispIndex));
            g = histeq(TEST,p);
            displacement(:,:,dispIndex) = g;
        end
        % Time axis.
        parfor imageDepth = 1:size(displacement,1)
            TEST = squeeze(displacement(imageDepth,:,:));
            Y = fft2(TEST);
            %         FFTFILTERBOT = imresize(FILT,[size(Y,1) size(Y,2)])
            %         FFTFILTERBOT(FFTFILTERBOT<0.1) = 0;
            %         FFTFILTERBOT(:,round(size(FFTFILTERBOT,2)/2):end) = 0;
            FFTFILTERRIGHT = ones([size(Y,1) size(Y,2)]);
            FFTFILTERRIGHT(1:round(size(FFTFILTERRIGHT,1)/2),1:round(size(FFTFILTERRIGHT,2)/2)) = 0;
            FFTFILTERRIGHT(round(size(FFTFILTERRIGHT,1)/2)+1:end,round(size(FFTFILTERRIGHT,2)/2)+1:end) = 0;
            FFTFILTERLEFT = fliplr(FFTFILTERRIGHT);
            NEWTOP = FFTFILTERLEFT .* Y;
            NEWBOT = FFTFILTERRIGHT .* Y;
            newleftifft = ifft2(ifftshift(NEWTOP));
            newrightifft = ifft2(ifftshift(NEWBOT));
            newDispLEFT(imageDepth,:,:) = newleftifft;
            newDispRIGHT(imageDepth,:,:) = newrightifft;
        end
        
        [Nz,Nx,Nt]= size(displacement);         % The dimensions of data in the z axis, the x axis and time.
        zaxis = linspace(0,(Nz-1)*delta_z,Nz)*1e3;                      %(mm) Aixial axis.
        xaxis = linspace(-(Nx-1)/2*delta_x,(Nx-1)/2*delta_x,Nx)*1e3;    %(mm) Lateral axis.
        taxis = linspace(0,(Nt-1)*delta_t,Nt);
        
        filename = [folders(folderIndex).folder,filesep,'threeViews.gif'];
        pauseAmount = 0.2;
        figure;
        for gifIndex = 1:size(displacement,3)
            subplot(3,1,1)
            imagesc(xaxis,zaxis,im2uint8((displacement(:,:,gifIndex))))
            subplot(3,1,2)
            imagesc(xaxis,zaxis,im2uint8(abs(newDispLEFT(:,:,gifIndex))))
            subplot(3,1,3)
            imagesc(xaxis,zaxis,im2uint8(abs(newDispRIGHT(:,:,gifIndex))))
            %     colormap(fireice)
            pause(0.2)
            if gifIndex == 1
                gif(filename,'DelayTime',pauseAmount)
            else
                gif
            end
        end
        
        filename = [folders(folderIndex).folder,filesep,'threeViewsFireIce.gif'];
        pauseAmount = 0.2;
        figure;
        for gifIndex = 1:size(displacement,3)
            subplot(3,1,1)
            imagesc(xaxis,zaxis,im2uint8((displacement(:,:,gifIndex))))
            subplot(3,1,2)
            imagesc(xaxis,zaxis,im2uint8(abs(newDispLEFT(:,:,gifIndex))))
            subplot(3,1,3)
            imagesc(xaxis,zaxis,im2uint8(abs(newDispRIGHT(:,:,gifIndex))))
            colormap(fireice)
            pause(0.2)
            if gifIndex == 1
                gif(filename,'DelayTime',pauseAmount)
            else
                gif
            end
        end
        save([folders(folderIndex).folder,filesep,'workspace.mat']);
 else
        % Do nothing, since we will load the file in the next step.
 end
end