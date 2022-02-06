function [particleVelocity,BScan] = calculateParticleVelocity(IQData,Parameters)
c = Parameters.c;
M = Parameters.M;
N = Parameters.N;
fc = Parameters.fc;
% Get dimensions
ROI_dim_z = size(IQData,1)-M;
ROI_dim_x = size(IQData,2);
ROI_dim_t = size(IQData,3)-N;

IQ = permute(IQData, [1 3 2] );    % For more convenient access!!!

tic
realIQ = real(IQ);
imagIQ = imag(IQ);
% BScan = abs(IQData(:,:,1));

% Savitzky-Golay filtering to smooth data
SGO = 3;
SGFL = 11;
for sgolayfiltindex = 1:size(realIQ,3)
    realIQ(:,:,sgolayfiltindex) = sgolayfilt(double(realIQ(:,:,sgolayfiltindex)),SGO,SGFL);
    imagIQ(:,:,sgolayfiltindex) = sgolayfilt(double(imagIQ(:,:,sgolayfiltindex)),SGO,SGFL);
end
particleVelocity = zeros([ROI_dim_z ROI_dim_x ROI_dim_t]);
%     parfor_progress(ROI_dim_z);


% Calculate particle velocity
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
        particleVelocity(ii_z,:,ii_t) = abs(c/(4*pi*fc)*atan(uu./ud)./(1+atan(du./dd)/(2*pi)));
        BScan(ii_z,:,ii_t) = abs(IQData( ii_z,:,1 ));
    end
    %         parfor_progress;
end
%     parfor_progress(0);
toc