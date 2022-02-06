
function [vector_complex_OCE_data,vec_phase_diff,particleVelocity] = ...
    vec_meth_Loupas(complex_OCE_data,sdl,vec_size)

close all force;
smooth_kernel = single(ones(vec_size).*(1/(vec_size^2)));
vec_phase_diff = single(zeros(size(complex_OCE_data)));

for t = 1:size(complex_OCE_data,3)
    
%     cdata1 = conv2(squeeze(complex_OCE_data(:,:,t-1)),smooth_kernel,'same');
%     cdata2 = conv2(squeeze(complex_OCE_data(:,:,t)),smooth_kernel,'same');
%     phase_diff_conj = conj(cdata1).*cdata2;
phase_diff_conj = complex_OCE_data(:,:,t);
    
    %% averaging vertically
    % 0.04s for loop, 0.3s parfor loop
    temp_1 = complex(single(zeros(size(phase_diff_conj,1),...
        size(phase_diff_conj,2),vec_size)));
    for x = 1:vec_size
        temp_1(:,:,x) = circshift(phase_diff_conj,x,2);
    end
    temp_2 = sum(temp_1,3)./sum(abs(temp_1),3);
    clear temp_1;
    
    %% averaging temporally
    temp_3 = complex(single(zeros(size(phase_diff_conj,1),...
        size(phase_diff_conj,2),vec_size)));
    for x = 1:vec_size
        temp_3(:,:,x) = circshift(temp_2,x,3);
    end
    temp_4 = sum(temp_3,3)./sum(abs(temp_3),3);
    clear temp_3;
    
%     vec_phase_diff(:,:,t) = angle(temp_4);

    vec_phase_diff(:,:,t) = angle(temp_4);
    
end
vec_phase_cumsum = cumsum(vec_phase_diff,3);
vector_complex_OCE_data = abs(complex_OCE_data).*exp(1i.*vec_phase_cumsum);

% for x = 1:10:size(vec_phase_diff,1)
%     
%     phase_profile_diff = unwrap(diff(angle(squeeze(complex_OCE_data(x,sdl(x),:)))));
%     phase_profile = cumsum(phase_profile_diff);
%     vec_profile_diff = unwrap(diff(angle(squeeze(vector_complex_OCE_data(x,sdl(x),:)))));
%     vec_profile = cumsum(vec_profile_diff);
%     
% %     phaseSNR = snr(phase_profile);
% %     vecSNR = snr(vec_profile);
%     figure(1);
%     set(gcf,'units','normalized','outerposition',[0.5,0.5,0.5,0.5]);
%     subplot(2,1,1);
%     plot(phase_profile,'LineWidth',2);
%     hold on;
%     plot(vec_profile,'LineWidth',2);
%     hold off; axis tight; title('Phase Profiles');
%     legend('Original','Vec. Meth.');
%     set(gca,'FontSize',20);
%     
%     subplot(2,1,2);
%     plot(phase_profile_diff,'LineWidth',2);
%     hold on;
%     plot(vec_profile_diff,'LineWidth',2);
%     hold off; axis tight; title('Phase Difference');
%     legend('Original','Vec. Meth.');
%     set(gca,'FontSize',20);
%     pause(0.001);
%     
% end
Parameters = evalin('base','Parameters')
IQData = vector_complex_OCE_data;
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
end