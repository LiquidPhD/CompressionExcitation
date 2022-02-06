%% calculates the phase difference using the vector method and adjusts the
% complex OCT data with vector method phase difference cumulative sum
% Matveyev et al., "Vector method for strain estimation in phase-sensitive
% optical coherence elastography." Laser Phys Lett, 2018
% Inputs:
% complex_OCE_data: complex OCT data matrix (x,z,t)
% sdl: vector of the surface depth levels
% vec_size: size of the vector kernel in pixels for z,t averaging
%
% Outputs:
% vector_complex_OCE_data: complex OCT data matrix after phase adjusted by
% vector method (x,z,t)
% vec_phase_diff: phase difference profiles (x,z,t)

function [vector_complex_OCE_data,vec_phase_diff] = ...
    vec_meth_snr(complex_OCE_data,sdl,vec_size)

close all force;
smooth_kernel = single(ones(vec_size).*(1/(vec_size^2)));
vec_phase_diff = single(zeros(size(complex_OCE_data)));

for t = 2:size(complex_OCE_data,3)
    
    cdata1 = conv2(squeeze(complex_OCE_data(:,:,t-1)),smooth_kernel,'same');
    cdata2 = conv2(squeeze(complex_OCE_data(:,:,t)),smooth_kernel,'same');
    phase_diff_conj = conj(cdata1).*cdata2;
    
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

close all force;
end