% converts x,z,t matrix of complex data to displacement by Loupas's
% algorithm
% INPUTS:
% complex_OCT_data: x,z,t matrix of complex OCT data (after FFT)
% sdl: vector of surface depth levels
% OCE_parameters: structure of OCE processing parameters
% Loupas_window_size: size of the window (along z) dimension
% OUTPUTS:
% Loupas_complex_OCT_data: x,z,t matrix of complex OCT data with Loupas
% phase shift cumsum
% Loupas_phase_shift: x,z,t matrix of phase differences/motion

function [Loupas_complex_OCE_data,Loupas_phase_shift] = ...
    Loupas_estimator(complex_OCE_data, sdl, OCE_parameters)

loup_x_size = OCE_parameters.vec_size(1);
loup_z_size = OCE_parameters.vec_size(2);
loup_t_size = OCE_parameters.vec_size(3);

Loupas_phase_shift = single(zeros(size(complex_OCE_data,1),size(complex_OCE_data,2),...
    size(complex_OCE_data,3)-1));
Loupas_complex_OCE_data = complex(single(zeros(size(complex_OCE_data,1),...
    size(complex_OCE_data,2)-loup_z_size,size(complex_OCE_data,3)-1)));

Loup_ts = (0:size(Loupas_phase_shift,3)-1)/30;
ts = (0:size(complex_OCE_data,3)-1)/30;
intensities = nanmean(20*log10(abs(complex_OCE_data)),3);
i_thresh_low = OCE_parameters.i_thresh_low;
close all force;
wait_bar = waitbar(0,'Processing Loupas phase shift');
for x = 1:size(complex_OCE_data,1)
    
    pos_slice = squeeze(complex_OCE_data(x,:,:));
    Loup_slice = single(zeros(size(complex_OCE_data,2)-loup_z_size,...
        length(Loup_ts)));
    
    %% Get the Loupas-estimated motion
    for t = floor(loup_t_size/2):1:size(complex_OCE_data,3)-loup_t_size
        
        for z1 = 1:OCE_parameters.how_many_depths-loup_z_size
            
            if z1 > sdl(x) - 1
                
                top_top_sum = 0; top_bot_sum = 0;
                for t1 = 1:loup_t_size-1
                    
                    for z2 = 1:loup_z_size
                        
                        if (intensities(x, z1 + z2 - 1) > i_thresh_low)
                            
                            imag_t0 = imag(pos_slice(z1+z2-1,t+t1-floor(loup_t_size/2)));
                            imag_t1 = imag(pos_slice(z1+z2-1,t+t1-floor(loup_t_size/2)+1));
                            
                            real_t0 = real(pos_slice(z1+z2-1,t+t1-floor(loup_t_size/2)));
                            real_t1 = real(pos_slice(z1+z2-1,t+t1-floor(loup_t_size/2)+1));
                            
                            top_top_sum = top_top_sum + ...
                                (imag_t0*real_t1 - real_t0*imag_t1);
                            top_bot_sum = top_bot_sum + ...
                                (real_t0*real_t1 + imag_t0*imag_t1);
                        end
                    end
                end
                
                
                bot_top_sum = 0; bot_bot_sum = 0;
                for t1 = 1:loup_t_size
                    
                    for z2 = 1:loup_z_size-1
                        
                        if (z1 >= sdl(x) - 5) && ...
                                (intensities(x, z1 + z2 - 1) > i_thresh_low)
                            
                            imag_z0 = imag(pos_slice(z1+z2-1,t+t1-floor(loup_t_size/2)+1));
                            imag_z1 = imag(pos_slice(z1+z2,t+t1-floor(loup_t_size/2)+1));
                            
                            real_z0 = real(pos_slice(z1+z2-1,t+t1-floor(loup_t_size/2)+1));
                            real_z1 = real(pos_slice(z1+z2,t+t1-floor(loup_t_size/2)+1));
                            
                            bot_top_sum = bot_top_sum + ...
                                (imag_z0*real_z1 - real_z0*imag_z1) + ...
                                (imag_z0*real_z1 - real_z0*imag_z1);
                            
                            bot_bot_sum = bot_bot_sum + ...
                                (real_z0*real_z1 + imag_z0*imag_z1) + ...
                                (real_z0*real_z1 + imag_z0*imag_z1);
                        end
                    end
                end
                
                Loupas_phase_shift(x,z1,t+1) = -1*atan(top_top_sum/top_bot_sum)/...
                    (1+(atan(bot_top_sum/bot_bot_sum))/(2*pi));
                Loup_slice(z1,t+1) = -1*atan(top_top_sum/top_bot_sum)/...
                    (1+(atan(bot_top_sum/bot_bot_sum))/(2*pi));
                
            end
            
        end
    end
    
    Loup_cumsum_slice = cumsum(Loup_slice,2);
    
    Loupas_complex_OCE_data(x,:,:) = ...
        squeeze(abs(complex_OCE_data(x,1:OCE_parameters.how_many_depths-loup_z_size,1:length(Loup_ts)))).*...
        exp(1i.*Loup_cumsum_slice);
    
    if mod(x,10) == 0
        
        phase_profile_diff = unwrap(diff(unwrap(angle(squeeze(complex_OCE_data(x,sdl(x),:))))));
        Loup_profile = (squeeze(Loupas_phase_shift(x,sdl(x),:)));
        phase_profile = (cumsum(phase_profile_diff));
        Loup_displacement_profile = (cumsum(Loup_profile));
        
        figure(1);
        set(gcf,'units','normalized','outerposition',[0.5,0.5,0.5,0.5]);
        subplot(2,1,1);
        plot(ts(1:length(phase_profile_diff)),phase_profile_diff,...
            ts(1:length(Loup_profile)),Loup_profile,'LineWidth',2);
        legend('FFT','Loupas');
        axis tight;
        title('Particle Velocity');
        set(gca,'FontSize',20);
        
        subplot(2,1,2);
        plot(ts(1:length(phase_profile)),phase_profile,...
            ts(1:length(Loup_displacement_profile)),Loup_displacement_profile,...
            'LineWidth',2);
        legend('FFT','Loupas');
        axis tight;
        title('Displacement');
        set(gca,'FontSize',20);
        pause(0.0001);
        
    end
    
    wait_bar = waitbar(x/size(complex_OCE_data,1),wait_bar,...
        ['Processing Loupas phase shift ',num2str(x/size(complex_OCE_data,1)*100),'%']);
end
close all force;
end