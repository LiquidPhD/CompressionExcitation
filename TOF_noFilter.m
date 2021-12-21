function [TOF_speed] = TOF_noFilter(particleVelocity,Parameters,folderIndex,position,N_radius)


delta_z = Parameters.PData.PDelta(3)*Parameters.Trans.lambda;     % (m)
delta_x = Parameters.PData.PDelta(1)*Parameters.Trans.lambda;     % (m)
delta_t = 1e-4           % (s)
disp = double(particleVelocity);
disp_z = disp(1:670,:,2:size(disp,3));  % The first frame is normally very noisy.
M = Parameters.M
% Basic parameters.
density = 1000;                   %(kg/m^3) Mass density of the medium.       [ADJUSTABLE]
[Nz,Nx,Nt]= size(disp_z);         % The dimensions of data in the z axis, the x axis and time.    
zaxis = linspace((M/2)*delta_z,(Nz-M/2)*delta_z,Nz)*1e3;        % (mm) Aixial axis.
xaxis = linspace(-(Nx-1)/2*delta_x,(Nx-1)/2*delta_x,Nx)*1e3;    % (mm) Lateral axis.
taxis = linspace(0,(Nt-1)*delta_t,Nt);                          % (s) Time axis.

% % Checking.
% disp_min = min(min(min(disp_z))); disp_max = max(max(max(disp_z)));
% figure
% for ii = 1:Nt
%     imagesc(disp_z(:,:,ii),[disp_min disp_max]/150)
%     title(ii)
%     pause(0.15)
% 	colormap('jet')
% end

%% I: PRE-PROCESSING 

% 1.1: Data smoothing using 'sgolay' filter.
disp_z_smooth = zeros(Nz,Nx,Nt);
tic();
for ii=1:Nt
    data = disp_z(:,:,ii);
    data_s = sgolayfilt(data,1,11,[],2);     % [ADJUSTABLE]
    disp_z_smooth(:,:,ii) = data_s;
end
toc();

% % Checking.
% disp_min = min(min(min(disp_z_smooth))); disp_max = max(max(max(disp_z_smooth)));
% figure
% for ii = 1:Nt
%     imagesc(disp_z_smooth(:,:,ii),[disp_min disp_max]/50)
%     title(ii)
%     pause(0.25)
% 	colormap('jet')
% end

% 1.2: Axial desampling. 
% Motivation: Normally, the axial resolution is much higher than the
% lateral resolution. For elasticity map, it is not nessary to have a high
% axial/lateral resolution ratio. 
% ratio_resolution = 1; % Define a ratio of axial/lateral resolution. [ADJUSTABLE]
% N_desamp = ceil((delta_x/delta_z)/ratio_resolution);
% Nz_desamp = ceil(Nz/N_desamp);
% disp_z_smooth_desamp = zeros(Nz_desamp,Nx,Nt);
% zaxis_desamp = zeros(1,Nz_desamp);
% tic();
% for ii = 1:Nz_desamp
%     index = (ii-1)*N_desamp+1;
%     disp_z_smooth_desamp(ii,:,:) = disp_z_smooth(index,:,:);
%     zaxis_desamp(ii) = zaxis(index);
% end
% toc();

% % Checking.
% disp_min = min(min(min(disp_z_smooth_desamp))); 
% disp_max = max(max(max(disp_z_smooth_desamp)));
% figure
% for ii = 1:Nt
%     imagesc(disp_z_smooth_desamp(:,:,ii),[disp_min disp_max]/10)
%     title(ii)
%     pause(0.1)
% 	colormap('jet')
% end

% 1.3: Data interpolation in time domain.
% Motivation: The 2018-07-04 signal processing project shows that temporal
%             interpolation can eliminate stairs of the TOF profiles.
N_interp = 5;       % Ratio of interpolation. [ADJUSTABLE]
taxis_interp = linspace(min(taxis),max(taxis),Nt*N_interp);
delta_t_interp = min(diff(taxis_interp));
tic();
disp_z_smooth_desamp = disp_z_smooth;
disp_z_smooth_desamp_interp = zeros([size(disp_z_smooth_desamp,1) size(disp_z_smooth_desamp,2) size(disp_z_smooth_desamp,3)*N_interp]);
for ii = 1:size(disp_z_smooth_desamp,1)
    signal = squeeze(disp_z_smooth_desamp(ii,:,:));
    signal_interp = interp1(taxis,signal',taxis_interp,'spline');
    disp_z_smooth_desamp_interp(ii,:,:) = signal_interp';
    num2str(ii)
end
toc();

% % Checking.
% disp_min = min(min(min(disp_z_smooth_desamp_interp))); 
% disp_max = max(max(max(disp_z_smooth_desamp_interp)));
% figure
% for ii = 1:length(taxis_interp)
%     imagesc(disp_z_smooth_desamp_interp(:,:,ii),[disp_min disp_max]/10)
%     title(ii)
%     pause(0.25)
% 	colormap('jet')
% end

%% II: PROCESSING

% % Checking
% Nz_interest = 450; % Choose axial index of interest.
% disp_xt = squeeze(disp_z_smooth_desamp_interp(Nz_interest,:,:));
% figure; imagesc(disp_xt');
% figure; surf(disp_xt');

TOF = zeros(size(disp_z_smooth_desamp_interp,1),Nx);
N_radius = 20; % Half distance for choosing reference waveform. [ADJUSTABLE]
tic();
for ii = 1:size(disp_z_smooth_desamp_interp,1)
    for jj = N_radius+1:Nx-N_radius
        waveform_left = squeeze(disp_z_smooth_desamp_interp(ii,jj-N_radius,:));
        waveform_right = squeeze(disp_z_smooth_desamp_interp(ii,jj+N_radius,:));
        [acor,lag] = xcorr(waveform_left,waveform_right);
        [~,I] = max(acor);
        TOF(ii,jj) = abs(lag(I)*delta_t_interp);
    end
end
toc();
TOF_speed = 2*N_radius*delta_x./TOF;
% elasticity = density*TOF_speed.^2/1000; 

% RestoreFigureToolBar;
% figure; imagesc(TOF_speed); colorbar; colormap('jet'); caxis([0 15])
end