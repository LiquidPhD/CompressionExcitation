function [displacement] = phase2displacement_comp(vector_phase,masks,kernel_size)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% lambda=0.84;
% refrac=1.376;
displacement=zeros(size(vector_phase));
% wait_bar=waitbar(0,'Unwrapping');

%masks_vector(masks_vector==0)=NaN;

for n=1:size(displacement,3)
%      wait_bar=waitbar(n/(size(displacement,3)),wait_bar,'unwrapping'); 
  
   phase_frame=(vector_phase(:,:,n).*masks(:,:,n));
%   surf_phase=median(phase_frame(1:5,:),1);
   for lat=1:size(displacement,2)
       if lat-kernel_size<1
           lat_start=1;
       else
           lat_start=lat-kernel_size;
       end
       if lat+kernel_size>size(phase_frame,2)
           lat_end=size(phase_frame,2);
       else
           lat_end=lat+kernel_size;
       end
       temp=(squeeze(phase_frame(:,lat_start:lat_end)));
        %surf_temp=median(temp(1:3,:),1);
        %temp=temp-surf_temp;
       med_phase=median(temp,2);
       displacement(:,lat,n)=unwrap(med_phase);
   end
   med_disp=median(displacement(1:10,:,n),1);
   displacement(:,:,n)=displacement(:,:,n)-med_disp;
%   displacement(:,:,n)=medfilt2(displacement(:,:,n),[7,7]);
%    figure;
%    imagesc(displacement(:,:,n).*masks(:,:,n)); colormap(jet); caxis([-3,3]);
end
% close(wait_bar);
% scaled_displacement=(displacement.*masks_vector).*(lambda/(4*pi*refrac));

end

