function [tr_matrix] = time_reversal_FFT(input_matrix,sdl,bdl,win_size)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
[z_num,x_num,t_num]=size(input_matrix);
corr_mapping=zeros(x_num,z_num);
Depth_window=bdl;
sdl_2=sdl;%round(max(top,[],2));
wait_bar=waitbar(0,'Reversing Time');



for x1=1:x_num
    wait_bar=waitbar(x1/x_num,wait_bar,'Reversing Time');

    for depth=sdl_2(x1)+win_size:Depth_window(x1)-win_size
        temp=squeeze(input_matrix(depth-win_size:depth+win_size,x1,:));
        %      temp=squeeze(disp_diff(x1,y1,:));
        corr_sig=fft2(temp);
        corr_sig_temp=ifft2(corr_sig.*conj(corr_sig));
        %corr_sig_temp=fft_shift(corr_sig_temp);
        
        corr_val=max(corr_sig_temp(:));
        corr_mapping(x1,depth)=corr_val;
    end
end


close(wait_bar);
tr_matrix=corr_mapping';
end

