function [tr_matrix] = time_reversal(input_matrix,sdl,bdl,win_size)
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
        temp=squeeze(input_matrix(depth,x1,:));
        corr_sig=xcorr(temp);
    
        corr_val=corr_sig(t_num);
        corr_mapping(x1,depth)=corr_val;
    end
end

close(wait_bar);

tr_matrix=corr_mapping';
end

