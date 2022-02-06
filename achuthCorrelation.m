% raw_data = images;
%     image1=imgaussfilt(mat2gray(20*log10(abs(raw_data(:,:,1))),[25,60]),3);
%     image_temp=imgaussfilt(mat2gray(20*log10(abs(raw_data(:,:,frame))),[25,60]),3);
    image1=imgaussfilt(images(1:650,:,1),3);
    image2=imgaussfilt(images(1:650,:,2),3);
    
%     mask1=single(image1>0.1 & image1<0.9);
%     mask2=single(image_temp>0.1 & image1<0.9);
%     
%     output=dftregistration(fft2(mask1),fft2(mask2),100);
%     
%     deltar = output(3);
%     deltac = output(4);
%     phase = 1;
%     [nr,nc]=size(image2);
%     Nr = ifftshift((-fix(nr/2):ceil(nr/2)-1));
%     Nc = ifftshift((-fix(nc/2):ceil(nc/2)-1));
%     [Nc,Nr] = meshgrid(Nc,Nr);
%     image2 = abs(ifft2(fft2(image_temp).*exp(1i*2*pi*(deltar*Nr/nr+deltac*Nc/nc))).*exp(-1i*phase));
    
    
%     figure; imagesc(image1); colormap(gray);
%     
%     figure; imagesc(image2); colormap(gray);
    
    
    p=3;
    q=p*3;
%       p=10;
%     q=p*3;
    corr_map=zeros(size(image1));
    dis_map_x=corr_map;
    dis_map_y=corr_map;
    
    for depth=1:size(image1,1)-q
        disp(depth)
        tic
        for pos=1:size(image1,2)-p
            for p_step=0:p
                for q_step=0:q
                    numer=(image1(depth+q_step,pos+p_step)-mean(image1(depth:depth+q,pos:pos+p),'all')).*(image2(depth+q_step,pos+p_step)-mean(image2(depth:depth+q,pos:pos+p),'all'));
                    denom=sqrt((image1(depth+q_step,pos+p_step)-mean(image1(depth:depth+q,pos:pos+p),'all'))^2 + (image2(depth+q_step,pos+p_step)-mean(image2(depth:depth+q,pos:pos+p),'all'))^2);
                    temp1(q_step+1)=numer/denom;
                end
                temp1_sum=nansum(temp1);
                temp2(p_step+1)=temp1_sum;
            end
            corr_map(depth,pos)=nansum(temp2);
            
            if abs(nansum(temp2))>0.1
                roi1=image1(depth:depth+q,pos:pos+p);
                roi2=image2(depth:depth+q,pos:pos+p);
                output=dftregistration(fft2(roi1),fft2(roi2),100);
                dis_map_x(depth,pos)=output(4);
                dis_map_y(depth,pos)=output(3);
            end
        end
        toc
    end
    
    corr_mapping(:,:,frame)=corr_map;
    dis_mapping_x(:,:,frame)=dis_map_x;
    dis_mapping_y(:,:,frame)=dis_map_y;
    