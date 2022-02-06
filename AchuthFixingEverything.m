figure; imagesc(20*log10(abs(IQData(:,:,2)))); colormap(gray);

vec_image=vec_phase_diff(:,:,2);

imagesc(vec_image); colormap(jet);

images=20*log10(abs(IQData));

for pos=1:size(mask4,2)
    sensor_sdl(pos)=find(mask4(:,pos,2)==1,1,'last');
end


for pos=1:size(mask7,2)
    sensor_bdl(pos)=find(mask7(:,pos)==1,1);
end

mask2=zeros(size(mask6(:,:,1)));
for pos=1:size(mask4,2);
    mask2(sensor_bdl(pos):end,pos)=1;
end

mask=zeros(size(mask6(:,:,1)));
for pos=1:size(mask4,2);
    mask(sensor_dl(pos):sensor_bdl(pos),pos)=1;
end


figure; imagesc(unwrapped_sample); 
figure; imagesc(abs(unwrapped_sample));
figure; imagesc(images(:,:,2)); colormap(gray);

vec_image_sensor=vec_image.*mask;
vec_image_sample=vec_image.*mask2;

unwrapped_sensor=phase_unwrap(vec_image_sensor);

for pos=1:size(vec_image_sensor,2)
    surf_dis(pos)=nanmedian(unwrapped_sensor(sensor_sdl(pos):sensor_sdl(pos)+10,pos));
    corr_sensor(:,pos)=(unwrapped_sensor(:,pos)-surf_dis(pos)).*mask(:,pos);
end


imagesc(corr_sensor);

unwrapped_sample=phase_unwrap(vec_image_sample);

for pos=1:size(vec_image_sample,2)
    surf_dis(pos)=nanmedian(unwrapped_sample(sensor_bdl(pos):sensor_bdl(pos)+10,pos));
    corr_sample(:,pos)=(unwrapped_sample(:,pos)-surf_dis(pos)).*mask2(:,pos);
end

strain_sensor=strain_calculation2(corr_sensor,zaxis(2),1,100);
figure;
imagesc(abs(strain_sensor)/1000); colormap(hot); colorbar; caxis([0,0.5]);

strain_sample=strain_calculation2(corr_sample,zaxis(2),1,100);
figure;
imagesc(abs(strain_sample)/1000); colormap(hot); colorbar; caxis([0,0.5]);

strain_whole=strain_sensor.*mask+strain_sample.*mask2;
figure;
imagesc(abs(strain_whole)/1000); colormap(hot); colorbar; caxis([0,0.5]);



mask3=modefilt(mask+mask2);

unwrapped_full=phase_unwrap(vec_image.*mask3);

for pos=1:size(vec_image_sample,2)
    surf_dis(pos)=median(unwrapped_full(sensor_dl(pos)+10,pos));
    corr_full(:,pos)=(unwrapped_full(:,pos)-surf_dis(pos)).*mask3(:,pos);
end

figure; imagesc(corr_full);
caxis([0,4]);

strain_full=strain_calculation2(corr_full,zaxis(2),1,100);

figure; imagesc(strain_full/1000); colormap(hot); caxis([0,0.5])