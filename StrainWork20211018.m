
for k = 1:size(vec_phase_diff,3)
    figure(1);
    imagesc(vec_phase_diff(:,:,k)')
    colormap(jet)
    title(num2str(k))
    pause(0.2)
end

for frame = 1:size(IQData,3)-1
    complex_phase_diff_comp(:,:,frame) = IQData(:,:,1).*conj(IQData(:,:,frame+1));
end

IQData = IQData(1:lowerBound,:,:);
sdl = ones([1 size(IQData,1)]);
VMIQ = permute(complex_phase_diff_comp,[2 1 3]);
[vector_complex_OCE_data,vec_phase_diff_2] = ...
    vec_meth_snr(VMIQ,sdl,20);
vec_phase_diff_2 = permute(vec_phase_diff_2,[2 1 3]);

figure;
for k = 1:size(vec_phase_diff_2,3)
    figure(1)
    imagesc(vec_phase_diff_2(:,:,k));
    colormap(jet)
    colorbar
    pause(0.2)
end
VMIQ = permute(vector_complex_OCE_data,[2 1 3]);

figure;
for k = 1:size(vec_phase_diff_2,3)
    figure(1)
plot(unwrap(squeeze(vec_phase_diff_2(:,150,k)),'b')); hold on;
plot(unwrap(squeeze(vec_phase_diff_2(:,400,k)),'r')); hold off;
    pause(0.2)
end

for k = 1:size(dis_cumsum,3)
  BScanStrainAllFast(:,:,k) = strain_2D_slash_no_sensor(baseFolder,Parameters.delta_z,...
dis_cumsum(:,:,k),25,75);
    
end


VPD_smoothed = permute(VPD_smoothed,[2 1 3]);

selectedSurface = 200;
k=20
% TEST = vec_phase_diff_2(selectedSurface:end,:,k)';
TEST = VPD_smoothed(selectedSurface:end,:,k);
TEST2 = phase_unwrap(TEST);
figure;
imagesc(TEST2);
figure;
plot(unwrap(squeeze(TEST2(:,100)),'b')); hold on;
plot(unwrap(squeeze(TEST2(:,450)),'r')); hold off;


figure; himage = imagesc(mat2gray(abs(IQData(selectedSurface:end,:,1)))); caxis([0 0.1])
colormap(gray);
  
h = impoly(gca);
PosTime = wait(h);

BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));
BW2(BW2>0) = 1;

figure; imagesc(BW2)

TESTLINE = TEST2;
TESTLINE(BW2~=1) = NaN;
figure; imagesc(TESTLINE)

TESTLINE = mean(TESTLINE,1,'omitnan');
figure; plot(TESTLINE)


%% Normalize phase diff to surface

figure; himage = imagesc(mat2gray(abs(IQData(selectedSurface:end,:,1)))); caxis([0 0.1])
colormap(gray);
[~,y] = ginput(1);

TESTABS = abs(TEST);

figure; plot(mean(TESTABS(round(y),:),1))
TESTFITTING = smooth(mean(TESTABS(round(y),:),1),33);
figure; plot(TESTFITTING)

XFITTING = 1:length(TESTFITTING);

TEST3 = phase_unwrap(TESTABS - TESTFITTING');
TEST3 = phase_unwrap(TESTABS - TESTLINE);
for lineIterator = round(y):size(TESTABS,1)
    newMatrix(lineIterator,:) = phase_unwrap(TESTABS(lineIterator,:)-TESTFITTING');
end

figure; imagesc(newMatrix)
figure; plot(TEST2(100,:)); hold on; plot(TEST3(100,:)); hold off; 

figure; imagesc(TEST2); figure; imagesc(TEST3)

% TEST3(1:selectedSurface,:) = [];
  TEST3STRAIN = strain_2D_robust_no_sensor2(Parameters.delta_z,...
    TEST3,100);  

% TEST3STRAINFAST = strain_2D_slash_no_sensor(baseFolder,Parameters.delta_z,...
% abs(TEST3),50,100);

% figure; imagesc(TEST3STRAINFAST)
figure; imagesc(mat2gray(TEST3STRAIN'));
caxis([0.5 0.8])

figure; imagesc(TEST3STRAIN');


%% YM stress calc
Instron = [58.81076 54.06557 56.22892 52.13362];
Instron = mean(Instron);
stiffInstron = [281.64155 277.05164 267.3551 273.22116];
stiffInstron = mean(stiffInstron);
softInstron = [140.13218 134.14949 138.51699 135.98369];
softInstron = mean(softInstron);
figure; himage = imagesc(mat2gray(abs(TEST3STRAIN'))); caxis([0 0.1])
colormap(gray);
  
h = impoly(gca);
PosTime = wait(h);

BW = createMask(h,himage);
BW2 = ceil(imgaussfilt(single(BW),3));
BW2(BW2>0) = 1;

figure; imagesc(BW2)

TESTLINE = TEST3STRAIN';
TESTLINE(BW2~=1) = NaN;
figure; imagesc(TESTLINE)

TESTLINE = mean(TESTLINE,1,'omitnan');
figure; plot(TESTLINE)

% E = epsilon / sigma
% YM = stress / strain

stressCalculated = Instron .* TESTLINE;
TEST3STRAIN = TEST3STRAIN';
stressCalculated2(1:250) = median(stressCalculated(1:250));
stressCalculated2(250:length(stressCalculated)) = median(stressCalculated(250:length(stressCalculated)));
figure; plot(stressCalculated2)
for k = 1:size( TEST3STRAIN,1)
YMCompression(k,:) = (stressCalculated2 ./ TEST3STRAIN(k,:));
end
figure; imagesc(YMCompression)
caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])

disp('Select bad region')
badRegion = roipoly;

leftBlock = roipoly;
rightBlock = roipoly;
YMCompression(badRegion==1)=NaN;
leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlock~=1) = NaN;
rightYM(rightBlock~=1)=NaN;
figure; imagesc(leftYM);
figure; imagesc(rightYM)

YMLeftPhantom = mean(rmoutliers(leftYM),'omitnan')
YMRightPhantom = mean(rmoutliers(rightYM),'omitnan')

YMLeftPhantom = mean(rmoutliers(abs(YMCompression(leftBlock))),'omitnan')
YMRightPhantom = mean(rmoutliers(abs(YMCompression(rightBlock))),'omitnan')

YMLNegRem = YMCompression;
YMLNegRem(YMLNegRem<0) = NaN;
YMRNegRem = YMLNegRem;

YMLeftPhantomNegRem = mean(rmoutliers(abs(YMLNegRem(leftBlock))),'omitnan')
YMRightPhantomNegRem = mean(rmoutliers(abs(YMRNegRem(rightBlock))),'omitnan')


figure; imagesc(YMCompression)
caxis([quantile(YMCompression(:),0.25) quantile(YMCompression(:),0.75)])



leftBlockTop = roipoly;
rightBlockTop = roipoly;

leftYM = YMCompression;
rightYM = YMCompression;
leftYM(leftBlockTop~=1) = NaN;
rightYM(rightBlockTop~=1)=NaN;
figure; imagesc(leftYM);
figure; imagesc(rightYM)

YMLeftPhantom = mean(rmoutliers(leftYM),'omitnan')
YMRightPhantom = mean(rmoutliers(rightYM),'omitnan')

YMLeftPhantom = mean(rmoutliers(abs(YMCompression(leftBlockTop))),'omitnan')
YMRightPhantom = mean(rmoutliers(abs(YMCompression(rightBlockTop))),'omitnan')

YMLNegRem = YMCompression;
YMLNegRem(YMLNegRem<0) = NaN;
YMRNegRem = YMLNegRem;

YMLeftPhantomNegRem = mean(rmoutliers(abs(YMLNegRem(leftBlockTop))),'omitnan')
YMRightPhantomNegRem = mean(rmoutliers(abs(YMRNegRem(rightBlockTop))),'omitnan')

