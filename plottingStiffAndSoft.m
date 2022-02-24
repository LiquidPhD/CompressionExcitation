figure; imagesc(mat2gray(abs(IQData(:,:,end))))

figure; plot(squeeze(vec_phase_diff(1:635,455,2)),'r'); hold on; plot(squeeze(vec_phase_diff(1:635,455,end)),'b'); hold off; 

figure; plot(squeeze(vec_phase_diff(1:635,130,2)),'r'); hold on; plot(squeeze(vec_phase_diff(1:635,130,end)),'b'); 
vline(270);
hold off; 

figure; plot(squeeze(vec_phase_diff(1:270,130,2)),'r'); hold on; plot(squeeze(vec_phase_diff(1:270,130,end)),'b'); 
figure; plot(squeeze(vec_phase_diff(271:635,130,2)),'r'); hold on; plot(squeeze(vec_phase_diff(271:635,130,end)),'b'); 

displacementData = cumsum(vec_phase_diff);
deltaZ = Parameters.delta_z;
figure;
h1 = plot(squeeze(displacementData(450,455,:))*deltaZ,'r'); 
hold on; 
h2 = plot(squeeze(displacementData(450,140,:))*deltaZ,'b')
h3 = plot(squeeze(displacementData(112,455,:))*deltaZ,'k'); 
h4 = plot(squeeze(displacementData(112,140,:))*deltaZ,'g')
 hold off;
 b = [h1 h2 h3 h4]
 legend(b,'Right','Left','Right Sensor','Left Sensor')
  legend(b,'Stiff','Soft','Stiff Sensor','Soft Sensor')
