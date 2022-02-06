YMSoft = 16.735;
YMStiff = 49.87333

for k = 1:size(leftYM,3)
   figure(1)
   subplot(3,1,1)
   plot(nonzeros(squeeze(leftYM(:,:,k))));
   ylim([0 100])
   title('Left')
   hold on;
   hline(YMSoft);
   hold off;
   subplot(3,1,2)
   plot(nonzeros(squeeze(rightYM(:,:,k))));
   ylim([0 100])
   title('Right')
   hold on;
   hline(YMStiff);
   hold off;
   subplot(3,1,3)
   plot(nonzeros(rightYM(:,:,k))./nonzeros(leftYM(:,:,k)));
   hold on;
   hline(YMStiff/YMSoft,'g');
   hold off; 
   pause
end

depth = 70
figure;
subplot(4,1,1)
plot(nonzeros(squeeze(leftYM(:,depth,:))))
title(['Left at depth = ',num2str(depth)])
ylim([0 100])
subplot(4,1,2)
plot(nonzeros(squeeze(rightYM(:,depth,:))))
title(['Right at depth = ', num2str(depth)]);
ylim([0 100])
subplot(4,1,3)
%    plot(nonzeros(rightYM(:,170,:))./nonzeros(leftYM(:,170,:)));
plot(nonzeros(squeeze(rightYM(:,depth,:))/squeeze(leftYM(:,depth,:))));
hold on;
  hold off; 
  subplot(4,1,4)
  plot(nonzeros(rightYM(:,depth,:))./nonzeros(leftYM(:,depth,:)));
  hold on
    hline(YMStiff/YMSoft,'g');
    hold off; 

for depth = 10:10:250
figure;
subplot(3,1,1)
plot(nonzeros(squeeze(leftYM(:,depth,:))))
title(['Left at depth = ',num2str(depth)])
ylim([0 100])
 hold on;
   hline(YMSoft);
   hold off;
subplot(3,1,2)
plot(nonzeros(squeeze(rightYM(:,depth,:))))
title(['Right at depth = ', num2str(depth)]);
ylim([0 100])
 hold on;
   hline(YMStiff);
   hold off;
subplot(3,1,3)
  plot(nonzeros(rightYM(:,depth,:))./nonzeros(leftYM(:,depth,:)));
  hold on
    hline(YMStiff/YMSoft,'g');
    hold off; 
    ylim([0 5])
    title('Array division')
    pause
end

for depth = 170
figure;
subplot(3,1,1)
plot(nonzeros(squeeze(leftYM(:,depth,:))))
title(['Left at depth = ',num2str(depth)])
ylim([0 100])
subplot(3,1,2)
plot(nonzeros(squeeze(rightYM(:,depth,:))))
title(['Right at depth = ', num2str(depth)]);
ylim([0 100])
subplot(3,1,3)
  plot(nonzeros(rightYM(:,depth,:))./nonzeros(leftYM(:,depth,:)));
  hold on
    hline(YMStiff/YMSoft,'g');
    hold off; 
    ylim([0 5])
    title('Array division')
    pause
end


%%09/16 test


for depth = 10:250
    
    LeftYM = nonzeros(squeeze(leftYM(:,depth,:)));
    RightYM = nonzeros(squeeze(rightYM(:,depth,:)));
    if (sum(YMSoft > LeftYM-0.15*LeftYM) > 5) && (sum(YMSoft < LeftYM + 0.15*LeftYM) > 5)
    else
        continue
    end
figure;
subplot(3,1,1)
plot(nonzeros(squeeze(leftYM(:,depth,:))))
title(['Left at depth = ',num2str(depth)])
ylim([0 100])
 hold on;
   hline(YMSoft);
   hold off;
subplot(3,1,2)
plot(nonzeros(squeeze(rightYM(:,depth,:))))
title(['Right at depth = ', num2str(depth)]);
ylim([0 100])
 hold on;
   hline(YMStiff);
   hold off;
subplot(3,1,3)
  plot(nonzeros(rightYM(:,depth,:))./nonzeros(leftYM(:,depth,:)));
  hold on
    hline(YMStiff/YMSoft,'g');
    hold off; 
    ylim([0 5])
    title('Array division')

    
end