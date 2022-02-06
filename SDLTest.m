
bottomCut = 550
for timePoint = 1:size(IQData,3)
    if timePoint >= size(IQData,3)/2;
        bottomCut = 450;
    end
    for counter = 1:size(IQData,2)
        Aline = abs(IQData(1:bottomCut,counter,timePoint));
%         figure(1);
%         plot(Aline);
        peakLoc = find(Aline == max(Aline));
%         hold on; vline(peakLoc);
%         hold off; 
        SDL(:,counter,timePoint) = peakLoc;
    end
    TEST(:,:,timePoint) = movmean(SDL(:,:,timePoint),10);
end

for k = 1:size(IQData,3)
    figure(1)
    imagesc(abs(IQData(:,:,k)));
    hold on;
    plot(squeeze(SDL(:,:,k)),'r')
    plot(squeeze(TEST(:,:,k)),'g');
    hold off;
    pause(0.1)
end

