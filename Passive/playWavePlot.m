function playWavePlot(displacement,xlocation,ylocation)

figure(999)
for k = 1:size(displacement,3)
   plot(squeeze(displacement(ylocation,xlocation,k)));
   pause(0.2)
end
end