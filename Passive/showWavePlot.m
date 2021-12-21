function showWavePlot(displacement,ylocation,xlocation)

figure(999)
   plot(squeeze(displacement(ylocation,xlocation,:)));
end