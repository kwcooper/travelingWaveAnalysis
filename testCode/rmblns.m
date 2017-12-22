%%%%%%%
figure; plot(root.b_lfp(1).signal(1:600));
hold on; plot(100*thetaPhs(1:600));
X=find(cycles(1:600)); plot(X,thetaPhs(X),'rx');
title("waveform")



for i=1:900, t=i*600:(i+2)*600;
plot(root.b_lfp(1).signal(t));
hold on; plot(100*thetaPhs(t));
X=find(cycles(t)); plot(X,100*thetaPhs(X),'rx');
pause;
end



CycleTs=root.b_lfp(1).ts(cycles);
Epochs = [CycleTs-0.125 CycleTs+0.125];
root.epoch=Epochs;



EegSnips=root.lfp.signal;
MeanThetaWave=nanmean(catPlus(3,EegSnips),3);


for I=1:16, 
root.active_lfp=I;
EegSnips=root.lfp.signal;
MeanThetaWave(I,:)=nanmean(catPlus(3,EegSnips),3);
end




plotLFP(rawWaves(8:2:16,:), Fs)

plotLFP(MeanThetaWave, Fs)