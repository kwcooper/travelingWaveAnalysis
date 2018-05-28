

function [thetaSpdMod_R] = thetaSpdAnalysis(root, chans)
% which of the channels should we analyse? Why not all of them! muh hahahah

for i = 1:chans
root.active_lfp = i;
[thetaSpdMod_amp, thetaSpdMod_frq, thetaSpdMod_R] = thetaSpdMod(root,[6 10],0);
end

root.active_lfp = 1;
[thetaSpdMod_amp, thetaSpdMod_frq, thetaSpdMod_R] = thetaSpdMod(root,[6 10],0);

t = thetaSpdMod_R;
figure; imagesc(t.vel,t.frq,t.thFrq_byVelHist); title('theta Freq by velocity');
figure; plot(t.vel,t.thAmp_byVel);
figure; plot(t.vel,t.thFrq_byVel);

end
