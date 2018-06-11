
% play with the theta & tracking data



function [thetaSpdMod_R] = thetaSpdAnalysis(root)
% which of the channels should we analyse? Why not all of them! muh hahahah
for i = 1:numChans
root.active_lfp = i;
[thetaSpdMod_amp, thetaSpdMod_frq, thetaSpdMod_R] = thetaSpdMod(root,[6 10],0);
end

root.active_lfp = 1;
[thetaSpdMod_amp, thetaSpdMod_frq, thetaSpdMod_R] = thetaSpdMod(root,[6 10],0);


figure; imagesc(t.vel,t.frq,t.thFrq_byVelHist); title('theta Freq by velocity');
figure; plot(t.vel,t.thAmp_byVel);
figure; plot(t.vel,t.thFrq_byVel);

end



%% theta mag stuff

tA = root.user_def.theta_amp;

tAs = spreadLFP(tA);

figure; plot(tAs(:,1:10000)');
figure; imagesc(tA(:,1:50000));


% % if we want avg amplitude...
% cycles = find(root.user_def.cycles(1,:));
% badCycles = root.user_def.cleanData_inds2cut(cycles);
% cycles(badCycles) = [];
% 
% % ud 180528: changed root.b_ts to lfp ts after addition of tracking data.  
% cycleTs=root.b_lfp(1).ts(cycles); %finds all theta cycles
% % finds bad theta cycles
% %epochSize = 0.100; %sets epoch size (moved to function argument)
% epochs = [cycleTs-epochSize cycleTs+epochSize]; % grabs epochs 
% root.epoch = epochs; 
% 
% root.b_myvar = tA'; % time across rows
% epchData = root.myvar; % this is then epoched with cmbhome
% % drop short epochs and trim long epochs until all are the same length
% nSamp = cellfun(@(c) length(c), epchData,'uni',1);
% epchLength = root.user_def.lfp_fs * 2 * epochSize;
% shortEpochs = nSamp < epchLength; epchData(shortEpochs) = [];
% epchData = cellfun(@(c) c(1:epchLength,:),epchData,'uni',0);  
% 
% % now compute the mean wave
% epchData = cat(3,epchData{:});
% avgThetaWave = mean(epchData,3);


lfp = root.user_def.lfp_origData;


pband = bandpower(lfp,500,[6 10])
ptot = bandpower(lfp,500,[0 100]);
per_power = 100*(pband/ptot)

figure; plot(pband)