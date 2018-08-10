function standingWaveXcorr(root)
% computes the cross corrolation 

% TD: add saving info for figs
% 


lfp = root.user_def.lfp_origData; % this will probably need to be updated for tracking
fs = root.user_def.lfp_fs;
freqBand = [5 9];

% filter signal for theta & grab amplitude
thetalfp = buttfilt(lfp',freqBand,fs,'bandpass',3)';
tAmp = root.user_def.theta_amp;



%% theta AMP
tAmpXC = nan(size(tAmp,1),(size(tAmp,2)*2)-1); 
% without allocation: 3.079304 seconds | with allocation: 2.860323 seconds 
for i = 1:size(tAmp,1)
  tAmpXC(i,:) = xcorr(tAmp(1,:), tAmp(i,:)); 
end

w = 125;
figure; plot(tAmpXC(2:8,size(tAmp,2)-w:size(tAmp,2)+w)');
title('Roble xcorr theta AMP');


%% theta AMP Normalized to CH1
tAmpNorm1 = tAmp / max(tAmp(1,:)')';
tAmpXC = nan(size(tAmp,1),(size(tAmp,2)*2)-1); 
%figure; plot(spreadLFP(tAmpXCNorm1(1:5,1:10000))')
% without allocation: 3.079304 seconds | with allocation: 2.860323 seconds 
for i = 1:size(tAmpNorm1,1)
  tAmpNorm1XC(i,:) = xcorr(tAmpNorm1(1,:), tAmpNorm1(i,:)); 
end

w = 125;
figure; plot(tAmpNorm1XC(2:8,size(tAmpNorm1,2)-w:size(tAmpNorm1,2)+w)');
title('Roble xcorr theta AMP Normalized to CH1');


%% theta AMP Normalized
tAmpNorm = tAmp ./ max(tAmp')';
tAmpNormXC = nan(size(tAmpNorm,1),(size(tAmpNorm,2)*2)-1); 
% without allocation: 3.079304 seconds | with allocation: 2.860323 seconds 
for i = 1:size(tAmpNorm,1)
  tAmpNormXC(i,:) = xcorr(tAmpNorm(1,:), tAmpNorm(i,:)); 
end

w = 125;
figure; plot(tAmpNormXC(2:8,size(tAmpNorm,2)-w:size(tAmpNorm,2)+w)');
title('Roble xcorr theta AMP');


%% Raw LFP
lfpXC = nan(size(tAmp,1),(size(lfp,2)*2)-1); 
for i = 1:size(lfp,1)
  lfpXC(i,:) = xcorr(lfp(1,:), lfp(i,:)); 
end

w = 250;
figure; plot(spreadLFP(lfpXC(2:8,size(lfp,2)-w:size(lfp,2)+w))');
title('Roble xcorr raw data ');


%% theta lfp
thetaXC = nan(size(thetalfp,1),(size(thetalfp,2)*2)-1); 
for i = 1:size(lfp,1)
  thetaXC(i,:) = xcorr(thetalfp(1,:), thetalfp(i,:)); 
end

w = 250;
figure; plot(spreadLFP(thetaXC(2:8,size(thetalfp,2)-w:size(thetalfp,2)+w))');
title('Roble xcorr theta');


%% theta lfp normalized
thetalfp_norm = thetalfp ./ max(thetalfp')';

thetaXCNorm = nan(size(thetalfp_norm,1),(size(thetalfp_norm,2)*2)-1); 
for i = 1:size(lfp,1)
  thetaXCNorm(i,:) = xcorr(thetalfp_norm(1,:), thetalfp_norm(i,:)); 
end

w = 250;
figure; plot(spreadLFP(thetaXCNorm(2:8,size(thetalfp_norm,2)-w:size(thetalfp_norm,2)+w))');
title('Roble xcorr theta Normalized 5-9');

%%

% save the data (tAmpNorm1XC)
saveName = 'xcorrData.mat';
savepath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave',saveName);
save(savepath, 'tAmpNorm1XC', '-v7.3');


end
%%

% a = sin(linspace(0,50,1000));
% b = sin((linspace(0,50,1000)) + pi/2);
% c = xcorr(a,a);
% d = xcorr(a,b);
% 
% figure; 
% subplot(4,1,1); plot(a);
% subplot(4,1,2); plot(b); 
% subplot(4,1,3); plot(c);% ylim([-200 200]);
% subplot(4,1,4); plot(d); %ylim([-200 200]);