function standingWaveXcorr(root)




lfp = root.user_def.lfp_origData; % this will probably need to be updated for tracking
thetalfp = buttfilt(lfp',[6 10],fs,'bandpass',3)';
tAmp = root.user_def.theta_amp;

% theta AMP
tAmpXC = nan(size(tAmp,1),(size(tAmp,2)*2)-1); 
% without allocation: 3.079304 seconds | with allocation: 2.860323 seconds 
for i = 1:size(tAmp,1)
  tAmpXC(i,:) = xcorr(tAmp(1,:), tAmp(i,:)); 
end

w = 125;
figure; plot(tAmpXC(2:8,size(tAmp,2)-w:size(tAmp,2)+w)');
title('Roble xcorr theta AMP');


% Raw LFP
lfpXC = nan(size(tAmp,1),(size(lfp,2)*2)-1); 
for i = 1:size(lfp,1)
  lfpXC(i,:) = xcorr(lfp(1,:), lfp(i,:)); 
end

w = 250;
figure; plot(spreadLFP(lfpXC(2:8,size(lfp,2)-w:size(lfp,2)+w))');
title('Roble xcorr raw data ');


% theta lfp
thetaXC = nan(size(thetalfp,1),(size(thetalfp,2)*2)-1); 
for i = 1:size(lfp,1)
  thetaXC(i,:) = xcorr(thetalfp(1,:), thetalfp(i,:)); 
end

w = 250;
figure; plot(spreadLFP(thetaXC(2:8,size(thetalfp,2)-w:size(thetalfp,2)+w))');
title('Roble xcorr theta');


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