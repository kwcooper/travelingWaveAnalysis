function standingWaveXcorr(ratStruct,type) % r
% Computes the cross corrolation; makes figs
%
% Amp types: tAmp tAmpNormCH1 tAmpNorm 
% LFP types: rawLFP thetaLFP thetaLFPNorm
%
% TWAVE analysis wrapper file. KWC 180811
% TD: add saving info for figs

%fs = root.user_def.lfp_fs;
fs = 500;
freqBand = [5 9];

%lfp = root.user_def.lfp_origData; % this will probably need to be updated for tracking
lfp = ratStruct.lfp;
thetaLFP = buttfilt(lfp',freqBand,fs,'bandpass',3)'; % filter signal for theta 
%tAmp = root.user_def.theta_amp; % grab amplitude
tAmp = ratStruct.thetaAmp;

savePath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', 'twImgDir');
%%
switch type
  
  case 'tAmp' % theta AMP
    tAmpXC = nan(size(tAmp,1),(size(tAmp,2)*2)-1);
    % without allocation: 3.079304 seconds | with allocation: 2.860323 seconds
    for i = 1:size(tAmp,1)
      tAmpXC(i,:) = xcorr(tAmp(1,:), tAmp(i,:));
    end
    w = 125; figure; plot(tAmpXC(2:size(ratStruct.lfp,1),size(tAmp,2)-w:size(tAmp,2)+w)');
    titleStr = [ratStruct.name, ' xcorr theta AMP']; title(titleStr);
    
    
  case 'tAmpNormCH1' % theta AMP Normalized to CH1
    tAmpNorm1 = tAmp / max(tAmp(1,:)')';
    tAmpNorm1XC = nan(size(tAmp,1),(size(tAmp,2)*2)-1);
    %figure; plot(spreadLFP(tAmpXCNorm1(1:5,1:10000))')
    % without allocation: 3.079304 seconds | with allocation: 2.860323 seconds
    for i = 1:size(tAmpNorm1,1)
      tAmpNorm1XC(i,:) = xcorr(tAmpNorm1(1,:), tAmpNorm1(i,:));
    end
    
    % normalize to channel one to fix axis shift
    tAmpNorm1XC = tAmpNorm1XC';
    w = 50; 
    tA = tAmpNorm1XC(((ceil((size(tAmpNorm1XC,1)/2))-w)):(ceil((size(tAmpNorm1XC,1)/2))+w),:); % 1322036:1322136 roble
    %tA = tAmpXC( max(tAmpXC(:,1))-50:max(tAmpXC(:,1))+50 , :);
    x = [min(tA,[],1);max(tA,[],1)];
    b = bsxfun(@minus,tA,x(1,:));
    b = bsxfun(@rdivide,b,diff(x,1,1));
%     for i = 1:size(tAmpNorm1XC)
%       b(i,:) = tAmpNorm1XC(i,:) ./ (max(tAmpNorm1XC(1,:)) - min(tAmpNorm1XC(1,:)));
%     end
    %w = 125; figure; plot(tAmpNorm1XC(2:size(ratStruct.lfp,1),size(tAmpNorm1,2)-w:size(tAmpNorm1,2)+w)');
    %figure; w = 50; plot(b(:,((size(b,2)/2)-w):(size(b,2)/2)+w)');
    figure; plot(spreadLFP(b')')
    titleStr = [ratStruct.name, ' xcorr theta AMP Normalized to CH1']; title(titleStr);
    
    
  case 'tAmpNorm' % theta AMP Normalized
    tAmpNorm = tAmp ./ max(tAmp')';
    tAmpNormXC = nan(size(tAmpNorm,1),(size(tAmpNorm,2)*2)-1);
    % without allocation: 3.079304 seconds | with allocation: 2.860323 seconds
    for i = 1:size(tAmpNorm,1)
      tAmpNormXC(i,:) = xcorr(tAmpNorm(1,:), tAmpNorm(i,:));
    end
    w = 125; figure; plot(tAmpNormXC(2:size(ratStruct.lfp,1),size(tAmpNorm,2)-w:size(tAmpNorm,2)+w)');
    titleStr = [ratStruct.name, ' xcorr theta AMP']; title(titleStr);

    
    
  case 'rawLFP' % Raw LFP
    lfpXC = nan(size(tAmp,1),(size(lfp,2)*2)-1);
    for i = 1:size(lfp,1)
      lfpXC(i,:) = xcorr(lfp(1,:), lfp(i,:));
    end
    w = 250; figure; plot(spreadLFP(lfpXC(2:size(ratStruct.lfp,1),size(lfp,2)-w:size(lfp,2)+w))');
    titleStr = [ratStruct.name, ' xcorr raw data ']; title(titleStr);

    
    
  case 'thetaLFP' % theta lfp
    thetaXC = nan(size(thetaLFP,1),(size(thetaLFP,2)*2)-1);
    for i = 1:size(lfp,1)
      thetaXC(i,:) = xcorr(thetaLFP(1,:), thetaLFP(i,:));
    end
    w = 250; figure; plot(spreadLFP(thetaXC(2:size(ratStruct.lfp,1),size(thetaLFP,2)-w:size(thetaLFP,2)+w))');
    titleStr = [ratStruct.name, ' xcorr theta']; title(titleStr);

    
    
  case 'thetaLFPNorm' % theta lfp normalized
    thetalfp_norm = thetaLFP ./ max(thetaLFP')';
    
    thetaXCNorm = nan(size(thetalfp_norm,1),(size(thetalfp_norm,2)*2)-1);
    for i = 1:size(lfp,1)
      thetaXCNorm(i,:) = xcorr(thetalfp_norm(1,:), thetalfp_norm(i,:));
    end
    w = 250; figure; plot(spreadLFP(thetaXCNorm(2:size(ratStruct.lfp,1),size(thetalfp_norm,2)-w:size(thetalfp_norm,2)+w))');
    titleStr = [ratStruct.name, ' xcorr theta Normalized 5-9']; title(titleStr);
    
  otherwise fprintf('StandingWaveXcorr: Please enter a type...\n')
%%

if 0
% save the data (tAmpNorm1XC)
saveName = 'xcorrData.mat';
savepath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave',saveName);
save(savepath, 'tAmpNorm1XC', '-v7.3');
end

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