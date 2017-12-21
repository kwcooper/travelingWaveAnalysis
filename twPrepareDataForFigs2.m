function root = twPrepareDataForFigs2(D,fs,refCh)

% create root object for use in epoching 
nSamp = size(D,2);
b_ts = linspace(0,nSamp/fs,nSamp+1); b_ts = b_ts(2:end); % cut off first to get rid of 0
root = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',b_ts,'fs_video',fs);
% store useful data in root
root.user_def.lfp_origData = D;
root.user_def.lfp_fs = fs;

% clean data and filter to epochs of good theta
%filtType = 'thetaDeltaRatio'; filtParams = 2; % works okay, some fine tuning may help further rule out low theta epochs and include 
filtType = 'thetaMag'; filtParams = [115 1000]; % works okay, some fine tuning may help further rule out low theta epochs and include 
[pctDataUsed,inds2cut] = cleanData_Intan(D(refCh,:),fs,filtType,filtParams); 
% store results in root
root.user_def.cleanData_inds2cut = inds2cut;
root.user_def.cleanData_pctDataUsed = pctDataUsed;

% extract needed variables
% Average waves needs to know cycle bounds & fs
[thetaPhs,thetaAmp,~] = extractThetaPhase(D(refCh,:),fs,'hilbert',[6 10]);
[cycles,~] = parseThetaCycles(thetaPhs,fs,[6 10]);
root.user_def.theta_phs = thetaPhs;
root.user_def.theta_amp = thetaAmp;
root.user_def.cycles = cycles;


