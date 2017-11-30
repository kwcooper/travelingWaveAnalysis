%this is keiland's version

% update these values or copy and paste them
Rat = 'Tio';
Session = '170703_1251_CircleTrack';
Recording = '2017-07-03_13-08-30';
tioChOrd = [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54]; chOrdTxt = 'Probe order';

% update this to match the name of the variable with your channel order
chOrd = tioChOrd;
dsFreq = 600; %what is this?

%% 
%ratLibPath = 'D:\Dropbox (NewmanLab)\ratsEphys';
workingDir = fullfile(ratLibPath,Rat,Session,Recording); cd(workingDir);
%disp(path)
nChan = length(chOrd);

% ??? Why is there a tmp root? k 
tmpRoot = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:0.01:21.0028*60],'fs_video',120);
tmpRoot.path_lfp = repmat({'experiment1_100.raw.kwd'},1,1);
tmpRoot = tmpRoot.LoadLFP(1,'downsample',dsFreq,'chOrd',chOrd(1));
tmpRoot.active_lfp = nChan;

fsVid = 120;
root = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:1/fsVid:max(tmpRoot.b_lfp(1).ts)],'fs_video',fsVid);
root.path_lfp = repmat({'experiment1_100.raw.kwd'},1,nChan);
%save experiment1.mat root
root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
root.active_lfp = nChan;

root = rmDataBlips(root);

%iterate through the channels
chans = [1:nChan];
dataDS = nan(length(chans),length(CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal)));
for i = 1:length(chans)
  root.active_lfp = chans(i);
  dataDS(i,:) = CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal);
end

%% Extract Theta

Fs = root.lfp.fs;
Wn_theta = [6/(Fs/2) 10/(Fs/2)];
[btheta,atheta] = butter(3,Wn_theta);


% extract theta with phase and power
theta_filt = nan(size(dataDS));
theta_phase =  nan(size(dataDS));
theta_amp =  nan(size(dataDS));
% why is this iterating through each data point, and not taking all of the data? k
for iD =  1:size(dataDS,1)
  theta_filt(iD,:) = filtfilt(btheta,atheta,dataDS(iD,:)); %filter the data
  theta_phase(iD,:) = atan2(imag(hilbert(theta_filt(iD,:))), theta_filt(iD,:));
  theta_amp(iD,:) = abs(hilbert(theta_filt(iD,:)));
end

%% 
% Focus on epochs of data with high theta amplitude
% high amplitude will be based on greater than 2 std above the mean
% so, compute session-wide mean and standard deviation of theta power
% using the last channel (update this if another channel makes more sense)

% I might need to change this k 
% Also, why not all theta? quick n dirty or the whole project?
ch = size(theta_amp,1);
meanAmp = mean(theta_amp(ch,:));
stdAmp = std(theta_amp(ch,:));
 
% now find the high theta 
%highTheta = find(theta_amp(ch,:)>(meanAmp+2*stdAmp));
%highTheta = find(theta_amp(end,:)>(meanAmp));
highTheta = find(theta_amp(end,:)>(meanAmp-stdAmp));
highThetaEp = mat2cell(highTheta, 1, diff([0 find([(diff(highTheta) > 1) 1])]));
lengthEp = cellfun(@length,highThetaEp);
inds_long = lengthEp>300;
highThetaEp_long = highThetaEp(inds_long);
 
clear thetaShift
for iT = 1:length(highThetaEp_long)
  for iC1 = 1:length(chans)
    for iC2 = 1:length(chans)
      thetaShift{iT}(iC1,iC2,:) = circDiff([theta_phase(iC1,highThetaEp_long{iT})', theta_phase(iC2,highThetaEp_long{iT})'],2,'rad'); % changed orientation of data
    end
  end
end

thetaShiftMat = cat(3,thetaShift{:}); % used cat instead of cell2mat
fprintf('Basing calculations off of %2.2f s of data\n',size(thetaShiftMat,3)/Fs);
[thetaShiftAngle,thetaShiftRbar] = circmean(thetaShiftMat,3);

%%
%Figures
[u,v] = pol2cart(thetaShiftAngle,thetaShiftRbar);
figure; quiver(u(2:2:end,2:2:end),v(2:2:end,2:2:end)); axis ij;  title([chOrdTxt ', Blank (1403)']);
circR = sqrt(u.^2 + v.^2); figure; imagesc(circR,[0 1]);

R = corr(theta_filt');
figure; imagesc(R);
FirstDegMCorr = mean(diag(R,1));
SecondDegMCorr = mean(diag(R,2)); 
title([chOrdTxt, ', Blank (1403), 1st vs 2nd neighbor corr ratio = ',num2str(FirstDegMCorr/SecondDegMCorr)]);
 
%thetaShiftMean = mean(thetaShiftMat,3);% adjusted dimension
%thetaShiftStd = std(thetaShiftMat,0,2);
 
% bins = [-pi:0.05:pi];
% figure; 
% for i = 1:size(thetaShiftMat,1), 
%   subplot(size(thetaShiftMat,1),1,i), hist(thetaShiftMat(i,:),bins); 
%   if i==1, title('8/17/2016 14:03 Recordings. Theta > 2 std. above mean power.'); end
%   ylabel([num2str(chans(i+1)) ' - ' num2str(chans(i))]);
% end
 
%highTheta = find(theta_amp(end,:)>(meanAmp+2*stdAmp));
highTheta = find(theta_amp(end,:)>(meanAmp));
highThetaEp = mat2cell(highTheta, 1, diff([0 find([(diff(highTheta) > 1) 1])]));
lengthEp = cellfun(@length,highThetaEp);
inds_long = lengthEp>300;
highThetaEp_long = highThetaEp(inds_long);
 
clear thetaShift
for iT = 1:length(highThetaEp_long),
  thetaShift{iT} = circDiff(theta_phase(:,highThetaEp_long{iT}),1,'rad');
end
thetaShiftMat = cell2mat(thetaShift);
 
bins = [-pi:0.05:pi];
figure; 
for i = 1:size(thetaShiftMat,1), 
  subplot(size(thetaShiftMat,1),1,i), hist(thetaShiftMat(i,:),bins); 
  if i==1, title('8/17/2016 14:03 Recordings. Theta > mean power.'); end
  ylabel([num2str(chans(i+1)) ' - ' num2str(chans(i))]);
end


%%
%new stuff

