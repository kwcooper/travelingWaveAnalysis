 function root = twPrepData(D,fs,metaData)
%FUNCTION TWPREPDATA %cleans and extracts useful properties of the data 
% as well as creates the root object with tracking data
%Part of the Traveling Wave Analysis package

%keyboard; 
% TD: 
%     need to look at timestamps

tracking = 0;
%% grab tracking data

% original script to align tracking, make root, and build fnamesets
%importDataWTrackingAuto(metaData)
if tracking
  % ... this should be changed? (I automated it)
  twalignDataTimebase(metaData.Rat, {metaData.Session}, metaData.Recording,[],[],metaData.trackRef,metaData.fType); % make position file
  
  %code adapted from invivoimport ~line 149 on 180516
  pos_file = fullfile(ratLibPath,metaData.Rat,metaData.Session,'alignedPositionData.csv');
  
  % dacq_pos_3p for 3 LEDs, dacq_pos for 2 LEDs
  % both files seem to return similar output data on Tio session 1
  [scaleFactor, vid_ts, vid_x, vid_y, vid_headdir, d_epoch] = dacq_pos_3p(pos_file);
  scaleFactor = 115/800; % S- CHANGE FOR PROJECT, make variable % k- these are also the same
  
  % make lfp timestamps
  b_ts = linspace(0,size(D,2)/fs,size(D,2)+1); b_ts = b_ts(2:end)'; % cut off first to get rid of 0
  
  % !! need to work on video timestamps...
  % calculate new stamps with mean of current, interpolate with the output stamps
  vid_ts_floor = floor(1/mean(diff(vid_ts)));
  
  
  root = CMBHOME.Session('name', metaData.Rat,...
    'epoch', [-inf inf],...
    'b_ts', vid_ts,... % as per invivoimport line 154
    'fs_video', 1/mean(diff(vid_ts)), ...
    'b_x', vid_x, ...
    'b_y', vid_y, ...
    'raw_pos', 1, ...
    'raw_headdir', 1, ...
    'b_headdir', vid_headdir, ...
    'date_created', now, ...
    'epoch', [-inf inf], ...
    'fs_video', 1/mean(diff(vid_ts)), ...
    'spatial_scale',scaleFactor);
  
  % import lfp time offset: invivoimport line 170
  ephys_align = fullfile(ratLibPath,metaData.Rat,metaData.Session,'alignedEphysData.mat');
  load(ephys_align,'timestamps','offset','starttime');
  
else
  b_ts = linspace(0,size(D,2)/fs,size(D,2)+1); b_ts = b_ts(2:end)';
  root = CMBHOME.Session('name', metaData.Rat, 'epoch', [-inf inf], 'b_ts', b_ts);
  
end
%%
% store useful data in root
root.user_def.lfp_origData = D;
root.user_def.lfp_fs = fs;
root.user_def.metaData = metaData;


for i = 1:size(D,1)
  root.b_lfp(i).signal = root.user_def.lfp_origData(i,:)'; 
  root.b_lfp(i).fs = root.user_def.lfp_fs; 
  root.b_lfp(i).ts = b_ts; 
end

% temporaly change D to the aligned LFP for correct sizing
%D = ;

% clean data and filter to epochs of good theta
%filtType = 'thetaDeltaRatio'; filtParams = 2; % works okay, some fine tuning may help further rule out low theta epochs and include 
filtType = 'thetaMag'; %filtParams = [115 1000]; % Update for each rat  
[pctDataUsed,inds2cut] = cleanData_Intan(D(metaData.ref,:),fs,filtType,metaData.filtParams); 
root.user_def.cleanData_inds2cut = inds2cut;
root.user_def.cleanData_pctDataUsed = pctDataUsed;


% Iterate through channels to extract needed variables
% Average waves needs to know cycle bounds & fs
thetaPhs = nan(size(D));
thetaAmp = nan(size(D));
cycles = nan(size(D));



for i = 1:size(D,1)
  [thetaPhs(i,:),thetaAmp(i,:),~] = extractThetaPhase(D(i,:),fs,'hilbert',[6 10]);
  [cycles(i,:),~] = parseThetaCycles(thetaPhs(i,:),fs,[6 10],0);         % grab peak data
  [cyclesTrough(i,:),~] = parseThetaCycles(thetaPhs(i,:),fs,[6 10],pi);  % grab trough data
end 

root.user_def.theta_phs = thetaPhs;
root.user_def.theta_amp = thetaAmp;
root.user_def.cycles = cycles;
%root.user_def.theta_phsTrough = thetaPhsTrough; % these are the same for peaks
%root.user_def.theta_ampTrough = thetaAmpTrough; % and for troughs, to be removed...
root.user_def.cyclesTrough = cyclesTrough;
root.user_def.thetaLFP =  buttfilt(D',[6 10],fs,'bandpass',3)'; % may need to D.T

end
  
  
  









