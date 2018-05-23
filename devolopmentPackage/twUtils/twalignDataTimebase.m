function twalignDataTimebase(ratName,prefixes,Recording,ephysFolder,bonsaiTime,ref,filetype)
% alignDataTimebase(ratName,prefixes,tet2save)
%
%   Aligns timestamps collected from Bonsai and Open Ephys to align the tracking data to the EEG.
%
%   INPUTS: 
%     ratName(string)  - folder of rat
%     prefixes(cell)   - session folders to edit
%     ephysFolder(cell)
%     bonsaiTime(cell)
%     ref(int)
%     filetype(string)
%
%   OUTPUTS:
%     alignedPositionData.csv - N X 5 file where N is the number of timestamps
%                               Columns contain, in order: [timestamp big-x big-y small-x small-y]
%     alignedEphysData.mat    - mat file with good inds and timestamps of ephys data
%                              

if ~exist('bonsaiTime','var') || isempty(bonsaiTime), bonsaiTime = cell(size(prefixes)); end
if ~exist('ephysFolder','var') || isempty(ephysFolder), ephysFolder = cell(size(prefixes)); end
if ~exist('filetype','var') || isempty(filetype), filetype = 'openephys'; end

for iP = 1:length(prefixes)
  basePath = fullfile(dropboxPath,'ratsEphys',ratName,prefixes{iP});
  
  %%%%BSC 170215 - if the final file resulting from this function already
  %exists, then the function has already been run on this session, so go
  %to the next session in prefixes
%   if exist([basePath filesep 'alignedEphysData.mat'],'file')
%     disp(['skipping ' prefixes{iP} ' - this session has already been aligned'])
%     continue
%   end
  %%%%
  
  % Prompt if multiple bonsai files are found
  if isempty(bonsaiTime{iP}), bonsaiTime{iP} = '*'; end
  bonsaiFile = dir([basePath filesep 'metadata' bonsaiTime{iP} '.csv']);
  ind = 1;
  if length(bonsaiFile)>1
    for iF = 1:length(bonsaiFile)
      disp([num2str(iF) ': ' bonsaiFile(iF).name]);
    end
    ind = input('Enter index of Bonsai file: ');   
  end
  tmp = extractBetween(bonsaiFile(ind).name,'metadata','.csv'); % for some reason this outputs a cell
  bonsaiTime{iP} = tmp{1};
  bonsaiFile = [basePath filesep bonsaiFile(ind).name];
  
  posFile = dir([basePath filesep 'tracking' bonsaiTime{iP} '.csv']); % ELN & BSC 170125
  posFile = [basePath filesep posFile.name];
  
  % SJV 170706 added option to read in open ephys for ECube. 
  % changed input parameters to alignVideo2LFP
  % Prompt if multiple ephys folders
  if isempty(ephysFolder{iP})
    folderList = CollectFolders(basePath);
    ind = 1;
    if length(folderList)>1
      % kwc add flag for prompting vs automation
      if 0
        for iF = 1:length(folderList)
          disp([num2str(iF) ': ' folderList{iF}]);
        end
        ind = input('Enter index of Open Ephys folder: ');
      else
        for iF = 1:length(folderList)
          if folderList{iF}(length(folderList{iF})-length(Recording)+1:end) == Recording
            ind = iF;
          end
        end
      end
      
    end
    ephysFolder{iP} = folderList{ind}(length(basePath)+2:end);
  end
    
  switch filetype
    case 'kwik'
      ephysFile = dir([basePath filesep ephysFolder{iP} filesep '*.kwd']);
      ephysFile = [ephysFile.folder filesep ephysFile.name];
      ephysRef = h5read(ephysFile,'/recordings/0/data',[ref 1],[1 inf]);
      fs = double(h5readatt(ephysFile,'/recordings/0/','sample_rate'));
      ephysTS = [1:length(ephysRef)]./fs; % KWC 180510
    case 'openephys'
      reference = dir([basePath filesep ephysFolder{iP} filesep '*_ADC' num2str(ref) '.continuous']); 
      if isempty(reference), reference = dir([basePath filesep ephysFolder{iP} filesep '*_PDI' num2str(ref) '.continuous']); end
      [ephysRef,ephysTS,info] = load_open_ephys_data([reference.folder filesep reference.name]); 
      fs = info.header.sampleRate;
    otherwise
      error('invalid file type');
  end
   
  bonsaiData = dlmread(bonsaiFile,' '); 
  bonsaiData(:,end) = [];
  
  fprintf('Computing time lag between recording and tracking... ');
  % 170824 SJV calculating a rough estimate of the time lag between the two
  % recordings using the file time stamps. This way the lags in the cross
  % correlation will be sure to include the time lag while also not having
  % to have some large defalt interval.
  e_ts = extractAfter(ephysFolder{iP},'_'); e_ts = split(e_ts,'-'); e_ts = cellfun(@str2double,e_ts);
  b_ts = extractAfter(bonsaiTime{iP},'T'); b_ts = split(b_ts,'_'); b_ts = cellfun(@str2double,b_ts);
  lags = e_ts - b_ts; lags = 3600*lags(1) + 60*lags(2) + lags(3);
  lags = [-(abs(lags)+2) (abs(lags)+2)]; % +2 second buffer
  
  [offset,ts_b,ts_e,badinds] = alignVideo2LFP(bonsaiData,ephysRef,ephysTS,fs,lags);
  
  cutinds_b = []; cutinds2_b = [];
  cutinds_e = []; cutinds2_e = [];
  
  % cut off beginnsave('rooing of the data from the program that started first
  if offset > 0 % bonsai started first
    ts_b = ts_b - offset;
  elseif offset < 0 % open ephys started first
    ts_e = ts_e + offset;
  end
  
  %%%%BSC 170214 - added to fix off-by-1 diffs; doesn't handle larger diffs
  posData = readSpottyTracking(posFile); % ELN & BSC 170125
  %posData = dlmread(posFile,',');
  
  if length(ts_b) == length(posData)+1, ts_b = ts_b(1:end-1); end
  if length(ts_b) == length(posData)-1, posData = posData(1:end-1); end
  %%%%
  
  % cut off the end of the data from the program that ended last
  if ts_b(end) > ts_e(end) % bonsai ended last
    cutinds_b = ts_b<0 | ts_b>ts_e(end);
    cutinds_e = ts_e<0;
  elseif ts_e(end) > ts_b(end) % open ephys ended last
    cutinds_b = ts_b<0;
    cutinds_e = ts_e<0 | ts_e>ts_b(end);
  end
  
  ts_b(cutinds_b) = [];
  ts_e(cutinds_e) = [];
    
  if ~isempty(badinds), posData(badinds,:) = []; end
  if any(cutinds_b), posData(cutinds_b,:) = []; end 
  
  posData = [ts_b posData]; % append timestamps
  dlmwrite([basePath filesep 'alignedPositionData.csv'],posData,'delimiter',',','precision',16); % save new position file
  
  goodInds = ~cutinds_e;
  timestamps = ts_e;
  
  % sjv 170824
  % the first timestamp given by open ephys is often not at 0. this is
  % already accounted for in the timestamps saved below; however when 
  % loading in spike waveforms and their corresponding timestamps 
  % (in dacq_tetrode), the "start time" offset is still present. Since we 
  % can't rely on the fact that there is a spike waveform at the first time 
  % point recorded by open ephys I'm saving the start time here so we can 
  % account for this offset later when loading inthe spike waveforms.
  starttime = ephysTS(1); 
  fprintf('Done!\nSaving work to data directory.\n');
  save([basePath filesep 'alignedEphysData.mat'],'goodInds','timestamps','offset','starttime');  
end

end 

function [offset,ts_sb,ts_se,badinds] = alignVideo2LFP(metadata,se,ts_se,fs,lags)

% extract Bonsai metadata
sb = metadata(:,1);
frame_sb = metadata(:,2);
t_sb = metadata(:,3);

% get bad frames
badinds = find(frame_sb == 0);
sb(badinds) = [];
frame_sb(badinds) = [];
t_sb(badinds) = [];

% convert timestamps to seconds
if all(t_sb==0)
  sec_sb(1) = 0; cycle_sb(1) = 0;
  frame_inc = diff(frame_sb);
  for i = 2:length(frame_sb)
    cycle_sb(i) = cycle_sb(i-1)+frame_inc(i-1)*8000/120-0.04;
    if cycle_sb(i)>= 8000
      sec_sb(i) = sec_sb(i-1)+1;
      cycle_sb(i) = cycle_sb(i)-8000;
    else
      sec_sb(i) = sec_sb(i-1);
    end
  end
else
  t_sb = dec2bin(t_sb,32);
  sec_sb = bin2dec(t_sb(:,1:7));
  cycle_sb = bin2dec(t_sb(:,8:20));
  offset_sb = bin2dec(t_sb(:,21:28)); % NOTE: last 4 bits (29-32) are inaccurate. Probably won't use this value in general.
end
  
% find restarts in time stamps and add to make continuous
wrapinds = find(diff(sec_sb)<0);
wrapinds = [0; wrapinds; length(sec_sb)];
ts_sb = zeros(length(sec_sb),1);
for iW = 1:length(wrapinds)-1
  for iN = wrapinds(iW)+1:wrapinds(iW+1)
    ts_sb(iN) = 128*(iW-1)+sec_sb(iN)+cycle_sb(iN)/8000;
  end
end
ts_sb = ts_sb-ts_sb(1);

sb_ = sb;
thresh = mean([min(sb) max(sb)]);
sb_(sb>thresh)=1;
sb_(sb<=thresh)=0;

% convert Open Ephys trigger signal to a logical with 0s and 1s
% Using the aligator clips to send the signal to the camera distorts the
% signal sent to open ephys, causeing the signal to not lie between 0 and
% one, so I'm using the mean value between the min and max to threshold
% the signal.
% I was also seeing a DC shift on the ADC signal from open ephys before the
% commutator was installed, the pulse seemed to be affected by movement
% artifacts of the rat. If this DC shift shows up in the trigger signal in
% open ephys, the following method will not threshold it correctly...
se_ = double(se);
thresh = mean([min(se) max(se)]);
se_(se>thresh) = 1;
se_(se<=thresh) = 0;

% create timestamps for Open Ephys
if ~exist('ts_se','var') || isempty(ts_se)
  ts_se = (0:length(se_)-1)'/fs;
else
  ts_se = ts_se - ts_se(1);
end

% get rising edge time
spk_sb = ts_sb(diff(sb_)>0);
spk_se =  ts_se(diff(se_)>0);

% get cross correlation of "spike" times
[r,lag] = CMBHOME.Spike.CrossCorr(spk_se,'ts2',spk_sb,'binsize',0.02,'lag',1.2*lags);
[maxVal,m] = max(r);
offset = lag(m);

if 0
  % plots each of the spikes along with the xcorr
  figure;
  subplot(3,1,1); plot(lag,r); hold on; plot(offset,maxVal,'*');
  subplot(3,1,2); i = 2;
  plot(ts_sb(ts_sb>=(i-1)+offset&ts_sb<(i-1)+10+offset),sb_(ts_sb>=(i-1)+offset&ts_sb<(i-1)+10+offset));
  axis([i-1+offset (i-1)+10+offset 0 1.5])
  subplot(3,1,3);
  plot(ts_se(ts_se>=(i-1)&ts_se<(i-1)+10),se_(ts_se>=(i-1)&ts_se<(i-1)+10));
  axis([i-1 (i-1)+10 0 1.5]);
end

if 0
  figure;
  for i=500:1000
    subplot(2,1,1);
    plot(ts_sb(ts_sb>=(i-1)&ts_sb<(i-1)+10),sb_(ts_sb>=(i-1)&ts_sb<(i-1)+10));
    axis([i-1 (i-1)+10 0 1.5])
    subplot(2,1,2);
    plot(ts_se(ts_se>=(i-1)&ts_se<(i-1)+10),se_(ts_se>=(i-1)&ts_se<(i-1)+10));
    axis([i-1 (i-1)+10 0 1.5]);
    pause;
  end
end

end
  
