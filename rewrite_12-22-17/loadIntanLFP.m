function [data,fs] = loadIntanLFP(chOrd,dataPath,fnames,dsRate,recordingNum)

  % get situated in the data folder
  curDir = pwd;
  if ~exist(dataPath,'dir'), error('Cannot find %s',dataPath); end
  cd(dataPath);
  
  % set fnames if not provided, quality check if provided
  if ~exist('fnames','var') || isempty(fnames)
    useKwd = true;
    fnames = dir('*.kwd');
    if length(fnames)>1, error('Found multiple kwd files, please specify which you want.'); end
    
    if isempty(fnames)
%       fnames = dir('*CH*.continuous');
%       fnames = cleanFnames(fnames);
      useKwd = false;
    end
  else
    if ~isstruct(fnames), error('fnames must be a struct, like the one made by dir().'); end
    for f = 1:size(fnames,1)
      if ~exist(fnames(f).name,'file') error('Did not find file %s\n',fnames(f).name); end
    end
    useKwd = strcmp(fnames(1).name(end-2:end),'kwd');
  end
  
  % initialize other variables as needed
  if ~exist('dsRate','var') || isempty(dsRate), dsRate = 500; end
  if ~exist('recordingNum','var') || isempty(recordingNum)
    if useKwd, recordingNum = 0;
    else, recordingNum = 100; end
  end  
  % load the data
  if useKwd
    ephysFile = fnames.name;
    data = h5read(ephysFile,['/recordings/',num2str(recordingNum),'/data']);
    fs = double(h5readatt(ephysFile,['/recordings/',num2str(recordingNum),'/'],'sample_rate'));
    data = data(chOrd,:);
    % downsample as needed
    dsStep = ceil(fs / dsRate);
    data = data(:,1:dsStep:end); fs = fs / dsStep;
    data = double(data);
  else
    ch = 1;
    fprintf('Loading(%i): ',length(chOrd));
    ephysFile = [num2str(recordingNum),'_CH',num2str(chOrd(ch)),'.continuous']; fprintf('%i ',ch);
    [tmpData, ~, info] = loadAndCorrectPhase(ephysFile, 1);
    %[data,ts] = load_open_ephys_data_faster(ephysFile); % using above line as it is presumably important to correct for the described phase errors
    fs = info.header.sampleRate;
    dsStep = ceil(fs / dsRate);
    tmpData = tmpData(1:dsStep:end); fs = fs/dsStep;
    % initialize data matrices
    data = nan(length(chOrd),length(tmpData));
    data(ch,:) = tmpData';
    % load the rest of the data
    for ch = 2:length(chOrd)
      ephysFile = [num2str(recordingNum),'_CH',num2str(chOrd(ch)),'.continuous'];fprintf('%i ',ch);
      tmpData = loadAndCorrectPhase(ephysFile, 1);
      %tmpData = load_open_ephys_data_faster(ephysFile); % using above line as it is presumably important to correct for the described phase errors
      tmpData = tmpData(1:dsStep:end);
      data(ch,:) = tmpData';
    end
    fprintf('\n');
  end
  
  cd(curDir);
    
end

function fnames = cleanFnames(fnames)
  
  % process filenames to extract channel and root numbers
  roots = cell(size(fnames));
  fnamesNums = nan(size(fnames));
  for i = 1:size(fnames,1)
    roots{i} = strtok(fnames(i).name,'_');
    fnamesNums(i) = str2num(fnames(i).name(strfind(fnames(i).name,'CH')+2:strfind(fnames(i).name,'.continuous')-1));
  end
  
  % confirm only one session in recording
  uniqRoots = unique(roots);
  if numel(uniqRoots)>1, error('Found multiple sessions, please specify which you want.'); end
  
  % sort channel numbers so that fnames is in numerical order instead of
  % alphabetical order
  [~,inds] = sort(fnamesNums);
  fnames = fnames(inds);
end