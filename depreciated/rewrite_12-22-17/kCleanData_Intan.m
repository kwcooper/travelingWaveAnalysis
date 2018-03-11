function [pctDataUsed,inds2cut,varargout] = cleanData_Intan(signal,fs,filtType,filtParams,varargin) %gammaamps,thetaangles)

% adapted from cleanData2.m on 171207 to handle Intan instead of Axona data

varargout = []; if ~exist('varargin','var'), varargin = []; end;
fs = double(fs);
buf = ceil(fs/2); % 500ms buffer around each artifact

% make sure signal is vector, whine if not
if ~isvector(signal), error('Signal must be a vector for a single channel'); end

%% compute specific filters if needed
if exist('filtType','var') && ~isempty(filtType)
  
  % NOTE: if you want scaled EEG, use scaledEEG.m to update the root object before submitting it to cleanData2;
  
  % find good amps
  switch(filtType)
    case 'theta' % relative amplitude of theta in percentage of it max value
      th_Amp = abs(hilbert(buttfilt(signal,[6 12],fs,'bandpass',4)));
      incInds = ~(th_Amp < max(th_Amp)*filtParams(1) | th_Amp > max(th_Amp)*filtParams(2));
      incEpchs = CMBHOME.Utils.OverThresholdDetect(incInds,0.5,ceil(fs/10),ceil(fs/2));
      excInds = true(size(incInds));
      for ep = 1:size(incEpchs,1)
        excInds(max(incEpchs(ep,1)-buf,1):min(incEpchs(ep,2)+buf,length(excInds))) = false;
      end

    case 'thetaMag' % absolute magnitude of theta in mV
      th_Amp = abs(hilbert(buttfilt(signal,[6 12],fs,'bandpass',4)));
      incInds = ~(th_Amp < filtParams(1) | th_Amp > filtParams(2));
      incEpchs = CMBHOME.Utils.OverThresholdDetect(incInds,0.5,ceil(fs/10),ceil(fs/2));
      excInds = true(size(incInds));
      for ep = 1:size(incEpchs,1)
        excInds(max(incEpchs(ep,1)-buf,1):min(incEpchs(ep,2)+buf,length(excInds))) = false;
      end
      
    case 'thetaDeltaRatio'
      th_Amp = abs(hilbert(buttfilt(signal,[6 12],fs,'bandpass',4)));
      dl_Amp = abs(hilbert(buttfilt(signal,[2 4],fs,'bandpass',4)));
      rat = th_Amp ./ dl_Amp;
      excInds = true(size(rat));
      incEpchs = CMBHOME.Utils.OverThresholdDetect(rat,filtParams(1),ceil(fs/4),ceil(fs/2));
      for ep = 1:size(incEpchs,1)
        excInds(max(incEpchs(ep,1)-buf,1):min(incEpchs(ep,2)+buf,length(excInds))) = false;
      end
      
    otherwise
      error('Unknown filtType');
  end
  excInds = reshape(excInds,1,[]);
else
  excInds = zeros(1,length(signal));
end

% record original data length
ODL = length(signal);

%% examine the signal quality itself to find noise artifacts
badInds = false(length(signal),1);

% drop saturated or flatlined data
t0 = [0-1e-4 0+1e-4]; % tresholds for no dV/dt
badInds_ = CMBHOME.Utils.ThresholdBandDetect(diff(signal),t0(1),t0(2),1,round(0.030 * fs)); % flatline data
for ep = 1:size(badInds_,1)
  badInds(max(badInds_(ep,1)-buf,0):min(badInds_(ep,2)+buf,length(badInds))) = true;
end

% drop high amplitude events
% high amplitude = 2.5 x RMS
lastBadInds = 0;
newBadIndsFound = true;
thresh = 5;
while newBadIndsFound
  sigRMS = norm(signal(~badInds))/sqrt(length(signal(~badInds))); % adapted from Hazem Saliba Baqaen. Hazem@brown.edu
  badInds_ = [];
  badInds_ = [badInds_; CMBHOME.Utils.OverThresholdDetect(signal,thresh*sigRMS,1,1)];% upward deflections
  badInds_ = [badInds_; CMBHOME.Utils.UnderThresholdDetect(signal,-thresh*sigRMS,1,1)];% downward deflections
  for ep = 1:size(badInds_,1)
    currInds = max(badInds_(ep,1)-buf,1):min(badInds_(ep,2)+buf,length(badInds));
    badInds(currInds) = true;
  end
  if lastBadInds == size(badInds_,1)
    newBadIndsFound = false;
  else
    lastBadInds = size(badInds_,1);
  end
end


sigRMS = norm(signal(~badInds))/sqrt(length(signal(~badInds))); % adapted from Hazem Saliba Baqaen. Hazem@brown.edu
badInds_ = [];
badInds_ = [badInds_; CMBHOME.Utils.OverThresholdDetect(signal,2.5*sigRMS,1,round(0.030 * fs))];% upward deflections
badInds_ = [badInds_; CMBHOME.Utils.UnderThresholdDetect(signal,-2.5*sigRMS,1,round(0.030 * fs))];% downward deflections
for ep = 1:size(badInds_,1)
  badInds(max(badInds_(ep,1)-buf,1):min(badInds_(ep,2)+buf,length(badInds))) = true;
end

%% review what is being cut 
if 1
  figure;
  plot(signal); hold on;
  plot([1 length(signal)],repmat(thresh*sigRMS,[2 1]),'-k');
  plot([1 length(signal)],repmat(-thresh*sigRMS,[2 1]),'-k');
  plot(find(badInds),signal(badInds),'oc')  % rejected by artifact detector
  plot(find(excInds),signal(excInds),'.r'); % doesn't meet filter
  keyboard
end


%% wrap up

% combine two exclusion vectors
inds2cut = badInds(:) | excInds(:);

% perform filter
for i = 1:length(varargin)
  if isvector(varargin{i})
    varargin{i}(inds2cut) = [];
  else
    varargin{i}(:,inds2cut) = [];
  end
end
varargout = varargin;

% document how much data remains
pctDataUsed = sum(~inds2cut)/ODL;
end


