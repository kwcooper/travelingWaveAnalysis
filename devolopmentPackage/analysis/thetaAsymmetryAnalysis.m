function [asymHist,asymWaves,peaksAndTrphs,asymScores] = thetaAsymmetryAnalysis(root)
% function [asymHist,asymWaves,peaksAndTrphs,asymScores] = thetaAsymmetryAnalysis(root)
%
% Asymetry is = log(risingDuration / fallingDuration)

%% Prep needed values
signal = root.b_lfp(root.active_lfp).signal;
fs = root.b_lfp(root.active_lfp).fs;

min_sep = fs / 20; % 50ms separation between peak and next trough
min_length = min_sep; % same for the separation between trough and next peak
max_dur = fs / 4;
min_dur = fs / 12;

%  method described by (Belluscio et al., 2012)
broadBand = [1 80];
signal = buttfilt(signal,broadBand,fs,'bandpass',3);

%% find filtered theta peaks and troughs
if exist('band','var') & ~isempty(band),	thetarange = minmax(band); 
else, thetarange = [4 12]; end
thetaPhs = angle(hilbert(buttfilt(signal,thetarange,fs,'bandpass',3))) + pi;

% find bounds of peaks as pi/2 < peaks < 3pi/2
peakBounds = CMBHOME.Utils.ThresholdBandDetect(thetaPhs, pi/2, 3*pi/2, min_sep, min_length);

% find bounds of troughs as  pi/2 < peaks < 3pi/2 after adding pi
trphBounds = CMBHOME.Utils.ThresholdBandDetect(mod(thetaPhs+pi,2*pi), pi/2, 3*pi/2, min_sep, min_length);

% find mean peak duration
meanPeakDur = ceil(mean(diff(peakBounds,[],2)));
peak_iM = ones(size(peakBounds,1),meanPeakDur); % prepare to pull 150 samples of data around each event
peak_iM(:,1) = peakBounds(:,1);
peak_iM = cumsum(peak_iM,2); % build indice matrix
peakBounds(max(peak_iM,[],2)>length(signal),:) = [];% drop cycle which extend past the end of signal
peak_iM(max(peak_iM,[],2)>length(signal),:) = []; % drop cycle which extend past the end of signal
peakSig = signal(peak_iM);
peakInd = peakBounds(:,1) + maxInd(peakSig,[],2) - 1;

% find mean trough duration
meanTrphDur = ceil(mean(diff(trphBounds,[],2)));
trph_iM = ones(size(trphBounds,1),meanTrphDur); % prepare to pull 150 samples of data around each event
trph_iM(:,1) = trphBounds(:,1);
trph_iM = cumsum(trph_iM,2); % build indice matrix
trphBounds(max(trph_iM,[],2)>length(signal),:) = [];% drop cycle which extend past the end of signal
trph_iM(max(trph_iM,[],2)>length(signal),:) = []; % drop cycle which extend past the end of signal
trphSig = signal(trph_iM);
trphInd = trphBounds(:,1) + minInd(trphSig,[],2) - 1;

%% assemble ordered list of peaks and troughs
peaksAndTrphs = [peakInd ones(size(peakInd)); trphInd -ones(size(trphInd))];
[~,srtInds] = sort(peaksAndTrphs(:,1));
peaksAndTrphs = peaksAndTrphs(srtInds,:);

% remove double-troughs and double-peaks
badInds = diff(peaksAndTrphs(:,2))==0;
peaksAndTrphs(badInds,:) = [];

% remove double-counted inds
badInds = find(diff(peaksAndTrphs(:,1))==0);
peaksAndTrphs(badInds,2) = nan; peaksAndTrphs(badInds+1,2) = nan;
peaksAndTrphs(isnan(peaksAndTrphs(:,2)),:) = [];

% start and end on a trough to get equal numbers of rise and fall times
firstTrph = find(peaksAndTrphs(:,2)==-1,1,'first');
lastTrph = find(peaksAndTrphs(:,2)==-1,1,'last');
peaksAndTrphs = peaksAndTrphs(firstTrph:lastTrph,:);

%% analyze rise and fall phases
% get rise and fall times
cyclDurs = diff(peaksAndTrphs);
riseDurs = cyclDurs(cyclDurs(:,2)==2,1);
fallDurs = cyclDurs(cyclDurs(:,2)==-2,1);
notCycles = ((riseDurs+fallDurs)>max_dur | (riseDurs+fallDurs)<min_dur);
riseDurs(notCycles) = []; fallDurs(notCycles) = [];

% put rise to fall ratio for corresponding peak in lfp.myvar in format [riseDur fallDur]
tmp = nan(size(signal,1),1);
peakInds = peaksAndTrphs(peaksAndTrphs(:,2)==1,1);
peakInds(notCycles) = [];
tmp(peakInds,:) = riseDurs./fallDurs;
root.b_lfp(root.active_lfp).b_myvar = tmp;

% get epoched scores
asymScores = root.lfp.myvar;
if iscell(asymScores), asymScores = cell2mat(asymScores); end
asymScores(isnan(asymScores)) = [];
asymScores = log(asymScores); % logarithmic scale

%% compile results 
% compute histogram of asymetries
[asymHist.counts,asymHist.edges,inds]= histcounts(asymScores,100);

% extract average waveforms
meanPeakDur = 63; buf = floor(meanPeakDur/2);
peak_iM = ones(sum(peakInds>buf),meanPeakDur); % prepare to data around each event
peak_iM(:,1) = peakInds(peakInds>buf)-buf;
peak_iM = cumsum(peak_iM,2); % build indice matrix
peakSig = signal(min(peak_iM,length(signal)));

% pull waveforms for each cycleType
allInds = unique(inds);
asymWaves = nan(meanPeakDur,length(allInds));
for i = 1:length(allInds)-1
  if sum(inds==i)/length(inds)>0.005
    asymWaves(:,i) = mean(peakSig(inds==i+1,:))';
  end
end
