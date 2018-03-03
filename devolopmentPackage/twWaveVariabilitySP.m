function [] = twWaveVariabilitySP(root, figData, scan)
%TWWAVEVARIBILITYSP stripped and designed to analyze the varibility 
% in the wave offset for slopes and peaks
%
% scan (1/0)  used to paroose
%
% Not the cleanest implementation, but allows you to use multiple variables
% look at twwavevaribility for a readable version

fprintf('\nExploring varibility in theta wave offset (slopes & peaks)...\n');
% set defaults and allocate mats
name = figData.ratInfo.name; 
if simu
  Rat = 'Test';
end

origData = root.user_def.lfp_origData; 
 % pi=toughs; 0=peaks;

thetaPhsPeaks = nan(size(origData));
cyclesTroughs = nan(size(origData));
thetaPhsPeaks = nan(size(origData));
cyclesTroughs = nan(size(origData));

phs = 0; fprintf(['Grabbing theta cycles (peaks)\n']);
for i = 1:size(origData,1)
  thetaPhsPeaks(i,:) = extractThetaPhase(origData(i,:),root.user_def.lfp_fs,'waveform',[6 10]);
  [cyclesPeaks(i,:),~] = parseThetaCycles(thetaPhsPeaks(i,:),root.user_def.lfp_fs,[6 10], phs); 
end 

phs = pi; fprintf(['Grabbing theta cycles (troughs)\n']);
for i = 1:size(origData,1)
  thetaPhsTroughs(i,:) = extractThetaPhase(origData(i,:),root.user_def.lfp_fs,'waveform',[6 10]);
  [cyclesTroughs(i,:),~] = parseThetaCycles(thetaPhsTroughs(i,:),root.user_def.lfp_fs,[6 10], phs); 
end 

% Let's clean up cycles a bit to remove the bad data
cyclePeakInds = find(cyclesPeaks(1,:)); cyclesSize = size(cyclesPeaks,2);
badCycles = root.user_def.cleanData_inds2cut(cyclePeakInds);
cyclePeakInds(badCycles) = []; 
fprintf(['Cleaned ', num2str(cyclesSize-size(cyclePeakInds,2)), ' peak data points.\n'])
cycleTroughInds = find(cyclesTroughs(1,:)); cyclesSize = size(cyclesTroughs,2);
badCycles = root.user_def.cleanData_inds2cut(cycleTroughInds);
cycleTroughInds(badCycles) = []; 
fprintf(['Cleaned ', num2str(cyclesSize-size(cycleTroughInds,2)), ' trough data points.\n'])

% epoch the cycles
epochSize = 0.100; cycleTs = root.b_ts(cyclePeakInds); 
epochs = [cycleTs-epochSize cycleTs+epochSize];  
root.epoch = epochs; root.b_myvar = cycles'; cyclePeaksEpoched = root.myvar;

epochSize = 0.100; cycleTs = root.b_ts(cycleTroughInds); 
epochs = [cycleTs-epochSize cycleTs+epochSize];
root.epoch = epochs; root.b_myvar = cycles'; cycleTroughsEpoched = root.myvar;

% to make uniform size: drop short epochs/trim long epochs, then make cycle mats 
nSampP = cellfun(@(c) length(c), cyclePeaksEpoched,'uni',1);
minSize = root.user_def.lfp_fs * 2 * epochSize; shortEpochs = nSampP < minSize; 
cyclePeaksEpoched(shortEpochs) = []; cycleTs(shortEpochs) = []; 
epochs = [cycleTs-epochSize cycleTs+epochSize]; root.epoch = epochs; 

for i = 1:size(cyclePeaksEpoched, 1)
  cyclePeaksEpoched{i} = cyclePeaksEpoched{i}(1:minSize, 1:size(cyclePeaksEpoched{i}, 2));
end

nSampT = cellfun(@(c) length(c), cycleTroughsEpoched,'uni',1);
minSize = root.user_def.lfp_fs * 2 * epochSize; shortEpochs = nSampT < minSize; 
cycleTroughsEpoched(shortEpochs) = []; cycleTs(shortEpochs) = []; 
epochs = [cycleTs-epochSize cycleTs+epochSize]; root.epoch = epochs; 

for i = 1:size(cycleTroughsEpoched, 1)
  cycleTroughsEpoched{i} = cycleTroughsEpoched{i}(1:minSize, 1:size(cycleTroughsEpoched{i}, 2));
end

cyclePeaksCycles = cell2mat(cyclePeaksEpoched);
cycleTroughsCycles = cell2mat(cycleTroughsEpoched);


%%

% Now let's grab the inds of each peak
dropCount = 0;
nanCount = 0;
for i = 1:size(cyclePeaksEpoched, 1)
  for j = 1:size(cyclePeaksEpoched{i}, 2)
    ind = find(cyclePeaksEpoched{i}(:,j));
    % if two inds found then keep the one closest to the midpoint of the epoch
    if size(ind, 1) > 1
      [meanDiff, k] = min(abs(ind-(size(cyclePeaksEpoched{i}, 1)/2)));
      ind = ind(k); dropCount = dropCount + 1;
    end
    if ~isempty(ind); cyclePeaksEpochedInds{i}(j,:) = ind;
    else nanCount = nanCount + 1;
    end
  end
end

dropCount = 0;
nanCount = 0;
for i = 1:size(cycleTroughsEpoched, 1)
  for j = 1:size(cycleTroughsEpoched{i}, 2)
    ind = find(cycleTroughsEpoched{i}(:,j));
    % if two inds found then keep the one closest to the midpoint of the epoch
    if size(ind, 1) > 1
      [meanDiff, k] = min(abs(ind-(size(cycleTroughsEpoched{i}, 1)/2)));
      ind = ind(k); dropCount = dropCount + 1;
    end
    if ~isempty(ind); cycleTroughsEpochedInds{i}(j,:) = ind;
    else nanCount = nanCount + 1;
    end
  end
end


% fit the points!
fprintf('Fitting the peak indicies...\n')
for i = 1:size(cyclePeaksEpochedInds, 2)
  if size(cyclePeaksEpochedInds{i},1) < 3
    b = [nan, nan]; stats.coeffcorr(2) = nan;
  else
    % figure; plot(1:size(cyclesEpochedInds{i}, 1), cyclesEpochedInds{i}','x')
    [b,stats] = robustfit(1:size(cyclePeaksEpochedInds{i}, 1), cyclePeaksEpochedInds{i}');
  end
  %p = polyfit(1:size(cyclesEpochedInds{i}, 1),cyclesEpochedInds{i}',1);
  cEPeaksSlopesRobust(i, 1) = b(2); 
  cEPeaksCorCoRobust(i, 1) = stats.coeffcorr(2); 
end

fprintf('Fitting the trough indicies...\n')
for i = 1:size(cycleTroughsEpochedInds, 2)
  if size(cycleTroughsEpochedInds{i},1) < 3
    b = [nan, nan]; stats.coeffcorr(2) = nan;
  else
    % figure; plot(1:size(cyclesEpochedInds{i}, 1), cyclesEpochedInds{i}','x')
    [b,stats] = robustfit(1:size(cycleTroughsEpochedInds{i}, 1), cycleTroughsEpochedInds{i}');
  end
  %p = polyfit(1:size(cyclesEpochedInds{i}, 1),cyclesEpochedInds{i}',1);
  cETroughsSlopesRobust(i, 1) = b(2); 
  cETroughsCorCoRobust(i, 1) = stats.coeffcorr(2); 
end



check = isnan(cETroughsSlopesRobust);
if sum(check) > 0
  ii = find(check)
  cETroughsSlopesRobust(ii) = 0;
end
%%
% let's take a look at each slope
figure; plot(1:size(cEPeaksSlopesRobust, 1), cEPeaksSlopesRobust, '.r');
hold on; plot(1:size(cETroughsSlopesRobust, 1), cETroughsSlopesRobust, '.b');
title([name, ' Robust ', phsTxt,'Slopes'])


ci = bootci(5000,{@(x) median(x),cEPeaksSlopesRobust},'type','per');
figure; histogram(cEPeaksSlopesRobust)
title([name, ' Robust slopes | Peak | median: [', num2str(ci(1)),' ',num2str(ci(2)),']'])
%title([name, ' Robust slopes | shift: ', num2str(phs)]) % statless version
%for debugging b 
xlabel('slope')
ylabel('num cycles')

% And a histogram of slopes, along with stats
ci = bootci(5000,{@(x) median(x),cETroughsSlopesRobust},'type','per');
figure; histogram(cETroughsSlopesRobust)
title([name, ' Robust slopes | Trough | median: [', num2str(ci(1)),' ',num2str(ci(2)),']'])
%title([name, ' Robust slopes | shift: ', num2str(phs)]) % statless version for debugging
xlabel('slope')
ylabel('num cycles')



% and a histogram of Coorolation coef (!) TD needs lookin' at
figure; histogram(cECorCoRobust)
title([name, ' ', phsTxt, ' Robust corCo'])
xlabel('Corrolation Coef')
ylabel('num cycles')

% TD: add code to save figs
% let's make the histogram's look good now
nBins = ceil(sqrt(size(cESlopesRobust,1))); % num bins ~ rounded sqrt size of data
figure;
h = histogram(cESlopesRobust, nBins)

%%
% quality control

% show number of dropped peaks and nans from inds grabbing
fprintf(['Dropped ', num2str(dropCount), ' peaks. (cycles with more than one)\n'])
fprintf(['There are ', num2str(nanCount), ' cycles without peaks.\n'])

% find the number of peaks counted
for i = 1:size(cyclesEpochedInds, 2)
  numPeaks(i, 1) = size(cyclesEpochedInds{i},1);
end
fprintf(['The average number of peaks is ', num2str(mean(numPeaks)), '.\n'])
fprintf(['The mean of your slopes is ', num2str(mean(cESlopesRobust)), '.\n'])


% curious what the mean data looks like?
if 0
  figure; plot(meanCycleInd')
  title(['Mean Indicies of Cycles'])
  ylabel('Mean'); xlabel('epochs')
end





%%
% Scan V2.0 ELN 220218
if scan
  root.b_myvar = origData'; epchLFP = root.myvar;
  root.b_myvar = thetaPhs'; epchThP = root.myvar;
  root.b_myvar = cycles';   epchCyc = root.myvar;
  figure; if ~exist('startInd','var'), startInd = 1; end
  for i = startInd:size(cyclesEpochedInds,2)
    cyclesEpochedTs = cyclesEpochedInds{i}; nChan = size(cyclesEpochedTs, 1);
    p = polyfit(1:nChan,cyclesEpochedTs',1);
    pv = polyval(p, 1:nChan); [b,stats] = robustfit(1:nChan,cyclesEpochedTs);
    b = fliplr(b'); bv = polyval(b, 1:nChan);
    tmpLFP = spreadLFP(epchLFP{i}')';
    hold off; plot(tmpLFP,'k');
    offsets = [0:nChan-1]*size(tmpLFP,1);
    tmpX = repmat(cyclesEpochedTs',2,1); tmpY = [1:nChan;tmpLFP(cyclesEpochedTs+offsets')'];
    hold on; plot(tmpX, tmpY, 'k:'); plot(cyclesEpochedTs', 1:nChan, 'o');
    pvX = repmat(pv,2,1); bvX = repmat(bv,2,1); chanY = reshape(repmat(0:nChan,2,1),[],1);
    plot(pv', 1:size(cyclesEpochedInds{i},1), 'b-.'); plot(bv', 1:size(cyclesEpochedInds{i},1), 'r');
    title(['set ', num2str(i), ' Slope: ', num2str(p(1)), ' RobustS: ', num2str(b(1)), 'RobustCorCo: ',num2str(stats.coeffcorr(2)) ]);
    xlabel('cycInd'); ylabel('offset');
    pause;
  end
end

end 






%%
% -=-=-=-=-=-=-=-= RIP =-=-=-=-=-=-=-=-=-



% % clc
% % clear
% % close all
% nbins=20;
% series1 = [10,25,90,35,16, 8, 25, 55, 55];
% series2 = [7,38,31,50,41,25,90,35,16];
% [series1,centers] = hist(series1,nbins);
% [series2] = hist(series2,centers);
% DataSum=series1+series2;
% figure
% width1 = 0.5;
% bar(centers,DataSum,width1,'FaceColor',[0.2,0.2,0.5],....
%                      'EdgeColor','none');
% hold on
% width2 = width1;
% bar(centers,series2,width2,'FaceColor',[0,0.7,0.7],...
%                      'EdgeColor',[0,0.7,0.7]);
% hold off
% legend('First Series','Second Series') % add legend







