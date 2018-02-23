%%
% Requires: Root Object

% td: simulated data
Rat = 'Test';

% lets set some defaults and allocate mats
name = Rat; 
origData = root.user_def.lfp_origData; 
epch = CTA.epDat;
atw = CTA.avgThetaWave;
phs = pi; % pi=toughs; 0=peaks;
%cycles = root.user_def.cycles; 

thetaPhs = nan(size(origData));
thetaAmp = nan(size(origData));
cycles = nan(size(origData));


% Let's grab the cycles and phase (root has one, but it used hilbert)
for i = 1:size(origData,1)
  thetaPhs(i,:) = extractThetaPhase(origData(i,:),root.user_def.lfp_fs,'waveform',[6 10]);
  [cycles(i,:),~] = parseThetaCycles(thetaPhs(i,:),root.user_def.lfp_fs,[6 10], phs); 
end 

%take a look at the cycles and the orgional data
if 0
  figure; plot(origData(:,1:3000)');
  hold on; plot(1000*cycles(:,1:3000)');
end

% Let's clean up cycles a bit to remove the bad data
%cycleInds = find(root.user_def.cycles(1,:)); cyclesSize = size(cycles,2); % Legacy code
cycleInds = find(cycles(1,:)); cyclesSize = size(cycles,2);
badCycles = root.user_def.cleanData_inds2cut(cycleInds);
cycleInds(badCycles) = []; 
fprintf(['Cleaned ', num2str(cyclesSize-size(cycleInds,2)), ' data points. (booted bad data)\n'])

%%
% now let's epoch the cycles, according to cycles!
epochSize = 0.100; 
cycleTs = root.b_ts(cycleInds); %finds all theta cycles
epochs = [cycleTs-epochSize cycleTs+epochSize]; % grabs epochs 
root.epoch = epochs; % change root epoch

root.b_myvar = cycles'; % add to root object to take advantage of the epoching
cyclesEpoched = root.myvar; % wha-la, we have epochs (a cell array)

% to make uniform size: drop short epochs/trim long epochs 
nSamp = cellfun(@(c) length(c), cyclesEpoched,'uni',1);
minSize = root.user_def.lfp_fs * 2 * epochSize;
shortEpochs = nSamp < minSize; 
cyclesEpoched(shortEpochs) = []; cycleTs(shortEpochs) = []; % drop short epochs
epochs = [cycleTs-epochSize cycleTs+epochSize]; % grabs epochs 
root.epoch = epochs; % change root epoch

% Resize cells
for i = 1:size(cyclesEpoched, 1)
  cyclesEpoched{i} = cyclesEpoched{i}(1:minSize, 1:size(cyclesEpoched{i}, 2));
end
%nSamp = cellfun(@(c) length(c), epchData,'uni',1); %for testing

cyclesCycles = cell2mat(cyclesEpoched);

% let's see what the data looks like
if 0
  fprintf('Press enter to see the next slide,\n press ctrl + C when finished!')
  figure;
  for i = 1:size(cyclesEpoched, 1)
    hold off;
    plot(cyclesEpoched{i})
    pause;
  end
end

% alt: we could check for epochs with only 4 peaks
%     sum(cyclesEpoched{i}, _) > 2 % or check other demention also for > 4  



cellSize = size(cyclesEpoched, 1);
numChan = size(cyclesEpoched{1}, 2);
testData = cell(cellSize, 1);
per = 0;
for i = 1:cellSize
  z = zeros(100,numChan);
  % Make perfect dataset
  if per
    nds(1,1) = 48;
    for j = 2:numChan;nds(1,j) = nds(1,j-1) + 101; end
    z(nds) = 1;
  % Make noisy dataset
  else 
    for j = 1:numChan
      ndI = random('Normal',48 + j,.5,1,1);
      z(round(ndI),j) = 1;
    end
  end
  testData{i} = z;
end
 

cyclesEpoched = testData;

% Now let's grab the inds of each peak
dropCount = 0;
nanCount = 0;
for i = 1:size(cyclesEpoched, 1)
  for j = 1:size(cyclesEpoched{i}, 2)
    ind = find(cyclesEpoched{i}(:,j));
    % (!) if two inds found then keep the one closest
    % to the midpoint of the epoch
    if size(ind, 1) > 1
      [meanDiff, k] = min(abs(ind-(size(cyclesEpoched{i}, 1)/2)));
      ind = ind(k);
      dropCount = dropCount + 1;
    end
    % handle if there arn't any peaks
    %     if isempty(ind)
    %       ind = nan;
    %       nanCount = nanCount + 1;
    %     end
    if ~isempty(ind)
      cyclesEpochedInds{i}(j,:) = ind;
    else
      nanCount = nanCount + 1;
    end
  end
end
fprintf(['Dropped ', num2str(dropCount), ' peaks. (cycles with multiple peaks)\n'])
fprintf(['Also, there are ', num2str(nanCount), ' nan''s. (cycles without peaks)\n'])

%%
%let's find the mean ind for all cycles per epoch
for i = 1:size(cyclesEpochedInds, 2)
  meanCycleInd(i, :) = mean(cyclesEpochedInds{i});
end

% (!) handle the case where there arn't cycles... (nan)
if find(isnan(meanCycleInd))
  meanCycleInd(isnan(meanCycleInd)) = 0;
end 

% curious what the mean data looks like?
if 0
  figure; plot(meanCycleInd')
  title(['Mean Indicies of Cycles'])
  ylabel('Mean'); xlabel('epochs')
end

% let's fit some points!
for i = 1:size(cyclesEpochedInds, 2)
  p = polyfit(1:size(cyclesEpochedInds{i}, 1),cyclesEpochedInds{i}',1);
  cyclesEpochedSlopes(i, 1) = p(1); 
end

% let's look a time series of slopes
%figure; plot(cyclesEpochedSlopes)
figure; plot(1:size(cyclesEpochedSlopes, 1), cyclesEpochedSlopes, '.')
% and now a histogram
ci = bootci(5000,{@(x) median(x),cyclesEpochedSlopes},'type','per');
figure; histogram(cyclesEpochedSlopes)
title([name, ' slopes | shift: ', num2str(phs),' | median: [',num2str(ci(1)),' ',num2str(ci(2)),']'])
xlabel('slope')
ylabel('num cycles')

%%

% find the number of peaks counted for analysis quality control
for i = 1:size(cyclesEpochedInds, 2)
  numPeaks(i, 1) = size(cyclesEpochedInds{i},1);
end
fprintf(['The average number of peaks is ', num2str(mean(numPeaks)), '.\n'])
%isnan(cyclesEpochedSlopes)
fprintf(['The mean of your slopes is ', num2str(mean(cyclesEpochedSlopes)), '.\n'])


% let's try it with robust fit
for i = 1:size(cyclesEpochedInds, 2)
  if size(cyclesEpochedInds{i},1) < 3
    b = [nan, nan];
    stats.coeffcorr(2) = nan;
  else
    % figure; plot(1:size(cyclesEpochedInds{i}, 1), cyclesEpochedInds{i}','x')
    [b,stats] = robustfit(1:size(cyclesEpochedInds{i}, 1), cyclesEpochedInds{i}'); % ELN fixed this to have it collect the stats
  end
  %p = polyfit(1:size(cyclesEpochedInds{i}, 1),cyclesEpochedInds{i}',1);
  cESlopesRobust(i, 1) = b(2); 
  %cECorCoRobust(i, 1) = b(1);  % KC code
  cECorCoRobust(i, 1) = stats.coeffcorr(2); % ELN fixed this to make it the correlation coefficient
end
% let's look at it
figure; plot(1:size(cESlopesRobust, 1), cESlopesRobust, '.')
title('Slopes')

% And a histogram of slopes
%ci = bootci(5000,{@(x) median(x),cESlopesRobust},'type','per');
figure; histogram(cESlopesRobust)
%title([name, ' Robust slopes | shift: ', num2str(phs),' | median: [',num2str(ci(1)),' ',num2str(ci(2)),']'])
title([name, ' Robust slopes | shift: ', num2str(phs)])
xlabel('slope')
ylabel('num cycles')

% and a histogram of Coorolation coef
figure; histogram(cECorCoRobust)
title([name, ' Robust corCo'])
xlabel('Corrolation Coef')
ylabel('num cycles')


% look at the differences between slope calculations
%size(cyclesEpochedSlopes)
%size(cESlopesRobust)
%figure; plot(cyclesEpochedSlopes - cESlopesRobust)
%figure; histogram(cyclesEpochedSlopes - cESlopesRobust)


% let's make the histogram's look good now
nBins = ceil(sqrt(size(cESlopesRobust,1))); % num bins ~ rounded sqrt size of data
figure;
h = histogram(cESlopesRobust, nBins)


% Let's look at what we plotted
fprintf('press enter to scan through data!\n');
fprintf('(Be sure to press ctrl + C when finished!)\n');
figure; if ~exist('startInd','var'), startInd = 1; end
for i = startInd:size(cyclesEpochedInds,2)
  hold off;
  p = polyfit(1:size(cyclesEpochedInds{i}, 1),cyclesEpochedInds{i}',1);
  pv = polyval(p, 1:size(cyclesEpochedInds{i}, 1));
  
  b = robustfit(1:size(cyclesEpochedInds{i}, 1),cyclesEpochedInds{i});
  b = fliplr(b');
  bv = polyval(b, 1:size(cyclesEpochedInds{i}, 1));
  plot(cyclesEpochedInds{i}', 1:size(cyclesEpochedInds{i},1), 'o');
  
  hold on;
  plot(pv', 1:size(cyclesEpochedInds{i},1), 'b-.');
  plot(bv', 1:size(cyclesEpochedInds{i},1), 'r');
  
  title(['set ', num2str(i), ' Slope: ', num2str(p(1)), ' RobustS: ', num2str(b(1)), 'RobustCorCo: ',num2str(b(2)) ]);
  xlabel('cycInd'); ylabel('offset');
  pause;
end

% Scan V2.0 ELN 220218
% Pass through root to epoch
  root.b_myvar = origData'; epchLFP = root.myvar;
  root.b_myvar = thetaPhs'; epchThP = root.myvar;
  root.b_myvar = cycles';   epchCyc = root.myvar;
  figure; if ~exist('startInd','var'), startInd = 1; end
for i = startInd:size(cyclesEpochedInds,2)
  cyclesEpochedTs = cyclesEpochedInds{i}; nChan = size(cyclesEpochedTs, 1);
  p = polyfit(1:nChan,cyclesEpochedTs',1);
  pv = polyval(p, 1:nChan);
  [b,stats] = robustfit(1:nChan,cyclesEpochedTs);
  b = fliplr(b');
  bv = polyval(b, 1:nChan);
  
  %plot the raw data
  tmpLFP = spreadLFP(epchLFP{i}')';
  hold off; plot(tmpLFP,'k');
  
  %plot the inds of each waveform
  offsets = [0:nChan-1]*size(tmpLFP,1);
  tmpX = repmat(cyclesEpochedTs',2,1); tmpY = [1:nChan;tmpLFP(cyclesEpochedTs+offsets')'];
  hold on; plot(tmpX, tmpY, 'k:');
  plot(cyclesEpochedTs', 1:nChan, 'o');
  %plot(spreadLFP(epchThP{i}',8)','c');
  %plot(spreadLFP(epchCyc{i}',2)','m:');
  %plot(spreadLFP(cyclesEpoched{i}',2)','m--');

  pvX = repmat(pv,2,1); bvX = repmat(bv,2,1); chanY = reshape(repmat(0:nChan,2,1),[],1);
  plot(pv', 1:size(cyclesEpochedInds{i},1), 'b-.');
  plot(bv', 1:size(cyclesEpochedInds{i},1), 'r');
  
  title(['set ', num2str(i), ' Slope: ', num2str(p(1)), ' RobustS: ', num2str(b(1)), 'RobustCorCo: ',num2str(stats.coeffcorr(2)) ]);
  xlabel('cycInd'); ylabel('offset');
  pause;
end

% Check how the tweedledorph looks over the raw data (Theta/data)
figure; plot(origData(1,1:1000)/500); hold on; plot(thetaPhs(1,1:1000))

%%
% Make simulated data: cyclesEpoched
cellSize = size(cyclesEpoched, 1);
numChan = size(cyclesEpoched{1}, 2);
testData = cell(cellSize, 1);
per = 0;
for i = 1:cellSize
  z = zeros(100,numChan);
  % Make perfect dataset
  if per
    nds(1,1) = 48;
    for j = 2:numChan;nds(1,j) = nds(1,j-1) + 101; end
    z(nds) = 1;
  % Make noisy dataset
  else 
    for j = 1:numChan
      ndI = random('Normal',48 + j,.5,1,1);
      z(round(ndI),j) = 1;
    end
  end
  testData{i} = z;
end
 

%%
% -=-=-=-=-=-=-=-= RIP =-=-=-=-=-=-=-=-=-
%          | Code Graveyard |

% Poetry is the theory of how life should be
% Poetry, like theory...

% % fun with robost fit
% x = (1:10)';
% y = 10 - 2*x + randn(10,1); y(10) = 0;
% bls = regress(y,[ones(10,1) x])
% brob = robustfit(x,y)
% figure; scatter(x,y);
% hold on
% plot(x,brob(1)+brob(2) * x,'r-', x, bls(1)+bls(2)*x,'m:')


% slopeVals= nan(1, size(indsMat, 2));
% for i = 1:size(indsMat, 2)
%   p = polyfit(1:size(indsMat(:,i), 1),indsMat(:,i)',1);
%   slopeVals(1, i) = p(1); 
% end

% 
% for i = 1:size(cyclesEpoched, 1)
%   for j = 1:size(cyclesEpoched{i},2)
%      chanInd = find(cyclesEpoched{i}(:, j));
%      % if more than one value, choose the one closest to the inds mean
%      if size(chanInd, 1) > 1
%        [c, index] = min(abs(chanInd-mean(meanCycleInd)));
%        chanInd = chanInd;
%      cyclesEpochedIndPC{i}(1,j) = chanInd;
%   end
%  end

% let's only look at the good ones

%




% scan = 1;
% if scan
%   for j = 1:size(cyclesEpoched,3)
%     for i = 1:size(cyclesEpoched,2)
%   hold off;
%   figure; hold on; plot(cyclesEpoched(:,:,j)); %plot(1:100,100 * cycles(:,i,j), 'k');
%     end
%   end
% end
% 
% 
% 
% 
% epchData = root.user_def.thetaEpoch; % calculated from the plotCycleTriggered average code
% 
% 
% for j = 1:size(epchData,3)
%   for i = 1:size(epchData,2)
%   thetaPhs(:,i,j) = extractThetaPhase(epchData(:,i,j),root.user_def.lfp_fs,'waveform',[6 10]);
%   [cycles(:,i,j),~] = parseThetaCycles(thetaPhs(:,i,j),root.user_def.lfp_fs,[6 10],0); % resetPh to pi instead of 0 for peaks!
%   end
% end
% 
% scan = 1;
% if scan
%   for j = 1:size(epchData,3)
%     for i = 1:size(epchData,2)
%   hold off;
%   figure; hold on; plot(epchData(:,:,j)); plot(1:100,100 * cycles(:,i,j), 'k');
%     end
%   end
% end
% 
% 
% 






