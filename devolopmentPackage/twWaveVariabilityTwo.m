
%remember to cut out the bad data!


% lets set some defaults and allocate mats
name = Rat;
origData = root.user_def.lfp_origData; 
epch = CTA.epDat;
atw = CTA.avgThetaWave;
%cycles = root.user_def.cycles; 

thetaPhs = nan(size(origData));
thetaAmp = nan(size(origData));
cycles = nan(size(origData));



% Let's grab the cycles and phase (root has one, but it used hilbert)
for i = 1:size(origData,1)
  thetaPhs(i,:) = extractThetaPhase(origData(i,:),root.user_def.lfp_fs,'waveform',[6 10]);
  [cycles(i,:),~] = parseThetaCycles(thetaPhs(i,:),root.user_def.lfp_fs,[6 10],0); % resetPh to pi instead of 0 for peaks!
end 

%take a look at the cycles and the orgional data
figure; plot(origData(:,1:3000)');
hold on; plot(1000*cycles(:,1:3000)');

% Let's clean up cycles a bit to remove the bad data
cycleInds = find(root.user_def.cycles(1,:)); cyclesSize = size(cycles,2);
badCycles = root.user_def.cleanData_inds2cut(cycleInds);
cycleInds(badCycles) = []; 
fprintf(['Cleaned ', num2str(cyclesSize-size(cycleInds,2)), ' data points.\n'])

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
cyclesEpoched(shortEpochs) = []; % drop short epochs

% Resize cells
for i = 1:size(cyclesEpoched, 1)
  cyclesEpoched{i} = cyclesEpoched{i}(1:minSize, size(cyclesEpoched{i}, 2));
end
%nSamp = cellfun(@(c) length(c), epchData,'uni',1); %for testing




%%%%% find out where the other channels are?!

figure; 
for i = 1:size(cyclesEpoched, 1)
  hold off;
  plot(size(cyclesEpoched{i}))
end















cyclesCycles = cell2mat(cyclesEpoched);


scan = 1;
if scan
  for j = 1:size(cyclesEpoched,3)
    for i = 1:size(cyclesEpoched,2)
  hold off;
  figure; hold on; plot(cyclesEpoched(:,:,j)); %plot(1:100,100 * cycles(:,i,j), 'k');
    end
  end
end




epchData = root.user_def.thetaEpoch; % calculated from the plotCycleTriggered average code


for j = 1:size(epchData,3)
  for i = 1:size(epchData,2)
  thetaPhs(:,i,j) = extractThetaPhase(epchData(:,i,j),root.user_def.lfp_fs,'waveform',[6 10]);
  [cycles(:,i,j),~] = parseThetaCycles(thetaPhs(:,i,j),root.user_def.lfp_fs,[6 10],0); % resetPh to pi instead of 0 for peaks!
  end
end

scan = 1;
if scan
  for j = 1:size(epchData,3)
    for i = 1:size(epchData,2)
  hold off;
  figure; hold on; plot(epchData(:,:,j)); plot(1:100,100 * cycles(:,i,j), 'k');
    end
  end
end









