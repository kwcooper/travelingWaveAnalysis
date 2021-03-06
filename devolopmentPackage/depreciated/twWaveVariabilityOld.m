
function twWaveVariability(root, CTA, name, scan, plt)

name = Rat;
origData = root.user_def.lfp_origData; %remember to cut out the bad data!
epch = CTA.epDat;
atw = CTA.avgThetaWave;
%cycles = root.user_def.cycles; 


thetaPhs = nan(size(origData));
thetaAmp = nan(size(origData));
cycles = nan(size(origData));

%look at the cycles and the orgional data
figure; plot(origData(:,1:3000)');
hold on; plot(1000*cycles(:,1:3000)');



% Let's grab the cycles and phase (root has one, but it used hilbert)
for i = 1:size(origData,1)
  thetaPhs(i,:) = extractThetaPhase(origData(i,:),root.user_def.lfp_fs,'waveform',[6 10]);
  [cycles(i,:),~] = parseThetaCycles(thetaPhs(i,:),root.user_def.lfp_fs,[6 10],0); % resetPh to pi instead of 0 for peaks!
end

%cycles = cycles(1:4, :); % k- let's just look at the first 4 channels
%cycles = cycles([3:5 7:8], :); % jesus's recomendations for romo

%Create a cell array to store the cycle indicies
for i = 1:size(cycles, 1)
  inds{i} = find(cycles(i,:));
end

% find cell bounds and resize all mats to match
[minSize, minIdx] = min(cellfun('size', inds, 2));
[maxSize, maxIdx] = max(cellfun('size', inds, 2));
fprintf(['Cutting ', num2str(maxSize - minSize), ' points.\n']);
for i = 1:size(inds, 2)
  inds{i} = inds{i}(1, 1:minSize);
end

% change orientation and convert to a mat
indsMat = cell2mat(inds');

%scaledI = (I-min(I(:))) ./ (max(I(:)-min(I(:))));


% slope values are caluclated porportionate to the indx values
% we can normalize them for consistancy
%v = max(indsMat'); % Grabs max value for each row indsMat2 = indsMat ./ v';
v = max(indsMat); % Grabs max value for each column
indsMat2 = indsMat ./ v;
indsMat = indsMat2; %temp

% now we can select which channels we want
goodChanRomo = [3:5 7:8];
indsMat = indsMat(goodChanRomo, :); %as per jesus's recomendation
indsMat = indsMat(:,1000:end);%the beginning is super messy... lets cut it out

scan = 1;
if scan
  % paruse the slopes
  figure;
  fprintf('press enter to scan through data!\n');
  for i = 1000:2000
    hold off;
    %plot(indsMat(:,i),1:size(indsMat,1), 'o');
    p = polyfit(1:size(indsMat(:,i),1),indsMat(:,i)',1);
    pv = polyval(p,1:size(indsMat(:,i),1));
    plot(indsMat(:,i)', 1:size(indsMat(:,i)',2), 'o');
    
    hold on;
    plot(pv', 1:size(indsMat(:,i),1));
    title(['set ', num2str(i), ' Slope: ', num2str(p(1))]);
    xlabel('cycInd'); ylabel('offset');
    
    pause;
  end
end

% need's work
figure;
subplot(5,5)
  for i = 1:1000
    if isequal(mod(i,100), 0)
      %plot(indsMat(:,i),1:size(indsMat,1), 'o');
      p = polyfit(1:size(indsMat(:,i),1),indsMat(:,i)',1);
      pv = polyval(p,1:size(indsMat(:,i),1));
      plot(indsMat(:,i)', 1:size(indsMat(:,i)',2), 'o');

      hold on;
      plot(pv', 1:size(indsMat(:,i),1));
      title(['set ', num2str(i), ' Slope: ', num2str(p(1))]);
      xlabel('channel'); ylabel('cycleInd');
    end
    pause;
  end

% find the slope for each column, store in slopeVals
slopeVals= nan(1, size(indsMat, 2));
for i = 1:size(indsMat, 2)
  p = polyfit(1:size(indsMat(:,i), 1),indsMat(:,i)',1);
  slopeVals(1, i) = p(1); 
end

% figure;
% hist(slopeVals)
figure; % this function is better
histogram(slopeVals)
%xlim([-0.001 0.015])
title([name," chan ",num2str(goodChanRomo),"(normalized/messy data cut)"])
xlabel('slope')
ylabel('num cycles')



























% dev 
% tstCyc = cycles(1,:);
% inds = find(tstCyc);
% 
% % finds the size of each interval
% minusVals = [];
% for i = 1:size(inds,2)-1
%   minusVals(:, i) = (inds(i+1) - inds(i));
% end
% epochMean = ceil(mean(minusVals));
% epochS = floor(size(tstCyc, 2) / epochMean);
% %trim data so it fits the reshape
% trim = epochMean * epochS;
% tstCyc = cycles(1,1:trim);
% ttcc = reshape(tstCyc, [epochS, epochMean]);
% ttcc = ttcc';
% old ramblings


% %td: combine these
% % Lets you scroll through the data to see what it looks like.
% figure;
% fprintf('press enter to scan through data!\n');
% for i = 1:100
%   hold off;
%   plot(indsMat(:,i),1:8, 'o');
%   title(num2str(i));
%   pause;
% end
% 
% % plot the first set of points with the fitted slope
% figure;
% i = 1;
% p = polyfit(1:size(indsMat(:,i), 1),indsMat(:,i)',1);
% pv = polyval(p,1:size(indsMat(:,i), 1));
% plot(1:size(indsMat(:,i), 1),indsMat(:,i)', 'o');
% title(['slope', num2str(p(1))]);
% hold on;
% plot(1:size(indsMat(:,i), 1),pv);
% 
% 
% 




% View the phase offset per cycle
%plot(cycles(4:8,1:500)')

%nnz(cycles) % counts number of nonzero elements

% l = [1:4; fliplr(1:4); 3:6]
% k = [3 2];
%
% for i = 1:size(l,1)
%   r = xcorr2(k,l);
% end
%
%
% r = xcorr2(atw, origData(:,1:1000));
% plot(r')
%
%
% %for i = 1:size(origData,1)
% clear r
% for i = 1
%   r(i,:) = xcorr(atw(i,:), origData(i,:));
% end
%
% for i = 1:size(epch, 3):
%
% end

end

