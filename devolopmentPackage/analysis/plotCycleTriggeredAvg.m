function [CTA] = plotCycleTriggeredAvg(root, epochSize, metaData, plt)

%!! Figure out which channel to do the calculations on (Reference?)
% TD add adjustments for the anatomical data. 

% grab cycles and clean bad ones
cycles = find(root.user_def.cycles(1,:));
badCycles = root.user_def.cleanData_inds2cut(cycles);
cycles(badCycles) = [];

% ud 180528: changed root.b_ts to lfp ts after addition of tracking data.  
cycleTs=root.b_lfp(1).ts(cycles); %finds all theta cycles
% finds bad theta cycles
%epochSize = 0.100; %sets epoch size (moved to function argument)
epochs = [cycleTs-epochSize cycleTs+epochSize]; % grabs epochs 
root.epoch = epochs; 

root.b_myvar = root.user_def.lfp_origData'; % time across rows
epchData = root.myvar; % this is then epoched with cmbhome
% drop short epochs and trim long epochs until all are the same length
nSamp = cellfun(@(c) length(c), epchData,'uni',1);
epchLength = root.user_def.lfp_fs * 2 * epochSize;
shortEpochs = nSamp < epchLength; epchData(shortEpochs) = [];
epchData = cellfun(@(c) c(1:epchLength,:),epchData,'uni',0);  

% now compute the mean wave
epchData = cat(3,epchData{:});
avgThetaWave = mean(epchData,3);


[tPts,nElecs] = size(avgThetaWave);
t = linspace(-epochSize,epochSize,epchLength);
lfp_ = avgThetaWave' / (-1 * 2.5 * rms(avgThetaWave(:)));
offsets = repmat([1:nElecs]',1,tPts);
lfp_ = lfp_ + offsets;

if plt % if breaks move h = figure back to top of chunck
  h = figure;
  plot(t,lfp_,'k'); axis ij;
  %Change axis to reflect proper channels
  
  xlabel('Time'); ylabel('Channel'); % need to update to distance
  title(['Averaged Waves: ', metaData.Rat])
  % ltr = text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top');
  % s = ltr.FontSize;
  % ltr.FontSize = 12;
  grid on
  
  
  plotName = [metaData.Recording '_' metaData.Rat];
  printFigure(gcf, [fullfile(metaData.savePath, 'cycTrigAvg',[plotName,'.',metaData.fig_type])],'imgType',metaData.fig_type);
  fprintf('Saved figure (CycTrig Avg)\n');
end

%Save and reorrient 
CTA.avgThetaWave = avgThetaWave';
CTA.lfp_ = lfp_;
CTA.t = t;
CTA.epDat = epchData;
root.user_def.atw = avgThetaWave';

%% for poster
% h = figure;
% plot(t,lfp_); axis ij;
% %Change axis to reflect proper channels
% ax = gca;
% ax.Visible = 'off';
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlabel('Time') %!! Is this correct? or should it be phase?
% ylabel('Channel') % !! What about the axis though...
% title(['Averaged Waves: ', figData.ratInfo.name])



end
