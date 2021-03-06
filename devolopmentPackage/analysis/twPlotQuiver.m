function [h,figData] = twPlotQuiver(root)
cycles = find(root.user_def.cycles);
badCycles = root.user_def.cleanData_inds2cut(cycles);
cycles(badCycles) = [];


cycleTs=root.b_ts(cycles); %finds all theta cycles
% finds bad theta cycles
epochSize = 0.100; %sets epoch size
epochs = [cycleTs-epochSize cycleTs+epochSize]; % grabs epochs 
root.epoch = epochs; 

root.b_myvar = root.user_def.lfp_origData'; % time across rows
epchData = root.myvar;

% % drop short epochs and trim long epochs until all are the same length
% nSamp = cellfun(@(c) length(c), epchData,'uni',1);
% epchLength = root.user_def.lfp_fs * 2 * epochSize;
% shortEpochs = nSamp < epchLength; epchData(shortEpochs) = [];
% epchData = cellfun(@(c) c(1:epchLength,:),epchData,'uni',0);

%extract theta phases

% Find the angle between vectors
%check orientation of ATW for this to work, also assumes ATW is in radians
%from each other
for i = 1:size(avgThetaWave, 1)-1
  phaseAngles(i,:) = acos(dot(avgThetaWave(i,:),avgThetaWave(i+1,:))/(norm(avgThetaWave(i,:)) * norm(avgThetaWave(i+1,:))));
end
phaseAngles = phaseAngles * 180/pi; 

%from the first channel
for i = 1:size(avgThetaWave, 1)-1
  phaseAngles(i,:) = acos(dot(avgThetaWave(1,:),avgThetaWave(i+1,:))/(norm(avgThetaWave(1,:)) * norm(avgThetaWave(i+1,:))));
end
phaseAnglesFromFirst = phaseAngles * 180/pi; 

% !! NEEDS ATTN Find the theta shift
thetaShift = {};

for iT = 1:size(epchData,3)
  for iC1 = 1:size(epchData, 2) - 1
    for iC2 = iC1 + 1:size(epchData,2) -1
      tmpmean = circmean(cirDiff([root.user_def.theta_phs(iC1,:),root.user_def.theta_phs(iC2,:)],2,'rad'));
      thetaShift{iT}(iC1,iC2) = tmpmean;
      tInfo.thetaShift{iT}(iC2,iC2) = tmpmean;
    end
  end
end
      
      
      tInfo.thetaShift{iT}(iC1,iC2) = tmpmean;
      tInfo.thetaShift{iT}(iC2,iC2) = tmpmean;


              
              

tmp = tInfo.thetaShiftAngle;
avgPkShft = nan(1,size(tmp,1));
for shft = 1:size(tmp,1)
    avgPkShft(shft) = mean(diag(tmp,-shft));
end



end
