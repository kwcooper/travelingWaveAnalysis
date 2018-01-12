function [ds] = twPlotPeakDiff(D, root)
keyboard;
%need to filter and run through hilbert... 
thetaPhase = root.user_def.theta_phs;
cycles = root.user_def.cycles;

%Finds the circular difference (currently just for the first chan to save time)
  for iC1 = 1:size(D(1,:),1)
    for iC2 = 1:size(D,1)
      dataShift(iC1,iC2,:) = circDiff([thetaPhase(iC1,:)', thetaPhase(iC2,:)'],2,'rad'); 
      % changed orientation of data
    end
  end

[dim, chan, len] = size(dataShift);
if dim == 1
  ds = reshape(dataShift, chan, len);
else
  warning('\n peakDiff, multiple dims found...\n');
  keyboard;
end


%lets try with avg theta wave
atw = figData.CTA.avgThetaWave;
for iC1 = 1:size(D(1,:),1)
  for iC2 = 1:size(D,1)
    dataShift(iC1,iC2,:) = circDiff([thetaPhase(iC1,:)', thetaPhase(iC2,:)'],2,'rad');
    % changed orientation of data
  end
end

% we should probably epoch it now...


% Then average the epochs... (output should be a single number?)
pocs = 

% Then   



end