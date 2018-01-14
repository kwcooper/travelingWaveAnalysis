function [ds] = twPlotPeakDiff(D, root)
keyboard;
%need to filter and run through hilbert... 
thetaPhase = root.user_def.theta_phs;

%Finds the circular difference (currently just for the first chan to save time)
  for iC1 = 1:size(thetaPhase(1,:),1)
    for iC2 = 1:size(thetaPhase,1)
      dataShift(iC1,iC2) = circmean(circDiff([thetaPhase(iC1,:)', thetaPhase(iC2,:)'],2,'rad'),2); 
      % changed orientation of data
    end
  end

[dim, chan, len] = size(dataShift);


dataShift = cumsum([0 circDiff(dataShift',1, 'rad')]);

%! make sure to adjust for all channels
x = 1:chan;
B = regress(dataShift', [x' ones(size(x))']);

%check result is reasonable.
figure;
plot(x,dataShift,'o');
hold on;
y = B(1)*x + B(2);
plot(x,y);



end