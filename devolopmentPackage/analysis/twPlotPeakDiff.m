function pd = twPlotPeakDiff(root, metaData, plt)
% Runs a linear regression through peak offset points; plots 

% TD: change axis to mm 

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
y = B(1)*x + B(2);

if plt % if brk move y= to after hold
  %check if result is reasonable.
  figure;
  plot(x,dataShift,'o');
  hold on;
  plot(x,y);
  title([metaData.Rat ' peak offset | Slope: ' num2str(B(1))]);
  %Axis: circular difrence (phase offset) and channels
  
  plotName = [metaData.Recording '_' metaData.Rat];
  printFigure(gcf, [fullfile(metaData.savePath, 'peakDiff',[plotName,'.',metaData.fig_type])],'imgType',metaData.fig_type);
  fprintf('Saved figure (peakDiff)\n');
end

pd.x = x; % elecNum
pd.dataShift = dataShift; %diffMat
pd.y = y; % phaseOffset
pd.B = B; % elecPhaseSlopeInt
end