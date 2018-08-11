function pd = twPlotPeakDiff(root, pt, metaData, plt)
% Runs a linear regression through peak offset points; plots 
keyboard;
% TD: change axis to mm 

% grab cycles according to pt
if pt == 0
  thetaPhase = root.user_def.theta_phs;
  phs = 'Peak';
elseif pt == pi
  thetaPhase = root.user_def.theta_phsTrough;
  phs = 'Trough';
end

%Finds the circular difference (currently just for the first chan to save time)
dataShift = nan(size(thetaPhase(1,:),1), size(thetaPhase,1));
for iC1 = 1:size(thetaPhase(1,:),1)
  for iC2 = 1:size(thetaPhase,1)
    dataShift(iC1,iC2) = circmean(circDiff([thetaPhase(iC1,:)', thetaPhase(iC2,:)'],2,'rad'),2);
  end
end

figure; sze = 1000;
subplot(3,1,1); plot(thetaPhase(1,1:sze));
subplot(3,1,2); plot(thetaPhase(2,1:sze));
subplot(3,1,3); plot(circmean(circDiff([thetaPhase(1,1:sze)', thetaPhase(1,1:sze)'],2,'rad'),2))


[dim, chan, len] = size(dataShift);
dataShift = cumsum([0 circDiff(dataShift',1, 'rad')]);

%! make sure to adjust for all channels
x = 1:chan;
B = regress(dataShift', [x' ones(size(x))']);
y = B(1)*x + B(2);

%% 

% make distance series
distAx(1) = 300; icmt = 560;
for i = 2:size(x,2)
  distAx(i) = distAx(i-1) + icmt;
end

if plt % if brk move y = to after hold
  %check if result is reasonable.
  figure;
  plot(x,dataShift,'o');
  hold on;
  plot(x,y);
  title([metaData.Rat, ' ', phs, ' offset | Slope: ' num2str(B(1))]);
  %Axis: circular difrence (phase offset) and channels
  
  plotName = [metaData.Recording '_' metaData.Rat '_' phs];
  printFigure(gcf, [fullfile(metaData.savePath, 'peakDiff',[plotName,'.',metaData.fig_type])],'imgType',metaData.fig_type);
  fprintf('Saved figure (peakDiff)\n');
end

pd.x = x; % elecNum
pd.dataShift = dataShift; %diffMat
pd.y = y; % phaseOffset
pd.B = B; % elecPhaseSlopeInt
end