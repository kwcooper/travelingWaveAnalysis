fprintf('\nWelcome to the Travleing Wave Project''s Multi-Rat Analysis Package!\n');
%To Do
%  Think of how to save the data from all rats to one struct... 
%      load it in from one location?
%  FLip tio's channel map
%  nov vs fam
clear all;


iteration = 1;
sess2run = [2 5 11 12 18 20 25]; %[2 5 11 12 13 14 18 20] best set without roble
ratLibPath = fullfile(ratLibPath,'ratsEphys'); % temp due to path errors
for i = sess2run %Select which sessions you would like to run
  
  sInd = i; % This selects the session you want to analyze TD: add user input
  force = 0; % maybe save the forcing to sessions?
  plt = 0; % would you like the code to generate figures?
  
  Rat =  sessions{sInd,1};
  Session =  sessions{sInd,2};
  Recording =  sessions{sInd,3};
  chTxt = sessions{sInd,4};
  chOrd = sessions{sInd,5};
  ref = sessions{sInd,6};
  badEnds = sessions{sInd,7}; % !! could update this to good ends?
  
  %% Path, fetching, preprocessing
  dataPath = fullfile(ratLibPath,Rat,Session,Recording); cd(dataPath);
  fprintf('> Changed directory\n');
  
  %dataPath = 'C:\Users\DarthMaul\Downloads\2017-12-19_20-04-50\2017-12-19_20-04-50'; fprintf('USING TEMP ROBLE PATH - FIX ME!\n');
  %fprintf(['__________________________________\n']);
  fprintf(['\n________{ Running session ', num2str(iteration), ' of ', num2str(length(sess2run)), ' }________\n'])
  fprintf(['Rat: ', Rat, '  Recording: ', Recording, '\n']);
  %fprintf(['Working dir: ', workingDir(:,34:end), '\n']);
  
  structName = [Session '_' Rat '_root.mat'];
  %Check if the data already exists
  if ~force && exist(fullfile(dataPath, structName), 'file')
    fprintf('> I found an existing root object, loading it in...\n');
    foundData = load(structName);
    
    % Check if the saved data is what we think it is
    fprintf('> > Checking session names...');
    if ~strcmp(foundData.Recording, sessions{sInd,3})
      warning('> > > Imported data fields do not match requested session!\n')
      force = 1;
      % TD add way to move into the forced data recalculation if fields
      % dont match.
      keyboard;
    else
      fprintf('Checks out!\n');
      root = foundData.root; %set
    end
  
  elseif  force || ~exist(fullfile(dataPath, structName), 'file')
    fprintf('> Processing data.\n')
    fnames = dir(fullfile(dataPath,'100_CH*.continuous'));
    fprintf('> > Fetching the lfp.\n')
    fprintf('(This may take some time)...\n')
    [D, fs] = kLoadIntanLFP(chOrd,dataPath,fnames);
    
    fprintf('> > Prepping the data.\n')
    fprintf('> > > (Cutting out bad epochs and grabbing theta cycles)\n')
    root = twPrepData(D,fs,ref);
    
    save(structName,'root','Recording','-v7.3');
    fprintf('> Saving precomputed data (root) so next time is a breeze!\n');
  end
 
%% Prepping
figData = struct;
figData.ratInfo.name = Rat;
figData.ratInfo.session = Session;
figData.ratInfo.recording = Recording;
figData.ratInfo.chOrdTxt = chTxt;
figData.ratInfo.ref = ref;

figData.saveFig = 1; %would you like the code to save the figures?
figData.figDir = 'twImgDir'; % where would you like them saved?
figData.savePath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', figData.figDir);
figData.fig_type = 'png'; % options = {'png','ps','pdf'} 

%Grab the intra rat struct
ratStruct = 'ratsComp.mat';
structPath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', ratStruct); 
if ~exist(structPath, 'file')
  % check if struct exists, if not, create it!
  fprintf('> No rat structure found, creating one... ');
  ratsComp = struct;
  save(structPath,'ratsComp','-v7.3'); 
  fprintf('saved.\n');
  foundStruct = load(structPath);
  ratsComp = foundStruct.ratsComp;
else
  % load it back in 
  foundStruct = load(structPath);
  fprintf('> Loading in struct\n');
  ratsComp = foundStruct.ratsComp;
end


% Add metadata to the struct
recordingPrime = matlab.lang.makeValidName(Recording); % Session names arn't compatable
ratsComp.ratRegress.(Rat).(recordingPrime) = struct;
ratsComp.ratRegress.(Rat).(recordingPrime).info = figData.ratInfo;

%% Plotting
% Makes the average wave plot  b
% td add argument for window width (two cycles)
epochSize = 0.100;

[CTA] = plotCycleTriggeredAvg(root, epochSize, figData, plt);
figData.CTA = CTA;

% Makes the cross corrolation plot %Need to update this
% twCrossCorr(root)
pd = twPlotPeakDiff(root, figData, plt);
figData.pd = pd;
ratsComp.ratRegress.(Rat).(recordingPrime).pd = figData.pd;
% Raw LFP: Grabs the raw data from the specified indicies
% td Set these to the specified ends

ind1 = 650;
ind2 = 750;
[rawWaves] = twGrabRawData(root.user_def.lfp_origData, ind1, ind2, plt);
figData.rawWaves = rawWaves;

%grab shifts (Degrees)
shiftsPerChanDeg = twPhaseShift(CTA.avgThetaWave);
ratsComp.shiftsPerChanDeg.(Rat) = shiftsPerChanDeg; 

% Histogram and slope mean
%interRatExploration(ratRegress, plt)
%%
% save the struct
save(structPath,'ratsComp','-v7.3'); 
iteration = iteration + 1;
end

fprintf('\nfin \n');