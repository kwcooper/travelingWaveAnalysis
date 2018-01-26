fprintf('\nWelcome to the Travleing Wave Project''s Multi-Rat Analysis Package!\n');
%To Do
%  Think of how to save the data from all rats to one struct... 
%      load it in from one location?
%  FLip tio's channel map
clear all;

%% Contains ephys file information, as well as channel mappings, and good epoch data
% Rat          Session                     Recording              selecte d           Channels?                                reference
sessions = {...
'RioNovelty',  '2017-07-27_16-00',          '2017-08-09_16-57-57',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8, []; ... % 1
'Rio',         '2017-08-10_CircleTrack',    '2017-08-10_19-14-01',   'sml',       [39 38 60 57],                                     1, []; ... % 2
'Rio',         '2017-08-10_CircleTrack',    '2017-08-10_19-14-01',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8, []; ... % 3
'Rio',         '2017-08-10_CircleTrack',    '2017-08-10_11-43-00',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8, []; ... % 4
'Rio',         '2017-08-10_CircleTrack',    '2017-08-10_11-43-00',   'sml',       [39 38 60 57],                                     1, []; ...
'Rio',         '2017-08-11_CircleTrack',    '2017-08-11_12-54-28',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8, []; ... % 6
'Rio',         '2017-08-20_CircleTrack',    '2017-08-20_12-41-36',   'all',       [11 12 14 13 8 7 5 6 27 28 26 25 20 19 21 22],     8, [];...
'Rio',         '2017-08-22_CircleTrack',    '2017-08-22_14-01-24',   'all',       [11 12 14 13 8 7 5 6 27 28 26 25 20 19 21 22],     8, [];... % 8
'Tio',         '170717_1824_CircleTrack',   '2017-07-17_18-30-47',   'sml',       [43 46 40 37 59 58 52 53],                         8, [];...
'Tio',         '170717_1824_CircleTrack',   '2017-07-17_18-30-47',   'sml',       [40 37 59 58],                                     1, [];... % 10
'Tio',         '170717_1824_CircleTrack',   '2017-07-17_18-30-47',   'smlFlp',    [53 52 58 59 37 40 46 43],                         8, [];...
'Tio',         '170717_1824_CircleTrack',   '2017-07-17_18-30-47',   'smlFlp',    [58 59 37 40],                                     1, [];... % 12
'Romo',        'CircleTrack_2017-11-21',    '2017-11-21_15-29-29',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-21',    '2017-11-21_16-49-48',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];... % 14
'Romo',        'CircleTrack_2017-11-22',    '2017-11-22_17-45-40',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-22',    '2017-11-22_18-03-50',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-24',    '2017-11-24_15-59-32',   'sleep',     [56 54 57 61 47 45 38 33],                         1, [];... %17
'Romo',        'CircleTrack_2017-11-26',    '2017-11-26_17-28-04',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-26',    '2017-11-26_17-49-55',   'sleep',     [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-27',    '2017-11-27_12-44-08',   'smlBst',    [56 54 57 61 47 45 38 33],                         1, [];... %20
'Romo',        'CircleTrack_2017-11-27',    '2017-11-27_13-10-45',   'sleep',     [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-28',    '2017-11-28_14-39-46',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-28',    '2017-11-28_14-39-46',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];... %23
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-19_20-04-50',   'all',       [28 32 4 8 11 15 19 23 26 30 2 6 9 13 17 21],      16, [650 750];...  
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-19_20-04-50',   'sml',       [28 32 4 8 11 15 19 23],                           1,  [650 750];...
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-20_17-31-46',   'all',       [28 32 4 8 11 15 19 23 26 30 2 6 9 13 17 21],      16, [650 750];... %26
};

iteration = 1;
sess2run = [2 5 11 12 18 20]; %[2 5 11 12 13 14 18 20] best set without roble
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
ratRegressStruct = 'ratRegress.mat';
structPath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', ratRegressStruct); 
if ~exist(structPath, 'file')
  % check if struct exists, if not, create it!
  fprintf('> No regression structure found, creating one... ');
  ratRegress = struct;
  save(structPath,'ratRegress','-v7.3'); 
  fprintf('saved.\n');
  foundStruct = load(structPath);
  ratRegress = foundStruct.ratRegress;
else
  % load it back in 
  foundStruct = load(structPath);
  fprintf('> Loading in struct\n');
  ratRegress = foundStruct.ratRegress;
end
% Add metadata to the struct
recordingPrime = matlab.lang.makeValidName(Recording); % Session names arn't compatable
ratRegress.(Rat).(recordingPrime) = struct;
ratRegress.(Rat).(recordingPrime).info = figData.ratInfo;

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
ratRegress.(Rat).(recordingPrime).pd = figData.pd;
% Raw LFP: Grabs the raw data from the specified indicies
% td Set these to the specified ends

ind1 = 650;
ind2 = 750;
[rawWaves] = twGrabRawData(root.user_def.lfp_origData, ind1, ind2, plt);
figData.rawWaves = rawWaves;

% Histogram and slope mean
%interRatExploration(ratRegress, plt)
%%
% save the struct
save(structPath,'ratRegress','-v7.3'); 
iteration = iteration + 1;
end

fprintf('\nfin \n');