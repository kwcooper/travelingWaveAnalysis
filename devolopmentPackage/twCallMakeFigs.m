fprintf('\nWelcome to the Travleing Wave Project''s Analysis Package!\n');
%To Do
%  Make subfigs; 
%  Save root and add force to avoid calculation and recalculation
%  Add raw lfp ends to the sessions structure

clear all;
sInd = 23; % This selects the session you want to analyze TD: add user input
force = 1; % forces a recalculation of the data 

%% Contains ephys file information, as well as channel mappings, and good epoch data
% Rat          Session                     Recording                notes        Channels                                     reference ends
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
'Romo',        'CircleTrack_2017-11-21',    '2017-11-21_15-29-29',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-21',    '2017-11-21_16-49-48',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];... % 12
'Romo',        'CircleTrack_2017-11-22',    '2017-11-22_17-45-40',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-22',    '2017-11-22_18-03-50',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-24',    '2017-11-24_15-59-32',   'sleep',     [56 54 57 61 47 45 38 33],                         1, [];... %15
'Romo',        'CircleTrack_2017-11-26',    '2017-11-26_17-28-04',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-26',    '2017-11-26_17-49-55',   'sleep',     [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-27',    '2017-11-27_12-44-08',   'smlBst',    [56 54 57 61 47 45 38 33],                         1, [];... %18
'Romo',        'CircleTrack_2017-11-27',    '2017-11-27_13-10-45',   'sleep',     [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-28',    '2017-11-28_14-39-46',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-28',    '2017-11-28_14-39-46',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];... %21
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-19_20-04-50',   'all',       [28 32 4 8 11 15 19 23 26 30 2 6 9 13 17 21],      16, [650 750];...  
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-19_20-04-50',   'sml',       [28 32 4 8 11 15 19 23],                           1,  [650 750];...
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-20_17-31-46',   'all',       [28 32 4 8 11 15 19 23 26 30 2 6 9 13 17 21],      16, [650 750];... %24
'Romo',        'CircleTrack_2017-11-27',    '2017-11-27_12-44-08',   'goodCh',    [57 61 47 38 33],                                   1, [];...
%'Regio'

};

Rat =  sessions{sInd,1};
Session =  sessions{sInd,2};
Recording =  sessions{sInd,3};
chTxt = sessions{sInd,4};  
chOrd = sessions{sInd,5}; 
ref = sessions{sInd,6};
badEnds = sessions{sInd,7}; % !! could update this to good ends?

%% Path, fetching, preprocessing
fprintf(['Grabbing data for ', Rat, '. Recording # ', Recording, '\n']);
chS = mat2str(chOrd); fprintf(['Channels: ',chS(2:end-1),'\n']);
rlb = fullfile(dropboxPath, 'ratsEphys');
dataPath = fullfile(rlb,Rat,Session,Recording); cd(dataPath);
%dataPath = 'C:\Users\DarthMaul\Downloads\2017-12-19_20-04-50\2017-12-19_20-04-50'; fprintf('USING TEMP ROBLE PATH - FIX ME!\n');

%Check if the data already exists
structName = [Session '_' Rat '_root.mat'];
%Check if the data already exists
if ~force && exist(fullfile(dataPath, structName), 'file')
  fprintf('I found an existing root object, loading it in...\n');
  foundData = load(structName);
  
  % Check if the saved data is what we think it is
  fprintf('Checking session names...');
  if ~strcmp(foundData.Recording, sessions{sInd,3})
    warning('Imported data fields do not match requested session!\n')
    force = 1;
    % TD add way to move into the forced data recalculation if fields
    % dont match.
    keyboard;
  else
    fprintf('Checks out!\n');
    root = foundData.root; %set
  end
end

if force
  fprintf('\nProcessing data.\n')
  fnames = dir(fullfile(dataPath,'100_CH*.continuous'));
  fprintf('\nFetching the lfp.\n')
  fprintf('This may take some time...\n')
  [D, fs] = kLoadIntanLFP(chOrd,dataPath,fnames);
  
  fprintf('\nPrepping the data.\n')
  fprintf('(Cutting out bad epochs and grabbing theta cycles)\n')
  root = twPrepData(D,fs,ref);
  
  save(structName,'root','Recording','-v7.3');
  fprintf('Saving precomputed data (root) so next time is a breeze!\n');
end


% root.user_def.name = Rat;
% root.user_def.name = Rat;
% root.user_def.session = Session;
% root.user_def.recording = Recording;
% root.user_def.chOrdTxt = chTxt;
% root.user_def.ref = ref;

%% Plotting
figData = struct;
figData.ratInfo.name = Rat;
figData.ratInfo.session = Session;
figData.ratInfo.recording = Recording;
figData.ratInfo.chOrdTxt = chTxt;
figData.ratInfo.ref = ref;

plt = 0;

figData.saveFig = 1;
figData.figDir = 'twImgDir';
figData.savePath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', figData.figDir);
figData.fig_type = 'png'; % options = {'png','ps','pdf'}


%Makes the average wave plot  b
% td add argument for window width (two cycles) 
epochSize = 0.140; %one hump 0.100;
[CTA] = plotCycleTriggeredAvg(root, epochSize, figData, plt);
figData.CTA = CTA;

%Makes the cross corrolation plot %Need to update this
%twCrossCorr(root)
pd = twPlotPeakDiff(root, figData, plt);
figData.pd = pd;

% Raw LFP: Grabs the raw data from the specified indicies
% td Set these to the specified ends (romo 2Hps: 650 - 815)
ind1 = 720;
ind2 = 1251; %866;
[rawWaves] = twGrabRawData(root.user_def.lfp_origData, ind1, ind2, plt);
figData.rawWaves = rawWaves;

% function to see variability in slope of data
twWaveVariability(root, figData, pi,  0, 0)
%% SubPlotting
%pause;
%twCombineFigs(figData, plt);




