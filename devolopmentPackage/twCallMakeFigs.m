fprintf('\nWelcome to the Travleing Wave Project''s Analysis Package!\n');
%To Do
%  Make subfigs; 
%  Save root and add force to avoid calculation and recalculation

sInd = 13; % This selects the session you want to analyze TD: add user input
force = 1; % forces a recalculation of the data 

%% Contains ephys file information, as well as channel mappings, and good epoch data
% Rat          Session                     Recording              selected           Channels?                                reference
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
'Romo',        'CircleTrack_2017-11-22',    '2017-11-22_17-45-40',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];...
'Romo',        'CircleTrack_2017-11-27',    '2017-11-27_12-44-08',   'sml',       [56 54 57 61 47 45 38 33],                         1, [];... % 12
'Roble',       'OpenFieldBigger_2017-12-19','2017-12-19_20-04-50',   'all',       [28 32 4 8 11 15 19 23 26 30 2 6 9 13 17 21],      16, [];... % 
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

dataPath = fullfile(ratLibPath,Rat,Session,Recording); 
%dataPath = 'C:\Users\DarthMaul\Downloads\2017-12-19_20-04-50\2017-12-19_20-04-50'; fprintf('USING TEMP ROBLE PATH - FIX ME!\n');
fnames = dir(fullfile(dataPath,'100_CH*.continuous'));

fprintf('\nFetching the lfp.\n')
fprintf('This may take some time...\n')
[D, fs] = kLoadIntanLFP(chOrd,dataPath,fnames);

fprintf('\nPrepping the data.\n')
fprintf('(Cutting out bad epochs and grabbing theta cycles)')
root = twPrepData(D,fs,ref);

%% Plotting 

%Makes the average wave plot
[h,figData] = plotCycleTriggeredAvg(root);
%figData.CTA.lfp_
%figData.CTA.t


%Makes the cross corrolation plot %Need to update this
%twCrossCorr(root)

%quiverplot and offsets (needs work)
%twPlotQuiver(root)

% Raw LFP


%% SubPlotting

%how to go about this, like the old analysis where we pass the data to it, or 
%calculate the figures in the function...



