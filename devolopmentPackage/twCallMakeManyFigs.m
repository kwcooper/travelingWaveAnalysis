 fprintf('\nTravleing Wave Project''s Multi-Rat Analysis Package\n');

 %To Do
% 1. tracking with multiple rats
% 2. double check root objects
% 3. fix/add figure saving 
% 4. confirm references/tracking
% 5. build slope vs velocity function
% 6. check timestamps
% 7. check slope gathering functions
% 8. check 

clear all;
dbstop if error;

iteration = 1;
%sess2run = [1 3 5 7 9]; % FAM
%sess2run = [2 4 6 8 10]; % NOV
%sess2run = [2 4 6 8 10 1 3 5 7 9]; %NOV & FAM 
%sess2run = [15 16]; % SCOP & SAL
%sess2run = [2 17 18 19 20]; % best data
sess2run = [20]; % Roble

% check for intraRat struct
ratStruct = 'ratsComp.mat';
structPath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', ratStruct); 

% check if struct exists, if not, check if we want to create it
if ~exist(structPath, 'file')
  makeStruct = 1;
else 
  makeStruct = 1; %input('Intrarat struct found, would you like to make a new one? (type 1 for yes, 0 for no): ');
  if makeStruct
    movefile(structPath, fullfile(structPath(1:length(structPath) - length(ratStruct)),'owStructs', ['ratsComp_', num2str(now), '.mat']));
  end
end
if makeStruct
  fprintf('> Creating intrarat struct... ');
  ratsComp = struct; ratsComp.data = {}; 
  save(structPath,'ratsComp','-v7.3'); fprintf('saved.\n');
  
else
  warning('> Using existing struct!');
end

tic
for i = sess2run %Select which sessions you would like to run
  
  sInd = i; % Selects the session to analyze 
  force = 1; %TD maybe save the forcing to sessions?
  plt = 0; % Generate figures?
  
  sessionsList; %import the session data  
  
  metaData = struct;
  metaData.Rat = sessions{sInd,1};
  metaData.Session = sessions{sInd,2};
  metaData.Recording = sessions{sInd,3};
  metaData.chTxt = sessions{sInd,4};
  metaData.chOrd = sessions{sInd,5};
  metaData.ref = sessions{sInd,6};
  metaData.trackRef = sessions{sInd,7};
  metaData.filtParams = sessions{sInd, 8};
  metaData.badEnds = sessions{sInd,9};
  
  metaData.saveFig = 1; 
  metaData.figDir = 'twImgDir'; 
  metaData.savePath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', metaData.figDir);
  metaData.fig_type = 'png'; % options = {'png','ps','pdf'}
  
  
  %% Path, fetching, preprocessing
  dataPath = fullfile(ratLibPath,metaData.Rat,metaData.Session,metaData.Recording); cd(dataPath);
  fprintf(['> Changed directory \n']);
   
  fprintf(['\n________{ Running session ', num2str(iteration), ' of ', num2str(length(sess2run)), ' }________\n'])
  fprintf(['Rat: ', metaData.Rat, '  Recording: ', metaData.Recording, ' ', metaData.chTxt, '\n']);
  %fprintf(['Working dir: ', workingDir(:,34:end), '\n']);
  
  %Check if the data already exists
  structName = [metaData.Session '_' metaData.Rat '_root.mat'];
  if ~force && exist(fullfile(dataPath, structName), 'file')
    fprintf('> I found an existing root object, loading it in...\n');
    foundData = load(structName);
    
    % Check if the saved data is what we think it is
    fprintf('   > Checking session names...');
    %disp(['\n Checking', foundData.Recording, sessions{sInd,3}])
    if ~strcmp(foundData.rcdng, sessions{sInd,3})
      warning('    > Imported data fields do not match requested session!\n')
      force = 1; % TD add way to move into the forced data recalculation if fields dont match.
      keyboard;
    else
      fprintf(' Checks out!\n');
      root = foundData.root; %set
    end
  
  elseif  force || ~exist(fullfile(dataPath, structName), 'file')
    fprintf('> Processing data.\n')
    fnames = dir(fullfile(dataPath,'100_CH*.continuous'));
    fprintf('  > Fetching the lfp. ');
    [D, fs, fType] = kLoadIntanLFP(metaData.chOrd,dataPath,fnames);
    metaData.fType = fType;
    
    % prepData: makes root, grabs tracking, and minor theta computations
    fprintf('\n  > Prepping data.\n')
    %fprintf('    > (Cutting out bad epochs and grabbing theta cycles)\n')
    root = twPrepData(D,fs,metaData); rcdng = metaData.Recording;
    save(structName,'root','rcdng','-v7.3');
    fprintf('> Saving precomputed data (root object) so next time is a breeze!\n');
  end
 
  
%% Grab the intra-rat struct
foundStruct = load(structPath); fprintf('> Loading in struct\n');
ratsComp = foundStruct.ratsComp;

% Add metadata to the struct
recordingPrime = matlab.lang.makeValidName(metaData.Recording); % Session names arn't compatable
%ratsComp.ratRegress.(metaData.Rat).meta = metaData;
%% Analysis
fprintf('> Running analyses \n')

% peak 0; trough pi;

% TD THIS BROKE WITH THE TRACKING DATA
% Average Theta Wave
epochSize = 0.100;
[CTA_P] = plotCycleTriggeredAvg(root, epochSize,  0, metaData, 0);
[CTA_T] = plotCycleTriggeredAvg(root, epochSize, pi, metaData, 0); 
root.user_def.CTA_P = CTA_P;
root.user_def.CTA_T = CTA_T;

% Make the cross corrolation plot %Need to update this
% twCrossCorr(root)

plt = 1;
% plots the slope of the average offset 
% pd is a struct that holds the line of fit info
pd = twPlotPeakDiff(root, pi, metaData, plt); 
root.user_def.pd = pd;

%grab shifts (Degrees)
shiftsPerChanDeg = twPhaseShift(CTA_P.avgThetaWave);

% plot PSD
% right now we are looking at only one of the channels
%plotPSD(root.lfp.signal(1,1:100000),root)

% Raw LFP: Grabs the raw data from the specified indicies
ind1 = 650; ind2 = 750; % td Set these to the specified ends from sessions
[rawWaves] = twGrabRawData(root.user_def.lfp_origData, ind1, ind2, plt);

% compute theta asym
[asmScores] = twComputeAsym(root);

% compute theta speed modulation
chans = root.user_def.metaData.chOrd;
%[thetaSpdMod_R] = thetaSpdAnalysis(root, chans);

% psd plot: need to pick a channel
%[psd] = plotPSD(,fs,tText);


%%
% update rat struct
ratsComp.idx = ['Rat  ', 'shiftsPerChanDeg  ', 'pd  ', 'asyemetry idx  ', 'metaData  '];
ratsComp.data = [ratsComp.data; {metaData.Rat, shiftsPerChanDeg, pd, asmScores, metaData}];


% Histogram and slope mean 
% td: this needs updated for the new struct version
%interRatExploration(ratRegress, plt)
%%
% save intra-struct & root
fprintf('> Saving root & intrastruct... ');
save(structPath,'ratsComp','-v7.3'); 
save(pwd, 'root', '-v7.3');
fprintf('done. \n');

iteration = iteration + 1;
end
endT = toc; fprintf(['Took ' num2str(endT/60) ' minutes to run.']);

% add multi rat computations here


fprintf('\nfin \n');