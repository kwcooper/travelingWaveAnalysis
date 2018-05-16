 fprintf('\nTravleing Wave Project''s Multi-Rat Analysis Package\n');
%To Do
%  store all in root and load in each root per rat...?
%  add intra sesson f() after loop
%  add time base alignment
%      Should probably verson this package before this
clear all;


iteration = 1;
%sess2run = [1 3 5 7 9]; % FAM
%sess2run = [2 4 6 8 10]; % NOV
%sess2run = [2 4 6 8 10 1 3 5 7 9]; %NOV & FAM 
sess2run = [15 16]; % SCOP & SAL
sess2run = [1];
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
  metaData.badEnds = sessions{sInd,8};
  
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
    fprintf('  > Fetching the lfp.\n')
    fprintf('(This may take some time)...\n')
    [D, fs, fType] = kLoadIntanLFP(metaData.chOrd,dataPath,fnames);
    metaData.fType = fType;
    
    % prepData: makes root, grabs tracking, and minor theta computations
    fprintf('  > Prepping data.\n')
    %fprintf('    > (Cutting out bad epochs and grabbing theta cycles)\n')
    root = twPrepData(D,fs,metaData); rcdng = metaData.Recording;
    save(structName,'root','rcdng','-v7.3');
    fprintf('> Saving precomputed data (root object) so next time is a breeze!\n');
  end
 
  
%%
%Grab the intra-rat struct
ratStruct = 'ratsComp.mat';
structPath = fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Keiland','Projects','travelingWave', ratStruct); 

% check if struct exists, if not, create it!
if ~exist(structPath, 'file')
  fprintf('> No rat structure found, creating one... ');
  ratsComp = struct; save(structPath,'ratsComp','-v7.3'); fprintf('saved.\n');
  foundStruct = load(structPath); ratsComp = foundStruct.ratsComp;
  ratsComp.data = {}; % add cell array
else % load it back in
  foundStruct = load(structPath); fprintf('> Loading in struct\n');
  ratsComp = foundStruct.ratsComp;
end

% Add metadata to the struct
recordingPrime = matlab.lang.makeValidName(metaData.Recording); % Session names arn't compatable
%ratsComp.ratRegress.(metaData.Rat).meta = metaData;


%% Plotting
fprintf('> Running analyses \n')
% Average Theta Wave
epochSize = 0.100;
[CTA] = plotCycleTriggeredAvg(root, epochSize, metaData, 1); 
root.user_def.CTA = CTA;

% Make the cross corrolation plot %Need to update this
% twCrossCorr(root)

% plots the slope of the average offset 
% pd is a struct that holds the line of fit info
pd = twPlotPeakDiff(root, metaData, plt); 
root.user_def.pd = pd;
%ratsComp.ratRegress.(metaData.Rat).(recordingPrime).pd = metaData.pd;


%grab shifts (Degrees)
shiftsPerChanDeg = twPhaseShift(CTA.avgThetaWave);
%ratsComp.shiftsPerChanDeg.(metaData.Rat) = shiftsPerChanDeg; 

% Raw LFP: Grabs the raw data from the specified indicies
ind1 = 650; ind2 = 750; % td Set these to the specified ends from sessions
[rawWaves] = twGrabRawData(root.user_def.lfp_origData, ind1, ind2, plt);
%metaData.rawWaves = rawWaves;

% compute theta asym
[asmScores] = twComputeAsym(root);

%%
% update rat struct
ratsComp.data = [ratsComp.data; {metaData.Rat, shiftsPerChanDeg, pd, asmScores, metaData}];


% Histogram and slope mean 
% td: this needs updated for the new struct version
%interRatExploration(ratRegress, plt)
%%
% save intra-struct & root
fprintf('> Saving root & intrastruct... ');
save(structPath,'ratsComp','-v7.3'); 
save(pwd, 'root', '-v7.3')
fprintf('done. \n');
iteration = iteration + 1;
end

% add post single rat computations here


fprintf('\nfin \n');