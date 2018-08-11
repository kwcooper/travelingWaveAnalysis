
% Grab data
cd(fullfile('D:','Dropbox (NewmanLab)','docs (1)','docs_Jesus','Projects','cleanDataUpdate','ratCleanLFP'))
files = dir('*.mat');
for i = 1:size(files,1)
  load(files(i).name)
end

%%
fs = 500;
freq = [5 9];

% concatenate ronaldo data
%ronaldo = [ronaldoF1, ronaldoF2, ronaldoF3, ronaldoF4, ronaldoF5];
ronaldo = cleanedData;

ratCleanedData = struct;
ratCleanedData(1).name = 'Regio';    ratCleanedData(1).lfp = cleanedData'; 
ratCleanedData(2).name = 'Rio';      ratCleanedData(2).lfp = rio';         
ratCleanedData(3).name = 'Romo';     ratCleanedData(3).lfp = romo';        
ratCleanedData(4).name = 'Ronaldo';  ratCleanedData(4).lfp = ronaldo'; 
ratCleanedData(5).name = 'Tio';      ratCleanedData(5).lfp = tio';   
ratCleanedData(6).name = 'Roble';    ratCleanedData(6).lfp = roble';   

%%
fprintf('Preprocessing Phase, Amplitude, and Cycles for each rat...\n')
for rat = 1:size(ratCleanedData,2)
  clear thetaPhs thetaAmp cycles % allocate these to save time and rm clear
  for i = 1:size(ratCleanedData(rat).lfp,1)
    [thetaPhs(i,:),thetaAmp(i,:),~] = extractThetaPhase(ratCleanedData(rat).lfp(i,:),fs,'hilbert',freq);
    [cycles(i,:),~] = parseThetaCycles(thetaPhs(i,:),fs,freq,0);         % grab peak data
    %[cyclesTrough(i,:),~] = parseThetaCycles(thetaPhs(i,:),fs,freq,pi);  % grab trough data
  end
  ratCleanedData(rat).thetaPhs = thetaPhs;
  ratCleanedData(rat).thetaAmp = thetaAmp;
  ratCleanedData(rat).cycles = cycles;
end

%% compute xcorr for all rats
fprintf('Computing Cross Correlation for LFP and Amplitude...\n')
for rat = 1:size(ratCleanedData,2)
  standingWaveXcorr(ratCleanedData(rat),'tAmpNormCH1')
  standingWaveXcorr(ratCleanedData(rat),'thetaLFPNorm')
end
















