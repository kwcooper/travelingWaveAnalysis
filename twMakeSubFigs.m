 function twMakeSubFigs(force)
% TWMAKESUBFIGS Analysis package for the traveling wave project!
% This is Keiland's version as of 8 Nov 17
% 
% Currently: working & in active devolopment



%2017-08-10_19-14-01 - weird one
%2017-08-10_11-43-00 - fine
%2017-08-11_12-54-28 - looks fine

% Contains ephys file information, as well as channel mappings
% Rat Session Recording selectedChannels? channels reference
% Split up by rat? New file?
sessions = {...
'RioNovelty',  '2017-07-27_16-00',         '2017-08-09_16-57-57',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_19-14-01',   'sml',       [39 38 60 57],                                     1; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_19-14-01',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_11-43-00',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_11-43-00',   'sml',       [39 38 60 57],                                     1; ...
'Rio',         '2017-08-11_CircleTrack',   '2017-08-11_12-54-28',   'all',       [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-20_CircleTrack',   '2017-08-20_12-41-36',   'all',       [11 12 14 13 8 7 5 6 27 28 26 25 20 19 21 22],     8;...
'Rio',         '2017-08-22_CircleTrack',   '2017-08-22_14-01-24',   'all',       [11 12 14 13 8 7 5 6 27 28 26 25 20 19 21 22],     8;...
'Tio',         '170717_1824_CircleTrack',  '2017-07-17_18-30-47',   'sml',       [43 46 40 37 59 58 52 53],                         8;...
'Tio',         '170717_1824_CircleTrack',  '2017-07-17_18-30-47',   'sml',       [40 37 59 58],                                     1;...
'Romo',        'CircleTrack_2017-11-22',   '2017-11-22_17-45-40',   'sml',       [56 54 57 61 47 45 38 33],                         1;...
'Romo',        'CircleTrack_2017-11-27',   '2017-11-27_12-44-08',   'sml',       [56 54 57 61 47 45 38 33],                         1;...
};

% This selects the session you want to analyze
sInd = 2;
 
Rat =  sessions{sInd,1};
Session =  sessions{sInd,2};
Recording =  sessions{sInd,3};
chTxt = sessions{sInd,4};
chOrd = sessions{sInd,5}; 
ref = sessions{sInd,6};

workingDir = fullfile(ratLibPath,Rat,Session,Recording); cd(workingDir);

% forces a recalculation of the data 
if ~exist('force') || isempty(force), force = 0; end
if force
    fprintf('Forcing data recalculation...\n')
end

% This is where all the data extraction is happening, the function is in a
% sepperate file. Definitely worth looking at!
% This returns two structs, root, and tInfo, which are used throughout the rest 
% of the script. Checking for the already computed data happens there too. 
[root,tInfo] = twPrepareDataForFigs(sInd,sessions,chOrd,force);


%%
% !! This needs to be moved or removed...
% new theta Extraction
disp('Extracting theta cycles... (new extraction)');
metho = 'hilbert';
disp('useing ' + metho)
root.epoch=[-inf,inf];
band = [6,10];
[thetaPhs,~,~] = extractThetaPhase(tInfo.signal(ref,:),tInfo.Fs,metho,band);
[cycles,~] = parseThetaCycles(thetaPhs,tInfo.Fs,band); 
inds = find(cycles);
%%


% Set up the figure info
figInfo = {};
figInfo.name = Rat;
figInfo.session = Session;
figInfo.recording = Recording;
figInfo.saveFig = 1;
figInfo.fnameRoot = fullfile('figs',[figInfo.name, '_', figInfo.session]);
figInfo.fig_type = 'png'; % options = {'png','ps','pdf'}
figInfo.chOrdTxt = chTxt;
figInfo.ref = ref;

% makes a folder called figs if it doesn't exist
if ~exist('figs', 'dir')
    mkdir figs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  THE BUSINESS END OF THE FUNCTION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

awData = processAvgWaves(root,tInfo.Fs, cycles, figInfo,chOrd);
pdData = processPeakDists(awData); 
quiverData = processQuiver(tInfo,chTxt, figInfo);
%  corrPlot(tInfo,chTxt, figInfo)
%  thetaGreaterMeanPower(tInfo,figInfo)
  %plotAvgWaveImg(root, cycles, ref, figInfo,chOrd)
% plotRawWaves(root,tInfo.Fs, figInfo)

  subplotOne(awData, pdData, quiverData, figInfo)

 keyboard

 end

% make avereaged waves 
function [awData] = processAvgWaves(root, Fs, cycles, figInfo, chOrd)

CycleTs=root.b_lfp(1).ts(cycles); %finds theta cycles
epochSize = 0.100; %sets epoch size
fprintf('calculating based off epochs of ', epochSize, '\n');
Epochs = [CycleTs-epochSize CycleTs+epochSize]; % grabs epochs %changed this from .125
root.epoch=Epochs; 

% Iterates over each channel 
for I=1:length(chOrd)
    root.active_lfp=I; % Set the active channel to one of the electrodes
    EegSnips=root.lfp.signal; % Snip the signal
    Snips(I,:) = EegSnips; 
    MeanThetaWave(I,:)=nanmean(catPlus(3,EegSnips),3); % Average the channel data
end
waveData = MeanThetaWave;

% plotLFP(data, Fs, 2.5,[],0);
% title([figInfo.name, 'Mean Theta Wave'])
% 
% if figInfo.saveFig
%     plotName = 'meanWaveTrace';
%     printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
% end

awData.waveData = waveData;
awData.Fs = Fs;
keyboard;

[nElecs,tPts,nSets] = size(awData.waveData);
nR = floor(sqrt(nSets));
nC = ceil(nSets/nR);
t = (1:tPts)/awData.Fs;
disp(2.5)
lfp_ = awData.waveData / (-1 * 2.5 * rms(awData.waveData(:)));
offsets = repmat([1:nElecs]',1,tPts,nSets);
lfp_ = lfp_ + offsets;
plot(t,lfp_,'k'); axis ij 
keyboard;
end

function pdData = processPeakDists(awData)
% 
% 
% M = randi(99, 4, 8);
% for k1 = 1:size(M,1)
%     [pks,loc] = findpeaks(M(k1,:));
%     P{k1} = [pks; loc];
% end
% 
% 
% plotLFP(data, Fs, 2.5,[],0);
% 
% title([figInfo.name, 'pk diff'])
% if figInfo.saveFig
%     plotName = 'meanWaveTrace';
%     printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
% end

%if we wanted to smooth it
% for I=1:16
%     root.active_lfp=I;
%     EegSnips=root.lfp.signal;
%     smoothMTW(I,:)=smoothts(nanmean(catPlus(3,EegSnips),3), 'g');
% end
% smoothMTW = MeanThetaWave(1:2:end,:);

%imagesc(tsa)
%figure; scatter(linspace(1,16,16),tInfo.thetaShiftAngle(1:16))
%figure; imagesc(diag(tInfo.thetaShiftAngle, 1));
%title('peak dists')


% tmp = tInfo.thetaShiftAngle;
% %tmp = tInfo.thetaShiftAngle(6:2:12,6:2:12);
% d2p = [mean(diag(tmp,0)) mean(diag(tmp,-1)) mean(diag(tmp,-2)) mean(diag(tmp,-3))];
% subplot(2,2,1); plot([0:3],rad2deg(d2p),'o-'); xlim([-0.5 3.5])


%looks for max val index in waves
for i = 1:size(awData.waveData,1)
    [~, pkInd(i)] = max(awData.waveData(i,:));
end
%figure; plot(pkInd(2:2:end))

%could subtract these to normalize across sessions
pdData.pkInd = pkInd; 
end


function [quiverData] = processQuiver(tInfo,chOrdTxt,figInfo)
% plots the polar coherence between electrodes
[u,v] = pol2cart(tInfo.thetaShiftAngle,tInfo.thetaShiftRbar);
%figure; quiver(u(2:2:end,2:2:end),v(2:2:end,2:2:end)); axis ij;
%title([figInfo.name, chOrdTxt, ', Blank (1403)']);

% if figInfo.saveFig
%     plotName = 'quiver';
%     printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
% end
quiverData.u = u;
quiverData.v = v;

% calculate the average peak shift 
tmp = tInfo.thetaShiftAngle;
avgPkShft = nan(1,size(tmp,1));
for shft = 1:size(tmp,1)
    avgPkShft(shft) = mean(diag(tmp,-shft));
end

% calculate the slope of the average peak shift
aps = rad2deg(avgPkShft(1,1:3)); % !! this is hardcoded to only pick first 3 to avoid the NaN's
x = linspace(0,size(aps,2),size(aps,2)); 
p = polyfit(x,aps,1);
pv = polyval(x,p);

B = regress(aps', [x' ones(size(x'))]);
keyboard;



quiverData.avgPkShft = rad2deg(avgPkShft); %convert to degrees for plotting
quiverData.p = p;
quiverData.pv = pv;
end

function corrPlot(tInfo, chOrdTxt, figInfo)
%plots the coherence between each cell
R = corr(tInfo.theta_filt');
figure; imagesc(R);
FirstDegMCorr = mean(diag(R,1));
SecondDegMCorr = mean(diag(R,2)); 
title([figInfo.name, chOrdTxt, ', Blank (1403), 1st vs 2nd neighbor corr ratio = ',num2str(FirstDegMCorr/SecondDegMCorr)]);
if figInfo.saveFig
    plotName = 'corLin';
    printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
end

[u,v] = pol2cart(tInfo.thetaShiftAngle,tInfo.thetaShiftRbar);
circR = sqrt(u.^2 + v.^2); figure; imagesc(circR,[0 1]);
title([figInfo.name, 'corr in radians']);
if figInfo.saveFig
    plotName = 'corCirc';
    printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
end
end

function thetaGreaterMeanPower(tInfo, figInfo)
%highTheta = find(theta_amp(end,:)>(meanAmp+2*stdAmp));
tInfo.highTheta = find(tInfo.theta_amp(end,:)>(tInfo.thetaMeanAmp));
tInfo.highThetaEp = mat2cell(tInfo.highTheta, 1, diff([0 find([(diff(tInfo.highTheta) > 1) 1])]));
lengthEp = cellfun(@length,tInfo.highThetaEp);
inds_long = lengthEp>300;
tInfo.highThetaEp_long = tInfo.highThetaEp(inds_long);
 
for iT = 1:length(tInfo.highThetaEp_long)
  tInfo.thetaShift_hiTheta{iT} = circDiff(tInfo.theta_phase(1:1:end,tInfo.highThetaEp_long{iT}),1,'rad');
end
tInfo.thetaShiftMat_hiTheta = cell2mat(tInfo.thetaShift_hiTheta);

figure;
imagesc(tInfo.thetaShiftMat_hiTheta);

bins = [-pi:0.05:pi];
figure; 
for i = 1:size(tInfo.thetaShiftMat_hiTheta,1)
  subplot(size(tInfo.thetaShiftMat_hiTheta,1),1,i), hist(tInfo.thetaShiftMat_hiTheta(i,:),bins); xlim([-pi pi]); grid on; ylim([0 5000]);
  if i==1, title([figInfo.name, '8/17/2016 14:03 Recordings. Theta > mean power.']); end
 % ylabel([num2str(chans(i+1)) ' - ' num2str(chans(i))]);
end

if figInfo.saveFig
    plotName = 'phaseDiffs';
    printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
end
end

function plotRawWaves(root, Fs, figInfo)
%make the raw waves from lfp signal
for i=1:length(root.b_lfp)
    rawWaves(i, :) = root.b_lfp(i).signal;
end
%plot raw eegWaves
kPlotLFP(rawWaves(8:2:16,:),1:1200,Fs)

%plotLFP(rawWaves(8:2:16,:), Fs);
%figInfo.name

if figInfo.saveFig
    plotName = 'egRawLFP';
    printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
end

end

function plotAvgWaveImg(root, cycles, ref, figInfo, chOrd)
CycleTs=root.b_lfp(ref).ts(cycles);
Epochs = [CycleTs-0.125 CycleTs+0.125];
root.epoch=Epochs;

%for each channel, fetch the eeg snips, then average them
for I=1:length(chOrd)
    root.active_lfp=I; 
    EegSnips=root.lfp.signal;
    MeanThetaWave(I,:)=nanmean(catPlus(3,EegSnips),3);
end

% figure; imagesc(MeanThetaWave);
% title([figInfo.name, 'Mean Theta Wave; ref:' + string(ref)])

% if figInfo.saveFig
%     plotName = 'meanWaveHeatmap';
%     printFigure(gcf, [figInfo.fnameRoot, '_',plotName,'.',figInfo.fig_type],'imgType',figInfo.fig_type);
% end

figure;
polarplot(MeanThetaWave(1:end))
%compass(MeanThetaWave(7:end))
end

function subplotOne(awData,pdData,quiverData, figInfo)
figure;

%% Raw Waves Pannel
subplot(2,2,1);
[nElecs,tPts,nSets] = size(awData.waveData);

for I=1:nElecs
    awData.waveDataSmooth(I,:)=smoothts(awData.waveData(I), 'g');
end

% x = linspace(0,size(awData.waveDataSmooth(1,1:end),2),size(awData.waveDataSmooth(1,1:end),2)); 
% p = polyfit(x,awData.waveDataSmooth(1,1:end,1));
% pv = polyval(x,p);
% 
% 
% nR = floor(sqrt(nSets));
% nC = ceil(nSets/nR);
% t = (1:tPts)/awData.Fs;
% disp(2.5)
% lfp_ = awData.waveData / (-1 * 2.5 * rms(awData.waveData(:)));
% offsets = repmat([1:nElecs]',1,tPts,nSets);
% lfp_ = lfp_ + offsets;
% plot(t,lfp_,'k'); axis ij 


title('Raw data');
t = text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = t.FontSize;
t.FontSize = 12;

%% avg quiver Pannel
subplot(2,2,2); 
%plot(pdData.pkInd(2:2:end)); hold on; %This looks at time between peaks
aps = quiverData.avgPkShft;
plot(aps) %Better; looks at average across phase
%text(0,0,['slope:', quiverData.p(1)]) %doesn't work well with subplot
xlabel('Electrode') %!! Is this correct? 
ylabel('Degree') % !! same here
title(['Average Peak Shift | Slope =', num2str(quiverData.p(1))]) %!! picking arbitrary slope
t = text(0.02,0.98,'B','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = t.FontSize;
t.FontSize = 12;
%axis should be electrodes? maybe?


%% Averaged Waves Pannel
subplot(2,2,3);
[nElecs,tPts,nSets] = size(awData.waveData);
nR = floor(sqrt(nSets));
nC = ceil(nSets/nR);
t = (1:tPts)/awData.Fs;
disp(2.5)
lfp_ = awData.waveData / (-1 * 2.5 * rms(awData.waveData(:)));
offsets = repmat([1:nElecs]',1,tPts,nSets);
lfp_ = lfp_ + offsets;
plot(t,lfp_,'k'); axis ij 
keyboard;

%add the find peaks function here? %Smooth?




%Change axis to reflect proper channels
xlabel('Time') %!! Is this correct? or should it be phase?
ylabel('Channel') % !! What about the axis though...
title('Averaged Waves')
t = text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = t.FontSize;
t.FontSize = 12;

%% Quiver Pannel
subplot(2,2,4);
quiver(quiverData.u(1:1:end,1:1:end),quiverData.v(1:1:end,1:1:end)); axis ij; 
xlabel('Channels') % !! yeah 
ylabel('Channels') % !! same
title([figInfo.name, ', ',figInfo.chOrdTxt, ', Blank (1403)'])
%title([figInfo.name, ', ',figInfo.chOrdTxt, ', Blank (1403)']);
%Axis should change to reflect the proper channels
t = text(0.02,0.98,'D','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = t.FontSize;
t.FontSize = 12;

suptitle([figInfo.name, ' Data: ', figInfo.session])

end
