 function twMakeSubFigs(force)
%figPath = fullfile(dropboxPath,SLpath,'..','figures');

%this is keiland's version

%update these values or copy and paste them


% Rat = 'Tio';
% Session = '170703_1251_CircleTrack';
% Recording = '2017-07-03_13-08-30';
% tioChOrd = [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54]; 
% chOrdTxt = 'Probe order';
% chOrd = tioChOrd;


%2017-08-10_19-14-01 - weird one
%2017-08-10_11-43-00 - fine
%2017-08-11_12-54-28 - looks fine

sessions = {...
'RioNovelty',  '2017-07-27_16-00',         '2017-08-09_16-57-57',   '',    [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_19-14-01',   'sml',    [39 38 60 57], 1; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_19-14-01',   'all',    [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-10_CircleTrack',   '2017-08-10_11-43-00',   '',    [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-11_CircleTrack',   '2017-08-11_12-54-28',   '',    [43 44 46 45 40 39 37 38 59 60 58 57 52 51 53 54], 8; ...
'Rio',         '2017-08-20_CircleTrack',   '2017-08-20_12-41-36',   '',    [11 12 14 13 8 7 5 6 27 28 26 25 20 19 21 22], 8;...
'Rio',         '2017-08-22_CircleTrack',   '2017-08-22_14-01-24',   '',    [11 12 14 13 8 7 5 6 27 28 26 25 20 19 21 22], 8;...
};
%s

sInd = 2;

Rat =  sessions{sInd,1};
Session =  sessions{sInd,2};
Recording =  sessions{sInd,3};
chTxt = sessions{sInd,4};
chOrd = sessions{sInd,5}; 
ref = sessions{sInd,6};

workingDir = fullfile(ratLibPath,Rat,Session,Recording); cd(workingDir);

if ~exist('force') || isempty(force), force = 0; end
if force
    fprintf('Forcing data recalculation...\n')
end
[root,tInfo] = prepareDataForFigs(sInd,sessions,chOrd,force);


%%
%new theta Extraction
disp("Extracting theta cycles...");
metho = 'hilbert';
disp("useing " + metho)
root.epoch=[-inf,inf];
band = [6,10];
[thetaPhs,~,~] = extractThetaPhase(tInfo.signal(ref,:),tInfo.Fs,metho,band);
[cycles,~] = parseThetaCycles(thetaPhs,tInfo.Fs,band);

inds = find(cycles);
%%
%Figures
%keyboard

%Set up the figure info
figInfo = {};
figInfo.name = Rat;
figInfo.session = Session;
figInfo.recording = Recording;
figInfo.saveFig = 1;
figInfo.fnameRoot = fullfile('figs',[figInfo.name, '_', figInfo.session]);
figInfo.fig_type = 'png'; % options = {'png','ps','pdf'}

figInfo.chOrdTxt = chTxt;
figInfo.ref = ref;

if ~exist('figs', 'dir')
    mkdir figs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  THE BUSINESS END OF THE SCRIPT %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 awData = processAvgWaves(root,tInfo.Fs, cycles, figInfo,chOrd);
 pdData = processPeakDists(awData);
  quiverData = processQuiver(tInfo,chTxt, figInfo);
%  corrPlot(tInfo,chTxt, figInfo)
%  thetaGreaterMeanPower(tInfo,figInfo)
  plotAvgWaveImg(root, cycles, ref, figInfo,chOrd)
% plotRawWaves(root,tInfo.Fs, figInfo)
  subplotOne(awData, pdData, quiverData, figInfo)

keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

tmp = tInfo.thetaShiftAngle;
avgPkShft = nan(1,size(tmp,1));
for shft = 1:size(tmp,1)
    avgPkShft(shft) = mean(diag(tmp,-shft));
end
quiverData.avgPkShft = avgPkShft;
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

function [awData] = processAvgWaves(root, Fs, cycles, figInfo, chOrd)

CycleTs=root.b_lfp(1).ts(cycles);
epochSize = .100;
Epochs = [CycleTs-0.100 CycleTs+0.100]; %changed this from .125
fprintf('calculating based of epochs of ', epochSize, '\n');
root.epoch=Epochs;

for I=1:length(chOrd)
    root.active_lfp=I;
    EegSnips=root.lfp.signal;
    MeanThetaWave(I,:)=nanmean(catPlus(3,EegSnips),3);
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
end

function subplotOne(awData,pdData ,quiverData, figInfo)
figure;
%Raw Waves Pannel
subplot(2,2,1);

% Pannel
subplot(2,2,2); 
%plot(pdData.pkInd(2:2:end)); hold on;
plot(rad2deg(quiverData.avgPkShft)) %This is a better way of computing it
title('Average Peak Shift')


%Averaged Waves Pannel
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

%Quiver Pannel
subplot(2,2,4);
quiver(quiverData.u(1:1:end,1:1:end),quiverData.v(1:1:end,1:1:end)); axis ij;  
title([figInfo.name, figInfo.chOrdTxt, ', Blank (1403)']);

end

function [root,tInfo] = prepareDataForFigs(sInd,sessions,chOrd,force)
% see if data has been pre-computed, if not, compute it!
if ~force && exist('twMakeFigs_workingData.mat','file')
    fprintf('Found precomputed data, loading it instead of recomputing it. . .');
    L = load('twMakeFigs_workingData.mat');
    root = L.root;
    tInfo = L.tInfo;
    fprintf('done.\n')
else
    fprintf('Computing basic data, this may take a min \nbut don''t worry, I''ll save it when I''m done\n');
     
     
    dsFreq = 600; %data sample frequency? what is this?
    nChan = length(chOrd);
    fsVid = 120;
 
    %%  Data Handeling
    if ~exist('root', 'var')
         
        %[data, timestamps, info] = load_open_ephys_data_faster(filename, varargin)
        % tmp to take advantage of epoching
        tmpRoot = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:0.01:21.0028*60],'fs_video',120);
        if exist('experiment1_100.raw.kwd','file')
            tmpRoot.path_lfp = repmat({'experiment1_100.raw.kwd'},1,length(chOrd));
        else
            for ch_ind = 1:length(sessions{sInd,5})
                path_lfp{ch_ind} = ['100_CH',num2str(sessions{sInd,5}(ch_ind)),'.continuous'];
                if ~exist(path_lfp{ch_ind},'file'), error('%s does not exist as lfp.\n',path_lfp{ch_ind}); end
            end
            tmpRoot.path_lfp = path_lfp;
        end
        tmpRoot = tmpRoot.LoadLFP(1,'downsample',dsFreq,'chOrd',chOrd(1));
         
         
        fprintf('Creating root object... \n')
        root = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:1/fsVid:max(tmpRoot.b_lfp(1).ts)],'fs_video',fsVid);
        root.path_lfp = tmpRoot.path_lfp;
        root.user_def.sessionInfo = sessions(sInd,:);
        %save experiment1.mat root
         
        fprintf('Loading lfp... \n')
        root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
        root.active_lfp = nChan;
    end
     
    %need to fix this somehow
    % %if isempty(root.b_lfp(1))
    %     fprintf('Loading lfp... \n')
    %     root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
    %     root.active_lfp = nChan;
    % %end
     
    root = rmDataBlips(root);
     
    %iterate through the channels
    chans = [1:nChan];
    dataDS = nan(length(chans),length(CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal)));
    for i = 1:length(chans)
        root.active_lfp = chans(i);
        dataDS(i,:) = CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal);
    end
     
    %% Extract Theta
     
    tInfo = {};
    tInfo.session = sessions(sInd,:);
    tInfo.Fs = root.lfp.fs;
    tInfo.Wn_theta = [6/(tInfo.Fs/2) 10/(tInfo.Fs/2)];
    [tInfo.btheta,tInfo.atheta] = butter(3,tInfo.Wn_theta);
    tInfo.signal = dataDS;
     
    % extract theta with phase and power
    fprintf('theta extraction \n')
     
    tInfo.theta_filt = nan(size(dataDS));
    tInfo.theta_phase =  nan(size(dataDS));
    tInfo.theta_amp =  nan(size(dataDS));
    % why is this iterating through each data point, and not taking all of the data? k
    for iD =  1:size(dataDS,1)
        tInfo.theta_filt(iD,:) = filtfilt(tInfo.btheta,tInfo.atheta,dataDS(iD,:)); %filter the data
        tInfo.theta_phase(iD,:) = atan2(imag(hilbert(tInfo.theta_filt(iD,:))), tInfo.theta_filt(iD,:));
        tInfo.theta_amp(iD,:) = abs(hilbert(tInfo.theta_filt(iD,:)));
    end
     
    %%
    % Focus on epochs of data with high theta amplitude
    % high amplitude will be based on greater than 2 std above the mean
    % so, compute session-wide mean and standard deviation of theta power
    % using the last channel (update this if another channel makes more sense)
     
    % I might need to change this k
    % Also, why not all theta? quick n dirty or the whole project?
    ch = size(tInfo.theta_amp,1);
    tInfo.thetaMeanAmp = mean(tInfo.theta_amp(ch,:));
    stdAmp = std(tInfo.theta_amp(ch,:));
     
    % now find the high theta
    %highTheta = find(theta_amp(ch,:)>(meanAmp+2*stdAmp));
    %highTheta = find(theta_amp(end,:)>(meanAmp));
    tInfo.highTheta = find(tInfo.theta_amp(end,:)>(tInfo.thetaMeanAmp-stdAmp));
    tInfo.highThetaEp = mat2cell(tInfo.highTheta, 1, diff([0 find([(diff(tInfo.highTheta) > 1) 1])]));
    lengthEp = cellfun(@length,tInfo.highThetaEp);
    inds_long = lengthEp>300;
    tInfo.highThetaEp_long = tInfo.highThetaEp(inds_long);
     
    %find the theta shift for high theta channels
    clear thetaShift
    for iT = 1:length(tInfo.highThetaEp_long)
        for iC1 = 1:length(chans) 
            for iC2 = 1:length(chans)
                tInfo.thetaShift{iT}(iC1,iC2,:) = circDiff([tInfo.theta_phase(iC1,tInfo.highThetaEp_long{iT})', ...
                    tInfo.theta_phase(iC2,tInfo.highThetaEp_long{iT})'],2,'rad'); % changed orientation of data
            end
        end
    end
     
    tInfo.thetaShiftMat = cat(3,tInfo.thetaShift{:}); % used cat instead of cell2mat
    
    fprintf('Basing calculations off of %2.2f s of data\n',size(tInfo.thetaShiftMat,3)/tInfo.Fs);
    [tInfo.thetaShiftAngle,tInfo.thetaShiftRbar] = circmean(tInfo.thetaShiftMat,3);
     
    % if it needed to be computed, save it out
    save('twMakeFigs_workingData.mat','root','tInfo','-v7.3');
    fprintf('Saving precomputed data to save time next time\n');
end
end