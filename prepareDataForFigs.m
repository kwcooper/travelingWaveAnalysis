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