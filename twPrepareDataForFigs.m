function [root,tInfo] = twPrepareDataForFigs(sInd,sessions,chOrd,force)
% see if data has been pre-computed, if not, compute it!
% !! TD: Should rename this to be specific for each session, loading it in for
% each, and display the name
if force
    fprintf('Forcing data recalculation...\n')
end

if ~force && exist('twMakeFigs_workingData.mat','file')
    fprintf('Found precomputed data, loading it instead of recomputing it. . .');
    L = load('twMakeFigs_workingData.mat');
    root = L.root;
    tInfo = L.tInfo;
    fprintf('done.\n')
else
    fprintf('Computing basic data, this may take a min \nbut don''t worry, I''ll save it when I''m done for next time\n');
    dsFreq = 600; 
    nChan = length(chOrd);
    fsVid = 120;
 
    %%  Data Handeling
    if ~exist('root', 'var') % !! what if it does exist...
        fprintf('Collecting the data...\n'); 
        %[data, timestamps, info] = load_open_ephys_data_faster(filename, varargin)
        % tmp to take advantage of epoching feature for ___________
        tmpRoot = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:0.01:21.0028*60],'fs_video',120);
        %check if there is a .kwd file
        if exist('experiment1_100.raw.kwd','file')
            fprintf('I found a .kwd file...\n');
            tmpRoot.path_lfp = repmat({'experiment1_100.raw.kwd'},1,length(chOrd));
        %if not, then there should be continuous files, let's iterate through them
        else
            fprintf('I found a .kwd file...\n');
            for ch_ind = 1:length(sessions{sInd,5})
                path_lfp{ch_ind} = ['100_CH',num2str(sessions{sInd,5}(ch_ind)),'.continuous'];
                if ~exist(path_lfp{ch_ind},'file'), error('%s does not exist as lfp.\n',path_lfp{ch_ind}); end
            end
            tmpRoot.path_lfp = path_lfp;
        end
        tmpRoot = tmpRoot.LoadLFP(1,'downsample',dsFreq,'chOrd',chOrd(1));
         
        % this is the actual CMB object we'll create, use, and return
        % !! TD I want to add notes with the session info to the object...
        fprintf('Creating root object... \n')
        root = CMBHOME.Session('name','experiment1','epoch',[-inf inf],'b_ts',[0:1/fsVid:max(tmpRoot.b_lfp(1).ts)],'fs_video',fsVid);
        root.path_lfp = tmpRoot.path_lfp;
        root.user_def.sessionInfo = sessions(sInd,:);
        %save experiment1.mat root
        
        fprintf('Loading lfp... \n')
        tic
        root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
        root.active_lfp = nChan;
        toc; fprintf('(To load the LFP) \n');
        
    end
     
    % need to fix this somehow -- 11/10 this probably isn't necessary anymore
    % %if isempty(root.b_lfp(1))
    %     fprintf('Loading lfp... \n')
    %     root = root.LoadLFP(1:nChan,'downsample',dsFreq,'chOrd',chOrd);
    %     root.active_lfp = nChan;
    % %end
     
    % looks for points 5 stdevs from the mean and then removes them
    % !! Needs some love (add the new function here?)
    root = rmDataBlips(root);
    
    
     
    %iterate through the channels and create a struct for the data
    % ContinuizeEpochs: Converts all cell arrays to arrays 
    chans = [1:nChan];
    dataDS = nan(length(chans),length(CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal)));
    for i = 1:length(chans)
        root.active_lfp = chans(i);
        dataDS(i,:) = CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal);
    end
     
    %% Extract Theta
    
    % We seem to have stopped using root at this time... 
    
    % create a tInfo struct which holds computed theta information 
    tInfo = {};
    tInfo.session = sessions(sInd,:);
    tInfo.Fs = root.lfp.fs;
    tInfo.Wn_theta = [6/(tInfo.Fs/2) 10/(tInfo.Fs/2)]; %filtering cutoff
    [tInfo.btheta,tInfo.atheta] = butter(3,tInfo.Wn_theta); %find filtering coef with butterworth filter
    tInfo.signal = dataDS;
    
    tic
    % extract theta with phase and power
    fprintf('theta extraction \n')
    tInfo.theta_filt = nan(size(dataDS));
    tInfo.theta_phase =  nan(size(dataDS));
    tInfo.theta_amp =  nan(size(dataDS));
    
    for iD =  1:size(dataDS,1)
        tInfo.theta_filt(iD,:) = filtfilt(tInfo.btheta,tInfo.atheta,dataDS(iD,:)); %filter the data
        tInfo.theta_phase(iD,:) = atan2(imag(hilbert(tInfo.theta_filt(iD,:))), tInfo.theta_filt(iD,:));
        tInfo.theta_amp(iD,:) = abs(hilbert(tInfo.theta_filt(iD,:)));
    end
    toc
    %% Focus on epochs with high theta amplitude
    
    % High theta currently defined as theta amplitude being greater than the mean amp - std amp
    
    % high amplitude will be based on greater than 2 std above the mean
    % so, compute session-wide mean and standard deviation of theta power
    % !! using the last channel (update this if another channel makes more sense)
    
     
    ch = size(tInfo.theta_amp,1);
    tInfo.thetaMeanAmp = mean(tInfo.theta_amp(ch,:));
    tInfo.stdAmp = std(tInfo.theta_amp(ch,:));
     
    % find the high theta
    %highTheta = find(theta_amp(ch,:)>(meanAmp+2*stdAmp));
    %highTheta = find(theta_amp(end,:)>(meanAmp));
    fprintf('finding high theta epochs\n');
    tInfo.highTheta = find(tInfo.theta_amp(end,:)>(tInfo.thetaMeanAmp-tInfo.stdAmp));
    tInfo.highThetaEp = mat2cell(tInfo.highTheta, 1, diff([0 find([(diff(tInfo.highTheta) > 1) 1])]));
    lengthEp = cellfun(@length,tInfo.highThetaEp);
    inds_long = lengthEp>300;
    tInfo.highThetaEp_long = tInfo.highThetaEp(inds_long);
    
    %%
    % find the theta shift for high theta channels
    fprintf('finding theta shift\n');
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
    
    
    %% new theta Extraction
    fprintf(' \n')
    fprintf('Extracting theta cycles... (new extraction)\n');
    metho = 'hilbert';
    fprintf(['useing '  metho '\n']);
    root.epoch=[-inf,inf];
    ref = sessions{sInd,6};
    band = [6,10];
    [thetaPhs,~,~] = extractThetaPhase(tInfo.signal(ref,:),tInfo.Fs,metho,band);
    [cycles,~] = parseThetaCycles(thetaPhs,tInfo.Fs,band); 
    %inds = find(cycles);
    
    tInfo.cycles = cycles;
    
    fprintf('Basing calculations off of %2.2f s of data\n',size(tInfo.thetaShiftMat,3)/tInfo.Fs);
    [tInfo.thetaShiftAngle,tInfo.thetaShiftRbar] = circmean(tInfo.thetaShiftMat,3);
     
    % if it needed to be computed, save it out
    % Should probably add the specific name to the session
    save('twMakeFigs_workingData.mat','root','tInfo','-v7.3');
    fprintf('Saving precomputed data so next time is a breeze\n');
end
end