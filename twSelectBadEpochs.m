function [xCords, yCords] = twSelectBadEpochs(root, chOrd)
% TWSELECTBADEPOCHS Use to find bad epochs
% Can use manual mode with GUI or automated detection
% Requires root object and channelOrder
%
% To Do: automated portion needs handling for multiple channel bad epoch
% detection

%%
manual = 0;
if manual
    one = 1;
    
    if one
        %This grabs all of the selected lfp
        fprintf("gathering all the lfp\n")
        for I=1:length(chOrd)
            root.active_lfp=I;
            signal = root.lfp.signal;
            lfp(I,:) = signal;
        end
        figure; plot(lfp(1:4,1:100:end)');
    else
        %this grabs one channel
        lfp = root.lfp.signal;
        figure; plot(lfp);
    end
    
    %now let's plot and select the epochs
    fprintf("\nSelect the bad epochs of the data\n");
    fprintf("Select only two points, one before the bad data and one after\n")
    fprintf("Press return when finished\n")
    
    %call the gui input
    [xCords,yCords]=ginput;
    
    %housekeeping
    fprintf([num2str(size(xCords,1)) " points selected."])
    
    
else
    %% Automated function
    % code adapted from cleanData_Intan from e 171208
    % td take the average
    signal = root.lfp.signal; % !! Needs handeling for multiple channels (loop over lfp array above?)
    fs = root.lfp.fs;
    
    
    % examine the signal quality itself to find noise artifacts
    badInds = false(length(signal),1);
    buf = ceil(fs/2); % 500ms buffer around each artifact
    
    % drop saturated or flatlined data
    t0 = [0-1e-4 0+1e-4]; % tresholds for no dV/dt
    badInds_ = CMBHOME.Utils.ThresholdBandDetect(diff(signal),t0(1),t0(2),1,round(0.030 * fs)); % flatline data
    for ep = 1:size(badInds_,1)
        badInds(max(badInds_(ep,1)-buf,0):min(badInds_(ep,2)+buf,length(badInds))) = true;
    end
    
    % drop high amplitude events
    % high amplitude = 2.5 x RMS
    sigRMS = norm(signal)/sqrt(length(signal)); % adapted from Hazem Saliba Baqaen. Hazem@brown.edu
    badInds_ = [];
    badInds_ = [badInds_; CMBHOME.Utils.OverThresholdDetect(signal,2.5*sigRMS,1,round(0.030 * fs))];% upward deflections
    badInds_ = [badInds_; CMBHOME.Utils.UnderThresholdDetect(signal,-2.5*sigRMS,1,round(0.030 * fs))];% downward deflections
    for ep = 1:size(badInds_,1)
        badInds(max(badInds_(ep,1)-buf,1):min(badInds_(ep,2)+buf,length(badInds))) = true;
    end
    
    %%
    % plot what will be removed
    figure;
    plot(signal); hold on;
    plot(find(badInds),signal(badInds), '.')
    
end
end