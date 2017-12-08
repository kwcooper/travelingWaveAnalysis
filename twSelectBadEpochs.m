function [xCords, yCords] = twSelectBadEpochs(root, chOrd)
% TWSELECTBADEPOCHS Gui ability to select the coordinates of the bad epochs
% in a recording session. 
% Requires root object.

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
    
    figure;
    plot(signal); hold on;
    plot(find(inds2cut),signal(inds2cut), '.')
    
end
end