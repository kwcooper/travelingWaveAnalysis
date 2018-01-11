function [h,figData] = twGrabRawData(D, ind1, ind2)

h = figure;
t = linspace(ind1, ind2, ind2-ind1);
lfp = D(:,ind1:ind2-1);
[nElecs, tPts] = size(lfp);
offsets = repmat([1:nElecs]',1,tPts);
lfp = (lfp/max(abs(lfp(:))));
lfpO = lfp + offsets;
plot(t,lfpO,'k'); axis ij;

xlabel('Time') %!! Is this correct? or should it be phase?
ylabel('Channel') % !! What about the axis though...
title('Raw Waves')
ltr = text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top');
s = ltr.FontSize;
ltr.FontSize = 12;
grid on

figData.rawWaves.lfpO = lfpO;
figData.rawWaves.t = t;

end