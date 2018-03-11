function plotOffsets(data)
% PLOTOFFSETS: plots the data in a more viewable manner
% pass the data as a chan X time matrix


[tPts,nElecs] = size(data);
% double check if the orientation is chill
if tPts < nElecs
  data = data';
  [tPts,nElecs] = size(data);
end

t = linspace(0,1,tPts);
lfp_ = data' / (-1 * 2.5 * rms(data(:))); % normalize
offsets = repmat([1:nElecs]',1,tPts);
lfp_ = lfp_ + offsets;

% (!) need to flip the data to match the slopes
h = figure;
plot(t,lfp_,'k'); axis ij;
%Change axis to reflect proper channels

grid on

end














% xlabel('Time') %!! Is this correct? or should it be phase?
% ylabel('Channel') % !! What about the axis though...
% title(['Averaged Waves: ', figData.ratInfo.name])
% % ltr = text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top');
% % s = ltr.FontSize;
% % ltr.FontSize = 12;