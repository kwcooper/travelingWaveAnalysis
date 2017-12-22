function twCrossCorr(root)

%assumes that the average theta wave was precalculated and saved to root

atw = root.user_def.atw;

% Calculate the lag from the first channel, for all channels
for i = 1:size(atw,1)-1
  [acor(i,:),lag(i,:)] = xcorr(atw(1,:),atw(i+1,:));
end

[~,I] = max(abs(acor'));

lagDiff = lag(I);
timeDiff = lagDiff/root.user_def.lfp_fs;
%figure; plot(timeDiff');

%plot them all together
figure; plot(lag',acor')
a3 = gca;
a3.XTick = sort([-3000:1000:3000 lagDiff(1)]);

%make a nice staggered plot
h = figure;
[nElecs,tPts] = size(lag);
acor_ = acor / (-1 * 2.5 * rms(acor(:)));
offsets = repmat([1:nElecs]',1,tPts);
acor_ = acor_ + offsets;
plot(lag',acor_','k'); axis ij;
a3 = gca;
a3.XTick = sort([-3000:1000:3000 lagDiff(1)]);

title('CrossCorr From 1st Channel');
xlabel('lag Difference');
ylabel('Channel');
end