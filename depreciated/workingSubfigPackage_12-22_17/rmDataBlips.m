function root = rmDataBlips(root)
% RMDATABLIPS function to find and remove the bad data points
%
% currently uses z score to find shitty inds... 5 std dev's ?
% (this was quickly written by K on 11/10, needs to be more accurately updated)

  badInds = abs(zscore(root.b_lfp(root.active_lfp).signal))>5;
  sum(badInds);
  badts = root.b_lfp(root.active_lfp).ts(badInds);
  badepchs = [badts-0.5,badts+0.5];
  root.epoch = badepchs;
  root = CMBHOME.Utils.MergeEpochs(root);
  badepchs = root.epoch';
  % goodepchs = [0 badepchs(1,1)];
  % goodepchs = [goodepchs; badepchs(1:end-1,1) badepchs(2:end,
  goodepchs = reshape([0 badepchs(:)' root.b_ts(end)],2,[])';
  goodepchs(goodepchs(:,2)-goodepchs(:,1)<=0,:) = [];
  root.epoch = goodepchs;
  %figure;plot(CMBHOME.Utils.ContinuizeEpochs(root.lfp.signal))
end