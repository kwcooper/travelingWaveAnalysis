
% code to compute slope figure from ratComp
slopes = [];
figure; hold on;
for i = 1:size(ratsComp.data, 1)
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.dataShift,'o');
  hold on;
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.y);
  slopes(i,:) = ratsComp.data{i,3}.B(1);
end

%%
% code to compute novel vs familair slopes
% assumes both in ratcomp
%TD need to make this distance
mNOV = []; mFAM = [];
figure;
for i = 1:size(ratsComp.data, 1)
  %plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.dataShift,'o');
  hold on;
  if i < 6
    c = 'r';
    mNOV(i,:) = ratsComp.data{i,3}.B(1);
  else
    c = 'b';
    mFAM(i-5,:) = ratsComp.data{i,3}.B(1);
  end
  
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.y, c);
end
title('rat theta shift: NOVEL (red) vs FAM (blue)')
% TD 

% plot individual slope vals
figure; hold on;
plot(mNOV', 'b')
plot(mFAM', 'r')
title('novel (blue) vs fam (red) slope vals')

% plot bar chart of the slope values
figure;
hb = bar([abs(mean(mNOV)),abs(mean(mFAM))]);
title('average theta slope vals')
set(gca,'xticklabel',{'novel','familiar'})
hb(1).FaceColor = 'b'; hb(1).FaceColor = 'r';



%%


% code to  compute SCP vs SAL
% assumes both in ratcomp
%TD need to make this distance
SCP = []; SAL = [];
figure;
for i = 1:size(ratsComp.data, 1)
  %plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.dataShift,'o');
  hold on;
  if i <= size(ratsComp.data, 1) / 2
    c = 'b';
    
    SCP(end+1) = ratsComp.data{i,3}.B(1);
  else
    c = 'r';
    SAL(end+1) = ratsComp.data{i,3}.B(1);
  end
  
  plot(ratsComp.data{i,3}.x, ratsComp.data{i,3}.dataShift,'o');
  plot(ratsComp.data{i,3}.x, ratsComp.data{i,3}.y, c);

end
title('rat theta shift: SCP (blue) vs SAL (red)')
xlabel('Channel'); % soon to be distance
ylabel('Offset from first channel');

% Need more data
% % plot individual slope vals
% figure; hold on;
% plot(SCP', 'b')
% plot(SAL', 'r')
% title('SCP (blue) vs SAL (red) slope vals')

% plot bar chart of the slope values
figure;
hb = bar([abs(mean(SCP)),abs(mean(SAL))]);
title('average theta slope vals')
set(gca,'xticklabel',{'SCP','SAL'})
hb(1).FaceColor = 'b'; %hb(2).FaceColor = 'r';






