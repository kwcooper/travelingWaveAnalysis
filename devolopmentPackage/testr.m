%sessionsList
%sessions{1,1}

sessTst = {...
'Tio',	       '170713_1325_CircleTrack_good',	  '2017-07-13_13-41-41',	'Familiar',   [53 52 58 59 37],           1, [];... %1
'Tio',	       '170717_1824_CircleTrack',	        '2017-07-17_18-30-47',  'Novel',      [53 52 58 59 37],           1, [];... 
'Rio',	       '2017-08-09_CircleTrack',	        '2017-08-09_13-33-51',	'Familiar',   [45 39 38 60 57],           1, [];...
};

data = [1,2,3,4,6,5,6,7,6,5,4,5,6,7,8,7,5,6,5,4,5,6,5,6,5,6,5,6,7,7,7,7,8,7,8];
data2 = sin(data);

a = {data, sessTst{2, 2}};
a = [a; {data2, sessTst{3, 2}}];

a = fliplr(a);

b = {};
for i = 1:size(sessTst, 1)
  b = [b; {sessTst{i,3}, i * data}];
end

b = [b; {data, data}];
b

% code to compute slope figure from ratComp
hold on;
slopes = [];
figure;
for i = 1:size(ratsComp.data, 1)
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.dataShift,'o');
  hold on;
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.y);
  slopes(i,:) = ratsComp.data{i,3}.B(1);
end

% code to compute novel vs familair slopes 
mNOV = [];
mFAM = [];
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
hb(1).FaceColor = 'b';
hb(1).FaceColor = 'r';





