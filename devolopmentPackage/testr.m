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

% make timestamps
cumsum(repmat(1/fs,1,length(ephysRef)))


