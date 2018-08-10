
% grab ronaldo's root object
load(fullfile('D:','Dropbox (NewmanLab)','ratsEphys','Ronaldo','WmazeNoTask_2018-06-22','2018-06-22_19-38-24','WmazeNoTask_2018-06-22_Ronaldo_root.mat'));

% grab the data out of it
for i =1:size(root.b_lfp,2)
  data(i,:) = root.b_lfp(i).signal;
end

% grab a sample
dsamp = data()




% "An oscilation is perdiction"


%autoregression? 






