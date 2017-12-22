%M = randi(99, 4, 30);

% M = MeanThetaWave;
% M = MeanThetaWave(1:2:end,:);
% %M = nanmean(M);
% for k1 = 1:size(M,1)
%     
%     [pks,loc] = findpeaks(nanmean(M(k1,:)));
%     P{k1} = [pks; loc];
% end
% 
% hold on;
% %plot(M(1,:));
% 
% for ii = 1:size(M,1)
%     
%     plot(P{ii}(2,:));
%     
% end


plot(diag(tsa,1))

