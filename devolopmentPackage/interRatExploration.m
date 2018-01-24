function interRatExploration(ratRegress)

% rats = ["Rio", "Tio", "Romo", "Roble"];
% for i = 1:size(rats, 2)
%   fprintf(rats(1,i))
%   ratRegress.(str2char(rats(1,i)));
% end


% cleaner way of doing it
%Gather the slopes from all the rats
fields = fieldnames(ratRegress);
slopes = nan(1, numel(fields));
for i = 1:numel(fields)
  slopes(1,i) = ratRegress.(fields{i}).pd.B(1);
end

fprintf(['The mean of the slopes is: ', num2str(mean(slopes)), ' (n = ', num2str(numel(fields)),')\n']);

histogram(slopes)
title(['Slope of wave offset, mean: ', num2str(mean(slopes)), ' (n = ', num2str(numel(fields)),')']);


end