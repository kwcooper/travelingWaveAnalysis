function interRatExploration(ratRegress, plt)

%to do:
% restructure the ratStruct by session for multiple per rat!

% rats = ["Rio", "Tio", "Romo", "Roble"];
% for i = 1:size(rats, 2)
%   fprintf(rats(1,i))
%   ratRegress.(str2char(rats(1,i)));
% end


% cleaner way of doing it
% Gather the slopes from all the rats
% assumes that each rat has the same number of recordings
namefields = fieldnames(ratRegress);
slopes = nan(1, numel(namefields));
for i = 1:numel(namefields)
  namefieldsfields(i,:) = fieldnames(ratRegress.(namefields{i}));
end


slopes(1,i) = ratRegress.(namefields{i}).pd.B(1);

fprintf(['The mean of the slopes is: ', num2str(mean(slopes)), ' (n = ', num2str(numel(namefields)),')\n']);

if plt
  histogram(slopes)
  title(['Slope of wave offset, mean: ', num2str(mean(slopes)), ' (n = ', num2str(numel(namefields)),')']);
end

end