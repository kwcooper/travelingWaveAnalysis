function interRatExploration(ratsComp, plt)


% make the slope offset figure
slopes = [];
figure; hold on;
for i = 1:size(ratsComp.data, 1)
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.dataShift,'o');
  hold on;
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.y);
  slopes(i,:) = ratsComp.data{i,3}.B(1);
end
title('Rat Slopes');
xlabel('Distance (Channel)');
ylabel('Offset from first channel (Radians)');


% make the slope offset figure with labels
slopes = [];
figure; hold on;
c = ['r', 'g', 'y', 'c', 'b', 'm'];
for i = 1:size(ratsComp.data, 1) % just look at the first half
  %plot points
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.dataShift,[c(i) '.']); hold on;
  %plot line
  plot(ratsComp.data{i,3}.x,ratsComp.data{i,3}.y, c(i));
  slopes(i,:) = ratsComp.data{i,3}.B(1);
  legendInfo{i} = [ratsComp.data{i,1}]; 
end
title('Rat Slopes');
xlabel('Distance (Offset from first Channel)'); ylabel('Offset from first channel (Radians)');

newLeg = cell(1,size(legendInfo,2)*2); j = 0;
for i = 1:size(newLeg,2)
  if mod(i,2) == 0; j = j + 1; newLeg{i} = legendInfo{j}; else newLeg{i} = ' '; end
end
    
legend(newLeg) %[findobj(gca,'Type','line')], [legendInfo legendInfo]

A = [6;2;2;1;1;0]
B=[A'; nan(1,numel(A))];
B=B(:)

figure; 
P = (0.1:0.1:4)';
for i = 2:11
    P(:,i) = P(:,1).^2 + (i-1);
    legendInfo{i-1} = ['X = ', ratsComp.data{i-1,1}]; 
end
plot(P(:,1),P(:,2:end))
legend(legendInfo)





% cleaner way of doing it
% Gather the slopes from all the rats
% assumes that each rat has the same number of recordings
namefields = fieldnames(ratsComp.ratRegress);
%slopes = nan(1, numel(namefields));
for i = 1:numel(namefields)
  namefieldsfields(i,:) = fieldnames(ratsComp.ratRegress.(namefields{i}))';
end

for k = 1:numel(namefields)
  for i = 1:size(namefieldsfields,2)
      slopes(k,i) = ratsComp.ratRegress.(namefields{k}).(char(namefieldsfields(k,i))).pd.B(1);
  end
end


%slopes(1,i) = ratRegress.(namefields{i}).pd.B(1); 
slopesPrime = reshape(slopes, [1,numel(slopes)]);

fprintf(['The mean of the slopes is: ', num2str(mean(slopesPrime)), ' (n = ', num2str(numel(namefields)),')\n']);

if plt
  histogram(slopesPrime)
  title(['Slope of wave offset, mean: ', num2str(mean(slopesPrime)), ' (n = ', num2str(numel(namefields)),')']);
end


% slope of offset scatter plot
y = [1 2 3 4];
%figure; scatter(slopesPrime, y)

numSess = 2;
numAnimal = 4;
slopesP = reshape(slopesPrime, [numSess,numAnimal])';
slopesP = slopesP * -1; % bring them back positive
sPAvg = mean(slopesP');
sPError = std(slopesP') / sqrt(numSess);

y = [1 2 3 4];
sze = []; c = [1,2,3,5];
figure; scatter(y,sPAvg,sze,c,'filled');
axis([0,4.5,0,.4]); hold on;
err = sPError;
errorbar(y, x, err, 'LineStyle', 'none');
set(gca,'xtick',0:4);
title('Average Theta Shift Slopes');
xlabel('Animal'); ylabel('Slope');
end


%code graveyard

% how to add a label to a line in a foor loop
% figure; 
% P = (0.1:0.1:4)';
% for i = 2:11
%     P(:,i) = P(:,1).^2 + (i-1);
%     legendInfo{i-1} = ['X = ', ratsComp.data{i-1,1}]; 
% end
% plot(P(:,1),P(:,2:end))
% legend(legendInfo)

